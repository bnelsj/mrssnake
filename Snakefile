import os
import sys
import pysam
from subprocess import CalledProcessError

import pandas as pd

SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)

shell.executable("/bin/bash")
shell.prefix("source %s/config.sh; set -euo pipefail; " % SNAKEMAKE_DIR)

if config == {}:
    configfile: "config.yaml"

MANIFEST = config["manifest"]
REFERENCE = config["reference"]
MASKED_REF = config[REFERENCE]["masked_ref"]
CONTIGS_FILE = config[REFERENCE]["contigs"]

BAM_PARTITIONS = config["bam_partitions"]
UNMAPPED_PARTITIONS = config["unmapped_partitions"]
if UNMAPPED_PARTITIONS == -1:
    UNMAPPED_PARTITIONS = max(BAM_PARTITIONS // 100, 1)
MAX_BP = config["max_bp_in_mem"]

MAX_EDIST = config["max_edist"]

TMPDIR = config["tmpdir"]
CLEAN_TEMP_FILES = config["clean_temp_files"]

if not os.path.exists("log"):
    os.makedirs("log")

CONTIGS = {}

with open(CONTIGS_FILE, "r") as reader:
    for line in reader:
        contig, size = line.rstrip().split()[0:2]
        CONTIGS[contig] = int(size)

SAMPLES = pd.read_table(MANIFEST)
SAMPLES.index = SAMPLES.sn

def get_sparse_matrices_from_sample(wildcards):
    return ["region_matrices/%s/%s.%d_%d" % (wildcards.sample, wildcards.sample, part, BAM_PARTITIONS) for part in range(BAM_PARTITIONS + UNMAPPED_PARTITIONS)]

localrules: all, get_headers, make_jobfile, clean, make_chunker

rule all:
    input:  expand("finished/{sample}.txt", sample = SAMPLES.sn)

rule clean:
    input: "mapping/{sample}/{sample}/wssd_out_file"
    output: touch("finished/{sample}.txt")
    priority: 50
    run:
        if CLEAN_TEMP_FILES:
            shell("rm region_matrices/{wildcards.sample}/*")

rule wssd_merge:
    input: wssd = expand("mapping/{{sample}}/{{sample}}/wssd_out_file.{contig}", contig = CONTIGS.keys()),
           shelve = expand("region_matrices/{{sample}}/{{sample}}.{part}_%d.h5" % BAM_PARTITIONS, part = range(BAM_PARTITIONS + UNMAPPED_PARTITIONS))
    output: "mapping/{sample}/{sample}/wssd_out_file"
    params: sge_opts="-l mfree=8G -l disk_free=10G -pe serial 1 -N merge_wssd -l h_rt=5:00:00 -soft -l gpfsstate=0"
    log: "log/wssd_merge/{sample}.txt"
    resources: mem=8
    priority: 40
    benchmark: "benchmarks/wssd_merge/{sample}.txt"
    run:
        tempfile = "%s/%s.wssd_out_file" % (TMPDIR, wildcards.sample)
        shell("python3 merger.py {tempfile} --infiles {input.wssd} --wssd_merge --contigs_file {CONTIGS_FILE}")
        shell("rsync {tempfile} {output}")

rule merge_sparse_matrices:
    input: expand("region_matrices/{{sample}}/{{sample}}.{part}_%d.h5" % BAM_PARTITIONS, part = range(BAM_PARTITIONS + UNMAPPED_PARTITIONS))
    output: temp("mapping/{sample}/{sample}/wssd_out_file.{contig}")
    params: sge_opts = "-l mfree=8G -l disk_free=10G -pe serial 1 -N merge_sample -l h_rt=5:00:00 -soft -l gpfsstate=0"
    log: "log/merge/{sample}.{contig}.txt"
    resources: mem=8
    priority: 30
    benchmark: "benchmarks/merge/{sample}.{contig}.txt"
    run:
        infile_glob = os.path.commonprefix(input) + "*"
        tempfile = "%s/%s.wssd_out_file.%s" % (TMPDIR, wildcards.sample, wildcards.contig)
        shell('python3 merger.py {tempfile} --infile_glob "{infile_glob}" --contig {wildcards.contig}')
        shell("rsync {tempfile} {output}")

rule map_and_count:
    input: bam = lambda wildcards: SAMPLES.ix[SAMPLES.sn == wildcards.sample, "bam"], index = lambda wildcards: SAMPLES.ix[SAMPLES.sn == wildcards.sample, "index"]
    output: temp("region_matrices/{sample}/{sample}.{part}_%d.h5" % (BAM_PARTITIONS))
    params: sge_opts = "-l mfree=10G -N map_count -l h_rt=10:00:00 -soft -l gpfsstate=0"
    benchmark: "benchmarks/counter/{sample}/{sample}.{part}.%d.txt" % BAM_PARTITIONS
    resources: mem=10
    priority: 20
    log: "log/map/{sample}/{part}_%s.txt" % BAM_PARTITIONS
    run:
        masked_ref_name = os.path.basename(MASKED_REF)
        fifo = "%s/mrsfast_fifo" % TMPDIR
        if TMPDIR != "":
            local_index = "%s/%s" % (TMPDIR, os.path.basename(input.index[0]))
        else:
            local_index = input.index[0]
        mrsfast_ref_path = "/var/tmp/mrsfast_index/%s" % masked_ref_name
        rsync_opts = """rsync {0}.index /var/tmp/mrsfast_index/ --bwlimit 10000 --copy-links -p;
                        rsync {2} {3} --bwlimit 10000 --copy-links -p; 
                        if [[ ! -e {1} ]]; then touch {1}; fi; 
                        echo Finished rsync from {0} to {1} >> /dev/stderr; 
                        echo Finished rsync from {2} to {3} >> /dev/stderr; """.format(MASKED_REF, mrsfast_ref_path, input.index[0], local_index)

        read_counter_args = "--max_basepairs_in_mem %d --max_edist %s" % (MAX_BP, MAX_EDIST)
        shell("hostname; echo part: {wildcards.part} nparts: {BAM_PARTITIONS} unmapped parts: {UNMAPPED_PARTITIONS}; mkfifo {fifo}; {rsync_opts}")
        shell("bin/bam_chunker_cascade -b {input.bam} -i {local_index} -p {wildcards.part} -n {BAM_PARTITIONS} -u {UNMAPPED_PARTITIONS} 2>> /dev/stderr | "
            "mrsfast --search {mrsfast_ref_path} -n 0 -e {MAX_EDIST} --crop 36 --seq /dev/stdin -o {fifo} --disable-nohit >> /dev/stderr | "
            "python3 read_counter.py {fifo} {output} {CONTIGS_FILE} {read_counter_args}"
            )

rule check_bam_files:
    input: [bam for bam in SAMPLES.bam]
    output: touch("BAMS_READABLE")
    params: sge_opts = "-l h_rt=1:0:0"
    priority: 50
    run:
        for bamfile in input:
            try:
                test = pysam.AlignmentFile(bamfile)
            except ValueError as e:
                print("Error: could not open %s as bam.\n%s\n" % (bamfile, str(e)), file = sys.stderr)
                sys.exit(1)

rule check_index:
    input: MASKED_REF
    output: touch("MRSFASTULTRA_INDEXED"), temp(".mrsfast_index_test_output.txt")
    params: sge_opts = "-l h_rt=1:0:0"
    run:
        try:
            shell("mrsfast --search {input[0]} --seq dummy.fq > {output[1]}")
        except CalledProcessError as e:
            sys.exit("Reference %s was not indexed with the current version of mrsfastULTRA. Please reindex." % input[0])

rule make_chunker:
    input: "src/chunker_cascade.cpp", "Makefile"
    output: "bin/bam_chunker_cascade"
    params: sge_opts = "-l h_rt=1:0:0"
    shell:
        "make"
