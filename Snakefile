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
    UNMAPPED_PARTITIONS = max(BAM_PARTITIONS // 500, 1)
MAX_BP = config["max_bp_in_mem"]

TMPDIR = config["tmpdir"]
LIVE_MERGE = config["live_merge"]
CLEAN_TEMP_FILES = config["clean_temp_files"]

if LIVE_MERGE:
    ruleorder: merge_sparse_matrices_live > merge_sparse_matrices
else:
    ruleorder: merge_sparse_matrices > merge_sparse_matrices_live

if not os.path.exists("log"):
    os.makedirs("log")

CONTIGS = {}

with open(CONTIGS_FILE, "r") as reader:
    for line in reader:
        contig, size = line.rstrip().split()
        CONTIGS[contig] = int(size)

SAMPLES = pd.read_table(MANIFEST)

def get_sparse_matrices_from_sample(wildcards):
    return ["region_matrices/%s/%s.%d_%d.pkl" % (wildcards.sample, wildcards.sample, part, BAM_PARTITIONS) for part in range(BAM_PARTITIONS + UNMAPPED_PARTITIONS)]

def get_multiple_contigs(sample, chr, num):
    names = SAMPLE_MAPPING_JOBS[sample]["%s.%s" % (chr, num)]
    return ["region_matrices/%s/%s.%s.pkl" % (sample, sample, region) for region in names]

localrules: all, get_headers, make_jobfile

rule all:
    input:  expand("mapping/{sample}/{sample}/wssd_out_file", sample = SAMPLES.sn)

rule wssd_merge:
    input: wssd = expand("mapping/{{sample}}/{{sample}}/wssd_out_file.{contig}", contig = CONTIGS.keys()),
           shelve = expand("region_matrices/{{sample}}/{{sample}}.{part}_%d.dat" % BAM_PARTITIONS, part = range(BAM_PARTITIONS + UNMAPPED_PARTITIONS))
    output: "mapping/{sample}/{sample}/wssd_out_file"
    params: sge_opts="-l mfree=8G -l data_scratch_ssd_disk_free=10G -pe serial 1 -N merge_wssd -l h_rt=5:00:00 -soft -l gpfsstate=0"
    log: "log/wssd_merge/{sample}.txt"
    resources: mem=8
    priority: 40
    benchmark: "benchmarks/wssd_merge/{sample}.txt"
    run:
        tempfile = "/data/scratch/ssd/%s.wssd_out_file" % wildcards.sample
        shell("python3 merger.py {tempfile} --wssd_merge {input.wssd}")
        shell("rsync {tempfile} {output}")
        shell("rm {tempfile}")

rule merge_sparse_matrices:
    input: expand("region_matrices/{{sample}}/{{sample}}.{part}_%d.dat" % BAM_PARTITIONS, part = range(BAM_PARTITIONS + UNMAPPED_PARTITIONS))
    output: temp("mapping/{sample}/{sample}/wssd_out_file.{contig}")
    params: sge_opts = "-l mfree=8G -l data_scratch_ssd_disk_free=10G -pe serial 1 -N merge_sample -l h_rt=5:00:00 -soft -l gpfsstate=0"
    log: "log/merge/{sample}.{contig}.txt"
    resources: mem=8
    priority: 30
    benchmark: "benchmarks/merge/{sample}.{contig}.txt"
    run:
        infile_glob = os.path.commonprefix(input) + "*"
        tempfile = "/data/scratch/ssd/%s.wssd_out_file.%s" % (wildcards.sample, wildcards.contig)
        shell('python3 merger.py {tempfile} --infile_glob "{infile_glob}" --contig {wildcards.contig} --live_merge')
        shell('python3 merger.py {tempfile} --infile_glob "{infile_glob}" --contig {wildcards.contig}')
        shell("rsync {tempfile} {output}")
        shell("rm {tempfile}")

rule merge_sparse_matrices_live:
    input: bam = lambda wildcards: SAMPLES.ix[SAMPLES.sn == wildcards.sample, "bam"], chunker = "bin/bam_chunker_cascade", bam_check = "BAMS_READABLE", index_check = "MRSFASTULTRA_INDEXED"
    output: temp("mapping/{sample}/{sample}/wssd_out_file.{contig}")
    params: sge_opts = "-l mfree=8G -l data_scratch_ssd_disk_free=1G -pe serial 1 -N merge_sample -l h_rt=48:00:00 -soft -l gpfsstate=0"
    log: "log/merge/{sample}.txt"
    resources: mem=8
    benchmark: "benchmarks/merger/{sample}.txt"
    priority: 10
    run:
        infile_glob = os.path.commonprefix(get_sparse_matrices_from_sample(wildcards)) + "*"
        tempfile = "/data/scratch/ssd/%s.wssd_out_file.%s" % (wildcards.sample, wildcards.contig)
        shell('python3 merger.py {tempfile} --infile_glob "{infile_glob}" --live_merge --contig {wildcards.contig}')
        shell("rsync {tempfile} {output}")
        shell("rm {tempfile}")

rule map_and_count:
    input: bam = lambda wildcards: SAMPLES.ix[SAMPLES.sn == wildcards.sample, "bam"], index = lambda wildcards: SAMPLES.ix[SAMPLES.sn == wildcards.sample, "index"], chunker = "bin/bam_chunker_cascade", readable = "BAMS_READABLE", mrsfast_indexed = "MRSFASTULTRA_INDEXED"
    output: [temp("region_matrices/{sample}/{sample}.{part}_%d.%s") % (BAM_PARTITIONS, ext) for ext in ["dat", "bak", "dir"]]
    params: sge_opts = "-l mfree=5G -N map_count -l h_rt=2:00:00 -soft -l gpfsstate=0"
    benchmark: "benchmarks/counter/{sample}/{sample}.{part}.%d.txt" % BAM_PARTITIONS
    priority: 20
    resources: mem=5
    log: "log/map/{sample}/{part}_%s.txt" % BAM_PARTITIONS
    run:
        masked_ref_name = os.path.basename(MASKED_REF)
        ofprefix = output[0].replace(".dat", "")
        fifo = "%s/mrsfast_fifo" % TMPDIR
        if TMPDIR != "":
            local_index = "%s/%s" % (TMPDIR, os.path.basename(input.index[0]))
        else:
            local_index = input.index[0]
        mrsfast_ref_path = "/var/tmp/mrsfast_index/%s" % masked_ref_name
        rsync_opts = """rsync {0}.index /var/tmp/mrsfast_index --bwlimit 10000 --copy-links;
                        rsync {2} {3} --bwlimit 10000 --copy-links; 
                        touch {1}; 
                        echo Finished rsync from {0} to {1} >> /dev/stderr; 
                        echo Finished rsync from {2} to {3} >> /dev/stderr; """.format(MASKED_REF, mrsfast_ref_path, input.index[0], local_index)

        read_counter_args = "--max_basepairs_in_mem %d" % MAX_BP
        shell("hostname; echo part: {wildcards.part} nparts: {BAM_PARTITIONS} unmapped parts: {UNMAPPED_PARTITIONS}; mkfifo {fifo}; {rsync_opts}")
        shell("{input.chunker} -b {input.bam} -i {local_index} -p {wildcards.part} -n {BAM_PARTITIONS} -u {UNMAPPED_PARTITIONS} 2>> /dev/stderr | "
            "mrsfast --search {mrsfast_ref_path} -n 0 -e 2 --crop 36 --seq /dev/stdin -o {fifo} --disable-nohit >> /dev/stderr | "
            "python3 read_counter.py {fifo} {ofprefix} {CONTIGS_FILE} {read_counter_args}"
            )

rule check_bam_files:
    input: [bam for bam in SAMPLES.bam]
    output: touch("BAMS_READABLE")
    params: sge_opts = ""
    priority: 50
    run:
        for bamfile in input:
            try:
                test = pysam.AlignmentFile(bamfile)
            except ValueError as e:
                print("Error: could not open %s as bam.\n%s\n" % (bamfile, str(e)), file = sys.stderr)
                sys.exit(1)

rule get_headers:
    input: expand("bam_headers/{sample}.txt", sample = SAMPLES.sn)

rule get_header:
    input: lambda wildcards: SAMPLES.ix[SAMPLES.sn == wildcards.sample, "bam"]
    output: "bam_headers/{sample}.txt"
    params: sge_opts = ""
    shell:
        "samtools view -H {input} > {output}"

rule check_index:
    input: MASKED_REF
    output: touch("MRSFASTULTRA_INDEXED"), temp(".mrsfast_index_test_output.txt")
    params: sge_opts = ""
    run:
        try:
            shell("mrsfast --search {input[0]} --seq dummy.fq > {output[1]}")
        except CalledProcessError as e:
            sys.exit("Reference %s was not indexed with the current version of mrsfastULTRA. Please reindex." % input[0])

rule make_chunker:
    input: "src/chunker_cascade.cpp", "Makefile"
    output: "bin/bam_chunker_cascade"
    params: sge_opts = ""
    shell:
        "make"
