import os
import sys
import pysam
import hashlib
import json
from MappingJob import *
from subprocess import CalledProcessError


if config == {}:
    configfile: "config.yaml"

MANIFEST = config["manifest"]
REFERENCE = config["reference"]
MASKED_REF = config[REFERENCE]["masked_ref"]
CONTIGS_FILE = config[REFERENCE]["contigs"]

BAM_PARTITIONS = config["bam_partitions"]
UNMAPPED_PARTITIONS = config["unmapped_partitions"]
AUTO_ASSIGN = config["auto_assign"]
MAX_BP = config["max_bp_in_mem"]

ARRAY_CONTIGS = config["array_contigs"]

AMAZON = config["amazon"]
BUCKET = config["bucket"]
TMPDIR = config["tmpdir"]
LIVE_MERGE = config["live_merge"]
CLEAN_TEMP_FILES = config["clean_temp_files"]

if LIVE_MERGE:
    ruleorder: merge_sparse_matrices_live > merge_sparse_matrices
else:
    ruleorder: merge_sparse_matrices > merge_sparse_matrices_live

if not AMAZON:
    shell.prefix("source config.sh; ")

if not os.path.exists("log"):
    os.makedirs("log")

CONTIGS = {}

with open(CONTIGS_FILE, "r") as reader:
    for line in reader:
        contig, size = line.rstrip().split()
        CONTIGS[contig] = int(size)

SAMPLES = {}

with open(MANIFEST, "r") as reader:
    for line in reader:
        sn, bam = line.rstrip().split()
        SAMPLES[sn] = bam

def get_sparse_matrices_from_sample(wildcards):
    return ["region_matrices/%s/%s.%d_%d.pkl" % (wildcards.sample, wildcards.sample, part, BAM_PARTITIONS) for part in range(BAM_PARTITIONS + UNMAPPED_PARTITIONS)]

def get_multiple_contigs(sample, chr, num):
    names = SAMPLE_MAPPING_JOBS[sample]["%s.%s" % (chr, num)]
    return ["region_matrices/%s/%s.%s.pkl" % (sample, sample, region) for region in names]

localrules: all, get_headers, make_jobfile

rule all:
    input:  expand("mapping/{sample}/{sample}/wssd_out_file", sample = SAMPLES.keys()),
            expand("finished/{sample}", sample = SAMPLES.keys())

rule clean:
    input:  to_remove = expand("region_matrices/{{sample}}/{{sample}}.{part}_%d.{ext}" % BAM_PARTITIONS, part = range(BAM_PARTITIONS + UNMAPPED_PARTITIONS), ext = ["dat", "bak", "dir"]),
            to_keep = "mapping/{sample}/{sample}/wssd_out_file"
    output: touch("finished/{sample}")
    params: sge_opts = "-l h_rt=01:00:00"
    run:
        if CLEAN_TEMP_FILES:
            remove_glob = os.path.commonprefix(input.to_remove) + "*"
            shell("rm {remove_glob}")

            if AMAZON:
                dirname = os.path.dirname(input.to_remove[0])
                shell('aws s3 rm s3://{BUCKET}/{dirname}/ --recursive --exclude "*wssd_out_file*"')

rule merge_sparse_matrices:
    input: expand("region_matrices/{{sample}}/{{sample}}.{part}_%d.dat" % BAM_PARTITIONS, part = range(BAM_PARTITIONS + UNMAPPED_PARTITIONS))
    output: "mapping/{sample}/{sample}/wssd_out_file"
    params: sge_opts = "-l mfree=40G -l data_scratch_ssd_disk_free=10G -pe serial 1 -N merge_sample -l h_rt=24:00:00"
    log: "log/merge/{sample}.txt"
    resources: mem=40
    benchmark: "benchmarks/merger/{sample}.json"
    run:
        infile_glob = os.path.commonprefix(input) + "*"
        if AMAZON:
            shell('python3 merger.py {output} --infile_glob "{infile_glob}"')
        else:
            shell("mkdir -p /data/scratch/ssd/{wildcards.sample}")
            shell('python3 merger.py /data/scratch/ssd/{wildcards.sample}/wssd_out_file --infile_glob "{infile_glob}"')
            shell("rsync /data/scratch/ssd/{wildcards.sample}/wssd_out_file {output}")
            shell("rm /data/scratch/ssd/{wildcards.sample}/wssd_out_file")
        if AMAZON:
            "aws s3 cp {output} s3://{BUCKET}/{output}"

rule merge_sparse_matrices_live:
    input: bam = lambda wildcards: SAMPLES[wildcards.sample], chunker = "bin/bam_chunker_cascade", bam_check = "BAMS_READABLE", index_check = "MRSFASTULTRA_INDEXED"
    output: "mapping/{sample}/{sample}/wssd_out_file"
    params: sge_opts = "-l mfree=40G -l data_scratch_ssd_disk_free=10G -pe serial 1 -N merge_sample -l h_rt=48:00:00"
    log: "log/merge/{sample}.txt"
    resources: mem=40
    benchmark: "benchmarks/merger/{sample}.json"
    priority: 20
    run:
        infile_glob = os.path.commonprefix(get_sparse_matrices_from_sample(wildcards)) + "*"
        if AMAZON:
            shell('python3 merger.py {output} --infile_glob "{infile_glob}" --live_merge')
        else:
            shell("mkdir -p /data/scratch/ssd/{wildcards.sample}")
            shell('python3 merger.py /data/scratch/ssd/{wildcards.sample}/wssd_out_file --infile_glob "{infile_glob}" --live_merge')
            shell("rsync /data/scratch/ssd/{wildcards.sample}/wssd_out_file {output}")
            shell("rm /data/scratch/ssd/{wildcards.sample}/wssd_out_file")
        if AMAZON:
            "aws s3 cp {output} s3://{BUCKET}/{output}"

rule map_and_count:
    input: lambda wildcards: SAMPLES[wildcards.sample], "bin/bam_chunker_cascade", "BAMS_READABLE", "MRSFASTULTRA_INDEXED"
    output: ["region_matrices/{sample}/{sample}.{part}_%d.%s" % (BAM_PARTITIONS, ext) for ext in ["dat", "bak", "dir"]]
    params: sge_opts = "-l mfree=4G -N map_count -l h_rt=10:00:00"
    benchmark: "benchmarks/counter/{sample}/{sample}.{part}.%d.json" % BAM_PARTITIONS
    priority: 20
    resources: mem=4
    log: "log/map/{sample}/{part}_%s.txt" % BAM_PARTITIONS
    shadow: AMAZON
    run:
        masked_ref_name = os.path.basename(MASKED_REF)
        ofprefix = os.path.commonprefix(output)
        if AMAZON:
            fifo = "mrsfast_fifo"
            rsync_opts = ""
            mrsfast_ref_path = MASKED_REF
        else:
            fifo = "%s/mrsfast_fifo" % TMPDIR
            rsync_opts = "rsync {0}.index /var/tmp/mrsfast_index; touch /var/tmp/mrsfast_index/{0}; echo Finished rsync from {0} to /var/tmp/mrsfast_index > /dev/stderr; ".format(MASKED_REF)
            mrsfast_ref_path = "/var/tmp/%s" % masked_ref_name

        if ARRAY_CONTIGS != [] and ARRAY_CONTIGS is not None:
            common_contigs = "--common_contigs " + " ".join(ARRAY_CONTIGS)
        else:
            common_contigs = ""

        if AUTO_ASSIGN:
            read_counter_args = "--max_basepairs_in_mem %d" % MAX_BP
        else:
            read_counter_args = ""

        shell(
            "hostname; "
            "mkfifo {fifo}; "
            "{rsync_opts}"
            "./bin/bam_chunker_cascade -b {input[0]} -p {wildcards.part} -n {BAM_PARTITIONS} -u {UNMAPPED_PARTITIONS} | "
            "mrsfast --search /var/tmp/mrsfast_index/{masked_ref_name} -n 0 -e 2 --crop 36 --seq /dev/stdin -o {fifo} --disable-nohit >> /dev/stderr | "
            "python3 read_counter.py {fifo} {ofprefix} {CONTIGS_FILE} {common_contigs} {read_counter_args}"
            )
        if AMAZON:
            shell("aws s3 cp {output} s3://{BUCKET}/{output}")

rule check_bam_files:
    input: [bam for key, bam in SAMPLES.items()]
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
    input: expand("bam_headers/{sample}.txt", sample = SAMPLES.keys())

rule get_header:
    input: lambda wildcards: SAMPLES[wildcards.sample]
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
