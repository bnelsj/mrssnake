import os
import sys
import pysam
import hashlib
import json
from MappingJob import *
from subprocess import CalledProcessError

shell.prefix("source config.sh; ")

if config == {}:
    configfile: "config.yaml"

MANIFEST = config["manifest"]
REFERENCE = config["reference"]
MASKED_REF = config[REFERENCE]["masked_ref"]
CONTIGS_FILE = config[REFERENCE]["contigs"]

BAM_PARTITIONS = config["bam_partitions"]
AUTO_ASSIGN = config["auto_assign"]
MAX_BP = config["max_bp_in_mem"]

USE_SOURCE_CONTIGS = config["use_source_contigs_for_array"]
ARRAY_CONTIGS = config["array_contigs"]

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
    names = SAMPLE_MAPPING_JOBS[wildcards.sample].keys()
    return ["region_matrices/%s/%s.%s.pkl" % (wildcards.sample, wildcards.sample, region) for region in names]

def get_multiple_contigs(sample, chr, num):
    names = SAMPLE_MAPPING_JOBS[sample]["%s.%s" % (chr, num)]
    return ["region_matrices/%s/%s.%s.pkl" % (sample, sample, region) for region in names]

localrules: all, get_headers, make_jobfile

rule all:
    input: expand("mapping/{sample}/{sample}/wssd_out_file", sample = SAMPLES.keys())

rule merge_sparse_matrices:
    input: expand("region_matrices/{sample}/{sample}.{part}_%d.pkl" % BAM_PARTITIONS, sample = SAMPLES.keys(), part = range(BAM_PARTITIONS))
    output: "mapping/{sample}/{sample}/wssd_out_file"
    params: sge_opts = "-l mfree=32G -l data_scratch_ssd_disk_free=10G -pe serial 1 -N merge_sample"
    benchmark: "benchmarks/merger/{sample}.json"
    run:
        shell("mkdir -p /data/scratch/ssd/{wildcards.sample}")
        shell("python3 merger.py /data/scratch/ssd/{wildcards.sample}/wssd_out_file --infiles {input}")
        shell("rsync /data/scratch/ssd/{wildcards.sample}/wssd_out_file {output}")
        shell("rm /data/scratch/ssd/{wildcards.sample}/wssd_out_file")

#rule merge_matrices_live:
#    input: lambda wildcards: SAMPLES[wildcards.sample]
#    output: "mapping/{sample}/{sample}/wssd_out_file_live"
#    params: sge_opts = "-l mfree=32G -N merge_{sample}", files = get_sparse_matrices_from_sample
#    priority: 10
#    shell:
#        "python3 live_merger.py {output} --infiles {params.files}"

rule map_and_count:
    input: lambda wildcards: SAMPLES[wildcards.sample]
    output: "region_matrices/{sample}/{sample}.{part}_%d.pkl" % BAM_PARTITIONS
    params: sge_opts = "-l mfree=4G -N map_count"
    benchmark: "benchmarks/counter/{sample}/{sample}.{part}.%d.json" % BAM_PARTITIONS
    priority: 20
    run:
        fifo = "$TMPDIR/mrsfast_fifo"
        masked_ref_name = os.path.basename(MASKED_REF)

        if ARRAY_CONTIGS != []:
            common_contigs = "--common_contigs " + " ".join(ARRAY_CONTIGS)
        else:
            common_contigs = ""

        if AUTO_ASSIGN:
            read_counter_args = "--max_basepairs_in_mem %d" % MAX_BP
        else:
            read_counter_args = ""

        shell(
            "mkfifo {fifo}; "
            "mkdir -p /var/tmp/mrsfast_index; "
            "rsync {MASKED_REF}.index /var/tmp/mrsfast_index; "
            "touch /var/tmp/mrsfast_index/{MASKED_REF}; "
            "echo Finished rsync from {MASKED_REF} to /var/tmp/mrsfast_index > /dev/stderr; "
            "./bin/bam_chunker {input[0]} {wildcards.part} {BAM_PARTITIONS} | "
            "mrsfast --search /var/tmp/mrsfast_index/{masked_ref_name} -n 0 -e 2 --crop 36 --seq /dev/stdin -o {fifo} --disable-nohit >> /dev/stderr | "
            "python3 read_counter.py {fifo} {output} {CONTIGS_FILE} {common_contigs} {read_counter_args}"
            )

rule check_bam_files:
    input: [bam for key, bam in SAMPLES.items()]
    output: touch("BAMS_READABLE")
    params: sge_opts = ""
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
    input: "src/chunker.cpp", "Makefile"
    output: "bin/bam_chunker"
    params: sge_opts = ""
    shell:
        "make"
