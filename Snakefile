import os

shell.prefix("source config.sh; ")


if config == {}:
    configfile: "config.json"

MANIFEST = config["manifest"]

MASKED_REF = config["masked_ref"]
CONTIGS_FILE = config["contigs"]
WINDOW_SIZE = config["window_size"]

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

def get_sparse_matrices_from_contig(wildcards):
    nums = get_region_nums_from_contig(wildcards)
    return ["region_matrices/%s.%s.%s.pkl" % (wildcards.sample, wildcards.chr, num) for num in nums]

def get_region_nums_from_contig(wildcards):
    contig_length = CONTIGS[wildcards.chr]
    nregions = contig_length // WINDOW_SIZE + 1
    return [str(x) for x in range(nregions)]

def get_region_from_contig_and_num(chr, num):
    contig_length = CONTIGS[chr]
    start = int(num) * WINDOW_SIZE
    end = min(start + WINDOW_SIZE, contig_length - 1)
    return "%s:%d-%d" % (chr, start, end)

localrules: all, run_tests

rule all:
    input: expand("mapping/{sample}/{sample}/wssd_out_file", sample = SAMPLES.keys())

rule run_tests:
    input: expand("fixed_output.pkl")

rule run_old:
    input: "fixed_output.sam"
    output: "fixed_output.pkl"
    params: sge_opts = "-l mfree=20G"
    benchmark: "benchmarks/mrsfast_counter.json"
    shell:
        "python mrsfast_simple_mapper.py {input} {output} {CONTIGS_FILE}"

rule get_wssd_out_file:
    input: expand("chr_matrices/{{sample}}.{chr}.pkl", chr = CONTIGS.keys())
    output: "mapping/{sample}/{sample}/wssd_out_file"
    params: sge_opts = "-l mfree=40G"
    shell:
        ""

rule merge_sparse_matrices:
    input: get_sparse_matrices_from_contig
    output: "chr_matrices/{sample}.{chr}.pkl"
    params: sge_opts = "-l mfree=8G"
    shell:
        ""

rule map_and_count:
    input: lambda wildcards: SAMPLES[wildcards.sample]
    output: "region_matrices/{sample}.{chr}.{num}.pkl"
    params: sge_opts = "-l mfree=40G"
    run:
        region_string = get_region_from_contig_and_num(wildcards.chr, wildcards.num)
        shell(
            "samtools view -h {region_string} {input} | "
            "samtools bam2fq -n - | "
            "mrsfast --search {MASKED_REF} -n 0 -e 2 --crop 36 --seq1 /dev/stdin -o /dev/stdout | "
            "python mrsfast_simple_mapper.py /dev/stdin {output} {CONTIGS_FILE}"
            )
