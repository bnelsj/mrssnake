import os

shell.prefix("source config.sh; ")


if config == {}:
    configfile: "config.json"

MANIFEST = config["manifest"]

MASKED_REF = config["masked_ref"]
CONTIGS_FILE = config["contigs"]
WINDOW_SIZE = config["window_size"]
TEMPLATE = config["template"]
TRIM_CONTIGS = config["trim_contigs"].upper().startswith("Y")

if not os.path.exists("log"):
    os.makedirs("log")

CONTIGS = {}

with open(CONTIGS_FILE, "r") as reader:
    for line in reader:
        contig, size = line.rstrip().split()
        CONTIGS[contig] = int(size)

REGIONS = []

for contig in CONTIGS.keys():
    nregions = CONTIGS[contig] // WINDOW_SIZE + 1
    REGIONS.extend(["%s.%d" % (contig, x) for x in range(nregions)])

SAMPLES = {}

with open(MANIFEST, "r") as reader:
    for line in reader:
        sn, bam = line.rstrip().split()
        SAMPLES[sn] = bam

def get_sparse_matrices_from_sample(wildcards):
    return ["region_matrices/%s.%s.pkl" % (wildcards.sample, region) for region in REGIONS]

def get_region_nums_from_contigs(wildcards):
    contig_length = CONTIGS[wildcards.chr]
    nregions = contig_length // WINDOW_SIZE + 1
    return [str(x) for x in range(nregions)]

def get_region_from_contig_and_num(chr, num):
    contig_length = CONTIGS[chr]
    start = int(num) * WINDOW_SIZE
    end = min(start + WINDOW_SIZE, contig_length - 1)
    return (chr, str(start), str(end))

localrules: all, run_tests

rule all:
    input: expand("mapping/{sample}/{sample}/wssd_out_file", sample = SAMPLES.keys())

rule merge_sparse_matrices:
    input: get_sparse_matrices_from_sample
    output: "mapping/{sample}/{sample}/wssd_out_file"
    params: sge_opts = "-l mfree=8G -l disk_free=20G"
    benchmark: "benchmarks/merger/{sample}.json"
    shell:
        "python3 merger.py $TMPDIR/wssd_out_file --infiles {input}; "
        "rsync TMPDIR/wssd_out_file {output}"

rule map_and_count:
    input: lambda wildcards: SAMPLES[wildcards.sample]
    output: "region_matrices/{sample}.{chr}.{num}.pkl"
    params: sge_opts = "-l mfree=6G"
    benchmark: "benchmarks/counter/{sample}.{chr}.{num}.json"
    run:
        chr, start, end = get_region_from_contig_and_num(wildcards.chr, wildcards.num)
        fifo = "$TMPDIR/mrsfast_fifo"
        if TRIM_CONTIGS:
            chr_trimmed = chr.replace("chr", "")
        else:
            chr_trimmed = chr
        shell(
            "mkfifo {fifo}; "
            "python3 chunker.py {input} {chr_trimmed} {start} {end} | "
            "mrsfast --search {MASKED_REF} -n 0 -e 2 --crop 36 --seq1 /dev/stdin -o {fifo} -u /dev/stdout | "
            "python3 mrsfast_parser.py {fifo} /dev/stdout {TEMPLATE} | "
            "python3 mrsfast_simple_mapper.py /dev/stdin {output} {CONTIGS_FILE} --common_contigs {chr}"
            )
