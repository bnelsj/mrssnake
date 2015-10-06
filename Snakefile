import os
import pysam

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

SAMPLES = {}

with open(MANIFEST, "r") as reader:
    for line in reader:
        sn, bam = line.rstrip().split()
        SAMPLES[sn] = bam

REGIONS = {}

for sample, bam in SAMPLES.items():
    REGIONS[sample] = []

    header = pysam.view("-H", bam)
    contigs = {line.split()[1].replace("SN:",""): int(line.split()[2].replace("LN:","")) for line in header if line.startswith("@SQ")}

    for contig, contig_length in contigs.items():
        nregions = contig_length // WINDOW_SIZE + 1
        coords = ["%d_%d" % (i * WINDOW_SIZE, (i+1) * WINDOW_SIZE) for i in range(nregions)]
        REGIONS[sample].extend(["%s.%s" % (contig, x) for x in coords])

def get_sparse_matrices_from_sample(wildcards):
    return ["region_matrices/%s/%s.%s.pkl" % (wildcards.sample, wildcards.sample, region) for region in REGIONS[wildcards.sample] + ["unmapped"]]

localrules: all

rule all:
    input: expand("mapping/{sample}/{sample}/wssd_out_file", sample = SAMPLES.keys())

rule merge_sparse_matrices:
    input: get_sparse_matrices_from_sample
    output: "mapping/{sample}/{sample}/wssd_out_file"
    params: sge_opts = "-l mfree=8G -l disk_free=20G"
    benchmark: "benchmarks/merger/{sample}.json"
    shell:
        "python3 merger.py $TMPDIR/wssd_out_file --infiles {input}; "
        "rsync $TMPDIR/wssd_out_file {output}"

rule map_and_count_unmapped:
    input: lambda wildcards: SAMPLES[wildcards.sample]
    output: "region_matrices/{sample}/{sample}.unmapped.pkl"
    params: sge_opts = "-l mfree=12G"
    benchmark: "benchmarks/counter/{sample}.unmapped.json"
    run:
        fifo = "$TMPDIR/mrsfast_fifo"
        shell(
            "mkfifo {fifo}; "
            'samtools view {input} "*" | '
            "python3 chunker.py {input} unmapped --fifo {fifo} | "
            "mrsfast --search {MASKED_REF} -n 0 -e 2 --crop 36 --seq1 /dev/stdin -o {fifo} -u /dev/stdout | "
            "python3 mrsfast_parser.py {fifo} /dev/stdout {TEMPLATE} | "
            "python3 read_counter.py /dev/stdin {output} {CONTIGS_FILE}"
            )

rule map_and_count:
    input: lambda wildcards: SAMPLES[wildcards.sample]
    output: "region_matrices/{sample}/{sample}.{chr}.{num}.pkl"
    params: sge_opts = "-l mfree=6G"
    benchmark: "benchmarks/counter/{sample}.{chr}.{num}.json"
    run:
        chr = wildcards.chr
        start, end = wildcards.num.split("_")
        fifo = "$TMPDIR/mrsfast_fifo"
        if TRIM_CONTIGS:
            chr_trimmed = chr.replace("chr", "")
        else:
            chr_trimmed = chr
        shell(
            "mkfifo {fifo}; "
            "python3 chunker.py {input} {chr_trimmed} --start {start} --end {end} --fifo {fifo} | "
            "mrsfast --search {MASKED_REF} -n 0 -e 2 --crop 36 --seq1 /dev/stdin -o {fifo} -u /dev/stdout | "
            "python3 mrsfast_parser.py {fifo} /dev/stdout {TEMPLATE} | "
            "python3 read_counter.py /dev/stdin {output} {CONTIGS_FILE} --common_contigs {chr}"
            )
