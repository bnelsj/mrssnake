import os
import sys
import pysam
import hashlib
import json

shell.prefix("source config.sh; ")


if config == {}:
    configfile: "config.yaml"

MANIFEST = config["manifest"]

MASKED_REF = config["masked_ref"]
CONTIGS_FILE = config["contigs"]
WINDOW_SIZE = config["window_size"]

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

JOBFILE = config["jobfile"]
# Load JOBFILE if it exists
if os.path.exists(JOBFILE):
    with open(JOBFILE, "r") as reader:
        SAMPLE_MAPPING_JOBS = json.load(reader)
else:
    print("""WARNING: no jobfile found. Run 'snakemake make_jobfile' to create one.""")

def create_jobfile(wildcards):
    # Combine short contigs for mapping
    SAMPLE_MAPPING_JOBS = {}
    for sn, bamfile in SAMPLES.items():
        bam = pysam.AlignmentFile(bamfile)
        contigs = {bam.references[i]: bam.lengths[i] for i in range(bam.nreferences)}
        bam.close()

        full_jobs = []
        jobs = []

        for contig, size in contigs.items():
            if size >= WINDOW_SIZE:
                nregions = size // WINDOW_SIZE + 1
                coords = ["%d_%d" % (i * WINDOW_SIZE, (i+1) * WINDOW_SIZE) for i in range(nregions)]
                full_jobs.extend(["%s.%s" % (contig, x) for x in coords])
               
            if size < WINDOW_SIZE:
                job_added = False
                for job in jobs:
                    if job[1] + size < WINDOW_SIZE:
                        job[0].append(contig)
                        job[1] += size
                        job_added = True
                        break
                if not job_added:
                    jobs.append([[contig], size])
        (contigs, size) = jobs[0]
        SAMPLE_MAPPING_JOBS[sn] = {"multiple.%s" % (hashlib.md5(" ".join(contigs).encode("utf-8")).hexdigest()): contigs for (contigs, size) in jobs}
        full_jobs_dict = {name: name.split(".") for name in full_jobs}
        SAMPLE_MAPPING_JOBS[sn].update(full_jobs_dict)
        SAMPLE_MAPPING_JOBS[sn].update({"unmapped": "unmapped"})
        return SAMPLE_MAPPING_JOBS

def get_sparse_matrices_from_sample(wildcards):
    names = SAMPLE_MAPPING_JOBS[wildcards.sample].keys()
    return ["region_matrices/%s/%s.%s.pkl" % (wildcards.sample, wildcards.sample, region) for region in names]

def get_multiple_contigs(sample, chr, num):
    names = SAMPLE_MAPPING_JOBS[sample]["%s.%s" % (chr, num)]
    return ["region_matrices/%s/%s.%s.pkl" % (sample, sample, region) for region in names]

localrules: all, get_headers, make_jobfile

rule all:
    input: expand("mapping/{sample}/{sample}/wssd_out_file", sample = SAMPLES.keys())#,
           #expand("region_matrices/list/{sample}.txt", sample = SAMPLES.keys())

rule list_sparse_matrices:
    input: get_sparse_matrices_from_sample
    output: "region_matrices/list/{sample}.txt"
    params: sge_opts = ""
    run:
        with open(output[0], "w") as of:
            for infile in input:
                of.write(infile + "\n")

rule merge_sparse_matrices:
    input: get_sparse_matrices_from_sample
    output: "mapping/{sample}/{sample}/wssd_out_file"
    params: sge_opts = "-l mfree=32G -l data_scratch_ssd_disk_free=10G -pe serial 1 -N merge_sample"
    benchmark: "benchmarks/merger/{sample}.json"
    run:
        shell("mkdir -p /data/scratch/ssd/{wildcards.sample}")
        shell("python3 merger.py /data/scratch/ssd/{wildcards.sample}/wssd_out_file --infiles {input}")
        shell("rsync /data/scratch/ssd/{wildcards.sample}/wssd_out_file {output}")
        shell("rm /data/scratch/ssd/{wildcards.sample}/wssd_out_file")

rule map_and_count_unmapped:
    input: lambda wildcards: SAMPLES[wildcards.sample], "BAMS_READABLE"
    output: "region_matrices/{sample}/{sample}.unmapped.pkl"
    params: sge_opts = "-pe orte 4 -l mfree=40G -N map_unmapped"
    benchmark: "benchmarks/counter/{sample}/{sample}.unmapped.json"
    priority: 50
    run:
        fifo = "$TMPDIR/mrsfast_fifo"
        masked_ref_name = os.path.basename(MASKED_REF)

        shell(
            "mkfifo {fifo}; "
            "mkdir -p /var/tmp/mrsfast_index; "
            "rsync {MASKED_REF}.index /var/tmp/mrsfast_index; "
            "touch /var/tmp/mrsfast_index/{MASKED_REF}; "
            "echo Finished rsync from {MASKED_REF} to /var/tmp/mrsfast_index > /dev/stderr; "
            "samtools view -h {input[0]} '*' | "
            "python3 chunker.py /dev/stdin --contig unmapped | "
            "mrsfast --search /var/tmp/mrsfast_index/{masked_ref_name} -n 0 -e 2 --crop 36 --seq /dev/stdin -o {fifo} --disable-nohit --threads 4 >> /dev/stderr | "
            "python3 read_counter.py {fifo} {output} {CONTIGS_FILE} --all_contigs"
            )

rule merge_matrices_live:
    input: lambda wildcards: SAMPLES[wildcards.sample]
    output: "mapping/{sample}/{sample}/wssd_out_file_live"
    params: sge_opts = "-l mfree=32G -N merge_{sample}", files = get_sparse_matrices_from_sample
    priority: 10
    shell:
        "python3 live_merger.py {output} --infiles {params.files}"

rule map_and_count:
    input: lambda wildcards: SAMPLES[wildcards.sample], "BAMS_READABLE"
    output: "region_matrices/{sample}/{sample}.{chr}.{num}.pkl"
    params: sge_opts = "-l mfree=4G -N map_count"
    benchmark: "benchmarks/counter/{sample}/{sample}.{chr}.{num}.json"
    priority: 20
    run:
        fifo = "$TMPDIR/mrsfast_fifo"
        masked_ref_name = os.path.basename(MASKED_REF)

        if (not USE_SOURCE_CONTIGS) and ARRAY_CONTIGS != []:
            common_contigs = ARRAY_CONTIGS
        else:
            common_contigs = wildcards.chr

        if wildcards.chr == "multiple":
            contigs = get_multiple_contigs(wildcards.sample, wildcards.chr, wildcards.num)
            chunker_args = "--contigs %s" % " ".join(contigs)
        else:
            start, end = wildcards.num.split("_")
            chunker_args = "--contig %s --start %s --end %s" % (wildcards.chr, start, end)

        shell(
            "mkfifo {fifo}; "
            "mkdir -p /var/tmp/mrsfast_index; "
            "rsync {MASKED_REF}.index /var/tmp/mrsfast_index; "
            "touch /var/tmp/mrsfast_index/{MASKED_REF}; "
            "echo Finished rsync from {MASKED_REF} to /var/tmp/mrsfast_index > /dev/stderr; "
            "python3 chunker.py {input[0]} {chunker_args} | "
            "mrsfast --search /var/tmp/mrsfast_index/{masked_ref_name} -n 0 -e 2 --crop 36 --seq /dev/stdin -o {fifo} --disable-nohit >> /dev/stderr | "
            "python3 read_counter.py {fifo} {output} {CONTIGS_FILE} --common_contigs {common_contigs}"
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

rule make_jobfile:
    input: MANIFEST
    output: JOBFILE
    params: sge_opts = ""
    run:
        sample_mapping_jobs = create_jobfile(wildcards)
        with open(output[0], "w") as writer:
            json.dump(sample_mapping_jobs, writer)

