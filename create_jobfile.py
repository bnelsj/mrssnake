import argparse
import hashlib
import pysam
import json

def create_jobdict(samples, window_size):
    # Combine short contigs for mapping
    SAMPLE_MAPPING_JOBS = {}
    unique_contig_dicts = {}
    for sn, bamfile in samples.items():
        try:
            bam = pysam.AlignmentFile(bamfile)
        except ValueError as e:
            print("Error: Can't open %s" % bamfile, file=sys.stderr)
        contigs = {bam.references[i]: bam.lengths[i] for i in range(bam.nreferences)}
        bam.close()

        full_jobs = []
        full_jobs_dict = {}
        jobs = []

        for contig, size in contigs.items():
            if size >= window_size:
                nregions = size // window_size + 1
                coords = ["%d_%d" % (i * window_size, (i+1) * window_size) for i in range(nregions)]
                full_jobs.extend(["%s.%s" % (contig, x) for x in coords])
               
            if size < window_size:
                job_added = False
                for job in jobs:
                    if job[1] + size < window_size:
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument()

    args = parser.parse_args()

elif snakemake is not None:
    samples = {}
    with open(snakemake.input[0], "r") as reader:
        for line in reader:
            sn, bam = line.rstrip().split()
            samples[sn] = bam

    jobdict = create_jobdict(samples, snakemake.params.window_size)
    with open(snakemake.output[0], "w") as writer:
        json.dump(jobdict, writer)

