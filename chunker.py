from __future__ import print_function
from __future__ import division

import sys
import argparse
import pysam
import numpy as np

def get_fetch_list(chr, start, end):
    if start is None or end is None:
        fetch_list = [chr]
    else:
        fetch_list = [chr, start, end]
    return fetch_list

def chunk_read(i, read, outfile, chunk_size, write_step = 1000000, line_buffer =1000):
    "Break a read into kmers"
    n_to_do = read.rlen // chunk_size
    for k in range(n_to_do):
        outfile.write(">0\n" + read.seq[k*chunk_size : k*chunk_size + chunk_size] + "\n")

    if i % line_buffer == 0:
        outfile.flush()

    if i % write_step == 0:
        sys.stderr.write("Chunker: %d reads chunked\n" % i)
        sys.stderr.flush()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("bamfile")
    parser.add_argument("chr", help="Input bam contig to chunk. Use 'unmapped' to get all unmapped reads.")
    parser.add_argument("--start", type=int)
    parser.add_argument("--end", type=int)
    parser.add_argument("--outfile", default = "/dev/stdout")
    parser.add_argument("--chunk_size", default = 36, type = int)

    args = parser.parse_args()

    outfile = open(args.outfile, "w")
    bamfile = pysam.AlignmentFile(args.bamfile)

    fetch_list = get_fetch_list(args.chr, args.start, args.end)

    msg = "Chunker: Chunking %s reads\n" % " ".join(map(str, fetch_list))
    sys.stderr.write(msg)

    try:
        if args.chr == "unmapped": # Assumes unmapped reads passed as sam by samtools view "*"
            for i, read in enumerate(bamfile.fetch(until_eof=True)):
                chunk_read(i, read, outfile, args.chunk_size)
        else:
            for i, read in enumerate(bamfile.fetch(*fetch_list, until_eof=True)):
                if args.start is not None and read.reference_start < args.start:
                    continue
                else:
                    chunk_read(i, read, outfile, args.chunk_size) 
        bamfile.close()
    except BrokenPipeError as e:
        sys.stdout.write(str(e))
        sys.stdout.flush()
        sys.exit(1)
    except ValueError as e:
        sys.stdout.write(str(e))
        sys.stdout.flush()
        sys.exit(1)
    finally:
        bamfile.close()

    #Handle regions where there are no reads
    try:
        sys.stderr.write("Chunker: finished chunking %d reads\n" % int(i + 1))
    except NameError:
        outfile.write("\n")
        sys.stderr.write("Chunker: found no reads to chunk\n")
    finally:
        outfile.close()

    sys.exit()
