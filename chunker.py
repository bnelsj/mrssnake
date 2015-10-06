from __future__ import print_function
from __future__ import division

import sys
import argparse
import pysam
import numpy as np

def get_fetch_list(chr, start, end):
    if start is None or end is None:
        fetch_list = chr
    else:
        fetch_list = [chr, start, end]
    return fetch_list

def chunk_read(read, outfile, chunk_size):
    n_to_do = read.rlen // chunk_size
    for k in range(n_to_do):
        outfile.write(">0\n" + read.seq[k*chunk_size : k*chunk_size + chunk_size] + "\n")

def chunk_reads(args, fetch_string):
    """This takes blocks of reads from a bam, then chunks them and sends them to mrsfast"""

    outfile = open(args.outfile, "w")
    if args.chr == "unmapped": # Assumes unmapped reads passed as sam by samtools view "*"
        with open(args.bamfile, "r") as input:
            for read in input:
                chunk_read(read, outfile, args.chunk_size)
    else:
        bamfile = pysam.AlignmentFile(args.bamfile)
        fetch_list = get_fetch_list(args.chr, args.start, args.end)
        for read in bamfile.fetch(*fetch_list, until_eof=True):
            chunk_read(read, outfile, args.chunk_size)

    #Handle regions where there are no reads
    try:
        l
    except NameError:
        if args.fifo is not None:
            with open(fifo, "w") as fifo_handle:
                fifo_handle.write("ERROR: no reads for %s\n" % fetch_string)
        outfile.write("\n")
        sys.stderr.write("No reads for %s\n" % fetch_string)
    finally:
        bamfile.close()
        outfile.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("bamfile")
    parser.add_argument("chr", help="Input bam contig to chunk. Use 'unmapped' to get all unmapped reads.")
    parser.add_argument("--start", type=int)
    parser.add_argument("--end", type=int)
    parser.add_argument("--outfile", default = "/dev/stdout")
    parser.add_argument("--chunk_size", default = 36, type = int)
    parser.add_argument("--fifo", help = "Path to fifo file")

    args = parser.parse_args()

    chunk_reads(args, fetch_string)
