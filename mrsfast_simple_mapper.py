from __future__ import print_function
from __future__ import division

try:
    import cPickle as pickle
except ImportError:
    import pickle

import sys
import argparse

import pysam
import numpy as np

from scipy.sparse import lil_matrix


def get_edist(read):
    return read.get_tag("NM")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", default="/dev/stdin")
    parser.add_argument("outfile", default="wssd_out_file")
    parser.add_argument("contig_lengths", help="tab-delimited file with contig names and lengths")
    parser.add_argument("--max_edist", default = 2, type = int, help = "Maximum edit distance of input reads")
    parser.add_argument("--common_contigs", default = [], nargs="+", help = "Create numpy array for common contigs (Much faster, more memory)")

    args = parser.parse_args()

    contigs = {}

    with open(args.contig_lengths, "r") as reader:
        for line in reader:
            contig, length = line.rstrip().split()
            #contig = contig.replace("chr", "")
            contigs[contig] = int(length)

    samfile = pysam.AlignmentFile(args.infile, "r", check_sq = False)

    read_dict = {}
    nrows = 2 * args.max_edist + 2
    nstart_rows = nrows // 2

    for i, read in enumerate(samfile):
        contig = samfile.getrname(read.reference_id)
        start = read.qstart
        end = read.qend

        edist = get_edist(read)
        if edist > args.max_edist:
            continue

        if contig not in read_dict:
            length = contigs[contig]
            if contig in args.common_contigs:
                read_dict[contig] = np.zeros((nrows, length), dtype=np.uint16)
                sys.stdout.write(contig + "(%s, %s) numpy array\n" % (nrows, length))
            else:
                read_dict[contig] = lil_matrix((nrows, length), dtype=np.uint16)
                sys.stdout.write(contig + " scipy lil_matrix\n")


        if contig in args.common_contigs:
            # Update read depth counts for numpy array
            read_dict[contig][edist + nstart_rows, start:end+1] += 1
        else:
            # Update read depth counts for sparse matrix
            slice = read_dict[contig][edist + nstart_rows, start:end+1].toarray()
            slice += 1
            read_dict[contig][:, start:end+1] = slice
 
        # Update read start counts
        read_dict[contig][edist, start] += 1

        if i % 100000 == 0:
            sys.stdout.write("%d reads processed\n" % i)
            sys.stdout.flush()

    samfile.close()

    for contig, array in read_dict.items():
        if contig in args.common_contigs:
            read_dict[contig] = lil_matrix(array)
            del(array)

    with open(args.outfile, "wb") as outfile:
        pickle.dump(read_dict, outfile, pickle.HIGHEST_PROTOCOL)
