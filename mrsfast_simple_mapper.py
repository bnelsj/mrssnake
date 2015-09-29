import sys
import argparse

import pysam
import numpy as np

from scipy.sparse import lil_matrix
import cPickle as pickle

def get_edist(read):
    return read.get_tag("NM")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", default="/dev/stdin")
    parser.add_argument("outfile", default="wssd_out_file")
    parser.add_argument("contig_lengths", help="tab-delimited file with contig names and lengths")
    parser.add_argument("--max_edist", default = 2, type = int, help = "Maximum edit distance of input reads")

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
    nstart_rows = nrows / 2

    for i, read in enumerate(samfile):
        contig = samfile.getrname(read.reference_id)
        start = read.qstart
        end = read.qend
        if contig not in read_dict:
            length = contigs[contig]
            read_dict[contig] = np.zeros((nrows, length), dtype=np.int)

        edist = get_edist(read)

        # Update read depth counts
        read_dict[contig][edist + nstart_rows, start:end+1] += 1
 
        # Update read start counts
        read_dict[contig][edist, start] += 1

        if i % 100000 == 0:
            print "%d reads processed" % i

    samfile.close()

    out_dict = {}

    for contig, array in read_dict.iteritems():
        out_dict[contig] = lil_matrix(array)
        del(array)

    with open(args.outfile, "wb") as outfile:
        pickle.dump(out_dict, outfile, pickle.HIGHEST_PROTOCOL)
