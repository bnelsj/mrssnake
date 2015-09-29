import argparse
import pysam
import numpy as np
from scipy.sparse import csc_matrix
import cPickle as pickle

def get_edist(read):
    tags = read.get_tags()
    for tag in tags:
        if tag[0] == "NM":
            return tag[1]

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", default="/dev/stdin")
    parser.add_argument("outfile", default="wssd_out_file")
    parser.add_argument("contig_lengths", help="tab-delimited file with contig names and lengths")
    parser.add_argument("--max_edist", default = 2, type = int, help = "Maximum edit distance of input reads")

    args = parser.parse_args()

    samfile = pysam.Samfile(args.infile)

    contigs = {}

    with open(args.contig_lengths, "r") as reader:
        for line in reader:
            contig, length = line.rstrip().split()
            contigs[contig] = int(length)

    read_dict = {}
    nrows = 2 * args.max_edist + 2
    nstart_rows = nrows / 2

    while True:
        try:
            read = samfile.next()
        except StopIteration:
            break

        contig = read.query_name
        start = read.qstart
        end = read.qend
        if contig not in read_dict:
            length = contigs[contig]
            read_dict[contig] = csc_matrix((nrows, length), dtype=np.int8)

        edist = get_edist(read)

        if edist > args.max_edist:
            print >> sys.stderr, "Max edit distance of %d exceeded: %s: %d" % (args.max_edist, read.query_name, edist)
            continue

        # Update read start counts
        read_dict[contig][edist, start] += 1

        # Update read depth counts
        read_dict[contig][edist + nstart_rows, start:end] += 1

    samfile.close()

    with open(args.outfile, "wb") as outfile:
        pickle.dump(read_dict, outfile, pickle.HIGHEST_PROTOCOL)
