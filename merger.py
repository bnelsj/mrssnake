from __future__ import print_function
from __future__ import division

try:
    import cPickle as pickle
except ImportError:
    import pickle

import argparse
from scipy.sparse import lil_matrix

def sum_counts(new_counts, old_counts):
    for contig, matrix in new_counts.items():
        if contig not in old_counts:
            old_counts[contig] = matrix
        else:
            old_counts[contig] += matrix
    return old_counts

def write_to_h5(counts, outfile):
    pass

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("outfile")
    parser.add_argument("--infiles", required=True, nargs="+")

    args = parser.parse_args()

    contigs = {}

    for infile in args.infiles:
        with open(infile, "r") as file:
            dat = pickle.load(file)
            old_counts = sum_counts(dat, contigs)
