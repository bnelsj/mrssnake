from __future__ import print_function
from __future__ import division

try:
    import cPickle as pickle
except ImportError:
    import pickle

import argparse
import tables
import numpy as np
from scipy.sparse import lil_matrix

def write_to_h5(counts, outfile):
    fout = tables.open_file(outfile, mode="w")
    group = fout.create_group(fout.root, "depthAndStarts_wssd")

    for contig, matrix in counts:
        nrows, ncols = matrix.shape
        nedists = nrows // 2
        carray_empty = tables.CArray(group, contig, tables.UInt32Atom(), (nedists, ncols), filters=tables.Filters(complevel=1, complib="lzo"))

        wssd_contig = [np.reshape(matrix[:, i], (2, nedists)) for i in range(ncols)]

        carray_empty[:, :, :] = wssd_contig
        del(wssd_contig)

    fout.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("outfile")
    parser.add_argument("--infiles", required=True, nargs="+")

    args = parser.parse_args()

    contigs = {}

    for infile in args.infiles:
        with open(infile, "r") as file:
            dat = pickle.load(file)

            for contig, matrix in dat.items():
                if contig not in contigs:
                    contigs[contig] = matrix
                else:
                    contigs[contig] += matrix
                del(matrix)
