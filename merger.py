from __future__ import print_function
from __future__ import division

try:
    import cPickle as pickle
except ImportError:
    import pickle

import sys
import argparse
import tables
import numpy as np
from scipy.sparse import lil_matrix

def write_to_h5(counts, fout):
    group = fout.create_group(fout.root, "depthAndStarts_wssd")

    for i, (contig, matrix) in enumerate(counts.items()):
        sys.stdout.write("Merger: %d Creating array for %s\n" %(i+1, contig))
        sys.stdout.flush()
        nrows, ncols = matrix.shape
        nedists = nrows // 2
        carray_empty = tables.CArray(group, contig, tables.UInt32Atom(), (ncols, nedists, 2), filters=tables.Filters(complevel=1, complib="lzo"))

        wssd_contig = matrix.toarray().T

        # Add depth counts
        carray_empty[:, :, 0] = wssd_contig[:, nedists:]

        # Add starts
        carray_empty[:, :, 1] = wssd_contig[:, 0:nedists]

        del(wssd_contig)

    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("outfile")
    parser.add_argument("--infiles", required=True, nargs="+")

    args = parser.parse_args()

    fout = tables.open_file(args.outfile, mode="w")
    sys.stdout.write("Successfully opened outfile %s\n" % args.outfile)
    contigs = {}

    ninfiles = len(args.infiles)

    for i, infile in enumerate(args.infiles):
        with open(infile, "rb") as file:
            sys.stdout.write("Loading pickle %d of %d: %s\n" % (i+1, ninfiles, infile))
            sys.stdout.flush()
            dat = pickle.load(file)

            for contig, matrix in dat.items():
                if contig not in contigs:
                    contigs[contig] = matrix.tocsr()
                else:
                    contigs[contig] += matrix.tocsr()
                del(matrix)

    sys.stdout.write("Finished loading pickles. Creating h5 file: %s\n" % args.outfile)
    sys.stdout.flush()

    for contig, array in contigs.items():  
        print("%s: %d" % (contig, np.count_nonzero(array.todense())))
    write_to_h5(contigs, fout)
    fout.close()
