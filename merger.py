from __future__ import print_function
from __future__ import division

try:
    import cPickle as pickle
except ImportError:
    import pickle

import time
import sys
import argparse
import tables
import numpy as np
from scipy.sparse import lil_matrix, csr_matrix

def write_to_h5(counts, fout):
    group = fout.create_group(fout.root, "depthAndStarts_wssd")

    for i, (contig, matrix) in enumerate(counts.items()):
        sys.stdout.write("Merger: %d Creating array for %s\n" %(i+1, contig))
        sys.stdout.flush()
        nrows, ncols = matrix.shape
        nedists = nrows // 2
        carray_empty = tables.CArray(group, contig, tables.UInt32Atom(), (ncols, nedists, 2), filters=tables.Filters(complevel=1, complib="lzo"))

        wssd_contig = matrix.T

        # Add depth counts
        carray_empty[:, :, 0] = wssd_contig[:, nedists:]

        # Add starts
        carray_empty[:, :, 1] = wssd_contig[:, 0:nedists]

        fout.flush()

        del(wssd_contig)

    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("outfile")
    parser.add_argument("--infiles", required=True, nargs="+")

    args = parser.parse_args()

    start_time = time.time()

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
                if isinstance(matrix, csr_matrix) or isinstance(matrix, lil_matrix):
                    if contig not in contigs:
                        contigs[contig] = matrix.toarray()
                    else:
                        contigs[contig] += matrix.toarray()
                elif isinstance(matrix, np.ndarray):
                    if contig not in contigs:
                        contigs[contig] = matrix
                    else:
                        contigs[contig] += matrix
                else:
                    print("Error: unrecognized data type for contig %s: %s" % (contig, matrix.__class__.__name__), file=sys.stderr, flush=True)
                    sys.exit(1)

                del(matrix)
    
    print("Finished loading pickles. Creating h5 file: %s" % args.outfile, file=sys.stdout, flush=True)

    write_to_h5(contigs, fout)
    finish_time = time.time()
    total_time = finish_time - start_time
    print("Finished writing wssd_out_file in %d seconds. Closing." % total_time, file=sys.stdout, flush=True)
    fout.close()
