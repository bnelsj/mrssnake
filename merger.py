from __future__ import print_function
from __future__ import division

import shelve

import time
import sys
import argparse
import glob
import tables
import numpy as np
from scipy.sparse import issparse

def add_contents_to_contigs(dat, contigs):
    """Take a dictionary-like object of matrices, add matrices to contigs dictionary.
       Converts matrices to np.ndarray automatically.
    """
    for contig, matrix in dat.items():
        if issparse(matrix):
            matrix = matrix.toarray()
        if isinstance(matrix, np.ndarray):
            if contig not in contigs:
                contigs[contig] = matrix
            else:
                contigs[contig] += matrix
        else:
            print("Error: unrecognized data type for contig %s: %s" % (contig, matrix.__class__.__name__), file=sys.stderr, flush=True)
            sys.exit(1)

        del(matrix)
    return(contigs)

def load_matrices_live(infiles, contigs):
    """Load infiles to contigs dictionary as they are finished."""
    infiles = set(infiles)
    total_infiles = len(infiles)
    processed_infiles = set()
    while len(infiles) > 0:
        for infile in infiles:
            # Check if infile exists and hasn't been modified in 5 minutes
            if os.path.isfile(infile) and time.time() - os.path.getmtime(infile) > 300:
                try:
                    dat = shelve.open(infile, flag="r")
                except error as e:
                    print("Error: %s: %s" % (infile, str(e)), file=sys.stderr, flush=True)
                    continue
                else:
                    contigs = add_contents_to_contigs(dat, contigs)
                    dat.close()
                    processed_infiles.add(infile)
                    print("Loaded pickle %d of %d: %s" % (len(processed_infiles), total_infiles, infile), file=sys.stdout, flush=True)
        infiles -= processed_infiles
        time.sleep(30)
    return contigs

def load_matrices_post(infiles, contigs):
    """Load infiles to contigs dictionary. Assumes all infiles are complete."""
    for i, infile in enumerate(infiles):
        with shelve.open(infile) as dat:
            print("Loading shelve %d of %d: %s" % (i+1, len(infiles), infile), file=sys.stdout, flush=True)

            contigs = add_contents_to_contigs(dat, contigs)
    return contigs

def write_to_h5(counts, fout):
    """Write counts (dictionary of contig matrices) to fout hdf5 file.
       Outfile is in wssd_out_file format.   
    """
    group = fout.create_group(fout.root, "depthAndStarts_wssd")

    for i, (contig, matrix) in enumerate(counts.items()):
        print("Merger: %d Creating array for %s" %(i+1, contig), file=sys.stdout, flush=True)
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
    parser.add_argument("outfile", help= "Path to output wssd_out_file")
    parser.add_argument("--infiles", nargs="+", default = None, help = "List of input shelves to merge")
    parser.add_argument("--infile_glob", default = None, help = "glob string for infiles")
    parser.add_argument("--live_merge", default = False, help = "Start merging infiles before they are all finished? (Default: %(default)s)")

    args = parser.parse_args()

    start_time = time.time()

    fout = tables.open_file(args.outfile, mode="w")
    print("Successfully opened outfile: %s" % args.outfile, file=sys.stdout, flush=True)
    contigs = {}

    infiles = []

    if args.infile_glob is not None:
        infiles.extend(glob.glob(args.infile_glob))

    if args.infiles is not None:
        infiles.extend(args.infiles)

    if args.live_merge:
        contigs = load_matrices_live(infiles, contigs)
    else:
        contigs = load_matrices_post(infiles, contigs)

    print("Finished loading shelves. Creating h5 file: %s" % args.outfile, file=sys.stdout, flush=True)

    write_to_h5(contigs, fout)
    finish_time = time.time()
    print("Finished writing wssd_out_file in %d seconds. Closing." % (finish_time - start_time), file=sys.stdout, flush=True)
    fout.close()
