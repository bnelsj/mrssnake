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

def convert_matrix(matrix):
    if issparse(matrix):
        matrix = matrix.toarray()
    if isinstance(matrix, np.ndarray):
        return matrix
    else:
        print("Error: unrecognized data type for contig %s: %s" % (contig, matrix.__class__.__name__), file=sys.stderr, flush=True)
        sys.exit(1)

def add_contents_to_contigs(dat, contigs):
    """Take a dictionary-like object of matrices, add matrices to contigs dictionary.
       Converts matrices to np.ndarray automatically.
    """
    for contig, matrix in dat.items():
        matrix = convert_matrix(matrix)
        if contig not in contigs:
            contigs[contig] = matrix
        else:
            contigs[contig] += matrix
    return contigs

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

def load_matrices_per_contig(infiles, contig):
    """Get counts from all infiles for a given contig dictionary.
    """
    contig_name = list(contig)[0]
    for i, infile in enumerate(infiles):
        with shelve.open(infile) as dat:
            print("Contig %s: loading shelve %d of %d: %s" % (contig_name, i+1, len(infiles), infile), file=sys.stdout, flush=True)
            
            if contig_name in dat:
                matrix = convert_matrix(dat["contig_name"])
                if contig[contig_name] is None:
                    contig[contig_name] = matrix
                else:
                    contig[contig_name] += matrix
    return contig

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
    parser.add_argument("--contigs_file", default = None, help="Tab-delimited table with contig names in the first column")
    parser.add_argument("--per_contig_merge", action="store_true", help="Merge matrices one contig at a time (low memory footprint)")

    args = parser.parse_args()

    if args.per_contig_merge:
        if args.contigs_file is None:
            print("Must specify --contigs file for per_contig_merge", file = sys.stderr)
            sys.exit(1)

    start_time = time.time()

    fout = tables.open_file(args.outfile, mode="w")
    print("Successfully opened outfile: %s" % args.outfile, file=sys.stdout, flush=True)

    contigs = {}
    contig_list = []

    if args.per_contig_merge:
        with open(args.contigs_file, "r") as contigs_file:
            for line in contigs_file:
                contig = line.rstrip().split()[0]
                contigs_list.append(contig)

    infiles = []

    if args.infile_glob is not None:
        infiles.extend(glob.glob(args.infile_glob))

    if args.infiles is not None:
        infiles.extend(args.infiles)

    if args.per_contig_merge:
        for contig_name in contigs_list:
            contig = {}
            contig[contig_name] = None
            contig = load_matrices_per_contig(infiles, contig)
            write_to_h5(contig, fout)
    else:
        if args.live_merge:
            contigs = load_matrices_live(infiles, contigs)
        else:
            contigs = load_matrices_post(infiles, contigs)
        print("Finished loading shelves. Creating h5 file: %s" % args.outfile, file=sys.stdout, flush=True)
        write_to_h5(contigs, fout)

    finish_time = time.time()
    print("Finished writing wssd_out_file in %d seconds. Closing." % (finish_time - start_time), file=sys.stdout, flush=True)
    fout.close()
