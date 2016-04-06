"""
This module contains functions for merging read count shelves into pytables wssd_out_files.
Merging can be done per contig and before all shelves are created.
"""

from __future__ import print_function
from __future__ import division

import shelve

import time
import os
import sys
import argparse
import glob

from dbm import error

import tables
from tables.exceptions import NoSuchNodeError
import numpy as np
from scipy.sparse import issparse

def convert_matrix(matrix):
    """
    Convert matrix to np.ndarray if it is sparse.
    """
    if issparse(matrix):
        matrix = matrix.toarray()
    if isinstance(matrix, np.ndarray):
        return matrix
    else:
        print("Error: unrecognized data type for matrix: %s" %
              (matrix.__class__.__name__),
              file=sys.stderr, flush=True)
        sys.exit(1)

def add_contents_to_contig(dat, contig):
    """Take a dictionary-like object of matrices, add matrices to contig dictionary.
       Converts matrices to np.ndarray automatically.
    """
    contig_string = list(contig)[0]
    if contig_string in dat:
        matrix = convert_matrix(dat[contig_string])
        if contig[contig_string] is None:
            contig[contig_string] = matrix
        else:
            contig[contig_string] += matrix
    return contig

def load_matrices_per_contig_live(matrices, contig):
    """Get counts from all matrices for a given contig dictionary as they are finished.
    """
    fileset = set(matrices)
    total_infiles = len(fileset)
    processed_infiles = set()
    while len(fileset) > 0:
        for infile in fileset:
            # Check if infile exists and hasn't been modified in 5 minutes
            # Adds .dat extension for shelve compatibility
            if os.path.isfile(infile + ".dat") and \
               time.time() - os.path.getmtime(infile + ".dat") > 300:
                try:
                    dat = shelve.open(infile, flag="r")
                except error as err:
                    print("Error: %s: %s" % (infile, str(err)), file=sys.stderr, flush=True)
                    continue
                else:
                    contig = add_contents_to_contig(dat, contig)
                    dat.close()
                    processed_infiles.add(infile)
                    print("Loaded shelve %d of %d: %s" %
                          (len(processed_infiles), total_infiles, infile),
                          file=sys.stdout, flush=True)
        fileset -= processed_infiles
        print("Checked all infiles. Sleeping 30s...", file=sys.stderr, flush=True)
        time.sleep(30)
    return contig

def load_matrices_per_contig(matrices, contig):
    """Get counts from all matrices for a given contig dictionary.
    """
    contig_string = list(contig)[0]
    for i, infile in enumerate(matrices):
        print("Contig %s: loading shelve %d of %d: %s" %
              (contig_string, i+1, len(matrices), infile),
              file=sys.stdout, flush=True)
        with shelve.open(infile, flag="r") as dat:
            contig = add_contents_to_contig(dat, contig)
    return contig

def write_to_h5(counts, fout_handle):
    """Write counts (dictionary of contig matrices) to fout hdf5 file.
       Outfile is in wssd_out_file format.
    """
    try:
        group = fout_handle.get_node(fout_handle.root, "depthAndStarts_wssd")
    except NoSuchNodeError:
        group = fout_handle.create_group(fout_handle.root, "depthAndStarts_wssd")
    finally:
        for i, (contig, matrix) in enumerate(counts.items()):
            print("Merger: %d Creating array for %s" %(i+1, contig), file=sys.stdout, flush=True)
            nrows, ncols = matrix.shape
            nedists = nrows // 2
            wssd_contig = matrix.T

            carray_empty = tables.CArray(group,
                                         contig,
                                         tables.UInt32Atom(),
                                         (ncols, nedists, 2),
                                         filters=tables.Filters(complevel=1, complib="lzo")
                                        )

            # Add depth counts
            carray_empty[:, :, 0] = wssd_contig[:, nedists:]

            # Add starts
            carray_empty[:, :, 1] = wssd_contig[:, 0:nedists]

            fout_handle.flush()

def write_wssd_to_h5(wssd_handle, fout_handle):
    """Append single contig wssd_out_file to fout hdf5 file.
       Outfile is in wssd_out_file format.
       Exits with error if multiple contigs are in a single wssd file.
    """
    try:
        group = fout_handle.get_node(fout_handle.root, "depthAndStarts_wssd")
    except NoSuchNodeError:
        group = fout_handle.create_group(fout_handle.root, "depthAndStarts_wssd")
    finally:
        nodes = wssd_handle.list_nodes("/depthAndStarts_wssd/")
        if len(nodes) > 0:
            for matrix in nodes:
                wssd_handle.copy_node(matrix, newparent=group)
                fout_handle.flush()
        else:
            print("Empty wssd file", file=sys.stderr)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("outfile", help="Path to output wssd_out_file")
    parser.add_argument("--infiles", nargs="+", default=None, help="List of input shelves to merge")
    parser.add_argument("--infile_glob", default=None, help="glob string for infiles")
    parser.add_argument("--live_merge", action="store_true",
                        help="Start merging infiles before they are all finished? \
                              (Default: %(default)s)")
    parser.add_argument("--contigs_file", default=None,
                        help="Tab-delimited table with contig names in the first column")
    parser.add_argument("--contig", default=None, help="Name of contig to merge")
    parser.add_argument("--wssd_merge", nargs="+", default=None,
                        help="Merge multiple wssd_out_files")

    args = parser.parse_args()

    contig_dict = {}
    contig_list = []

    if not args.wssd_merge:
        if args.contigs_file is None and args.contig is None:
            print("Error: Must specify --contigs_file or --contig for shelve merge",
                  file=sys.stderr)
            sys.exit(1)
        else:
            if args.contig is not None:
                contig_list.append(args.contig)
            if args.contigs_file is not None:
                with open(args.contigs_file, "r") as contigs_file:
                    for line in contigs_file:
                        contig_name = line.rstrip().split()[0]
                        contig_list.append(contig_name)
    else:
        if args.infile_glob is not None or args.infiles is not None:
            print("Error: cannot specify infile_glob or infiles with wssd_merge", file=sys.stderr)
            sys.exit(1)


    if args.live_merge:
        load_func = load_matrices_per_contig_live
    else:
        load_func = load_matrices_per_contig_live

    start_time = time.time()

    infiles = []

    if args.infile_glob is not None:
        infiles.extend(glob.glob(args.infile_glob))

    if args.infiles is not None:
        infiles.extend(args.infiles)

    # Remove extensions and get unique shelves
    if infiles != []:
        infiles = [x.replace(".dat", "").replace(".bak", "").replace(".dir", "") for x in infiles]
        infiles = list(set(infiles))

    if args.infile_glob is not None or args.infiles is not None:
        if infiles == []:
            print("No infiles found. Exiting...", file=sys.stderr)
            sys.exit(1)

    if args.wssd_merge is None:
        with tables.open_file(args.outfile, mode="w") as fout:
            print("Successfully opened outfile: %s" % args.outfile, file=sys.stdout, flush=True)
            for contig_name in contig_list:
                contig_dict = {contig_name: None}
                contig_dict = load_func(infiles, contig_dict)
                write_to_h5(contig_dict, fout)

    else:
        # Merge wssd files into single wssd_out_file
        with tables.open_file(args.outfile, mode="w") as fout:
            print("Successfully opened outfile: %s" % args.outfile, file=sys.stdout, flush=True)
            for wssd_file in args.wssd_merge:
                print("Reading wssd_file: %s" % wssd_file, file=sys.stdout, flush=True)
                with tables.open_file(wssd_file, mode="r") as wssd:
                    write_wssd_to_h5(wssd, fout)

    finish_time = time.time()
    print("Finished writing wssd_out_file in %d seconds. Closing." %
          (finish_time - start_time),
          file=sys.stdout,
          flush=True
         )
