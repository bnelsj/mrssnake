"""
This module contains functions for merging read count hdf5 files into wssd_out_files using pytables.
Merging can be done per contig and before all files are created.
"""

from __future__ import print_function
from __future__ import division

import time
import sys
import argparse
import glob

import tables
#import numpy as np
from scipy.sparse import csr_matrix

def merge_contigs_to_wssd(infile_list, fout_handle):
    """
        Create a wssd_out_file by copying nodes from a list of single contig
        wssd files.
    """
    group = fout_handle.create_group(fout_handle.root, "depthAndStarts_wssd")
    dest_contigs = []
    print("Successfully opened outfile: %s" % args.outfile, file=sys.stdout, flush=True)
    for file in infile_list:
        with tables.open_file(file, mode="r") as h5_handle:
            nodes = [node for node in h5_handle.list_nodes]
            if len(nodes) == 0:
                print("Error: no contigs in infile: %s" % (file), file=sys.stderr)
                sys.exit(1)
            if len(nodes) > 1:
                print("Error: more than one contig (%d contigs) in infile: %s" % (len(nodes), file), file=sys.stderr)
                sys.exit(1)
            contig = nodes[0]
            if contig.name in dest_contigs:
                print("Error: contig %s in multiple infiles including: %s" % (contig.name, file), file=sys.stderr)
                sys.exit(1)
            h5_handle.copy_node(contig, newparent=group)
            fout_handle.flush()

def load_sparse_matrix(name, root):
    """
        Adapted from:
        http://stackoverflow.com/questions/11129429/storing-numpy-sparse-matrix-in-hdf5-pytables
    """
    pars = []
    for par in ('data', 'indices', 'indptr', 'shape'):
        pars.append(getattr(root, '%s_%s' % (name, par)).read())
    mat = csr_matrix(tuple(pars[:3]), shape=pars[3])
    return mat

def merge_sparse_h5_to_wssd(infile_list, contig_list, fout_handle):
    """
        Create a wssd_out_file with depth and counts for each contig in contig list
        for each infile in infile_list. Assumes the h5 files contain sparse csr_matrices.
    """

    out_group = fout_handle.create_group(fout_handle.root, "depthAndStarts_wssd")
    for contig in contig_list:
        contig_array = None
        for h5file in infile_list:
            print("Reading contig %s from h5file: %s" % (contig, h5file),
                  file=sys.stdout,
                  flush=True)
            with tables.open_file(h5file, mode="r") as fin:
                in_group = fin.get_node(fin.root, "depthAndStarts_wssd")
                job_matrix = load_sparse_matrix(contig, in_group)
                if contig_array is None:
                    contig_array = job_matrix.toarray()
                else:
                    contig_array += load_sparse_matrix(contig, in_group)

        if contig_array is None:
            print("Contig %s not found in infiles" % contig)
        print("Writing contig %s to h5file" % (contig))
        nrows, ncols = contig_array.shape
        nedists = nrows // 2
        wssd_contig = contig_array.T
        carray_empty = tables.CArray(out_group,
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
        del carray_empty
        del wssd_contig
        del contig_array

def merge_h5_to_wssd(infile_list, contig_list, fout_handle):
    """
        Create a wssd_out_file with depth and counts for each contig in contig list
        for each infile in infile_list.
    """
    group = fout_handle.create_group(fout_handle.root, "depthAndStarts_wssd")
    print("Successfully opened outfile: %s" % args.outfile, file=sys.stdout, flush=True)
    for contig in contig_list:
        contig_array = None
        for h5file in infile_list:
            print("Reading contig %s from h5file: %s" % (contig, h5file),
                  file=sys.stdout,
                  flush=True)
            with tables.open_file(h5file, mode="r") as h5_handle:
                source_nodes = [node.name for node in h5_handle.list_nodes("/depthAndStarts_wssd")]
                if contig in source_nodes:
                    if contig_array is None:
                        contig_array = h5_handle.get_node("/depthAndStarts_wssd", name=contig)[:]
                    else:
                        contig_array += h5_handle.get_node("/depthAndStarts_wssd", name=contig)[:]

        if contig_array is not None:
            carray = fout_handle.create_carray(group, contig, obj=contig_array,
                                               filters=tables.Filters(complevel=1, complib="lzo"))
            fout_handle.flush()
        else:
            print("Contig %s not found in infiles" % contig)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("outfile", help="Path to output wssd_out_file")
    parser.add_argument("--infiles", nargs="+", default=None,\
                        help="List of input hdf5 files to merge")
    parser.add_argument("--infile_glob", default=None, help="glob string for infiles")
    parser.add_argument("--live_merge", action="store_true",
                        help="Start merging infiles before they are all finished? \
                              (Default: %(default)s)")
    parser.add_argument("--contigs_file", default=None,
                        help="Tab-delimited table with contig names in the first column")
    parser.add_argument("--contig", default=None, help="Name of contig to merge")
    parser.add_argument("--wssd_merge", action="store_true", help="Merge wssd contigs to wssd_out_file")

    args = parser.parse_args()

    contig_dict = {}
    contigs = []

    if args.contigs_file is None and args.contig is None:
        print("Error: Must specify --contigs_file or --contig for merge",
              file=sys.stderr)
        sys.exit(1)
    
    if args.contig is not None:
        contigs.append(args.contig)
    if args.contigs_file is not None:
        with open(args.contigs_file, "r") as contigs_file:
            for line in contigs_file:
                contig_name = line.rstrip().split()[0]
                contigs.append(contig_name)

    start_time = time.time()

    infiles = []

    if args.infile_glob is not None:
        infiles.extend(glob.glob(args.infile_glob))

    if args.infiles is not None:
        infiles.extend(args.infiles)

    if args.infile_glob is not None or args.infiles is not None:
        if infiles == []:
            print("No infiles found. Exiting...", file=sys.stderr)
            sys.exit(1)

    if args.wssd_merge:
        with tables.open_file(args.outfile, mode="w") as fout:
            merge_contigs_to_wssd(infiles, fout)
    else:
        # Merge hdf5 files into single wssd_out_file
        with tables.open_file(args.outfile, mode="w") as fout:
            merge_sparse_h5_to_wssd(infiles, contigs, fout)

    finish_time = time.time()
    print("Finished writing wssd_out_file in %d seconds. Closing." %
          (finish_time - start_time),
          file=sys.stdout,
          flush=True
         )
