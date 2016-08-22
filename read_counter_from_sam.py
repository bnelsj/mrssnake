# cython: infer_types=True
# cython: boundscheck=False
# cython: wraparound=False
"""
Convert mrsfast sam output to matrix of read depth and read start depth values.
Uses sparse matrices to reduce memory footprint.
"""

from __future__ import print_function
from __future__ import division

import sys
import argparse

from functools import total_ordering
import time

import re

import tables
import numpy as np
import pandas as pd
from tables import NoSuchNodeError

def create_array(contig, contigs_file, max_edist):
    """Create global numpy array 'matrix' for given contig.
    shape = (contig length, nedists, 2)
    third dimension is for depth and starts, respectively.
    """

    contig_dat = pd.read_table(contigs_file, header=None, names=["contig", "length"])
    contig_dat.index = contig_dat.contig

    contig_length = contig_dat.ix[contig_dat.contig == contig, "length"]
    global matrix
    matrix = np.ndarray((contig_length, max_edist+1, 2), dtype=np.uint32)

def add_to_array(pos, edist, rlen=36):
    """Add read entry to contig. Assumes pos is 0-based.
    """

    end = pos + rlen
    matrix[pos:end, edist, 0] += 1
    matrix[pos, edist, 1] += 1

def process_samfile(samfile, contig, max_edist, rlen=36):
    regex_full = re.compile("[^ @\t]+\t[0-9]+\t(%s)\t([0-9]+)\t.+NM:i:([0-9]+)" % contig)
    nhits = 0

    with open(samfile, "r") as sam:
        for line in sam:
            match = regex_full.match(line)
            if match is not None:
                contig, pos, edist = match.group(1,2,3)
                pos, edist = int(pos), int(edist)
                add_to_array(pos, edist, rlen=36)
                nhits += 1

    return nhits

def write_to_h5(contig, fout_handle, chunksize=1000000):
    """Write counts (dictionary of contig matrices) to fout hdf5 file
    in increments of chunksize bases. Outfile is in wssd_out_file format.
    """
    try:
        group = fout_handle.get_node(fout_handle.root, "depthAndStarts_wssd")
    except NoSuchNodeError:
        group = fout_handle.create_group(fout_handle.root, "depthAndStarts_wssd")
    finally:
        carray_empty = tables.CArray(group,
                                     contig,
                                     tables.UInt32Atom(),
                                     matrix.shape,
                                     filters=tables.Filters(complevel=1, complib="lzo")
                                    )

        contig_len = matrix.shape[0]
        nchunks = contig_len // chunksize
        if nchunks * chunksize < contig_len:
            nchunks += 1

        for i in range(nchunks):
            s = i * chunksize
            e = s + chunksize
            if e > contig_len:
                e = contig_len
            carray_empty[s:e, :, :] = matrix[s:e, :, :]

        fout_handle.flush()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("samfile", help="Sam-format input file")
    parser.add_argument("outfile", help="HDF5 file with read depths and starts")
    parser.add_argument("contig", help="Name of contig")
    parser.add_argument("contig_lengths", help="tab-delimited file with contig names and lengths")
    parser.add_argument("--max_edist",
                        default=2,
                        type=int,
                        help="Maximum edit distance of input reads"
                       )
    parser.add_argument("--read_length", default=36, type=int, help="Length of input reads")
    parser.add_argument("--log", default=sys.stderr, help="Path to log file. Default: sys.stderr")

    args = parser.parse_args()

    run_start = time.time()

    if args.log is not sys.stderr:
        logfile = open(args.log, "w")
    else:
        logfile = sys.stderr

    create_array(args.contig, args.contig_lengths, args.max_edist)

    nhits = process_samfile(args.samfile, args.contig, args.max_edist)

    mp_end = time.time() - run_start

    print("Finished reading samfile in %d seconds." % mp_end, file=logfile)
    print("Contig %s had %d hits" % (args.contig, nhits), file=logfile)
    print("Writing to hdf5: %s" % args.outfile, file=logfile)
    
    with tables.open_file(args.outfile, mode="w") as h5file:
        write_to_h5(args.contig, h5file)
        del matrix

    run_end = time.time()
    total_runtime = run_end - run_start

    print("Counter: finished writing matrices to hdf5: %s" % args.outfile, file=logfile)
    print("Counter: total runtime: %d seconds" % total_runtime, file=logfile, flush=True)
    if logfile is not sys.stderr:
        logfile.close()
    sys.exit()
