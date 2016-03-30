"""
Convert mrsfast sam output to matrix of read depth and read start depth values.
Uses sparse matrices to reduce memory footprint.
"""

from __future__ import print_function
from __future__ import division

from pickle import HIGHEST_PROTOCOL

import shelve

import sys
import argparse

from functools import total_ordering
import time

import numpy as np

from scipy.sparse import lil_matrix, bsr_matrix
import pysam

@total_ordering
class Contig:
    """
    Stores name, size, reads, and total reads for a contig
    and implements ordering based on number of reads
    """
    def __init__(self, contig_name, contig_size):
        self.name = contig_name
        self.size = contig_size
        self.total_reads = 0
        self.reads = 0
    def __eq__(self, other):
        return self.reads / self.size == other.reads / other.size

    def __lt__(self, other):
        return self.reads / self.size < other.reads / other.size

class ContigManager:
    """
    Stores list of contigs and manages which contigs use numpy arrays
    and which ones use sparse matrices.

    Contains methods for switching contigs between numpy arrays
    and sparse matrices based on recent read density.
    """
    def __init__(self, max_bases, contigs_seen=None, array_contigs=None):
        self.max_bases = max_bases
        if contigs_seen is not None:
            self.contigs_seen = contigs_seen
        else:
            self.contigs_seen = {}
        if array_contigs is not None:
            self.array_contigs = array_contigs
        else:
            self.array_contigs = []
        self.used_bases = sum([contigs_seen[contig].size for contig in self.array_contigs])

    def add_contig(self, contig, contig_size=None):
        """
        Add contig to contig manager
        """
        if contig_size is not None:
            contig = Contig(contig, contig_size)
        if contig.name not in self.contigs_seen.keys():
            self.contigs_seen[contig.name] = contig
        if self.used_bases + contig.size <= self.max_bases and \
        contig.name not in self.array_contigs:
            self.array_contigs.append(contig.name)
            self.used_bases += contig.size

    def add_contig_to_array_contigs(self, contig, contig_size=None):
        """
        Switch contig read matrix to numpy array
        """
        if contig_size is not None:
            contig = Contig(contig, contig_size)
        self.add_contig(contig)
        if contig.name not in self.array_contigs:
            self.array_contigs.append(contig.name)
            self.used_bases += contig.size

    def reset_read_counts(self):
        """
        Reset the read counts for all contigs in ContigManager
        """
        for name, contig in self.contigs_seen.items():
            contig.total_reads += contig.reads
            contig.reads = 0

    def rebalance(self):
        """
        Use recent read density to reassign contigs to numpy
        or sparse matrices
        """
        new_used_bases = 0
        new_array_contigs = []
        for name, contig in sorted(self.contigs_seen.items(), key=lambda x: x[1], reverse=True):
            contig = self.contigs_seen[name]
            if new_used_bases + contig.size <= self.max_bases:
                new_array_contigs.append(contig.name)
                new_used_bases += contig.size

        self.array_contigs = new_array_contigs
        self.used_bases = new_used_bases
        self.reset_read_counts()


def get_array_contigs(contigs, args):
    """
    Get names of contigs to store in ContigManager
    """
    array_contigs = []
    if args.all_contigs:
        array_contigs = contigs.keys()

    else:
        if args.common_contigs is not None:
            array_contigs.extend(args.common_contigs)
        if args.noncanonical_contigs:
            canonical = ["chr%s" % str(x) for x in list(range(1, 25)) + ["X", "Y", "M"]]
            array_contigs.extend([contig for contig in contigs.keys() if contig not in canonical])
    return array_contigs

def update_read_depth_and_start(matrix, edist, start, end, nedists=3):
    """
    Update read matrix based on start and edist of new read
    """
    if isinstance(matrix, np.ndarray):
        matrix[edist + nedists, start:end] += 1
    else:
        mat_slice = matrix[edist + nedists, start:end].toarray()
        mat_slice += 1
        matrix[edist + nedists, start:end] = mat_slice

    matrix[edist, start] += 1

    return matrix

#def count_reads_sans_pysam(input, contigs, array_contigs, args):
#    read_dict = {}
#    nrows = 2 * args.max_edist + 2
#    nstart_rows = nrows // 2
#
#    with open(input, "r") as infile:
#        for line in infile:
#            if line.startswith("@"):
#                continue
#            contig, start, cigar, edist = line.split()[2,3,5,11]
#
#            if edist > args.max_edist:
#                continue
#
#            start = start - 1
#            rlen = int(cigar[:-1])
#            end = start + rlen
#            edist = int(edist[4:])
#
#            if contig not in read_dict:
#                length = contigs[contig]
#                if contig in array_contigs:
#                    read_dict[contig] = np.zeros((nrows, length), dtype=np.uint16)
#                    sys.stderr.write("Counter: %s (%d, %d) numpy array\n" %
#                                     (contig, nrows, length))
#                    sys.stderr.flush()
#                else:
#                    read_dict[contig] = lil_matrix((nrows, length), dtype=np.uint16)
#                    sys.stderr.write("Counter: %s scipy lil_matrix\n" % contig)
#                    sys.stderr.flush()
#
#            if end > contigs[contig]:
#                end = contigs[contig]
#
#            if contig in array_contigs:
#                # Update read depth counts for numpy array
#                read_dict[contig][edist + nstart_rows, start:end] += 1
#            else:
#                # Update read depth counts for sparse matrix
#                slice = read_dict[contig][edist + nstart_rows, start:end].toarray()
#                slice += 1
#                read_dict[contig][edist + nstart_rows, start:end] = slice
#
#            # Update read start counts
#            read_dict[contig][edist, start] += 1
#
#            if i % 1000000 == 0:
#                sys.stderr.write("Counter: %d reads processed\n" % i)
#                sys.stderr.flush()
#
#    return read_dict

def rebalance_contigs(contig_manager, read_dict):
    """
    Reassign contigs to numpy or sparse matrices based on recent read depth
    """
    if args.max_basepairs_in_mem > 0:
        print("Counter: rebalancing array contigs...", file=logfile)
        start_time = time.time()
        old_contigs = contig_manager.array_contigs
        contig_manager.rebalance()
        print("Contigs removed:", " ".join(
            [contig for contig in old_contigs
             if contig not in contig_manager.array_contigs]
            ),
              sep=" ",
              file=logfile
             )
        print("Contigs added:  ", " ".join(
            [contig for contig in contig_manager.array_contigs
             if contig not in old_contigs]), sep=" ", file=logfile, flush=True
             )
        for contig, array in read_dict.items():
            if contig not in contig_manager.array_contigs:
                if not isinstance(array, lil_matrix):
                    read_dict[contig] = lil_matrix(array)
                    del array
            else:
                if not isinstance(array, np.ndarray):
                    read_dict[contig] = array.toarray()
                    del array
        finish_time = time.time()
        total_time = finish_time - start_time
        print("Counter: rebalancing finished in %d sec..." %
              total_time,
              file=logfile,
              flush=True
             )


def count_reads(samfile, contig_manager, args):
    """
    Update matrices for each contig based on mrsfast read mappings
    """
    read_dict = {}
    nrows = 2 * args.max_edist + 2
    nstart_rows = nrows // 2

    for i, read in enumerate(samfile):
        contig = samfile.getrname(read.reference_id)
        start = read.reference_start
        end = read.reference_end

        edist = read.get_tag("NM")
        if edist > args.max_edist:
            continue

        if contig not in read_dict:
            length = contig_manager.contigs_seen[contig].size
            if contig in contig_manager.array_contigs:
                read_dict[contig] = np.zeros((nrows, length), dtype=np.uint16)
                print("Counter: %s (%d, %d) numpy array" %
                      (contig, nrows, length),
                      file=logfile,
                      flush=True
                     )
            else:
                read_dict[contig] = lil_matrix((nrows, length), dtype=np.uint16)
                print("Counter: %s scipy lil_matrix" % contig, file=logfile, flush=True)


        if isinstance(read_dict[contig], np.ndarray):
            # Update read depth counts for numpy array
            read_dict[contig][edist + nstart_rows, start:end] += 1
        else:
            # Update read depth counts for sparse matrix
            mat_slice = read_dict[contig][edist + nstart_rows, start:end].toarray()
            mat_slice += 1
            read_dict[contig][edist + nstart_rows, start:end] = mat_slice

        # Update read start counts
        read_dict[contig][edist, start] += 1
        contig_manager.contigs_seen[contig].reads += 1

        if i % args.rebalance_freq == 0 and i > 0:
            print("Counter: %d reads processed" % i, file=logfile, flush=True)
            rebalance_contigs(contig_manager, read_dict)
           
    return read_dict


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", default="/dev/stdin")
    parser.add_argument("outfile", default="wssd_out_file")
    parser.add_argument("contig_lengths", help="tab-delimited file with contig names and lengths")
    parser.add_argument("--max_edist",
                        default=2,
                        type=int,
                        help="Maximum edit distance of input reads"
                       )
    parser.add_argument("--common_contigs",
                        default=[],
                        nargs="+",
                        help="Create numpy array for common contigs (Much faster, more memory)")
    parser.add_argument("--all_contigs",
                        action="store_true",
                        help="Create numpy array for all contigs (Fast, high mem requirement)")
    parser.add_argument("--noncanonical_contigs",
                        action="store_true",
                        help="Create numpy array for all noncanonical contigs \
                        (Fast, high mem requirement).")
    parser.add_argument("--max_basepairs_in_mem",
                        type=int,
                        default=0,
                        help="Number of reference bp to keep in memory")
    parser.add_argument("--rebalance_freq",
                        type=int,
                        default=100000,
                        help="Reads between contig array rebalancing"
                       )
    parser.add_argument("--log", default=sys.stderr, help="Path to log file. Default: sys.stderr")

    args = parser.parse_args()

    run_start = time.time()

    if args.log is not sys.stderr:
        logfile = open(args.log, "w")
    else:
        logfile = sys.stderr

    contigs = {}

    with open(args.contig_lengths, "r") as reader:
        for line in reader:
            contig, length = line.rstrip().split()
            contigs[contig] = int(length)

    numpy_contigs = get_array_contigs(contigs, args)

    contig_manager = ContigManager(args.max_basepairs_in_mem)

    # Add array contigs first
    for name in numpy_contigs:
        if name not in contigs:
            continue
        size = contigs[name]
        contig_manager.add_contig_to_array_contigs(name, size)

    # Add remaining contigs to array_contigs if they fit
    for name, size in contigs.items():
        if name in numpy_contigs:
            continue
        else:
            contig_manager.add_contig(name, size)

    print("Max bases: ", contig_manager.max_bases, "Used bases: ", contig_manager.used_bases,
          file=logfile
         )

    samfile = pysam.AlignmentFile(args.infile, "r", check_sq=False)
    print("Counter: got samfile header", file=logfile, flush=True)

    try:
        read_dict = count_reads(samfile, contig_manager, args)
    except OSError as e:
        print(str(e), file=logfile, flush=True)
        sys.exit(1)
    finally:
        samfile.close()

    total_reads = sum([contig.total_reads for name, contig in contig_manager.contigs_seen.items()])
    print("Counter: printing read counts", file=logfile)
    for name, contig in sorted(contig_manager.contigs_seen.items(),
                               key=lambda x: x[1],
                               reverse=True):
        print(name, contig.total_reads, sep=" ", file=logfile)

    print("Counter: finished counting %d reads" % total_reads, file=logfile, flush=True)

    for contig, array in read_dict.items():
        read_dict[contig] = bsr_matrix(array)
        del array

    print("Counter: finished converting numpy arrays and lil matrices to bsr_matrix",
          file=logfile, flush=True)

    with shelve.open(args.outfile, protocol=HIGHEST_PROTOCOL) as outfile:
        outfile.update(read_dict)

    run_end = time.time()
    total_runtime = run_end - run_start

    print("Counter: finished pickling matrices: %s" % args.outfile, file=logfile)
    print("Counter: total runtime: %d" % total_runtime, file=logfile, flush=True)
    if logfile is not sys.stderr:
        logfile.close()
    sys.exit()
