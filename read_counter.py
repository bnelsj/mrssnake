from __future__ import print_function
from __future__ import division

try:
    import cPickle as pickle
except ImportError:
    import pickle

import sys
import argparse

import pysam
import numpy as np

from scipy.sparse import lil_matrix
from scipy.sparse import csr_matrix
from functools import total_ordering
import time

@total_ordering
class Contig:
    def __init__(self, name, size):
        self.name = name
        self.size = size
        self.reads = 0
    def __eq__(self, other):
        return self.reads / self.size == other.reads / other.size

    def __lt__(self, other):
        return self.reads / self.size == other.reads / other.size

class ContigManager:
    def __init__(self, max_bases, contigs_seen = {}, array_contigs = []):
        self.max_bases = max_bases
        self.contigs_seen = contigs_seen
        self.array_contigs = array_contigs
        self.used_bases = sum([contigs_seen[contig].size for contig in self.array_contigs])

    def add_contig(self, contig, size=None):
        if size is not None:
            contig = Contig(contig, size)
        if contig.name not in self.contigs_seen.keys():
            self.contigs_seen[contig.name] = contig
        if self.used_bases + contig.size <= self.max_bases and contig.name not in self.array_contigs:
            self.array_contigs.append(contig.name)
            self.used_bases += contig.size

    def add_contig_to_array_contigs(self, contig, size = None):
        if size is not None:
            contig = Contig(contig, size)
        self.add_contig(contig)
        if contig.name not in self.array_contigs:
            self.array_contigs.append(contig.name)
            self.used_bases += contig.size

    def rebalance(self):
        new_used_bases = 0
        new_array_contigs = []
        for name, contig in sorted(self.contigs_seen.items()):
            if new_used_bases + contig.size <= self.max_bases:
                new_array_contigs.append(contig.name)
                new_used_bases += contig.size

        self.array_contigs = new_array_contigs
        self.used_bases = new_used_bases


def get_array_contigs(contigs, args):

    array_contigs = []
    if args.all_contigs:
        array_contigs = contigs.keys()

    else:
        if args.common_contigs is not None:
            array_contigs.extend(args.common_contigs)
        if args.noncanonical_contigs:
            canonical = ["chr%s" % str(x) for x in list(range(1,25)) + ["X", "Y", "M"]]
            array_contigs.extend([contig for contig in contigs.keys() if contig not in canonical])
    return array_contigs

def update_read_depth_and_start(matrix, edist, start, end, nedists=3):
    if isinstance(matrix, np.ndarray):
        matrix[edist + nedists, start:end] += 1
    else:
        slice = matrix[edist + nedists, start:end].toarray()
        slice += 1
        matrix[edist + nedists, start:end] = slice

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
#                    sys.stderr.write("Counter: %s (%d, %d) numpy array\n" % (contig, nrows, length))
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


def count_reads(samfile, contig_manager, args):
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
                sys.stderr.write("Counter: %s (%d, %d) numpy array\n" % (contig, nrows, length))
                sys.stderr.flush()
            else:
                read_dict[contig] = lil_matrix((nrows, length), dtype=np.uint16)
                sys.stderr.write("Counter: %s scipy lil_matrix\n" % contig)
                sys.stderr.flush()


        if isinstance(read_dict[contig], np.ndarray):
            # Update read depth counts for numpy array
            read_dict[contig][edist + nstart_rows, start:end] += 1
        else:
            # Update read depth counts for sparse matrix
            slice = read_dict[contig][edist + nstart_rows, start:end].toarray()
            slice += 1
            read_dict[contig][edist + nstart_rows, start:end] = slice
 
        # Update read start counts
        read_dict[contig][edist, start] += 1
        contig_manager.contigs_seen[contig].reads += 1

        if i % 1000000 == 0 and i > 0:
            sys.stderr.write("Counter: %d reads processed\n" % i)
            if args.max_basepairs_in_mem > 0:
                sys.stderr.write("Counter: rebalancing array contigs...\n")
                start_time = time.time()
                print(contig_manager.array_contigs, sep=" ")
                contig_manager.rebalance()
                print(contig_manager.array_contigs, sep=" ")
                for contig, array in read_dict.items():
                    if contig not in contig_manager.array_contigs:
                        if not isinstance(array, lil_matrix):
                            read_dict[contig] = lil_matrix(array)
                    else:
                        if not isinstance(array, np.ndarray):
                            read_dict[contig] = array.toarray()
                finish_time = time.time()
                total_time = finish_time - start_time
                sys.stderr.write("Counter: rebalancing finished in %d sec...\n" % total_time)
            sys.stderr.flush()

    return read_dict


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", default="/dev/stdin")
    parser.add_argument("outfile", default="wssd_out_file")
    parser.add_argument("contig_lengths", help="tab-delimited file with contig names and lengths")
    parser.add_argument("--max_edist", default = 2, type = int, help = "Maximum edit distance of input reads")
    parser.add_argument("--common_contigs", default = [], nargs="+", help = "Create numpy array for common contigs (Much faster, more memory)")
    parser.add_argument("--all_contigs", action="store_true", help="Create numpy array for all contigs (Fast, high mem requirement)")
    parser.add_argument("--noncanonical_contigs", action="store_true", help="Create numpy array for all noncanonical contigs (Fast, high mem requirement).")
    parser.add_argument("--max_basepairs_in_mem", type=int, default = 0, help="Number of reference bp to keep in memory")

    args = parser.parse_args()

    contigs = {}

    with open(args.contig_lengths, "r") as reader:
        for line in reader:
            contig, length = line.rstrip().split()
            contigs[contig] = int(length)

    array_contigs = get_array_contigs(contigs, args)

    contig_manager = ContigManager(args.max_basepairs_in_mem)

    # Add array contigs first
    for name in array_contigs:
        if name not in contigs:
            continue
        size = contigs[name]
        contig_manager.add_contig_to_array_contigs(name, size)

    # Add remaining contigs to array_contigs if they fit
    for name, size in contigs.items():
        if name in array_contigs:
            continue
        else:
            contig_manager.add_contig(name, size)

    print("Max bases: ", contig_manager.max_bases, "Used bases: ", contig_manager.used_bases)

    samfile = pysam.AlignmentFile(args.infile, "r", check_sq = False)
    sys.stderr.write("Counter: got samfile header\n")
    sys.stdout.flush()

    try:
        read_dict = count_reads(samfile, contig_manager, args)
    except OSError as e:
        sys.stderr.write(str(e))
        sys.stderr.flush()
        sys.exit(1)
    finally:
        samfile.close()

    sys.stderr.write("Counter: finished counting reads\n")
    sys.stderr.flush()

    for contig, array in read_dict.items():
        read_dict[contig] = csr_matrix(array)
        del(array)

    sys.stderr.write("Counter: finished converting numpy arrays and lil matrices to csr_matrix\n")
    sys.stderr.flush()

    with open(args.outfile, "wb") as outfile:
        pickle.dump(read_dict, outfile, pickle.HIGHEST_PROTOCOL)

    sys.stderr.write("Counter: finished pickling matrices\n")
    sys.stderr.flush()
    sys.exit()
