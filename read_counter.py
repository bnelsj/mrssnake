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

def update_read_depth_and_start(matrix, edist, nedists=3, start, end):
    if isinstance(matrix, np.ndarray):
        matrix[edist + nedists, start:end] += 1
    else:
        slice = matrix[edist + nedists, start:end].toarray()
        slice += 1
        matrix[edist + nedists, start:end] = slice

    matrix[edist, start] += 1

    return matrix

def count_reads_sans_pysam(input, contigs, array_contigs, args):
    read_dict = {}
    nrows = 2 * args.max_edist + 2
    nstart_rows = nrows // 2

    with open(input, "r") as infile:
        for line in infile:
            if line.startswith("@"):
                continue
            contig, start, cigar, edist = line.split()[2,3,5,11]

            if edist > args.max_edist:
                continue

            start = start - 1
            rlen = int(cigar[:-1])
            end = start + rlen
            edist = int(edist[4:])

            if contig not in read_dict:
                length = contigs[contig]
                if contig in array_contigs:
                    read_dict[contig] = np.zeros((nrows, length), dtype=np.uint16)
                    sys.stderr.write("Counter: %s (%d, %d) numpy array\n" % (contig, nrows, length))
                    sys.stderr.flush()
                else:
                    read_dict[contig] = lil_matrix((nrows, length), dtype=np.uint16)
                    sys.stderr.write("Counter: %s scipy lil_matrix\n" % contig)
                    sys.stderr.flush()

            if end > contigs[contig]:
                end = contigs[contig]

            if contig in array_contigs:
                # Update read depth counts for numpy array
                read_dict[contig][edist + nstart_rows, start:end] += 1
            else:
                # Update read depth counts for sparse matrix
                slice = read_dict[contig][edist + nstart_rows, start:end].toarray()
                slice += 1
                read_dict[contig][edist + nstart_rows, start:end] = slice
     
            # Update read start counts
            read_dict[contig][edist, start] += 1

            if i % 1000000 == 0:
                sys.stderr.write("Counter: %d reads processed\n" % i)
                sys.stderr.flush()

    return read_dict


def count_reads(samfile, contigs, array_contigs, args):
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
            length = contigs[contig]
            if contig in array_contigs:
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

        if i % 1000000 == 0:
            sys.stderr.write("Counter: %d reads processed\n" % i)
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

    args = parser.parse_args()

    contigs = {}

    with open(args.contig_lengths, "r") as reader:
        for line in reader:
            contig, length = line.rstrip().split()
            contigs[contig] = int(length)

    array_contigs = get_array_contigs(contigs, args)

    samfile = pysam.AlignmentFile(args.infile, "r", check_sq = False)
    sys.stderr.write("Counter: got samfile header\n")
    sys.stdout.flush()

    try:
        read_dict = count_reads(samfile, contigs, array_contigs, args)
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
