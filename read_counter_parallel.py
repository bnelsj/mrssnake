# cython: infer_types=True
# cython: boundscheck=False
# cython: wraparound=False
"""
Convert mrsfast sam output to matrix of read depth and read start depth values.
Uses sparse matrices to reduce memory footprint.
"""

from __future__ import print_function
from __future__ import division

import os, sys
import argparse

from functools import total_ordering
import time

import re

import tables
import numpy as np
import pandas as pd
from tables import NoSuchNodeError

import multiprocessing
import signal
import create_depth_array

def process_block(block, contig_hits, regex_full, nhits):
    for (contig, pos, edist) in regex_full.findall(block):
        pos = int(pos) - 1
        contig_hits[contig].append("{}\t{}".format(pos, edist))
        nhits += 1
    return contig_hits, nhits

def get_contig_dict(contigs, of_prefix, proc_num):
    contig_dict = {}
    contig_hits = {}
    for contig in contigs:
        contig_hits[contig] = []
        contig_dict[contig] = open("{}.{}.{}.txt".format(of_prefix, contig, proc_num), "w")
    return contig_dict, contig_hits

def write_output(contig_handles, contig_hits):
    for contig in contig_handles.keys():
        if contig_hits[contig] != []:
            print("\n".join(contig_hits[contig]), file=contig_handles[contig])

def reset_hits(contig_hits):
    for contig in contig_hits.keys():
        contig_hits[contig] = []
    return contig_hits

def close_contig_handles(contig_handles):
    for contig in contig_handles.keys():
        contig_handles[contig].close()

def worker(contigs, of_prefix, proc_num, block_size=1024*1024):
    contig_handles, contig_hits = get_contig_dict(contigs, of_prefix, proc_num)
    contigs_string = "|".join(contigs)
    regex_full = re.compile("[^ @\t]+\t[0-9]+\t(%s)\t([0-9]+)\t.+NM:i:([0-9]+)" % contigs_string)
    nhits = 0
    batch = 0

    dat = None

    while dat != "":
        lock.acquire()
        dat = sam_handle.read(block_size)
        dat += sam_handle.readline()
        lock.release()

        contig_hits, nhits = process_block(dat, contig_hits, regex_full, nhits)
        write_output(contig_handles, contig_hits)
        contig_hits = reset_hits(contig_hits)
        if batch % 100 == 0:
            print("Process {} finished batch {}".format(proc_num, batch))
        batch += 1
    close_contig_handles(contig_handles)
    return proc_num


def setup(sh, l):
    global sam_handle
    global lock
    sam_handle = sh
    lock = l

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("samfile", help="Sam-format input file")
    parser.add_argument("of_prefix", help="prefix for output files")
    parser.add_argument("--contigs_file", required=True,
                        help="tab-delimited file with contig names and lengths")
    parser.add_argument("--log", default=sys.stderr, help="Path to log file. Default: sys.stderr")
    parser.add_argument("--threads", default=1, type=int, help="Number of threads to use (Default: %(default)s)")
    parser.add_argument("--block_size", type=int, default=1024*1024, help="Size of block to read from disk (Default: %(default)s)")

    args = parser.parse_args()

    run_start = time.time()

    if args.log is not sys.stderr:
        logfile = open(args.log, "w")
    else:
        logfile = sys.stderr

    contigs = pd.read_table(args.contigs_file, header=None, names=["contigs", "size"])["contigs"].tolist()

    l = multiprocessing.Lock()
    sh = open(args.samfile, mode="r")

    setup(sh, l)
    threads = []

    # Save reference to original SIGINT signal
    default_handler = signal.getsignal(signal.SIGINT)

    # Set SIGINT to ignore mode
    # To prevent sending it to child processes
    signal.signal(signal.SIGINT, signal.SIG_IGN)

    dirname = os.path.dirname(args.of_prefix)
    if not os.path.exists(dirname):
        os.makedirs(dirname)

    for i in range(args.threads):
        p = multiprocessing.Process(target=worker, args=(contigs, args.of_prefix, i, args.block_size))
        threads.append(p)
        p.start()

    signal.signal(signal.SIGINT, default_handler)

    try:
        for thread in threads:
            thread.join()

    except KeyboardInterrupt:
        print("Caught keyboard interrupt. Shutting down.")
        for thread in threads:
            thread.terminate()
        sys.exit(1)

    except:
        print("Unexpected error:", sys.exec_info()[0])
        for thread in threads:
            thread.terminate()
        sys.exit(1)

    mp_end = time.time() - run_start

    run_end = time.time()
    total_runtime = run_end - run_start

    print("Counter: total runtime: %d seconds" % total_runtime, file=logfile, flush=True)
    if logfile is not sys.stderr:
        logfile.close()
    sys.exit()
