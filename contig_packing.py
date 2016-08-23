from __future__ import print_function
from __future__ import division

import argparse
import attr
from attr.validators import instance_of
import bisect
from pandas import read_table

@attr.s
class Contig(object):
    """Class with contig name and length attributes.
    Comparisons between instances of this class are based on length.
    """
    name = attr.ib(validator=instance_of(str), cmp=False)
    length = attr.ib(validator=instance_of(int))

@attr.s
class ContigBatch(object):
    """Class that contains a batch of Contig objects.
    Comparisons between instances of this class are based on total_length.
    """
    contigs = attr.ib(default=attr.Factory(list), init=False, cmp=False)
    total_length = attr.ib(default=0, init=False)

    def add_contig(self, contig):
        """Add Contig contig to set of contigs and update total size"""
        self.contigs.append(contig)
        self.total_length += contig.length

    def merge(self, other):
        """Merge two batches.
        """
        self.contigs.extend(other.contigs)
        self.total_length += other.total_length

    def print(self):
        contigs = ",".join([contig.name for contig in self.contigs])
        return "{}\t{}".format(contigs, self.total_length)

class ContigGrouping(object):
    """Class that contains batches of contigs.
    Contains methods for combining batches based on total_length.
    """
    
    def __init__(self, contig_list):
        self.batches = []
        for contig in contig_list:
            batch = ContigBatch()
            batch.add_contig(contig)
            self.add_batch(batch)

    def add_batch(self, batch):
        bisect.insort_left(self.batches, batch)

    def merge_batches(self):
        """Create result batch by adding smallest batch to second-smallest,
        insert result into batches and remove smallest batches.
        Maintains sorted order.
        """
        first = self.batches.pop(0)
        second = self.batches.pop(0)
        first.merge(second)

        bisect.insort_left(self.batches, first)

    def write(self, outfile):
        with open(outfile, "w") as of:
            for i, batch in enumerate(self.batches):
                print(i, batch.print(), file=of, sep="\t")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("contig_list", help="Tab-delimited list with contig and length columns")
    parser.add_argument("outfile", help="Tab-delimited output file with batch number \
                                         and comma-separated list of contigs")
    parser.add_argument("--nbatches", type=int, default=24, help="Number of batches \
                                                                  to group contigs into")

    args = parser.parse_args()

    contig_df = read_table(args.contig_list, header=None, names=["contig", "size"])

    contigs = []
    for i, entry in contig_df.iterrows():
        contigs.append(Contig(entry.contig, entry["size"]))

    grouping = ContigGrouping(contigs)

    print(grouping.batches)

    while len(grouping.batches) > args.nbatches:
        grouping.merge_batches()

    grouping.write(args.outfile)
