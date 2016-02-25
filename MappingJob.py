import argparse
import bisect
from functools import total_ordering

@total_ordering
class MappingJob:
    def __init__(self, contig, size):
        self.contigs = [contig]
        self.size = size

    def add_contig(self, contig, size):
        self.contigs.append(contig)
        self.size += size

    def __eq__(self, other):
        return self.size == other.size

    def __lt__(self, other):
        return self.size < other.size
        
class MappingJoblist:
    """Class for assigning small contigs to separate jobs"""
    def __init__(self, max_size):
        self.jobsets = []
        self.max_size = max_size

    def insert_job(self, job):
        """Create a new jobset for job while maintaining sorted order"""
        bisect.insort_left(self.jobsets, job)

    def add_job(self, contig, size):
        """Add job to a jobset with space for job or insert into joblist"""
        index = bisect.bisect_right(self.jobsets, MappingJob("", self.max_size - size))
        if index < len(self.jobsets):
            self.jobsets[index].add_contig(contig, size)
        else:
            self.insert_job(job)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument()

    args = parser.parse_args()

