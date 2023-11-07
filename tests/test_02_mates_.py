from Bio.Seq import Seq
import os
import re
import sys
import unittest

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src'))

from Biogrinder import Biogrinder
from SimulatedRead import SimulatedRead


class Test_02_Mates(unittest.TestCase):
    def setUp(self):
        self.factory = Biogrinder('-rf', 'data/shotgun_database_extended.fa',
                                  '-tr', '100',
                                '-rd', '48',
                             '-id', '250')


    def test_paired_reads(self):
        self.assertTrue(self.factory.next_lib(), 'Mate pairs')
        self.nof_reads = 0
        while True:
            read = self.factory.next_read()
            if not read:
                break
            self.nof_reads += 1
            ok_mate(read, None, self.nof_reads)

        self.assertEqual(self.nof_reads, 100)


def ok_mate(read, req_strand, nof_reads):
    assert isinstance(read, SimulatedRead), "read is not an instance of SimulatedRead"
    source = read.reference_id
    strand = read.strand
    if req_strand is None:
        req_strand = strand
    else:
        assert strand == req_strand, "Strand does not match required strand"

    letters = {'seq1': 'a',
               'seq2': 'c',
               'seq3': 'g',
               'seq4': 't',
               'seq5': 'atg',
               }.get(source, '')
    if req_strand == -1:  # Take the reverse complement
        letters = str(Seq(letters).reverse_complement())

    assert re.match(f"[{letters}]+", str(read.seq)), "Sequence does not match the expected pattern"
    id = f"{(nof_reads + 1) // 2}/{'1' if nof_reads % 2 else '2'}"
    assert read.id == id, "Read ID does not match the number of reads"


if __name__ == '__main__':
    unittest.main()