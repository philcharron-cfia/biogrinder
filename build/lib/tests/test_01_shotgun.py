from Bio.Seq import Seq
import os
import re
import sys
import unittest

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src'))

from Biogrinder import Biogrinder
from SimulatedRead import SimulatedRead


class Test_01_Shotgun(unittest.TestCase):
    def test_short_arguments(self):
        factory = Biogrinder('-rf', 'data/shotgun_database_extended.fa',
                             '-tr', '10')
        self.assertIsNotNone(factory)
        self.assertTrue(factory.next_lib())

    def test_long_arguments(self):
        factory = Biogrinder('--reference_file', 'data/shotgun_database_extended.fa',
                             '--read_dist', '48',
                             '--total_reads', '100')
        self.assertIsNotNone(factory)
        self.assertTrue(factory.next_lib())
        
    def test_read_count(self):
        factory = Biogrinder('--reference_file', 'data/shotgun_database_extended.fa',
                                  '--read_dist', '48',
                                  '--total_reads', '100',
                                  '--insert_dis', '0')
    
        nof_reads = 0
        factory.next_lib()
        read = factory.next_read()
        while read:
            nof_reads += 1
            ok_read(read, None, nof_reads)
            read = factory.next_read()
        self.assertEqual(nof_reads, 100)

def ok_read(read, req_strand, nof_reads):
    assert isinstance(read, SimulatedRead), "read is not an instance of SimulatedRead"
    source = read.reference_id
    strand = read.strand
    if req_strand is None:
        req_strand = strand
    else:
        assert strand == req_strand, "Strand does not match required strand"

    letters = {'seq1': 'A',
               'seq2': 'C',
               'seq3': 'G',
               'seq4': 'T',
               'seq5': 'ATG',
               }.get(source, '')
    if req_strand == -1:  # Take the reverse complement
        letters = str(Seq(letters).reverse_complement())

    assert re.match(f"[{letters}]+", str(read.seq)), "Sequence does not match the expected pattern"
    assert read.id == str(nof_reads), "Read ID does not match the number of reads"

if __name__ == '__main__':
    unittest.main()