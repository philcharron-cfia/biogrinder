from Bio.Seq import Seq
import os
import re
import sys
import unittest

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src'))
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'tests'))
current_dir = os.path.dirname(os.path.abspath(__file__))

from functions_test import *
from Biogrinder import Biogrinder
from SimulatedRead import SimulatedRead

class Test_17_Libraries(unittest.TestCase):
    def test_multiple_shotgun_libraries(self):
        factory = Biogrinder('-rf', current_dir + '/data/shotgun_database.fa',
                             '-tr', '99',
                             '-rd', '48',
                             '-nl', '4',
                             '-id', '0')
        nof_libs = 0
        lib = factory.next_lib()
        while lib:
            nof_libs += 1
            nof_reads = 0
            read = factory.next_read()
            while read:
                nof_reads += 1
                ok_read(read, None, nof_reads, nof_libs)
                read = factory.next_read()
            self.assertEqual(nof_reads, 99)
            lib = factory.next_lib()
        self.assertEqual(nof_libs, 4)
                
    def test_multiple_paired_end_shotgun_libraries(self):
        factory = Biogrinder('-rf', current_dir + '/data/shotgun_database.fa',
                             '-tr', '100',
                             '-rd', '48',
                             '-nl', '4',
                             '-id', '250')
        nof_libs = 0
        lib = factory.next_lib()
        while lib:
            nof_libs += 1
            nof_reads = 0
            read = factory.next_read()
            while read:
                nof_reads += 1
                ok_mate(read, None, nof_reads, nof_libs)
                read = factory.next_read()
            self.assertEqual(nof_reads, 100)
            lib = factory.next_lib()
        self.assertEqual(nof_libs, 4)
                                              
def ok_read(read, req_strand, nof_reads, nof_libs):
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
    assert read.id == str(nof_libs) + "_" + str(nof_reads), "Read ID does not match the number of reads"

def ok_mate(read, req_strand, nof_reads, nof_libs):
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
    read_id = f"{(nof_reads + 1) // 2}/{'1' if nof_reads % 2 else '2'}"
    id = str(nof_libs) + "_" + read_id
    assert read.id == id, "Read ID does not match the number of reads"
    
if __name__ == '__main__':
    unittest.main()