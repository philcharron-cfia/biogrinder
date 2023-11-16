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


class Test_03_Amplicon(unittest.TestCase):
    def test_forward_reverse_primers_all_sequncing(self):
        factory = Biogrinder('-rf', current_dir + '/data/amplicon_database.fa',
                             '-fr', current_dir + '/data/forward_reverse_primers.fa',
                             '-lb', '0',
                             '-rd', '48',
                             '-tr', '100',
                             '-id', '0')
        self.assertTrue(factory.next_lib())
        nof_reads = 0
        while True:
            read = factory.next_read()
            if not read:
                break
            nof_reads += 1
            ok_read(read, None, nof_reads)
        self.assertEqual(nof_reads, 100)

    def test_forward_reverse_primers_forward_sequncing(self):
        factory = Biogrinder('-rf', current_dir + '/data/amplicon_database.fa',
                             '-fr', current_dir + '/data/forward_reverse_primers.fa',
                             '-lb', '0',
                             '-rd', '48',
                             '-tr', '100',
                             '-un', '1',
                             '-id', '0')
        self.assertTrue(factory.next_lib())
        nof_reads = 0
        while True:
            read = factory.next_read()
            if not read:
                break
            nof_reads += 1
            ok_read(read, 1, nof_reads)
        self.assertEqual(nof_reads, 100)

    def test_forward_reverse_primers_reverse_sequncing(self):
        factory = Biogrinder('-rf', current_dir + '/data/amplicon_database.fa',
                             '-fr', current_dir + '/data/forward_reverse_primers.fa',
                             '-lb', '0',
                             '-rd', '48',
                             '-tr', '100',
                             '-un', '-1',
                             '-id', '0')
        self.assertTrue(factory.next_lib())
        nof_reads = 0
        while True:
            read = factory.next_read()
            if not read:
                break
            nof_reads += 1
            ok_read(read, -1, nof_reads)
        self.assertEqual(nof_reads, 100)

    def test_reverse_forward_primers_all_sequncing(self):
        factory = Biogrinder('-rf', current_dir + '/data/amplicon_database.fa',
                             '-fr', current_dir + '/data/reverse_forward_primers.fa',
                             '-lb', '0',
                             '-rd', '48',
                             '-tr', '100',
                             '-id', '0')
        self.assertTrue(factory.next_lib())
        nof_reads = 0
        while True:
            read = factory.next_read()
            if not read:
                break
            nof_reads += 1
            ok_read(read, None, nof_reads)
        self.assertEqual(nof_reads, 100)

    def test_reverse_forward_primers_forward_sequncing(self):
        factory = Biogrinder('-rf', current_dir + '/data/amplicon_database.fa',
                             '-fr', current_dir + '/data/reverse_forward_primers.fa',
                             '-lb', '0',
                             '-rd', '48',
                             '-tr', '100',
                             '-un', '1',
                             '-id', '0')
        self.assertTrue(factory.next_lib())
        nof_reads = 0
        while True:
            read = factory.next_read()
            if not read:
                break
            nof_reads += 1
            ok_read(read, 1, nof_reads)
        self.assertEqual(nof_reads, 100)

    def test_reverse_forward_primers_reverse_sequncing(self):
        factory = Biogrinder('-rf', current_dir + '/data/amplicon_database.fa',
                             '-fr', current_dir + '/data/reverse_forward_primers.fa',
                             '-lb', '0',
                             '-rd', '48',
                             '-tr', '100',
                             '-un', '-1',
                             '-id', '0')
        self.assertTrue(factory.next_lib())
        nof_reads = 0
        while True:
            read = factory.next_read()
            if not read:
                break
            nof_reads += 1
            ok_read(read, -1, nof_reads)
        self.assertEqual(nof_reads, 100)
    
def ok_read(read, req_strand, nof_reads):
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
  
    if req_strand == -1:
        letters = str(Seq(letters).reverse_complement())
    
    source_key = next((key for key in letters if source.startswith(key)), None)
    if source_key:
        letter_set = letters[source_key]
        if req_strand == -1:
            letter_set = str(Seq(letter_set).reverse_complement())
        assert re.match(f"[{letter_set}]+", str(read.seq)), "Sequence does not match the expected pattern"
    
    assert read.id == str(nof_reads), "Read ID does not match the number of reads"

if __name__ == '__main__':
    unittest.main()