from Bio.Seq import Seq
import os
import sys
import unittest
from functions_test import *

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src'))

from Biogrinder import Biogrinder
from SimulatedRead import SimulatedRead

class Test_28_ReverseComplementAmplicon(unittest.TestCase):
    def test_forward_reverse_primers_forward_sequencing(self):
        factory = Biogrinder('-rf', 'data/revcom_amplicon_database.fa',
                             '-fr', 'data/forward_reverse_primers.fa',
                             '-lb', '0',
                             '-un', '1',
                             '-rd', '48',
                             '-tr', '100',
                             '-id', '0')
        factory.next_lib()
        nof_reads = 0
        while True:
            read = factory.next_read()
            if not read:
                break
            nof_reads += 1
            self.ok_read(read, None, nof_reads)
        self.assertEqual(nof_reads, 100)

    def test_reverse_forward_primers_reverse_sequencing(self):
        factory = Biogrinder('-rf', 'data/revcom_amplicon_database.fa',
                             '-fr', 'data/forward_reverse_primers.fa',
                             '-lb', '0',
                             '-un', '-1',
                             '-rd', '48',
                             '-tr', '100',
                             '-id', '0')
        factory.next_lib()
        nof_reads = 0
        while True:
            read = factory.next_read()
            if not read:
                break
            nof_reads += 1
            self.ok_read(read, None, nof_reads)
        self.assertEqual(nof_reads, 100)

    def ok_read(self, read, req_strand, nof_reads):
        self.assertIsInstance(read, SimulatedRead)
        source = read.reference_id
        strand = read.strand
        req_strand = strand if req_strand is None else req_strand
        self.assertEqual(strand, req_strand)

        letters = 'T' if re.match(r'^seq1', source) else ''
        if req_strand == -1:
            letters = str(Seq(letters).reverse_complement())

        self.assertRegex(str(read.seq), f"[{letters}]+")
        self.assertEqual(read.id, str(nof_reads))
    
if __name__ == '__main__':
    unittest.main()

