from Bio.Seq import Seq
import os
import sys
import unittest
from io import StringIO

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src'))
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'tests'))
current_dir = os.path.dirname(os.path.abspath(__file__))

from functions_test import *
from Biogrinder import Biogrinder
from SimulatedRead import SimulatedRead

class Test_27_Stdin(unittest.TestCase):
    def setUp(self):
        # Simulating the __DATA__ section in Perl
        data = """>seq1 this is the first sequence
aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
>seq2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
>seq3
gggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggg
gggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggg
>seq4
tttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt
>seq5 last sequence, last comment
aaaaaaaaaattttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttgggggggggg"""
        
        self.stdin_backup = sys.stdin
        sys.stdin = StringIO(data)

    def tearDown(self):
        sys.stdin = self.stdin_backup

    def test_input_from_stdin(self):

        factory = Biogrinder('-rf', '-',
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

        letters = {'seq1': 'A', 'seq2': 'C', 'seq3': 'G', 'seq4': 'T', 'seq5': 'ATG'}[source]
        if req_strand == -1:
            letters = str(Seq(letters).reverse_complement())

        self.assertRegex(str(read.seq), f"[{letters}]+")
        self.assertEqual(read.id, str(nof_reads))
    
if __name__ == '__main__':
    unittest.main()

