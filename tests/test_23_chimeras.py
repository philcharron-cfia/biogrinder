import os
import re
import sys
import unittest

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src'))
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'tests'))
current_dir = os.path.dirname(os.path.abspath(__file__))

from functions_test import *
from Biogrinder import Biogrinder



class Test_23_Chimeras(unittest.TestCase):
    def test_no_chimeras(self):
        factory = Biogrinder('-rf', current_dir + '/data/amplicon_database.fa',
                             '-fr', current_dir + '/data/forward_reverse_primers.fa',     
                             '-lb', '0',
                             '-un', '1',
                             '-cp', '0',
                             '-cd', '1',
                             '-ck', '0',
                             '-tr', '100',
                             '-id', '0')
        factory.next_lib()
        while True:
            read = factory.next_read()
            if not read:
                break
            self.assertEqual(len(get_references(read)), 1)
            seq = remove_primers(str(read.seq), 'AACT.AAA.GAATTG.CGG', 'G.ACACACCGCCCGT')
            self.assertRegex(seq, r'^(A+|C+|G+|T+)+$')
   
    def test_50_percent_chimeras(self):
        factory = Biogrinder('-rf', current_dir + '/data/amplicon_database.fa',
                             '-fr', current_dir + '/data/forward_reverse_primers.fa',     
                             '-lb', '0',
                             '-un', '1',
                             '-cp', '50',
                             '-cd', '1',
                             '-ck', '0',
                             '-tr', '100',
                             '-id', '0')
        factory.next_lib()
        nof_chimeras = 0
        nof_regulars = 0
        while True:
            read = factory.next_read()
            if not read:
                break
            nof_chimeras += len(get_references(read))
            nof_regulars += len(get_references(read))
        self.assertAlmostEqual(nof_chimeras / nof_regulars, 1, delta=0.1)

    def test_100_percent_bimeras_chimeras(self):
        factory = Biogrinder('-rf', current_dir + '/data/amplicon_database.fa',
                             '-fr', current_dir + '/data/forward_reverse_primers.fa',     
                             '-lb', '0',
                             '-un', '1',
                             '-cp', '100',
                             '-cd', '1',
                             '-ck', '0',
                             '-tr', '100',
                             '-id', '0')
        factory.next_lib()
        while True:
            read = factory.next_read()
            if not read:
                break
            self.assertEqual(len(get_references(read)), 2)
    
    def test_100_percent_trimeras_chimeras(self):
        factory = Biogrinder('-rf', current_dir + '/data/amplicon_database.fa',
                             '-fr', current_dir + '/data/forward_reverse_primers.fa',     
                             '-lb', '0',
                             '-un', '1',
                             '-cp', '100',
                             '-cd', '0', '1',
                             '-ck', '0',
                             '-tr', '100',
                             '-id', '0')
        factory.next_lib()
        while True:
            read = factory.next_read()
            if not read:
                break
            self.assertEqual(len(get_references(read)), 3)
    
    def test_100_percent_quadrameras_chimeras(self):
        factory = Biogrinder('-rf', current_dir + '/data/amplicon_database.fa',
                             '-fr', current_dir + '/data/forward_reverse_primers.fa',     
                             '-lb', '0',
                             '-un', '1',
                             '-cp', '100',
                             '-cd', '0', '0', '1',
                             '-ck', '0',
                             '-tr', '100',
                             '-id', '0')
        factory.next_lib()
        while True:
            read = factory.next_read()
            if not read:
                break
            self.assertEqual(len(get_references(read)), 4)
    
    def test_100_percent_all_chimeras(self):
        factory = Biogrinder('-rf', current_dir + '/data/amplicon_database.fa',
                             '-fr', current_dir + '/data/forward_reverse_primers.fa',     
                             '-lb', '0',
                             '-un', '1',
                             '-cp', '100',
                             '-cd', '1', '1', '1',
                             '-ck', '0',
                             '-tr', '1000',
                             '-id', '0')
        factory.next_lib()
        chim_sizes = {}
        while True:
            read = factory.next_read()
            if not read:
                break
            nof_refs = len(get_references(read))
            chim_sizes[nof_refs] = chim_sizes.get(nof_refs, 0) + 1
            self.assertIn(nof_refs, range(2, 5))
        for size in [2, 3, 4]:
            self.assertAlmostEqual(chim_sizes.get(size, 0), 333.3, delta=333.3*0.1)

def remove_primers(seq, forward_re, reverse_re):
    """
    Remove forward and reverse primers from the sequence.
    """
    seq = re.sub(forward_re, '', seq, flags=re.IGNORECASE)
    seq = re.sub(reverse_re, '', seq, flags=re.IGNORECASE)
    return seq

if __name__ == '__main__':
    unittest.main()

