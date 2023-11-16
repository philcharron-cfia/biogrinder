import os
import sys
import unittest

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src'))
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'tests'))
current_dir = os.path.dirname(os.path.abspath(__file__))

from functions_test import *
from Biogrinder import Biogrinder

class Test_15_Multiplex(unittest.TestCase):
    blockPrint()
    def test_single_mid_shotgun(self):
        factory = Biogrinder('-rf', current_dir + '/data/shotgun_database.fa',
                             '-mi', current_dir + '/data/mids.fa',
                             '-nl', '1',                             
                             '-tr', '10',
                             '-rd', '52',
                             '-rs', '1234',
                             '-id', '0')
        factory.next_lib()
        while True:
            read = factory.next_read()
            if not read:
                break
            self.assertEqual(len(read.seq), 52)
            self.assertEqual(read.seq[:4], 'ACGT')

    def test_two_mid_shotgun(self):
        factory = Biogrinder('-rf', current_dir + '/data/shotgun_database.fa',
                             '-mi', current_dir + '/data/mids.fa',
                             '-nl', '2',                             
                             '-tr', '10',
                             '-rd', '52',
                             '-rs', '1234',
                             '-id', '0')
        factory.next_lib()
        while True:
            read = factory.next_read()
            if not read:
                break
            self.assertEqual(len(read.seq), 52)
            self.assertEqual(read.seq[:4], 'ACGT')
        
        factory.next_lib()
        while True:
            read = factory.next_read()
            if not read:
                break
            self.assertEqual(len(read.seq), 52)
            self.assertEqual(read.seq[:8], 'AAAATTTT')

    def test_single_mid_amplicon(self):
        factory = Biogrinder('-rf', current_dir + '/data/single_amplicon_database.fa',
                             '-mi', current_dir + '/data/mids.fa',
                             '-fr', current_dir + '/data/forward_reverse_primers.fa',
                             '-nl', '1',                             
                             '-tr', '10',
                             '-rd', '70',
                             '-un', '1',
                             '-st', '1',
                             '-id', '0')
        factory.next_lib()
        while True:
            read = factory.next_read()
            if not read:
                break
            self.assertEqual(str(read.seq), 'ACGTAAACTTAAAGGAATTGACGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGTACACACCGC')
                                                 
    def test_single_mid_amplicon_too_long(self):
        factory = Biogrinder('-rf', current_dir + '/data/single_amplicon_database.fa',
                             '-mi', current_dir + '/data/mids.fa',
                             '-fr', current_dir + '/data/forward_reverse_primers.fa',
                             '-nl', '1',                             
                             '-tr', '10',
                             '-rd', '80',
                             '-un', '1',
                             '-st', '1',
                             '-id', '0')
        factory.next_lib()
        while True:
            read = factory.next_read()
            if not read:
                break
            self.assertEqual(str(read.seq), 'ACGTAAACTTAAAGGAATTGACGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGTACACACCGCCCGT')

    def test_two_mid_amplicon(self):
        factory = Biogrinder('-rf', current_dir + '/data/single_amplicon_database.fa',
                             '-mi', current_dir + '/data/mids.fa',
                             '-fr', current_dir + '/data/forward_reverse_primers.fa',
                             '-nl', '2',
                             '-sp', '100',                             
                             '-tr', '10',
                             '-rd', '74',
                             '-un', '1',
                             '-st', '1',
                             '-id', '0')
        factory.next_lib()
        while True:
            read = factory.next_read()
            if not read:
                break
            self.assertEqual(str(read.seq), 'ACGTAAACTTAAAGGAATTGACGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGTACACACCGCCCGT') 
if __name__ == '__main__':
    unittest.main()