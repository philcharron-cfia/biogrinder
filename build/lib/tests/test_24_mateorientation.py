import os
import sys
import unittest
from functions_test import *

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src'))

from Biogrinder import Biogrinder



class Test_24_MateOrientation(unittest.TestCase):
    def test_FR_oriented_mates_uni_for(self):
        factory = Biogrinder('-rf', 'data/oriented_database.fa',
                             '-rd', '80',
                             '-un', '1',
                             '-tr', '100',
                             '-rs', '12345',
                             '-id', '240',
                             '-mo', 'FR')
        factory.next_lib()
        nof_reads = 0
        while True:
            read1 = factory.next_read()
            if not read1:
                break
            read2 = factory.next_read()
            nof_reads += 2
            seqtype = type(read1, read2)
            self.assertEqual(seqtype, 'FR')
        self.assertEqual(nof_reads, factory.total_reads)
    
    def test_FR_oriented_mates_uni_rev(self):
        factory = Biogrinder('-rf', 'data/oriented_database.fa',
                             '-rd', '80',
                             '-un', '-1',
                             '-tr', '100',
                             '-rs', '12345',
                             '-id', '240',
                             '-mo', 'FR')
        factory.next_lib()
        nof_reads = 0
        while True:
            read1 = factory.next_read()
            if not read1:
                break
            read2 = factory.next_read()
            nof_reads += 2
            seqtype = type(read1, read2)
            self.assertEqual(seqtype, 'FR')
        self.assertEqual(nof_reads, factory.total_reads)
    
    def test_FR_oriented_mates(self):
        factory = Biogrinder('-rf', 'data/oriented_database.fa',
                             '-rd', '80',
                             '-un', '0',
                             '-tr', '100',
                             '-rs', '12345',
                             '-id', '240',
                             '-mo', 'FR')
        factory.next_lib()
        nof_reads = 0
        while True:
            read1 = factory.next_read()
            if not read1:
                break
            read2 = factory.next_read()
            nof_reads += 2
            seqtype = type(read1, read2)
            self.assertEqual(seqtype, 'FR')
        self.assertEqual(nof_reads, factory.total_reads)
    
    def test_FF_oriented_mates_uni_for(self):
        factory = Biogrinder('-rf', 'data/oriented_database.fa',
                             '-rd', '80',
                             '-un', '1',
                             '-tr', '100',
                             '-rs', '12345',
                             '-id', '240',
                             '-mo', 'FF')
        factory.next_lib()
        nof_reads = 0
        while True:
            read1 = factory.next_read()
            if not read1:
                break
            read2 = factory.next_read()
            nof_reads += 2
            seqtype = type(read1, read2)
            self.assertEqual(seqtype, 'FF')
        self.assertEqual(nof_reads, factory.total_reads)
    
    def test_FF_oriented_mates_uni_rev(self):
        factory = Biogrinder('-rf', 'data/oriented_database.fa',
                             '-rd', '80',
                             '-un', '-1',
                             '-tr', '100',
                             '-rs', '12345',
                             '-id', '240',
                             '-mo', 'FF')
        factory.next_lib()
        nof_reads = 0
        while True:
            read1 = factory.next_read()
            if not read1:
                break
            read2 = factory.next_read()
            nof_reads += 2
            seqtype = type(read1, read2)
            self.assertEqual(seqtype, 'RR')
        self.assertEqual(nof_reads, factory.total_reads)
    
    def test_FF_oriented_mates(self):
        factory = Biogrinder('-rf', 'data/oriented_database.fa',
                             '-rd', '80',
                             '-un', '0',
                             '-tr', '100',
                             '-rs', '12345',
                             '-id', '240',
                             '-mo', 'FF')
        factory.next_lib()
        nof_reads = 0
        while True:
            read1 = factory.next_read()
            if not read1:
                break
            read2 = factory.next_read()
            nof_reads += 2
            seqtype = type(read1, read2)
            self.assertRegex(seqtype, r'(FF|RR)')

        self.assertEqual(nof_reads, factory.total_reads)
    
    def test_RF_oriented_mates_uni_for(self):
        factory = Biogrinder('-rf', 'data/oriented_database.fa',
                             '-rd', '80',
                             '-un', '1',
                             '-tr', '100',
                             '-rs', '12345',
                             '-id', '240',
                             '-mo', 'RF')
        factory.next_lib()
        nof_reads = 0
        while True:
            read1 = factory.next_read()
            if not read1:
                break
            read2 = factory.next_read()
            nof_reads += 2
            seqtype = type(read1, read2)
            self.assertEqual(seqtype, 'RF')
        self.assertEqual(nof_reads, factory.total_reads)
    
    def test_RF_oriented_mates_uni_rev(self):
        factory = Biogrinder('-rf', 'data/oriented_database.fa',
                             '-rd', '80',
                             '-un', '-1',
                             '-tr', '100',
                             '-rs', '12345',
                             '-id', '240',
                             '-mo', 'RF')
        factory.next_lib()
        nof_reads = 0
        while True:
            read1 = factory.next_read()
            if not read1:
                break
            read2 = factory.next_read()
            nof_reads += 2
            seqtype = type(read1, read2)
            self.assertEqual(seqtype, 'RF')
        self.assertEqual(nof_reads, factory.total_reads)
    
    def test_RF_oriented_mates(self):
        factory = Biogrinder('-rf', 'data/oriented_database.fa',
                             '-rd', '80',
                             '-un', '0',
                             '-tr', '100',
                             '-rs', '12345',
                             '-id', '240',
                             '-mo', 'RF')
        factory.next_lib()
        nof_reads = 0
        while True:
            read1 = factory.next_read()
            if not read1:
                break
            read2 = factory.next_read()
            nof_reads += 2
            seqtype = type(read1, read2)
            self.assertEqual(seqtype, 'RF')
        self.assertEqual(nof_reads, factory.total_reads)

    def test_RR_oriented_mates_uni_for(self):
        factory = Biogrinder('-rf', 'data/oriented_database.fa',
                             '-rd', '80',
                             '-un', '1',
                             '-tr', '100',
                             '-rs', '12345',
                             '-id', '240',
                             '-mo', 'RR')
        factory.next_lib()
        nof_reads = 0
        while True:
            read1 = factory.next_read()
            if not read1:
                break
            read2 = factory.next_read()
            nof_reads += 2
            seqtype = type(read1, read2)
            self.assertEqual(seqtype, 'RR')
        self.assertEqual(nof_reads, factory.total_reads)
    
    def test_RR_oriented_mates_uni_rev(self):
        factory = Biogrinder('-rf', 'data/oriented_database.fa',
                             '-rd', '80',
                             '-un', '-1',
                             '-tr', '100',
                             '-rs', '12345',
                             '-id', '240',
                             '-mo', 'RR')
        factory.next_lib()
        nof_reads = 0
        while True:
            read1 = factory.next_read()
            if not read1:
                break
            read2 = factory.next_read()
            nof_reads += 2
            seqtype = type(read1, read2)
            self.assertEqual(seqtype, 'FF')
        self.assertEqual(nof_reads, factory.total_reads)
    
    def test_RR_oriented_mates(self):
        factory = Biogrinder('-rf', 'data/oriented_database.fa',
                             '-rd', '80',
                             '-un', '0',
                             '-tr', '100',
                             '-rs', '12345',
                             '-id', '240',
                             '-mo', 'RR')
        factory.next_lib()
        nof_reads = 0
        while True:
            read1 = factory.next_read()
            if not read1:
                break
            read2 = factory.next_read()
            nof_reads += 2
            seqtype = type(read1, read2)
            self.assertRegex(seqtype, r'(FF|RR)')

        self.assertEqual(nof_reads, factory.total_reads)

    
def type(read1, read2):
        read1_t = {
            'CCCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA': 'F',
            'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTGGG': 'R',
        }
        read2_t = {
            'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTT': 'F',
            'AAATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT': 'R',
        }
        type_ = ''
        if str(read1.seq) in read1_t and str(read2.seq) in read2_t:
            type_ = read1_t[str(read1.seq)] + read2_t[str(read2.seq)]
        elif str(read2.seq) in read1_t and str(read1.seq) in read2_t:
            type_ = read1_t[str(read2.seq)] + read2_t[str(read1.seq)]
        return type_
   
if __name__ == '__main__':
    unittest.main()

