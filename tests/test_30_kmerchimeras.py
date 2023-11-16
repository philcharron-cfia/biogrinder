import os
import sys
import unittest
from functions_test import *

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src'))

from Biogrinder import Biogrinder

class Test_30_KmerChimeras(unittest.TestCase):
    delta = 0.2
    def test_kmer_chimeras_bimeras(self):
        factory = Biogrinder('-rf', 'data/kmers.fa',
                             '-lb', '0',
                             '-un', '1',
                             '-tr', '100',
                             '-cp', '100',
                             '-cd', '1',
                             '-ck', '8',
                             '-id', '0')
        factory.next_lib()
        refs = {}
        while True:
            read = factory.next_read()
            if not read:
                break
            references = get_references(read)
            self.assertEqual(len(references), 2)
            for ref in references:
                refs[ref] = refs.get(ref, 0) + 1
        self.assertIn('seq1', refs)
        self.assertIn('seq2', refs)
        self.assertIn('seq3', refs)
        self.assertIn('seq4', refs)
        self.assertNotIn('seq5', refs)
    
    def test_kmer_chimeras_trimeras(self):
        factory = Biogrinder('-rf', 'data/kmers.fa',
                             '-lb', '0',
                             '-un', '1',
                             '-tr', '300',
                             '-cp', '100',
                             '-cd', '0', '1',
                             '-ck', '8',
                             '-id', '0')
        factory.next_lib()
        refs = {}
        while True:
            read = factory.next_read()
            if not read:
                break
            references = get_references(read)
            self.assertEqual(len(references), 3)
            for ref in references:
                refs[ref] = refs.get(ref, 0) + 1
        self.assertIn('seq1', refs)
        self.assertIn('seq2', refs)
        self.assertIn('seq3', refs)
        self.assertIn('seq4', refs)
        self.assertNotIn('seq5', refs)
        
    def test_kmer_chimeras_quadrameras(self):
        factory = Biogrinder('-rf', 'data/kmers.fa',
                             '-lb', '0',
                             '-un', '1',
                             '-tr', '300',
                             '-cp', '100',
                             '-cd', '0', '0', '1',
                             '-ck', '8',
                             '-id', '0')
        factory.next_lib()
        refs = {}
        while True:
            read = factory.next_read()
            if not read:
                break
            references = get_references(read)
            self.assertEqual(len(references), 4)
            for ref in references:
                refs[ref] = refs.get(ref, 0) + 1
        self.assertIn('seq1', refs)
        self.assertIn('seq2', refs)
        self.assertIn('seq3', refs)
        self.assertIn('seq4', refs)
        self.assertNotIn('seq5', refs)
    
    def test_kmer_chimeras_all(self):
        factory = Biogrinder('-rf', 'data/kmers.fa',
                             '-lb', '0',
                             '-un', '1',
                             '-tr', '1000',
                             '-cp', '100',
                             '-cd', '1', '1', '1',
                             '-ck', '8',
                             '-id', '0')
        factory.next_lib()
        refs = {}
        chim_sizes = {}
        while True:
            read = factory.next_read()
            if not read:
                break
            references = get_references(read)
            nof_refs = len(references)
            chim_sizes[nof_refs] = chim_sizes.get(nof_refs, 0) + 1
            self.assertTrue(2 <= nof_refs <= 4)
            for ref in references:
                refs[ref] = refs.get(ref, 0) + 1
        
        # Check the distribution of chimeras
        self.assertTrue(333.3 * (1 - self.delta) <= chim_sizes[2] <= 333.3 * (1 + self.delta))
        self.assertTrue(333.3 * (1 - self.delta) <= chim_sizes[3] <= 333.3 * (1 + self.delta))
        self.assertTrue(333.3 * (1 - self.delta) <= chim_sizes[4] <= 333.3 * (1 + self.delta))
        self.assertIn('seq1', refs)
        self.assertIn('seq2', refs)
        self.assertIn('seq3', refs)
        self.assertIn('seq4', refs)
        self.assertNotIn('seq5', refs)
    
    def test_kmer_chimeras_same_abundance(self):
        factory = Biogrinder('-rf', 'data/kmers2.fa',
                             '-lb', '0',
                             '-un', '1',
                             '-tr', '1000',
                             '-cp', '100',
                             '-cd', '0', '1',
                             '-ck', '8',
                             '-id', '0')
        factory.next_lib()
        refs = {}
        while True:
            read = factory.next_read()
            if not read:
                break
            references = get_references(read)
            self.assertEqual(len(references), 3)
            for ref in references:
                refs[ref] = refs.get(ref, 0) + 1
        expected = {'seq1': 800, 'seq2': 1000, 'seq3': 1200}
        self.assertTrue(expected['seq1'] * (1 - self.delta) <= refs['seq1'] <= expected['seq1'] * (1 + self.delta))
        self.assertTrue(expected['seq2'] * (1 - self.delta) <= refs['seq2'] <= expected['seq2'] * (1 + self.delta))
        self.assertTrue(expected['seq3'] * (1 - self.delta) <= refs['seq3'] <= expected['seq3'] * (1 + self.delta))

    def test_kmer_chimeras_different_abundance(self):
        factory = Biogrinder('-rf', 'data/kmers2.fa',
                             '-af', 'data/abundance_kmers.txt',
                             '-lb', '0',
                             '-un', '1',
                             '-tr', '1000',
                             '-cp', '100',
                             '-cd', '0', '1',
                             '-ck', '8',
                             '-id', '0')
        factory.next_lib()
        refs = {}
        while True:
            read = factory.next_read()
            if not read:
                break
            references = get_references(read)
            self.assertEqual(len(references), 3)
            for ref in references:
                refs[ref] = refs.get(ref, 0) + 1
        expected = {'seq1': 1400, 'seq2': 450, 'seq3': 1150}

        self.assertTrue(expected['seq1'] * (1 - self.delta) <= refs['seq1'] <= expected['seq1'] * (1 + self.delta))
        self.assertTrue(expected['seq2'] * (1 - self.delta) <= refs['seq2'] <= expected['seq2'] * (1 + self.delta))
        self.assertTrue(expected['seq3'] * (1 - self.delta) <= refs['seq3'] <= expected['seq3'] * (1 + self.delta))


if __name__ == '__main__':
    unittest.main()

