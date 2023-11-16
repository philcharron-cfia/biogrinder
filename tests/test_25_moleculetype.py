import os
import sys
import unittest
from functions_test import *

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src'))

from Biogrinder import Biogrinder



class Test_25_MoleculeType(unittest.TestCase):
    def test_dna_database(self):
        dna_want_chars = set('ACGT')
        factory = Biogrinder('-rf', 'data/database_dna.fa',
                             '-rd', '240',
                             '-tr', '100',
                             '-mr', '100', '0',
                             '-id', '0',
                             '-md', 'uniform', '20')
        factory.next_lib()
        self.assertEqual(factory.alphabet, 'dna')

        while True:
            read = factory.next_read()
            if not read:
                break
            got_chars = set(str(read.seq))
            self.assertEqual(got_chars, dna_want_chars)

    def test_rna_database(self):
        rna_want_chars  = set('ACGU')
        factory = Biogrinder('-rf', 'data/database_rna.fa',
                             '-rd', '240',
                             '-tr', '100',
                             '-mr', '100', '0',
                             '-id', '0',
                             '-md', 'uniform', '20')
        factory.next_lib()
        self.assertEqual(factory.alphabet, 'rna')

        while True:
            read = factory.next_read()
            if not read:
                break
            got_chars = set(str(read.seq))
            self.assertEqual(got_chars, rna_want_chars)

    def test_protein_database(self):
        protein_want_chars = set('ARNDCEQGHILKMFPSTWYV')

        factory = Biogrinder('-rf', 'data/database_protein.fa',
                             '-rd', '240',
                             '-tr', '100',
                             '-mr', '100', '0',
                             '-id', '0',
                             '-un', '1',
                             '-md', 'uniform', '20')
        factory.next_lib()
        self.assertEqual(factory.alphabet, 'protein')
        while True:
            read = factory.next_read()
            if not read:
                break
            got_chars = set(str(read.seq))
            self.assertEqual(got_chars, protein_want_chars)
    
    def test_mixed_database(self):

        with self.assertRaises(Exception) as context:
            factory = Biogrinder('-rf', 'data/database_mixed.fa',
                                '-tr', '100',
                                '-id', '0',
                                '-un', '1')
            factory.next_lib()
            
        self.assertEqual(str(context.exception), "Error: Cannot determine what type of molecules the reference sequences are. Got 1 sequences of type 'dna' and 2 others.")
   
if __name__ == '__main__':
    unittest.main()

