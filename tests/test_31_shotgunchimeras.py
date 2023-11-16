import os
import sys
import unittest

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src'))
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'tests'))
current_dir = os.path.dirname(os.path.abspath(__file__))

from functions_test import *
from Biogrinder import Biogrinder

class Test_31_ShotgunChimeras(unittest.TestCase):
    def test_shotgun_chimeras(self):
        factory = Biogrinder('-rf', current_dir + '/data/shotgun_database_shared_kmers.fa',
                             '-di', '5',
                             '-tr', '300',
                             '-cp', '100',
                             '-cd', '1', '1', '1',
                             '-ck', '10',
                             '-id', '0')
        factory.next_lib()
        refs = {}
        while True:
            read = factory.next_read()
            if not read:
                break
            references = get_references(read)
            self.assertTrue(2 <= len(references) <= 4)
            for ref in references:
                refs[ref] = refs.get(ref, 0) + 1
        
        self.assertIsNone(factory.next_lib())
    
    def test_partial_reference_sequences(self):
        factory = Biogrinder('-rf', current_dir + '/data/shotgun_database_shared_kmers.fa',
                             '-di', '3',
                             '-tr', '300',
                             '-cp', '100',
                             '-cd', '1',
                             '-ck', '2',
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
        self.assertEqual(len(refs.keys()), 3)
        self.assertIsNone(factory.next_lib())

if __name__ == '__main__':
    unittest.main()

