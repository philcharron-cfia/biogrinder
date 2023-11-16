import os
import sys
import unittest

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src'))

from Biogrinder import Biogrinder

class Test_07_Diversity(unittest.TestCase):
    
    def setUp(self):
        self.sources = {}

    def test_single_library_single_diversity(self):
        factory = Biogrinder('-rf', 'data/shotgun_database.fa',
                             '-rs', '1234',
                             '-tr', '100',
                             '-di', '2',
                             '-id', '0')

        factory.next_lib()   
        while True:
            read = factory.next_read()
            if not read:
                break
            source = read.reference_id
            self.sources[source] = None

        self.assertEqual(len(self.sources), 2)
        self.sources.clear()

    def test_two_libraries_single_diversity(self):
        factory = Biogrinder('-rf', 'data/shotgun_database.fa',
                             '-rs', '1234',
                             '-tr', '100',
                             '-di', '2',
                             '-nl', '2',
                             '-id', '0')
        factory.next_lib()
        while True:
            read = factory.next_read()
            if not read:
                break
            source = read.reference_id
            self.sources[source] = None

        self.assertEqual(len(self.sources), 2)
        self.sources.clear()

        factory.next_lib()
        while True:
            read = factory.next_read()
            if not read:
                break
            source = read.reference_id
            self.sources[source] = None

        self.assertEqual(len(self.sources), 2)
        self.sources.clear()


    def test_two_libraries_two_diversities(self):
        factory = Biogrinder('-rf', 'data/shotgun_database.fa',
                             '-rs', '1234',
                             '-tr', '100',
                             '-di', '2', '3',
                             '-nl', '2',
                             '-id', '0')
        factory.next_lib()
        while True:
            read = factory.next_read()
            if not read:
                break
            source = read.reference_id
            self.sources[source] = None

        self.assertEqual(len(self.sources), 2)
        self.sources.clear()

        factory.next_lib()
        while True:
            read = factory.next_read()
            if not read:
                break
            source = read.reference_id
            self.sources[source] = None


        self.assertEqual(len(self.sources), 3)
        self.sources.clear()

if __name__ == '__main__':
    unittest.main()
