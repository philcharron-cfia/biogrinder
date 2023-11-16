import os
import sys
import unittest

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src'))

from Biogrinder import Biogrinder

class Test_08_Shared(unittest.TestCase):
    
    def setUp(self):
        self.sources = {}
        self.shared = {}

    def check_shared_species(self, lib_num, read, shared_threshold):
        source = read.reference_id
        self.sources.setdefault(lib_num, {})[source] = None
        # Check if this source genome is shared
        if lib_num == shared_threshold and all(source in self.sources[i] for i in range(1, shared_threshold)):
            self.shared[source] = None

    def test_no_species_shared(self):
        factory = Biogrinder('-rf', 'data/shotgun_database.fa',
                        '-rs', '1234',
                        '-am', 'uniform',
                        '-lb', '0',
                        '-tr', '100',
                        '-nl', '3',
                        '-sp', '0',
                        '-id', '0')
        for lib_num in range(1, 4):
            factory.next_lib()
            while True:
                read = factory.next_read()
                if not read:
                    break
                self.check_shared_species(lib_num, read, 3)

        self.assertEqual(len(self.sources), 3)
        for lib in self.sources.values():
            self.assertEqual(len(lib), 1)
        self.assertEqual(len(self.shared), 0)
        self.sources.clear()
        self.shared.clear()

    def test_50_species_shared(self):
        factory = Biogrinder('-rf', 'data/shotgun_database.fa',
                        '-rs', '1234',
                        '-am', 'uniform',
                        '-lb', '0',
                        '-tr', '100',
                        '-nl', '3',
                        '-sp', '50',
                        '-id', '0')
        for lib_num in range(1, 4):
            factory.next_lib()
            while True:
                read = factory.next_read()
                if not read:
                    break
                self.check_shared_species(lib_num, read, 3)

        self.assertEqual(len(self.sources), 3)
        for lib in self.sources.values():
            self.assertEqual(len(lib), 2)
        self.assertEqual(len(self.shared), 1)
        self.sources.clear()
        self.shared.clear()

    def test_67_species_shared(self):
        factory = Biogrinder('-rf', 'data/shotgun_database.fa',
                        '-rs', '1234',
                        '-am', 'uniform',
                        '-lb', '0',
                        '-tr', '100',
                        '-nl', '3',
                        '-sp', '67',
                        '-id', '0')
        for lib_num in range(1, 4):
            factory.next_lib()
            while True:
                read = factory.next_read()
                if not read:
                    break
                self.check_shared_species(lib_num, read, 3)

        self.assertEqual(len(self.sources), 3)
        for lib in self.sources.values():
            self.assertEqual(len(lib), 3)
        self.assertEqual(len(self.shared), 2)
        self.sources.clear()
        self.shared.clear()  
    
    def test_all_species_shared(self):
        factory = Biogrinder('-rf', 'data/shotgun_database.fa',
                        '-rs', '1234',
                        '-am', 'uniform',
                        '-lb', '0',
                        '-tr', '100',
                        '-nl', '3',
                        '-sp', '100',
                        '-id', '0')
        for lib_num in range(1, 4):
            factory.next_lib()
            while True:
                read = factory.next_read()
                if not read:
                    break
                self.check_shared_species(lib_num, read, 3)

        self.assertEqual(len(self.sources), 3)
        for lib in self.sources.values():
            self.assertEqual(len(lib), 5)
        self.assertEqual(len(self.shared), 5)
        self.sources.clear()
        self.shared.clear()  
        
    def test_inequal_richness(self):
        factory = Biogrinder('-rf', 'data/shotgun_database.fa',
                        '-rs', '1234',
                        '-am', 'uniform',
                        '-lb', '0',
                        '-tr', '100',
                        '-nl', '2',
                        '-di', '3', '5',
                        '-sp', '100',
                        '-id', '0')
        for lib_num in range(1, 4):
            factory.next_lib()
            while True:
                read = factory.next_read()
                if not read:
                    break
                self.check_shared_species(lib_num, read, 2)

        self.assertEqual(len(self.sources), 2)
        self.assertEqual(len(list(self.sources.values())[0]), 3)
        self.assertEqual(len(list(self.sources.values())[1]), 5)
        self.assertEqual(len(self.shared), 3)
        self.sources.clear()
        self.shared.clear() 


if __name__ == '__main__':
    unittest.main()
