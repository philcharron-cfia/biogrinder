import os
import sys
import unittest

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src'))

from Biogrinder import Biogrinder

class Test_06_Seed(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.factory = Biogrinder('-rf', 'data/shotgun_database_extended.fa',
                                      '-tr', '10',
                                      '-id', '0')
        cls.seed2 = cls.factory.random_seed
        cls.factory.next_lib()        
        cls.dataset1 = []
        while True:
            read = cls.factory.next_read()
            if not read:
                break
            cls.dataset1.append(read)

    def test_set_the_seed(self):
        factory = Biogrinder('-rf', 'data/shotgun_database_extended.fa',
                             '-rs', '1234',
                             '-tr', '10',
                             '-id', '0')
        seed1 = factory.random_seed
        self.assertEqual(seed1, 1234, "Set the seed")
 
    def test_get_seed_automatically(self):
        self.assertGreater(Test_06_Seed.seed2, 0, "Get a seed automatically")

    def test_specify_the_same_seed(self):
        factory = Biogrinder('-rf', 'data/shotgun_database_extended.fa',
                             '-tr', '10',
                             '-rs', str(Test_06_Seed.seed2),
                             '-id', '0')
        seed3 = factory.random_seed
        self.assertEqual(seed3, Test_06_Seed.seed2, "Specify the same seed")
        factory.next_lib()   
        dataset2 = []
        while True:
            read = factory.next_read()
            if not read:
                break
            dataset2.append(read)

        for read1, read2 in zip(Test_06_Seed.dataset1, dataset2):
            self.assertEqual(read1.seq, read2.seq, "Sequences generated with the same seed should be equal")
            self.assertEqual(read1.id, read2.id, "IDs generated with the same seed should be equal")

if __name__ == '__main__':
    unittest.main()
