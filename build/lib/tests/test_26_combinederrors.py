import os
import sys
import unittest
from functions_test import *

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src'))

from Biogrinder import Biogrinder
from SimulatedRead import SimulatedRead

class Test_26_CombinedErrors(unittest.TestCase):
    def test_combined_errors_uniform(self):
        factory = Biogrinder('-rf', 'data/shotgun_database_extended.fa',
                             '-un', '1',
                             '-rd', '48',
                             '-tr', '1000',
                             '-hd', 'balzer',
                             '-mr', '100', '0',
                             '-md', 'uniform', '10',
                             '-cp', '10',
                             '-cd', '100',
                             '-ck', '0',
                             '-id', '0',)
        factory.next_lib()
        nof_reads = 0
    

        while True:
            read = factory.next_read()
            if not read:
                break
            nof_reads += 1
            self.assertIsInstance(read, SimulatedRead)
            self.assertEqual(read.id, str(nof_reads))
        self.assertEqual(nof_reads, 1000)
    
    def test_combined_errors_linear(self):
        factory = Biogrinder('-rf', 'data/shotgun_database_extended.fa',
                             '-un', '1',
                             '-rd', '20', 'normal', '10',
                             '-tr', '1000',
                             '-hd', 'balzer',
                             '-mr', '85', '15',
                             '-md', 'linear', '2', '2',
                             '-cp', '10',
                             '-cd', '100',
                             '-ck', '0',
                             '-id', '0',)
        factory.next_lib()
        nof_reads = 0
    

        while True:
            read = factory.next_read()
            if not read:
                break
            nof_reads += 1
            self.assertIsInstance(read, SimulatedRead)
            self.assertEqual(read.id, str(nof_reads))
        self.assertEqual(nof_reads, 1000)
        


    
if __name__ == '__main__':
    unittest.main()

