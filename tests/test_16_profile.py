import os
import sys
import unittest

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src'))
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'tests'))
current_dir = os.path.dirname(os.path.abspath(__file__))

from functions_test import *
from Biogrinder import Biogrinder

class Test_16_Profile(unittest.TestCase):
    def test_no_profile(self):
        factory = Biogrinder('-rf', current_dir + '/data/single_seq_database.fa',
                             '-tr', '100',
                             '-rd', '50',
                             '-un', '1',
                             '-rs', '1234',
                             '-id', '0')
        factory.next_lib()
        while True:
            read = factory.next_read()
            if not read:
                break
            self.assertEqual(str(read.seq), 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA')
    
    def test_profile(self):
        factory = Biogrinder('-pf', current_dir + '/data/profile.txt')
        factory.next_lib()
        while True:
            read = factory.next_read()
            if not read:
                break
            self.assertEqual(str(read.seq), 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA')

    def test_mix_profile(self):
        factory = Biogrinder('-pf', current_dir + '/data/profile.txt',
                             '-dt', '0',
                             '-nl', '2',
                             '-mi', current_dir + '/data/mids.fa',
                             '-sp', '100')
        factory.next_lib()
        while True:
            read = factory.next_read()
            if not read:
                break
            self.assertEqual(str(read.seq), 'ACGTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA')
            self.assertEqual(read.description, '<unknown description>', 'Tracking should not be defined')  
                                              
if __name__ == '__main__':
    unittest.main()