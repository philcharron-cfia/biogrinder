import os
import sys
import unittest

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src'))
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'tests'))
current_dir = os.path.dirname(os.path.abspath(__file__))

from functions_test import *
from Biogrinder import Biogrinder

class Test_05_Forbidden(unittest.TestCase):

    def test_with_dubious_chars(self):
        factory = Biogrinder('-rf', current_dir + '/data/dirty_database.fa',
                             '-rd', '80',
                             '-rs', '10',
                             '-tr', '10',
                             '-id', '0')
        factory.next_lib()
        read = factory.next_read()
        self.assertRegex(str(read.seq), r'[N-]', 'With dubious chars')


    def test_exclude_chars(self):
        factory = Biogrinder('-rf', current_dir + '/data/dirty_database.fa',
                             '-ec', 'n-',
                             '-rd', '25',
                             '-rs', '10',
                             '-tr', '10',
                             '-id', '0')
        factory.next_lib()
        while True:
            read = factory.next_read()
            if not read:
                break
            self.assertNotRegex(str(read.seq), r'[N-]', 'Exclude chars')

    def test_cannot_generate_read(self):
        factory = Biogrinder('-rf', current_dir + '/data/dirty_database.fa',
                             '-ec', 'N-',
                             '-rs', '1234',
                             '-rd', '71',
                             '-rs', '10',
                             '-tr', '10',
                             '-id', '0')
        factory.next_lib()
        error_message = r'Error: Could not take a random shotgun read without forbidden characters from reference sequence seq2 \(10 attempts made\)\.'
        with self.assertRaisesRegex(Exception, error_message):
            factory.next_read()

    def test_delete_chars(self):
        factory = Biogrinder('-rf', current_dir + '/data/dirty_database.fa',
                             '-dc', 'N-',
                             '-rd', '60',
                             '-rs', '10',
                             '-tr', '10',
                             '-id', '0')
        factory.next_lib()
        while True:
            read = factory.next_read()
            if not read:
                break
            self.assertNotRegex(str(read.seq), r'[N-]', 'Exclude chars')

if __name__ == '__main__':
    unittest.main()
