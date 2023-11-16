import os
import re
import sys
import unittest

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src'))
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'tests'))
current_dir = os.path.dirname(os.path.abspath(__file__))

from functions_test import *
from Biogrinder import Biogrinder

class Test_11_Tracking(unittest.TestCase):
    def test_bidirectional_shotgun_tracking(self):
        factory = Biogrinder('-rf', current_dir + '/data/shotgun_database_extended.fa',
                            '-tr', '10',
                            '-un', '0',
                            '-dt', '1',
                            '-id', '0')
        factory.next_lib()
        while True:
            read = factory.next_read()
            if not read:
                break
            self.assertIsNotNone(re.search(r'reference=.*position=\(?(complement\()?\d+\.\.\d+(\))?\)?', read.description))

    def test_forward_shotgun_tracking(self):
        factory = Biogrinder('-rf', current_dir + '/data/shotgun_database_extended.fa',
                            '-tr', '10',
                            '-un', '1',
                            '-dt', '1',
                            '-id', '0')
        factory.next_lib()
        while True:
            read = factory.next_read()
            if not read:
                break
            self.assertIsNotNone(re.search(r'reference=.*position=\(?\d+\.\.\d+\)?', read.description))
   
    def test_reverse_shotgun_tracking(self):
        factory = Biogrinder('-rf', current_dir + '/data/shotgun_database_extended.fa',
                            '-tr', '10',
                            '-un', '-1',
                            '-dt', '1',
                            '-id', '0')
        factory.next_lib()
        while True:
            read = factory.next_read()
            if not read:
                break
            self.assertIsNotNone(re.search(r'reference=.*position=(complement\()?\d+\.\.\d+(\))?', read.description))
    
    def test_forward_amplicon_tracking(self):
        factory = Biogrinder('-rf', current_dir + '/data/amplicon_database.fa',
                             '-fr', current_dir + '/data/forward_reverse_primers.fa',
                             '-lb', '0',
                            '-tr', '10',
                            '-un', '1',
                            '-dt', '1',
                            '-id', '0')
        factory.next_lib()
        while True:
            read = factory.next_read()
            if not read:
                break
            self.assertIsNotNone(re.search(r'reference=.*\_\d+F\+\d+R.*position=\(?\d+\.\.\d+\)?', read.description))

    def test_reverse_amplicon_tracking(self):
        factory = Biogrinder('-rf', current_dir + '/data/amplicon_database.fa',
                             '-fr', current_dir + '/data/forward_reverse_primers.fa',
                             '-lb', '0',
                            '-tr', '10',
                            '-un', '-1',
                            '-dt', '1',
                            '-id', '0')
        factory.next_lib()
        while True:
            read = factory.next_read()
            if not read:
                break
            self.assertIsNotNone(re.search(r'reference=.*\_\d+F\+\d+R.*position=(complement\()?\d+\.\.\d+(\))?', read.description))

    def test_bimeric_amplicon_tracking(self):
        factory = Biogrinder('-rf', current_dir + '/data/amplicon_database.fa',
                             '-fr', current_dir + '/data/forward_reverse_primers.fa',
                             '-lb', '0',
                            '-tr', '10',
                            '-un', '1',
                            '-dt', '1',
                            '-cp', '100',
                            '-cd', '1',
                            '-ck', '0',
                            '-id', '0')
        factory.next_lib()
        while True:
            read = factory.next_read()
            if not read:
                break
            self.assertIsNotNone(re.search(r'reference=.*\_\d+F\+\d+R.*,.*\_\d+F\+\d+R.*position=\(?\d+\.\.\d+\)?', read.description))

    def test_no_tracking(self):
        factory = Biogrinder('-rf', current_dir + '/data/shotgun_database.fa',
                            '-tr', '10',
                            '-dt', '0',
                            '-id', '0')
        factory.next_lib()
        read = factory.next_read()
        self.assertEqual(read.description, '<unknown description>', 'Tracking should not be defined')

    def test_default_tracking(self):
        factory = Biogrinder('-rf', current_dir + '/data/shotgun_database_extended.fa',
                            '-tr', '10',
                            '-id', '0')
        factory.next_lib()
        while True:
            read = factory.next_read()
            if not read:
                break
            self.assertIsNotNone(re.search(r'reference=.*position=\(?(complement\()?\d+\.\.\d+(\))?\)?', read.description))

if __name__ == '__main__':
    unittest.main()

