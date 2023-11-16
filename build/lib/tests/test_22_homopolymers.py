import os
import re
import statistics
import sys
import unittest
from functions_test import *

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src'))

from Biogrinder import Biogrinder



class Test_22_Homopolymers(unittest.TestCase):
    nof_refs = 5
    max_refs = 10

    def test_no_errors(self):
        factory = Biogrinder('-rf', 'data/single_seq_database.fa',
                             '-rd', '50',
                             '-tr', '100',
                             '-un', '1',
                             '-rs', '123',
                             '-id', '0')
        factory.next_lib()
        while True:
            read = factory.next_read()
            if not read:
                break
            self.assertEqual(str(read.seq), 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA')
            self.assertNotIn('errors', read.description)
        
    def test_substitutions(self):
        factory = Biogrinder('-rf', 'data/single_seq_database.fa',
                             '-rd', '50',
                             '-tr', '1000',
                             '-un', '1',
                             '-mr', '100', '0',
                             '-md', 'uniform', '10',
                             '-rs', '123',
                             '-id', '0')
        factory.next_lib()
        while True:
            read = factory.next_read()
            if not read:
                break
            match = re.search(r'errors=(\S+)', read.description)
            if match:
                error_str = match.group(1)
                self.assertIn('%', error_str)
                self.assertNotIn('+', error_str)
                self.assertNotIn('-', error_str)
            else:
                self.assertTrue(True)

    def test_indels(self):
        factory = Biogrinder('-rf', 'data/single_seq_database.fa',
                             '-rd', '50',
                             '-tr', '1000',
                             '-un', '1',
                             '-mr', '0', '100',
                             '-md', 'uniform', '10',
                             '-rs', '123',
                             '-id', '0')
        factory.next_lib()
        while True:
            read = factory.next_read()
            if not read:
                break
            match = re.search(r'errors=(\S+)', read.description)
            if match:
                error_str = match.group(1)
                self.assertNotIn('%', error_str)
                self.assertRegex(error_str, r'[-+]')
            else:
                self.assertTrue(True)
                        
    def test_indels_substitutions(self):
        factory = Biogrinder('-rf', 'data/single_seq_database.fa',
                             '-rd', '50',
                             '-tr', '1000',
                             '-un', '1',
                             '-mr', '50', '50',
                             '-md', 'uniform', '10',
                             '-rs', '123',
                             '-id', '0')
        factory.next_lib()
        nof_indels = 0
        nof_subs = 0

        while True:
            read = factory.next_read()
            if not read:
                break
            match = re.search(r'errors=(\S+)', read.description)
            if match:
                error_str = match.group(1)
                self.assertRegex(error_str, r'[-+%]')
                nof_indels +=  error_str.count('-') + error_str.count('+')
                nof_subs +=  error_str.count('%')
            else:
                self.assertTrue(True)
        self.assertTrue(0.92 <= nof_indels/nof_subs <= 1.08)
    
    def test_uniform_frequent(self):
        factory = Biogrinder('-rf', 'data/single_seq_database.fa',
                             '-rd', '50',
                             '-tr', '1000',
                             '-un', '1',
                             '-mr', '50', '50',
                             '-md', 'uniform', '10',
                             '-rs', '123456',
                             '-id', '0')
        factory.next_lib()
        epositions = []
    
        while True:
            read = factory.next_read()
            if not read:
                break
            positions = error_positions(read)
            if positions:
                epositions.extend(positions)
        prof = hist(epositions, 1, 50)
        mean_value = statistics.mean(prof)

        self.assertTrue(70 <= prof[0] <= 130)
        self.assertTrue(70 <= prof[24] <= 130)
        self.assertTrue(70 <= prof[-1] <= 130)
        self.assertTrue(97 <= mean_value <= 103)

    def test_uniform_rare(self):
        factory = Biogrinder('-rf', 'data/single_seq_database.fa',
                             '-rd', '50',
                             '-tr', '10000',
                             '-un', '1',
                             '-mr', '50', '50',
                             '-md', 'uniform', '0.1',
                             '-rs', '1234',
                             '-id', '0')
        factory.next_lib()
        epositions = []
    
        while True:
            read = factory.next_read()
            if not read:
                break
            positions = error_positions(read)
            if positions:
                epositions.extend(positions)
        prof = hist(epositions, 1, 50)
        mean_value = statistics.mean(prof)

        self.assertTrue(7 <= prof[0] <= 13)
        self.assertTrue(7 <= prof[24] <= 13)
        self.assertTrue(7 <= prof[-1] <= 13)
        
        self.assertTrue(9 <= mean_value <= 11)

    def test_linear(self):
        factory = Biogrinder('-rf', 'data/single_seq_database.fa',
                             '-rd', '50',
                             '-tr', '1000',
                             '-un', '1',
                             '-mr', '100', '0',
                             '-md', 'linear', '5', '15',
                             '-rs', '123456',
                             '-id', '0')
        factory.next_lib()
        epositions = []
    
        while True:
            read = factory.next_read()
            if not read:
                break
            positions = error_positions(read)
            if positions:
                epositions.extend(positions)
        prof = hist(epositions, 1, 50)
        mean_value = statistics.mean(prof)

        self.assertTrue(30 <= prof[0] <= 70)
        self.assertTrue(70 <= prof[24] <= 130)
        self.assertTrue(120 <= prof[-1] <= 180)
        self.assertTrue(97 <= mean_value <= 103)
        
    def test_polynomial(self):
        factory = Biogrinder('-rf', 'data/single_seq_database.fa',
                             '-rd', '100',
                             '-tr', '1000',
                             '-un', '1',
                             '-mr', '100', '0',
                             '-md', 'poly4', '1', '4.4e-7',
                             '-rs', '123',
                             '-id', '0')
        factory.next_lib()
        epositions = []
    
        while True:
            read = factory.next_read()
            if not read:
                break
            positions = error_positions(read)
            if positions:
                epositions.extend(positions)
        prof = hist(epositions, 1, 100)
        mean_value = statistics.mean(prof)
        
        self.assertTrue(1 <= prof[0] <= 27)
        self.assertTrue(7 <= prof[24] <= 67)
        self.assertTrue(410 <= prof[-1] <= 488)
        self.assertTrue(90 <= mean_value <= 110)

if __name__ == '__main__':
    unittest.main()

