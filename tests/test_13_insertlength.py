import os
import statistics
import sys
import unittest
from functions_test import *
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src'))

from Biogrinder import Biogrinder

def insert_length(mate1, mate2):
    if mate1.end > mate2.end:
        mate1, mate2 = mate2, mate1
    length = mate2.end - mate1.start + 1
    return length

class Test_13_InsertLength(unittest.TestCase):
    def setUp(self):
        self.lengths = []

    def test_same_size_inserts(self):
        factory = Biogrinder('-rf', 'data/single_seq_database.fa',
                            '-tr', '1000',
                            '-rs', '123456',
                            '-rd', '50',
                            '-id', '150')
        factory.next_lib()
        while True:
            read1 = factory.next_read()
            read2 = factory.next_read()
            if not read1 or not read2:
                break
            self.lengths.append(insert_length(read1, read2))
        
        min_value = min(self.lengths)
        max_value = max(self.lengths)
        mean_value = statistics.mean(self.lengths)
        stdev_value = statistics.stdev(self.lengths)
        self.assertEqual(min_value, 150, 'Minimum value not 150')
        self.assertEqual(max_value, 150, 'Maximum value not 150')
        self.assertEqual(mean_value, 150, 'Mean value not 150')
        self.assertEqual(stdev_value, 0, 'Standard deviation value not 0')

        self.lengths.clear()
        
    def test_uniform_distribution_size_inserts(self):
        factory = Biogrinder('-rf', 'data/single_seq_database.fa',
                            '-tr', '1000',
                            '-rs', '123456',
                            '-rd', '50',
                            '-id', '150', 'uniform', '15')
        factory.next_lib()
        while True:
            read1 = factory.next_read()
            read2 = factory.next_read()
            if not read1 or not read2:
                break
            self.lengths.append(insert_length(read1, read2))

        min_value = min(self.lengths)
        max_value = max(self.lengths)
        mean_value = statistics.mean(self.lengths)
        stdev_value = statistics.stdev(self.lengths)
        self.assertTrue(min_value >= 135)
        self.assertTrue(max_value <= 165)
        self.assertTrue(148 <= mean_value <= 152)
        self.assertTrue(7 <= stdev_value <= 10)

        histogram = hist(self.lengths, 50, 250)
        ehistogram = uniform(50, 250, 135, 165, 1000)
        coeff = corr_coeff(histogram, ehistogram, mean_value)
        self.assertTrue(coeff > 0.99)
        
        self.lengths.clear()
        
    def test_normal_distribution_length_reads(self):
        factory = Biogrinder('-rf', 'data/single_seq_database.fa',
                            '-tr', '1000',
                            '-rs', '123456',
                            '-rd', '50',
                            '-id', '150', 'normal', '10')
        factory.next_lib()
        while True:
            read1 = factory.next_read()
            read2 = factory.next_read()
            if not read1 or not read2:
                break
            self.lengths.append(insert_length(read1, read2))
        
        mean_value = statistics.mean(self.lengths)
        stdev_value = statistics.stdev(self.lengths)
        self.assertTrue(149 <= mean_value <= 151)
        self.assertTrue(9 <= stdev_value <= 11)

        histogram = hist(self.lengths, 50, 250)
        ehistogram = normal(50, 250, mean_value, stdev_value**2, 1000)
        coeff = corr_coeff(histogram, ehistogram, mean_value)
        self.assertTrue(coeff > 0.99)

        self.lengths.clear()        
        
if __name__ == '__main__':
    unittest.main()