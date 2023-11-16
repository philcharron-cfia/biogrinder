import os
import statistics
import sys
import unittest
from functions_test import *
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src'))

from Biogrinder import Biogrinder

class Test_12_ReadLength(unittest.TestCase):
    def setUp(self):
        self.lengths = []

    def test_same_length_reads(self):
        factory = Biogrinder('-rf', 'data/shotgun_database.fa',
                            '-tr', '100',
                            '-rs', '123456',
                            '-rd', '25',
                            '-id', '0')
        factory.next_lib()
        while True:
            read = factory.next_read()
            if not read:
                break
            self.lengths.append(len(read))
        
        min_value = min(self.lengths)
        max_value = max(self.lengths)
        mean_value = statistics.mean(self.lengths)
        stdev_value = statistics.stdev(self.lengths)
        self.assertEqual(min_value, 25, 'Minimum value not 25')
        self.assertEqual(max_value, 25, 'Maximum value not 25')
        self.assertEqual(mean_value, 25, 'Mean value not 25')
        self.assertEqual(stdev_value, 0, 'Standard deviation value not 0')

        self.lengths.clear()
 
    def test_uniform_distribution_length_reads(self):
        factory = Biogrinder('-rf', 'data/shotgun_database.fa',
                            '-tr', '1000',
                            '-rs', '123456',
                            '-rd', '50', 'uniform', '10',
                            '-id', '0')
        factory.next_lib()
        while True:
            read = factory.next_read()
            if not read:
                break
            self.lengths.append(len(read))

        min_value = min(self.lengths)
        max_value = max(self.lengths)
        mean_value = round(statistics.mean(self.lengths))
        stdev_value = statistics.stdev(self.lengths)
        self.assertEqual(min_value, 40, 'Minimum value not 40')
        self.assertEqual(max_value, 60, 'Maximum value not 60')
        self.assertEqual(mean_value, 50, 'Mean value not 50')
        self.assertTrue(5.3 <= stdev_value <= 6.3)

        histogram = hist(self.lengths, 1, 100)
        ehistogram = uniform(1, 100, 40, 60, 1000)
        coeff = corr_coeff(histogram, ehistogram, mean_value)
        self.assertTrue(coeff > 0.99)
        
        self.lengths.clear()

    def test_normal_distribution_length_reads(self):
        factory = Biogrinder('-rf', 'data/shotgun_database.fa',
                            '-tr', '1000',
                            '-rs', '12345',
                            '-rd', '50', 'normal', '5',
                            '-id', '0')
        factory.next_lib()
        while True:
            read = factory.next_read()
            if not read:
                break
            self.lengths.append(len(read))
        
        mean_value = round(statistics.mean(self.lengths))
        stdev_value = statistics.stdev(self.lengths)
        self.assertTrue(49 <= mean_value <= 51)
        self.assertTrue(4.5 <= stdev_value <= 5.5)

        histogram = hist(self.lengths, 1, 100)
        ehistogram = normal(1, 100, mean_value, stdev_value**2, 1000)
        coeff = corr_coeff(histogram, ehistogram, mean_value)
        self.assertTrue(coeff > 0.99)

        self.lengths.clear()        

if __name__ == '__main__':
    unittest.main()