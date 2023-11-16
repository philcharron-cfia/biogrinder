import os
import statistics
import sys
import unittest
from functions_test import *
from math import exp, log

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src'))

from Biogrinder import Biogrinder



class Test_20_CommunityStructure(unittest.TestCase):
    nof_refs = 5
    max_refs = 10

    def test_uniform_community_structure(self):
        factory = Biogrinder('-rf', 'data/shotgun_database_extended.fa',
                             '-rd', '48',
                             '-tr', '1000',
                             '-lb', '0',
                             '-am', 'uniform', '0',
                             '-id', '0')
        factory.next_lib()
        reads = []
        while True:
            read = factory.next_read()
            if not read:
                break
            reads.append(read.reference_id)
        ra = rank_abundance(reads, self.max_refs)
        mean_value = statistics.mean(ra)
        era = uniform_cstruct(self.max_refs, self.nof_refs, 1000)
        coeff = corr_coeff(ra, era, mean_value)
        self.assertTrue(coeff > 0.97)

    def test_linear_community_structure(self):
        factory = Biogrinder('-rf', 'data/shotgun_database_extended.fa',
                             '-rd', '48',
                             '-tr', '1000',
                             '-lb', '0',
                             '-am', 'linear', '0',
                             '-id', '0')
        factory.next_lib()
        reads = []
        while True:
            read = factory.next_read()
            if not read:
                break
            reads.append(read.reference_id)
        ra = rank_abundance(reads, self.max_refs)
        mean_value = statistics.mean(ra)
        era = linear_cstruct(self.max_refs, self.nof_refs, 1000)
        coeff = corr_coeff(ra, era, mean_value)
        self.assertTrue(coeff > 0.97)

    def test_powerlaw_community_structure(self):
        factory = Biogrinder('-rf', 'data/shotgun_database_extended.fa',
                             '-rd', '48',
                             '-tr', '1000',
                             '-lb', '0',
                             '-am', 'powerlaw', '0.5',
                             '-id', '0')
        factory.next_lib()
        reads = []
        while True:
            read = factory.next_read()
            if not read:
                break
            reads.append(read.reference_id)
        ra = rank_abundance(reads, self.max_refs)
        mean_value = statistics.mean(ra)
        era = powerlaw_cstruct(self.max_refs, self.nof_refs, 0.5, 1000)
        coeff = corr_coeff(ra, era, mean_value)
        self.assertTrue(coeff > 0.97)
    
    def test_logarithmic_community_structure(self):
        factory = Biogrinder('-rf', 'data/shotgun_database_extended.fa',
                             '-rd', '48',
                             '-tr', '1000',
                             '-lb', '0',
                             '-am', 'logarithmic', '0.5',
                             '-id', '0')
        factory.next_lib()
        reads = []
        while True:
            read = factory.next_read()
            if not read:
                break
            reads.append(read.reference_id)
        ra = rank_abundance(reads, self.max_refs)
        mean_value = statistics.mean(ra)
        era = logarithmic_cstruct(self.max_refs, self.nof_refs, 0.5, 1000)
        coeff = corr_coeff(ra, era, mean_value)
        self.assertTrue(coeff > 0.97)
    
    def test_exponential_community_structure(self):
        factory = Biogrinder('-rf', 'data/shotgun_database_extended.fa',
                             '-rd', '48',
                             '-tr', '1000',
                             '-lb', '0',
                             '-am', 'exponential', '0.5',
                             '-id', '0')
        factory.next_lib()
        reads = []
        while True:
            read = factory.next_read()
            if not read:
                break
            reads.append(read.reference_id)
        ra = rank_abundance(reads, self.max_refs)
        mean_value = statistics.mean(ra)
        era = exponential_cstruct(self.max_refs, self.nof_refs, 0.5, 1000)
        coeff = corr_coeff(ra, era, mean_value)
        self.assertTrue(coeff > 0.97)

def exponential_cstruct(x_max, max, param, num):
    ys = []
    total = 0
    for x in range(1, max + 1):
        y = exp(-x * param)
        total += y
        ys.append(y)

    for x in range(max):
        ys[x] *= num / total

    ys.extend([0] * (x_max - len(ys)))
    return ys

def linear_cstruct(x_max, max, num):
    ys = []
    total = 0
    for x in range(max, 0, -1):
        y = x
        total += y
        ys.append(y)

    for x in range(max):
        ys[x] *= num / total

    ys.extend([0] * (x_max - len(ys)))
    return ys

def logarithmic_cstruct(x_max, max, param, num):
    ys = []
    total = 0
    for x in range(1, max + 1):
        y = (log(x + 1)) ** (-param)
        total += y
        ys.append(y)

    for x in range(max):
        ys[x] *= num / total

    ys.extend([0] * (x_max - len(ys)))
    return ys

def powerlaw_cstruct(x_max, max, param, num):
    ys = []
    total = 0
    for x in range(1, max + 1):
        y = x ** (-param)
        total += y
        ys.append(y)

    for x in range(max):
        ys[x] *= num / total

    ys.extend([0] * (x_max - len(ys)))
    return ys

def rank_abundance(data, max):
    hash_map = {}
    for val in data:
        hash_map[val] = hash_map.get(val, 0) + 1

    y_data = sorted(hash_map.values(), reverse=True)
    y_data.extend([0] * (max - len(y_data)))
    return y_data

def uniform_cstruct(x_max, max, num):
    ys = []
    width = max
    for x in range(1, x_max + 1):
        y = num / width if x <= max else 0
        ys.append(y)
    return ys
if __name__ == '__main__':
    unittest.main()