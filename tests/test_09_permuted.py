import os
import sys
import unittest

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src'))
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'tests'))
current_dir = os.path.dirname(os.path.abspath(__file__))

from functions_test import *
from Biogrinder import Biogrinder

def compare_ranks(ranks1, ranks2, rank1_perm):
    # Top genomes that should be permuted
    perm_ids = ranks1[:rank1_perm] if rank1_perm > 0 else []

    # Copy ranks to dictionaries with genome ID as key and rank as value
    refs1 = {genome_id: rank for rank, genome_id in enumerate(ranks1, 1)}
    refs2 = {genome_id: rank for rank, genome_id in enumerate(ranks2, 1)}

    # Test that permuted genomes have a different rank
    for perm_id in perm_ids:
        rank1 = refs1.get(perm_id)
        rank2 = refs2.get(perm_id)
        assert rank1 != rank2, f"Rank for genome {perm_id} should differ between sets"
        del refs1[perm_id]
        del refs2[perm_id]

    # Remaining genomes should have identical ranks (because it is 100% shared)
    sorted_refs1 = sorted(refs1, key=refs1.get)
    sorted_refs2 = sorted(refs2, key=refs2.get)

    # Length of the smallest array (number of species shared is relative to this)
    min_arr_len = min(len(sorted_refs1), len(sorted_refs2))
    for i in range(min_arr_len):
        id1 = sorted_refs1[i]
        id2 = sorted_refs2[i]
        assert id1 == id2, f"Genome {id1} and {id2} should have the same rank"

    return True


def get_ranks(factory):
    sources = {}
    while True:
        read = factory.next_read()
        if read is None:
            break
        source = read.reference_id
        if source not in sources:
            sources[source] = 1
        else:
            sources[source] += 1
    # Sort the sources by count in descending order and return the list of sorted source IDs
    ranks = sorted(sources, key=sources.get, reverse=True)
    return ranks

class Test_09_Permuted(unittest.TestCase):
    
    def test_no_species_permuted(self):
        factory = Biogrinder('-rf', current_dir + '/data/shotgun_database.fa',
                            '-rs', '1234',
                            '-am', 'powerlaw', '1.8',
                            '-lb', '0',
                            '-tr', '1000',
                            '-nl', '2',
                            '-sp', '100',
                            '-pp', '0',
                            '-id', '0')

        self.assertTrue(factory.next_lib())
        ranks1 = get_ranks(factory)
        self.assertEqual(len(ranks1), 5)
        self.assertTrue(factory.next_lib())
        ranks2 = get_ranks(factory)
        self.assertEqual(len(ranks2), 5)

        rank1_perm = 0
        self.assertTrue(compare_ranks(ranks1, ranks2, rank1_perm))

    def test_40_species_permuted(self):
        factory = Biogrinder('-rf', current_dir + '/data/shotgun_database.fa',
                            '-rs', '1234',
                            '-am', 'powerlaw', '1.8',
                            '-lb', '0',
                            '-tr', '1000',
                            '-nl', '2',
                            '-sp', '100',
                            '-pp', '40',
                            '-id', '0')

        self.assertTrue(factory.next_lib())
        ranks1 = get_ranks(factory)
        self.assertEqual(len(ranks1), 5)
        self.assertTrue(factory.next_lib())
        ranks2 = get_ranks(factory)
        self.assertEqual(len(ranks2), 5)
        
        rank1_perm = 2
        self.assertTrue(compare_ranks(ranks1, ranks2, rank1_perm))

    def test_60_species_permuted(self):
        factory = Biogrinder('-rf', current_dir + '/data/shotgun_database.fa',
                            '-rs', '1234',
                            '-am', 'powerlaw', '1.8',
                            '-lb', '0',
                            '-tr', '1000',
                            '-nl', '2',
                            '-sp', '100',
                            '-pp', '60',
                            '-id', '0')

        self.assertTrue(factory.next_lib())
        ranks1 = get_ranks(factory)
        self.assertEqual(len(ranks1), 5)
        self.assertTrue(factory.next_lib())
        ranks2 = get_ranks(factory)
        self.assertEqual(len(ranks2), 5)
        
        rank1_perm = 3
        self.assertTrue(compare_ranks(ranks1, ranks2, rank1_perm))
    
    def test_80_species_permuted(self):
        factory = Biogrinder('-rf', current_dir + '/data/shotgun_database.fa',
                            '-rs', '1234',
                            '-am', 'powerlaw', '1.8',
                            '-lb', '0',
                            '-tr', '1000',
                            '-nl', '2',
                            '-sp', '100',
                            '-pp', '80',
                            '-id', '0')

        self.assertTrue(factory.next_lib())
        ranks1 = get_ranks(factory)
        self.assertEqual(len(ranks1), 5)
        self.assertTrue(factory.next_lib())
        ranks2 = get_ranks(factory)
        self.assertEqual(len(ranks2), 5)
        
        rank1_perm = 4
        self.assertTrue(compare_ranks(ranks1, ranks2, rank1_perm))

    def test_all_species_permuted(self):
        factory = Biogrinder('-rf', current_dir + '/data/shotgun_database.fa',
                            '-rs', '123456',
                            '-am', 'powerlaw', '1.8',
                            '-lb', '0',
                            '-tr', '1000',
                            '-nl', '2',
                            '-sp', '100',
                            '-pp', '100',
                            '-id', '0')

        self.assertTrue(factory.next_lib())
        ranks1 = get_ranks(factory)
        self.assertEqual(len(ranks1), 5)
        self.assertTrue(factory.next_lib())
        ranks2 = get_ranks(factory)
        self.assertEqual(len(ranks2), 5)
        
        rank1_perm = 5
        self.assertTrue(compare_ranks(ranks1, ranks2, rank1_perm))

    def test_inequal_richness(self):
        factory = Biogrinder('-rf', current_dir + '/data/shotgun_database.fa',
                            '-rs', '123456',
                            '-am', 'powerlaw', '1.8',
                            '-lb', '0',
                            '-tr', '1000',
                            '-nl', '2',
                            '-di', '3', '5',
                            '-sp', '100',
                            '-pp', '100',
                            '-id', '0')
        self.assertTrue(factory.next_lib())
        ranks1 = get_ranks(factory)
        self.assertEqual(len(ranks1), 3)
        self.assertTrue(factory.next_lib())
        ranks2 = get_ranks(factory)
        self.assertEqual(len(ranks2), 5)

        rank1_perm = 3
        self.assertTrue(compare_ranks(ranks1, ranks2, rank1_perm))
if __name__ == '__main__':
    unittest.main()

