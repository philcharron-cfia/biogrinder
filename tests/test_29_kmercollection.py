from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import sys
import unittest
from functions_test import *

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src'))

from KmerCollection import KmerCollection


class Test_29_KmerCollection(unittest.TestCase):
    def setUp(self):
        self.seq1 = SeqRecord(Seq('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'), id='seq1')
        self.seq2 = SeqRecord(Seq('CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCGGGGGGGGCCCCCCCC'), id='seq4')
        self.col = KmerCollection(k=8)
    
    def test_kmer_initialization(self):
        self.assertIsInstance(self.col, KmerCollection)
        self.assertEqual(self.col.k, 8)

    def test_add_seqs(self):
        self.col.add_seqs([self.seq1])
        self.col.add_seqs([self.seq2])
        
        # Check if kmers are added correctly
        self.assertIn('seq1', self.col.collection_by_kmer['AAAAAAAA'])
        self.assertIn('seq4', self.col.collection_by_kmer['AAAAAAAA'])
        self.assertIn('seq4', self.col.collection_by_kmer['CCCCCCCC'])
        self.assertIn('seq4', self.col.collection_by_kmer['CCCCGGGG'])
        self.assertIn('seq4', self.col.collection_by_kmer['ACCCCCCC'])

        self.assertIn('AAAAAAAA', self.col.collection_by_seq['seq1'])
        self.assertIn('AAAAAAAA', self.col.collection_by_seq['seq4'])
        self.assertIn('CCCCCCCC', self.col.collection_by_seq['seq4'])
        self.assertIn('CCCCGGGG', self.col.collection_by_seq['seq4'])
        self.assertIn('ACCCCCCC', self.col.collection_by_seq['seq4'])

    def test_filter_rare(self):
        self.col.add_seqs([self.seq1, self.seq2])
        self.col.filter_rare(2)

        # Check if filtering works correctly
        self.assertIn('seq1', self.col.collection_by_kmer['AAAAAAAA'])
        self.assertIn('seq4', self.col.collection_by_kmer['AAAAAAAA'])
        self.assertIn('seq4', self.col.collection_by_kmer['CCCCCCCC'])
        self.assertNotIn('CCCCGGGG', self.col.collection_by_kmer)
        self.assertNotIn('ACCCCCCC', self.col.collection_by_kmer)

        self.assertIn('AAAAAAAA', self.col.collection_by_seq['seq1'])
        self.assertIn('AAAAAAAA', self.col.collection_by_seq['seq4'])
        self.assertIn('CCCCCCCC', self.col.collection_by_seq['seq4'])
        self.assertNotIn('CCCCGGGG', self.col.collection_by_seq['seq4'])
        self.assertNotIn('ACCCCCCC', self.col.collection_by_seq['seq4'])
    
    def test_counts(self):
        self.col.add_seqs([self.seq1, self.seq2])

        # Count of all kmers
        kmers_counts = self.col.counts()
        kmers = sorted(kmers_counts[0])
        counts = sorted(kmers_counts[1])
        kmer_expected = ['AAAAAAAA','AAAAAAAC','AAAAAACC','AAAAACCC','AAAACCCC',
                         'AAACCCCC','AACCCCCC','ACCCCCCC','CAAAAAAA','CCAAAAAA',
                         'CCCAAAAA','CCCCAAAA','CCCCCAAA','CCCCCCAA','CCCCCCCA',
                         'CCCCCCCC','CCCCCCCG','CCCCCCGG','CCCCCGGG','CCCCGGGG',
                         'CCCGGGGG','CCGGGGGG','CGGGGGGG','GCCCCCCC','GGCCCCCC',
                         'GGGCCCCC','GGGGCCCC','GGGGGCCC','GGGGGGCC','GGGGGGGC',
                         'GGGGGGGG']
        counts_expected = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                           1,1,1,43,74]

        assert kmers == kmer_expected
        assert counts == counts_expected

        # Frequency of kmers from position >= 40
        kmers_freqs = self.col.counts(start=40, freq=1)
        kmers = sorted(kmers_freqs[0])
        freqs = sorted(kmers_freqs[1])
        freqs = [round(num, 5) for num in freqs]
        kmer_expected = ['AAAAAAAA','AACCCCCC','ACCCCCCC','CCCCCCCC','CCCCCCCG',
                         'CCCCCCGG','CCCCCGGG','CCCCGGGG','CCCGGGGG','CCGGGGGG',
                         'CGGGGGGG','GCCCCCCC','GGCCCCCC','GGGCCCCC','GGGGCCCC',
                         'GGGGGCCC','GGGGGGCC','GGGGGGGC','GGGGGGGG']
        freqs_expected = [0.01471,0.01471,0.01471,0.01471,0.01471,0.01471,
                          0.01471,0.01471,0.01471,0.01471,0.01471,0.01471,
                          0.01471,0.01471,0.01471,0.01471,0.01471,0.25,0.5]
        assert kmers == kmer_expected
        assert freqs == freqs_expected

        # Frequency of kmers from position >= 40 in seq1
        kmers_freqs = self.col.counts(id="seq1", start=40, freq=1)
        kmers = sorted(kmers_freqs[0])
        freqs = sorted(kmers_freqs[1])
        kmer_expected = ['AAAAAAAA']
        freqs_expected = [1]
        assert kmers == kmer_expected
        assert freqs == freqs_expected

        # Frequency of kmers from position >= 40 in seq4
        kmers_freqs = self.col.counts(id="seq4", start=40, freq=1)
        kmers = sorted(kmers_freqs[0])
        freqs = sorted(kmers_freqs[1])
        freqs = [round(num, 5) for num in freqs]
        kmer_expected = ['AACCCCCC','ACCCCCCC','CCCCCCCC','CCCCCCCG','CCCCCCGG',
                         'CCCCCGGG','CCCCGGGG','CCCGGGGG','CCGGGGGG','CGGGGGGG',
                         'GCCCCCCC','GGCCCCCC','GGGCCCCC','GGGGCCCC','GGGGGCCC',
                         'GGGGGGCC','GGGGGGGC','GGGGGGGG']
        freqs_expected = [0.02941,0.02941,0.02941,0.02941,0.02941,0.02941,
                          0.02941,0.02941,0.02941,0.02941,0.02941,0.02941,
                          0.02941,0.02941,0.02941,0.02941,0.02941,0.5]
        assert kmers == kmer_expected
        assert freqs == freqs_expected

    def test_counts_filtered(self):
        self.col.add_seqs([self.seq1, self.seq2])
        self.col.filter_shared(2)
        assert isinstance(self.col, KmerCollection)

        # Get counts of all kmers
        kmers, counts = self.col.counts()
        assert sorted(kmers) == ['AAAAAAAA']
        assert sorted(counts) == [74]

        # Check collection by kmer
        by_kmer = self.col.collection_by_kmer
        assert 'AAAAAAAA' in by_kmer
        assert 'seq1' in by_kmer['AAAAAAAA']
        assert 'seq4' in by_kmer['AAAAAAAA']
        assert 'CCCCCCCC' not in by_kmer
        assert 'CCCCGGGG' not in by_kmer
        assert 'ACCCCCCC' not in by_kmer

        # Check collection by sequence
        by_seq = self.col.collection_by_seq
        assert 'AAAAAAAA' in by_seq['seq1']
        assert 'AAAAAAAA' in by_seq['seq4']
        assert 'CCCCCCCC' not in by_seq['seq4']
        assert 'CCCCGGGG' not in by_seq['seq4']
        assert 'ACCCCCCC' not in by_seq['seq4']

        # Check sources for a specific kmer
        sources, counts = self.col.sources('AAAAAAAA')
        assert sources == ['seq1', 'seq4']
        assert counts == [73, 1]

        # Check sources for a specific kmer in a specific sequence
        sources, counts = self.col.sources('AAAAAAAA', 'seq1')
        assert sources == ['seq4']
        assert counts == [1]

        # Check sources for a non-existent kmer
        sources, counts = self.col.sources('ZZZZZZZZ')
        assert sources == []
        assert counts == []

        # Check kmers for specific sequences
        kmers, counts = self.col.kmers('seq1')
        assert kmers == ['AAAAAAAA']
        assert counts == [73]

        kmers, counts = self.col.kmers('seq4')
        assert kmers == ['AAAAAAAA']
        assert counts == [1]

        kmers, counts = self.col.kmers('asdf')
        assert kmers == []
        assert counts == []

        # Check positions of a kmer in a sequence
        positions = self.col.positions('AAAAAAAA', 'seq1')
        assert positions == list(range(1, 74))

        positions = self.col.positions('AAAAAAAA', 'seq4')
        assert positions == [34]

        positions = self.col.positions('CCCCGGGG', 'seq4')
        assert positions == []

        positions = self.col.positions('AAAAAAAA', 'seq3')
        assert positions == []

    def test_kmer_ids(self):
        self.col.add_seqs([self.seq1], ['abc'])
        self.col.add_seqs([self.seq2], ['123'])
        self.col.filter_rare(2)

        assert isinstance(self.col, KmerCollection)

        # Check collection by kmer
        by_kmer = self.col.collection_by_kmer
        assert 'AAAAAAAA' in by_kmer and 'abc' in by_kmer['AAAAAAAA']
        assert 'AAAAAAAA' in by_kmer and '123' in by_kmer['AAAAAAAA']

        # Check collection by sequence
        by_seq = self.col.collection_by_seq
        assert 'AAAAAAAA' in by_seq['abc']
        assert 'AAAAAAAA' in by_seq['123']

        # Check sources for a specific kmer
        sources, counts = self.col.sources('AAAAAAAA')
        assert sources == ['abc', '123']
        assert counts == [73, 1]

        # Check sources for a specific kmer in a specific sequence
        sources, counts = self.col.sources('AAAAAAAA', 'abc')
        assert sources == ['123']
        assert counts == [1]

    def test_kmer_weights(self):
        # Using weights
        self.col.add_seqs([self.seq1, self.seq2])
        self.col.filter_shared(2)
        _weights = {'seq1': 10, 'seq4': 0.1}
        self.col.weights = _weights

        # Check sources and counts with weights
        sources, counts = self.col.sources('AAAAAAAA')
        assert sources == ['seq1', 'seq4']
        assert counts == [730, 0.1]

        # Check kmers and counts
        kmers, counts = self.col.counts()
        assert kmers == ['AAAAAAAA']
        assert counts == [730.1]

        # Check kmers and counts for specific sequences
        kmers, counts = self.col.kmers('seq1')
        assert kmers == ['AAAAAAAA']
        assert counts == [730]

        kmers, counts = self.col.kmers('seq4')
        assert kmers == ['AAAAAAAA']
        assert counts == [0.1]

        # Reset weights
        self.col.weights = {}

    def test_kmer_from_file(self):
        # Read from file (assuming the file contains sequence data)
        file = 'data/kmers.fa'  # Replace with the actual file path
        col = KmerCollection(k=8, file=file)
        assert isinstance(col, KmerCollection)
           
if __name__ == '__main__':
    unittest.main()

