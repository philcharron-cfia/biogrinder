import os
import re
import sys
import unittest

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src'))

from Biogrinder import Biogrinder

class Test_19_GeneCopyBias(unittest.TestCase):
    def test_abundance_single_library_no_copy_bias(self):
        factory = Biogrinder('-rf', 'data/multiple_amplicon_database.fa',
                             '-fr', 'data/forward_reverse_primers.fa',
                             '-af', 'data/abundances2.txt',
                             '-cb', '0',
                             '-rd', '48',
                             '-tr', '1000',
                             '-un', '1',
                             '-rs', '12345',
                             '-id', '0')
        factory.next_lib()
        sources = {}
        while True:
            read = factory.next_read()
            if not read:
                break
            source = read.reference_id
            source = re.search(r'seq\d+', source).group()

            sources[source] = sources.get(source, 0) + 1
        
        self.assertIn('seq1', sources)
        self.assertIn('seq2', sources)
        self.assertIn('seq3', sources)
        self.assertNotIn('seq4', sources)
        
        self.assertTrue(580 <= sources['seq1'] <= 620)
        self.assertTrue(280 <= sources['seq2'] <= 320)
        self.assertTrue(80 <= sources['seq3'] <= 120) 

    def test_abundance_single_library_copy_bias(self):
        factory = Biogrinder('-rf', 'data/multiple_amplicon_database.fa',
                             '-fr', 'data/forward_reverse_primers.fa',
                             '-af', 'data/abundances2.txt',
                             '-cb', '1',
                             '-rd', '48',
                             '-tr', '1000',
                             '-un', '1',
                             '-rs', '12345',
                             '-id', '0')
        factory.next_lib()
        sources = {}
        while True:
            read = factory.next_read()
            if not read:
                break
            source = read.reference_id
            source = re.search(r'seq\d+', source).group()

            sources[source] = sources.get(source, 0) + 1
        
        self.assertIn('seq1', sources)
        self.assertIn('seq2', sources)
        self.assertIn('seq3', sources)
        self.assertNotIn('seq4', sources)
        
        self.assertTrue(367 <= sources['seq1'] <= 407)
        self.assertTrue(560 <= sources['seq2'] <= 600)
        self.assertTrue(12 <= sources['seq3'] <= 52)     

if __name__ == '__main__':
    unittest.main()