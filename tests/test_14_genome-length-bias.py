import os
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

class Test_14_GenomeLengthBias(unittest.TestCase):

    def test_single_library_genome_abundance(self):
        factory = Biogrinder('-rf', 'data/shotgun_database.fa',
                             '-af', 'data/abundances.txt',
                             '-lb', '1',                             
                            '-tr', '1000',
                            '-rs', '123456',
                            '-id', '0')
        factory.next_lib()
        sources = {}
        while True:
            read = factory.next_read()
            if not read:
                break
            source = read.reference_id
            sources[source] = sources.get(source, 0) + 1
        
        self.assertIn('seq1', sources)
        self.assertIn('seq2', sources)
        self.assertNotIn('seq3', sources)
        self.assertIn('seq4', sources)
        self.assertIn('seq5', sources)
        
        self.assertTrue(414 <= sources['seq1'] <= 477)
        self.assertTrue(303 <= sources['seq2'] <= 363)
        self.assertTrue(81 <= sources['seq4'] <= 141)
        self.assertTrue(81 <= sources['seq5'] <= 141)      
        
if __name__ == '__main__':
    unittest.main()