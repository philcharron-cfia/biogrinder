import os
import sys
import unittest

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src'))

from Biogrinder import Biogrinder

class Test_10_Quality(unittest.TestCase):
    
    def test_with_quality_scores(self):
        factory = Biogrinder('-rf', 'data/shotgun_database_extended.fa',
                            '-rd', '52',
                            '-tr', '10',
                            '-ql', '30', '10',
                            '-id', '0')
        
        self.assertIsNotNone(factory, 'No quality scores')
        factory.next_lib()
        read = factory.next_read()
        qual = read.letter_annotations["phred_quality"]
        self.assertEqual(len(qual), 52, 'Quality score list should be of length 52')
        self.assertEqual(qual, [30] * 52, 'Quality scores should match expected values')
       
if __name__ == '__main__':
    unittest.main()

