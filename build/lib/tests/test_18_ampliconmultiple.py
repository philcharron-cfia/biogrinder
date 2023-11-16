import os
import re
import sys
import unittest

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src'))

from Biogrinder import Biogrinder
from SimulatedRead import SimulatedRead


class Test_18_AmpliconMultiple(unittest.TestCase):
    def test_forward_reverse_primers_multiple_amplicons(self):
        factory = Biogrinder('-rf', 'data/multiple_amplicon_database.fa',
                             '-fr', 'data/forward_reverse_primers.fa',
                             '-lb', '0',
                             '-rd', '100',
                             '-tr', '100',
                             '-un', '1',
                             '-id', '0')
        factory.next_lib()
        nof_reads = 0
        got_amplicons = {}
        while True:
            read = factory.next_read()
            if not read:
                break
            amp_desc = re.search(r'seq\d+_\d+F\+\d+R(_rev)?_LEN\d+_\d+', read.description).group()
            nof_reads += 1
            got_amplicons[amp_desc] = None
            ok_read(read, 1, nof_reads)
        self.assertEqual(nof_reads, 100)
        expected_amplicons = {
            'seq1_926F+1392R_LEN95_0': None,
            'seq1_926F+1392R_LEN335_1': None,
            'seq1_926F+1392R_LEN175_2': None,
            'seq2_926F+1392R_LEN95_0': None,
            'seq2_926F+1392R_LEN195_1': None, 
            'seq2_926F+1392R_LEN355_2': None,
            'seq2_926F+1392R_LEN455_3': None, 
            'seq2_926F+1392R_LEN615_4': None,
            'seq2_926F+1392R_LEN95_5': None,
            'seq2_926F+1392R_LEN195_6': None,
            'seq2_926F+1392R_LEN355_7': None, 
            'seq2_926F+1392R_LEN95_8': None,
            'seq3_926F+1392R_LEN95_0': None,
            'seq4_926F+1392R_rev_LEN95_0': None, 
            'seq4_926F+1392R_LEN95_0': None,
        }
        assert got_amplicons == expected_amplicons
    
def ok_read(read, req_strand, nof_reads):
    assert isinstance(read, SimulatedRead), "Read is not an instance of SimulatedRead" 
    assert re.match(r'^seq\d+.*', read.reference_id), "Reference ID does not match pattern"
    strand = read.strand
    if req_strand is None:
        req_strand = strand
    else:
        assert strand == req_strand, "Strand does not match required strand" 
      
    assert read.id == str(nof_reads), "Read ID does not match the number of reads"

if __name__ == '__main__':
    unittest.main()