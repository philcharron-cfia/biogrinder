from Bio.Seq import Seq
import os
import sys
import unittest

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src'))
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'tests'))
current_dir = os.path.dirname(os.path.abspath(__file__))

from functions_test import *
from Biogrinder import Biogrinder

class Test_32_Database(unittest.TestCase):
    def test_create_database(self):
        factory = Biogrinder('-rf', current_dir + '/data/shotgun_database.fa')
        db = factory.database
        ids = sorted(list(db['ids'].keys()))
        self.assertListEqual(ids, ['seq1', 'seq2', 'seq3', 'seq4', 'seq5'])
  
    def test_retrieving_sequences(self):
        factory = Biogrinder('-rf', current_dir + '/data/shotgun_database.fa')
        db = factory.database
        
        self.assertIsNone(get_seq(db, 'zzz'))

        seq_id, seq = get_seq(db, 'seq5')
        self.assertEqual(seq_id, 'seq5')
        self.assertEqual(seq, 'AAAAAAAAAATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTGGGGGGGGGG')
        
        seq_id, seq = get_seq(db, 'seq5:2..11')
        self.assertEqual(seq_id, 'seq5')
        self.assertEqual(seq, 'AAAAAAAAAT')
        
        seq_id, seq = get_seq(db, 'seq5/-1')
        self.assertEqual(seq_id, 'seq5')
        self.assertEqual(seq, 'CCCCCCCCCCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTTTTTTTTT')
        
        seq_id, seq = get_seq(db, 'seq5:2..11/-1')
        self.assertEqual(seq_id, 'seq5')
        self.assertEqual(seq, 'ATTTTTTTTT')      
 
def get_seq(db, id):
    parts = id.split('/')
    id_strand = parts[0]
    strand = 1 if len(parts) == 1 else int(parts[1])
    
    parts = id_strand.split(':')
    id = parts[0]
    start_stop = None if len(parts) == 1 else parts[1]
    start = stop = None
    if start_stop:
        start, stop = map(int, start_stop.split('..'))
    if id not in db['ids']:
        return None

    seqs = db['db']
    for seq_record in seqs:
        if seq_record.id == id:
            seq_id = seq_record.id
            
            if start_stop:
                seq = str(seq_record.seq)[start-1:stop]         
                
            else:
                seq = str(seq_record.seq)
                
    if (strand < 0):
            seq = str(Seq(seq).reverse_complement())

    return seq_id, seq


if __name__ == '__main__':
    unittest.main()

