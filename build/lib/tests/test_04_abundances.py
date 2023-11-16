import os
import sys
import unittest

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src'))

from Biogrinder import Biogrinder

class Test_04_Abundances(unittest.TestCase):

    def test_single_shotgun_library_abundance(self):
        factory = Biogrinder('-rf', 'data/shotgun_database_extended.fa',
                             '-af', 'data/abundances.txt',
                             '-lb', '0',
                             '-rs', '10',
                             '-tr', '1000',
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

        self.assertTrue(220 <= sources['seq1'] <= 280)
        self.assertTrue(220 <= sources['seq2'] <= 280)
        self.assertTrue(220 <= sources['seq4'] <= 280)
        self.assertTrue(220 <= sources['seq5'] <= 280)

        self.assertIsNone(factory.next_lib())
        sources.clear()

    def test_single_amplicon_library_abundance(self):
        factory = Biogrinder('-rf', 'data/amplicon_database.fa',
                             '-af', 'data/abundances2.txt',
                             '-fr', 'data/forward_reverse_primers.fa',
                             '-lb', '0',
                             '-cb', '0',
                             '-un', '1',
                             '-rd', '25',
                             '-rs', '1234',
                             '-tr', '1000',
                             '-st', '1',
                             '-id', '0')
        factory.next_lib()
        sources = {}
        while True:
            read = factory.next_read()
            
            if not read:
                break
            source = read.reference_id
            for prefix in ['seq1', 'seq2', 'seq3']:
                if source.startswith(prefix):
                    sources[prefix] = sources.get(prefix, 0) + 1
                    break

        self.assertIn('seq1', sources)
        self.assertIn('seq2', sources)
        self.assertIn('seq3', sources)
        # Use assert methods to replace between_ok
        self.assertTrue(570 <= sources['seq1'] <= 630)
        self.assertTrue(270 <= sources['seq2'] <= 330)
        self.assertTrue(70 <= sources['seq3'] <= 130)

        self.assertIsNone(factory.next_lib())
        sources.clear()
    def test_multiple_shotgun_libraries_abundance(self):
        factory = Biogrinder('-rf', 'data/shotgun_database_extended.fa',
                             '-af', 'data/abundances_multiple.txt',
                             '-lb', '0',
                             '-rs', '10',
                             '-tr', '1000',
                             '-nl', '3',
                             '-id', '0')
        sources = get_reads(factory)
        sources = get_reads(factory)
        sources = get_reads(factory)

        self.assertNotIn('seq1', sources)
        self.assertNotIn('seq2', sources)
        self.assertIn('seq3', sources)
        self.assertNotIn('seq4', sources)
        self.assertNotIn('seq5', sources)

        self.assertEqual(sources.get('seq3', 0), 1000)

        self.assertIsNone(factory.next_lib())

def get_reads(factory):
    factory.next_lib()
    sources = {}
    nof_reads = 0
    while True:
        read = factory.next_read()
        if not read:
            break
        nof_reads += 1
        source = read.reference_id
        for prefix in ['seq1', 'seq2', 'seq3', 'seq4', 'seq5']:
            if source.startswith(prefix):
                sources[prefix] = sources.get(prefix, 0) + 1
                break
    assert nof_reads == 1000, "Number of reads does not match"

    return sources

if __name__ == '__main__':
    unittest.main()
