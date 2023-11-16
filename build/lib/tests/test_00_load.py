import sys
import os
import unittest

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src'))

import main
import AmpliconSearch
import arguments
import Biogrinder
import KmerCollection
import misc
import reports
import SimulatedRead

class Test_00_Load(unittest.TestCase):
    def test_imports(self):
        self.assertIsNotNone(main)
        self.assertIsNotNone(AmpliconSearch)
        self.assertIsNotNone(arguments)
        self.assertIsNotNone(Biogrinder)
        self.assertIsNotNone(KmerCollection)
        self.assertIsNotNone(misc)
        self.assertIsNotNone(reports)
        self.assertIsNotNone(SimulatedRead)        

    def test_version(self):
        self.assertEqual(main.__version__, '0.1.0')

if __name__ == '__main__':
    unittest.main()