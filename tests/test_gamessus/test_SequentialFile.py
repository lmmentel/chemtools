
import unittest
from chemtools.gamessus import SequentialFile
import os

class TestSequentialFile(unittest.TestCase):

    def setUp(self):
        seqfile = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data/he_mini_hf.F08')
        self.daf = SequentialFile(seqfile)

    def tearDown(self):
        pass

    def test_read_irect_he(self):
        pass

if __name__ == "__main__":
    unittest.main()
