
import unittest
from chemtools.gamessus import GamessReader, readseq


class TestGamessReader(unittest.TestCase):

    def setUp(self):
        pass

    def test_he_mini_hf_2electron_ints(self):
        ints = readseq('data/he_mini_hf.F08')
        self.assertAlmostEqual(ints[0], 0.12856711538353)


    def tearDown(self):
        pass
