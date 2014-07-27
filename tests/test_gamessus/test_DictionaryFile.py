
import unittest
from chemtools.gamessus import DictionaryFile
import os

class TestDictionaryFile(unittest.TestCase):

    def setUp(self):
        dictfile = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data/he_mini_hf.F10')
        self.daf = DictionaryFile(dictfile)

    def tearDown(self):
        pass

    def test_read_irect_he(self):
        pass
    def test_read_ioda_he(self):
        pass
    def test_read_ifilen_he(self):
        pass
    def test_read_is_he(self):
        pass
    def test_read_ipk_he(self):
        pass
    def test_geometry(self):
        pass

if __name__ == "__main__":
    unittest.main()
