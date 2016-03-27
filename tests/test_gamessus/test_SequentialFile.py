
import unittest
from chemtools.calculators.gamessreader import SequentialFile
import os
import numpy as np

class TestSequentialFileHeMini(unittest.TestCase):

    def setUp(self):
        seqfile = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data/he/he_mini_hf.F08')
        self.seq = SequentialFile(seqfile)

    def tearDown(self):
        self.seq.file.close()

    def test_read_irect_he(self):
        ints = np.array([0.12856711538353], dtype=float)
        self.assertAlmostEqual(self.seq.readseq(), ints)

class TestSequentialFileNeDZAO(unittest.TestCase):

    def setUp(self):
        seqfile = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data/ne/ne_dz_guga.F08')
        self.seq = SequentialFile(seqfile)

    def tearDown(self):
        self.seq.file.close()

    def test_read_ao_integrals_ne(self):
        ints = np.load(os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data/ne/ne_dz_guga_aoints.npy'))
        self.assertTrue(np.allclose(self.seq.readseq(), ints))

class TestSequentialFileNeDZMO(unittest.TestCase):

    def setUp(self):
        seqfile = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data/ne/ne_dz_guga.F09')
        self.seq = SequentialFile(seqfile)

    def tearDown(self):
        self.seq.file.close()

    def test_read_mo_one_electron_integrals_ne(self):
        #ints = np.load(os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data/ne/ne_dz_guga_moints.npy'))
        #self.assertTrue(np.allclose(self.seq.readseq(mos=True), ints))
        pass

    def test_read_mo_two_electron_integrals_ne(self):
        ints = np.load(os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data/ne/ne_dz_guga_moints.npy'))
        self.assertTrue(np.allclose(self.seq.readseq(mos=True, skip_first=True), ints))

class TestSequentialFileNeDZ2RDM(unittest.TestCase):

    def setUp(self):
        seqfile = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data/ne/ne_dz_guga.F15')
        self.seq = SequentialFile(seqfile)

    def tearDown(self):
        self.seq.file.close()

    def test_read_two_particle_density_matrix_ne(self):
        ints = np.load(os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data/ne/ne_dz_guga_rdm2.npy'))
        self.assertTrue(np.allclose(self.seq.readseq(mos=True, skip_first=False), ints))

if __name__ == "__main__":
    unittest.main()
