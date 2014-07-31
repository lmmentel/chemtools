import unittest
from chemtools import basisset

class TestBasisSetModule(unittest.TestCase):

    def test_eventemp_if_nf_is_integer(self):
        self.assertRaises(TypeError, basisset.eventemp, (0.0, (1.0, 2.0)))

    def test_eventemp_if_number_of_params_is_2(self):
        self.assertRaises(ValueError, basisset.eventemp, 2, (1.0, 2.0, 3.0))

    def test_eventemp_zero_nf(self):
        self.assertListEqual(basisset.eventemp(0, (1.0, 2.0)), [])

    def test_eventemp_small(self):
        self.assertListEqual(basisset.eventemp(1, (2.0, 3.0)), [2.0])

    def test_eventemp_medium(self):
        self.assertListEqual(basisset.eventemp(4, (0.5, 2.0)), [0.5, 1.0, 2.0, 4.0])

    def test_welltemp_if_nf_is_integer(self):
        self.assertRaises(TypeError, basisset.welltemp, (0.0, (1.0, 2.0, 3.0, 4.0)))

    def test_welltemp_if_number_of_params_is_4(self):
        self.assertRaises(ValueError, basisset.welltemp, 2, (1.0, 2.0, 3.0, 4.0, 5.0))

    def test_welltemp_zero_nf(self):
        self.assertListEqual(basisset.welltemp(0, (1.0, 2.0, 3.0, 4.0)), [])

    def test_welltemp_small(self):
        self.assertListEqual(basisset.welltemp(1, (1.0, 2.0, 3.0, 4.0)), [4.0])

    def test_welltemp_medium(self):
        self.assertListEqual(basisset.welltemp(4, (0.5, 2.0, 3.0, 4.0)), [0.5, 1.0, 2.0, 16.0])

    def test_legendre_if_nf_is_integer(self):
        self.assertRaises(TypeError, basisset.legendre, (0.0, (1.0, 2.0, 3.0, 4.0)))

    def test_legendre_if_number_of_params_is_at_least_1(self):
        self.assertRaises(ValueError, basisset.legendre, 2, ())

    def test_legendre_zero_nf(self):
        self.assertListEqual(basisset.legendre(0, (1.0, 2.0, 3.0, 4.0)), [])

    def test_legendre_small(self):
        self.assertListEqual(basisset.legendre(2, (1.0, 1.0)), [1.0, 7.3890560989306504])

    def test_legendre_medium(self):
        self.assertListEqual(basisset.legendre(4, (1.0, 1.2, 1.5)), [3.6692966676192444, 1.1051709180756477, 2.459603111156949, 40.447304360067399])

if __name__ == "__main__":
    unittest.main()
