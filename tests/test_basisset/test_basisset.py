import os
import numpy as np
import unittest
from chemtools import basisset
from chemtools.basisset import BasisSet

class TestSequenceFunctions(unittest.TestCase):

    def test_eventemp_if_nf_is_integer(self):
        self.assertRaises(TypeError, basisset.eventemp, (0.0, (1.0, 2.0)))

    def test_eventemp_if_number_of_params_is_2(self):
        self.assertRaises(ValueError, basisset.eventemp, 2, (1.0, 2.0, 3.0))

    def test_eventemp_zero_nf(self):
        self.assertTrue(np.allclose(basisset.eventemp(0, (1.0, 2.0)),
            np.array([])))

    def test_eventemp_small(self):
        self.assertTrue(np.allclose(basisset.eventemp(1, (2.0, 3.0)),
            np.array([2.0])))

    def test_eventemp_medium(self):
        self.assertTrue(np.allclose(basisset.eventemp(4, (0.5, 2.0)),
            np.array([4.0, 2.0, 1.0, 0.5])))

    def test_welltemp_if_nf_is_integer(self):
        self.assertRaises(TypeError, basisset.welltemp, (0.0, (1.0, 2.0, 3.0, 4.0)))

    def test_welltemp_if_number_of_params_is_4(self):
        self.assertRaises(ValueError, basisset.welltemp, 2, (1.0, 2.0, 3.0, 4.0, 5.0))

    def test_welltemp_zero_nf(self):
        self.assertTrue(np.allclose(basisset.welltemp(0, (1.0, 2.0, 3.0, 4.0)),
            np.array([])))

    def test_welltemp_small(self):
        self.assertTrue(np.allclose(basisset.welltemp(1, (1.0, 2.0, 3.0, 4.0)),
            np.array([4.0])))

    def test_welltemp_medium(self):
        self.assertTrue(np.allclose(basisset.welltemp(4, (0.5, 2.0, 3.0, 4.0)),
            np.array([16.0, 3.8984375, 1.1875, 0.505859375])))

    def test_legendre_if_nf_is_integer(self):
        self.assertRaises(TypeError, basisset.legendre, (0.0, (1.0, 2.0, 3.0, 4.0)))

    def test_legendre_if_number_of_params_is_at_least_1(self):
        self.assertRaises(ValueError, basisset.legendre, 2, ())

    def test_legendre_zero_nf(self):
        self.assertTrue(np.allclose(basisset.legendre(0, (1.0, 2.0, 3.0, 4.0)),
            np.array([])))

    def test_legendre_small(self):
        self.assertTrue(np.allclose(basisset.legendre(2, (1.0, 1.0)),
            np.array([7.3890560989306504, 1.0])))

    def test_legendre_medium(self):
        self.assertTrue(np.allclose(basisset.legendre(4, (2.0, 1.5, 0.5)),
            np.array([54.59815003,  10.3122585 ,   3.79366789,   2.71828183])))

class TestBasisSetInitialization(unittest.TestCase):

    def test_init_eventemp_s(self):

        f = [('s', 'et', 5, (0.5, 2.0))]
        bs = BasisSet.from_sequence(functs=f, element='H', name='small et')
        self.assertIsInstance(bs, BasisSet)
        self.assertEqual(bs.name, 'small et')
        self.assertEqual(bs.element, 'H')
        self.assertTrue(np.allclose(bs.functions['s']['e'],
            np.array([8.0, 4.0, 2.0, 1.0, 0.5])))

    def test_init_eventemp_sp(self):

        f = [('s', 'et', 5, (0.5, 2.0)), ('p', 'et', 4, (1.0, 3.0))]
        bs = BasisSet.from_sequence(functs=f, element='H', name='small et')
        self.assertIsInstance(bs, BasisSet)
        self.assertTrue(np.allclose(bs.functions['s']['e'],
            np.array([8.0, 4.0, 2.0, 1.0, 0.5])))
        self.assertTrue(np.allclose(bs.functions['p']['e'],
            np.array([27.0, 9.0, 3.0, 1.0])))

    def test_init_eventemp_df(self):

        f = [('d', 'et', 4, (0.5, 2.0)), ('f', 'et', 3, (1.0, 3.0))]
        bs = BasisSet.from_sequence(functs=f, element='H', name='small et')
        self.assertIsInstance(bs, BasisSet)
        self.assertTrue(np.allclose(bs.functions['d']['e'],
            np.array([4.0, 2.0, 1.0, 0.5])))
        self.assertTrue(np.allclose(bs.functions['f']['e'],
            np.array([9.0, 3.0, 1.0])))

    def test_init_welltemp_s(self):

        f = [('s', 'wt', 4, (0.5, 2.0, 3.0, 4.0))]
        bs = BasisSet.from_sequence(functs=f, element='H', name='small wt')
        self.assertIsInstance(bs, BasisSet)
        self.assertTrue(np.allclose(bs.functions['s']['e'],
            np.array([16.000, 3.8984375, 1.1875, 0.50585938])))

    def test_init_welltemp_sp(self):

        wt_s = np.array([30.4,  13.78515502,   6.21305899,   2.78349551,
                         1.24082247,   0.55241203])
        wt_p = np.array([72.9, 21.58054659, 6.35310115, 1.86463574, 0.54936271])
        f = [('s', 'wt', 6, (0.5, 2.0, 0.9, 1.2)), ('p', 'wt', 5, (0.5, 3.0, 0.8, 1.3))]
        bs = BasisSet.from_sequence(functs=f, element='H', name='small wt')
        self.assertIsInstance(bs, BasisSet)
        self.assertTrue(np.allclose(bs.functions['s']['e'], wt_s))
        self.assertTrue(np.allclose(bs.functions['p']['e'], wt_p))

    def test_init_legendre_s(self):

        le_s = np.array([ 20.08553692,   5.29449005,   1.39561243,   0.36787944])

        f = [('s', 'le', 4, (1.0, 2.0))]
        bs = BasisSet.from_sequence(functs=f, element='H', name='small et')
        self.assertIsInstance(bs, BasisSet)
        self.assertTrue(np.allclose(bs.functions['s']['e'], le_s))

    def test_init_legendre_spd(self):

        leg_s = np.array([1.48029993e+03,   7.57611094e+00,   1.22140276e+00,
                          6.53769785e-01,   1.22456428e-01])
        leg_p = np.array([330.29955991,    5.21663132,    4.70275331,    1.8221188])
        leg_d = np.array([36.59823444,   1.8221188,   1.56831219])
        f = [('s', 'le', 5, (1.0, 3.5, 1.6, 1.2)), ('p', 'le', 4, (2.0, 1.5, 1.2, 1.1)), ('d', 'le', 3, (1.0, 1.5, 1.1))]
        bs = BasisSet.from_sequence(functs=f, element='H', name='small et')
        self.assertIsInstance(bs, BasisSet)
        self.assertTrue(np.allclose(bs.functions['s']['e'], leg_s))
        self.assertTrue(np.allclose(bs.functions['p']['e'], leg_p))
        self.assertTrue(np.allclose(bs.functions['d']['e'], leg_d))

    def test_init_mixed_formula_spd(self):

        et_s = np.array([8.0, 4.0, 2.0, 1.0, 0.5])
        wt_p = np.array([16.000, 3.8984375, 1.1875, 0.50585938])
        leg_d = np.array([36.59823444,   1.8221188,   1.56831219])
        f = [('s', 'et', 5, (0.5, 2.0)), ('p', 'wt', 4, (0.5, 2.0, 3.0, 4.0)), ('d', 'le', 3, (1.0, 1.5, 1.1))]
        bs = BasisSet.from_sequence(functs=f, element='H', name='small et')
        self.assertIsInstance(bs, BasisSet)
        self.assertTrue(np.allclose(bs.functions['s']['e'], et_s))
        self.assertTrue(np.allclose(bs.functions['p']['e'], wt_p))
        self.assertTrue(np.allclose(bs.functions['d']['e'], leg_d))

class TestBasisSetPrintingFormatsUncontracted(unittest.TestCase):

    def setUp(self):

        fs = [('s', 'et', 4, (0.5, 2.0)), ('p', 'et', 3, (1.0, 3.0)), ('d', 'et', 2, (1.5, 4.0))]
        self.bs = BasisSet.from_sequence(name='test', element='H', functs=fs)
        self.maxDiff = None

    def tearDown(self):
        del self.bs

    def test_printing_to_gamessus(self):

        bsstr = '''S  1
  1        4.0000000000     1.00000000
S  1
  1        2.0000000000     1.00000000
S  1
  1        1.0000000000     1.00000000
S  1
  1        0.5000000000     1.00000000
P  1
  1        9.0000000000     1.00000000
P  1
  1        3.0000000000     1.00000000
P  1
  1        1.0000000000     1.00000000
D  1
  1        6.0000000000     1.00000000
D  1
  1        1.5000000000     1.00000000

'''
        self.assertEqual(self.bs.to_gamessus(), bsstr)

    def test_printing_to_molpro(self):

        bsstr = '''s, H, 4.0000000000, 2.0000000000, 1.0000000000, 0.5000000000
c, 1.1, 1.00000000
c, 2.2, 1.00000000
c, 3.3, 1.00000000
c, 4.4, 1.00000000
p, H, 9.0000000000, 3.0000000000, 1.0000000000
c, 1.1, 1.00000000
c, 2.2, 1.00000000
c, 3.3, 1.00000000
d, H, 6.0000000000, 1.5000000000
c, 1.1, 1.00000000
c, 2.2, 1.00000000
'''
        self.assertEqual(self.bs.to_molpro(), bsstr)


    def test_printing_to_cfour(self):

        bsstr = '''
H:test


  3
    0    1    2
    4    3    2
    4    3    2

     4.00000000     2.00000000     1.00000000     0.50000000

     1.00000000     0.00000000     0.00000000     0.00000000
     0.00000000     1.00000000     0.00000000     0.00000000
     0.00000000     0.00000000     1.00000000     0.00000000
     0.00000000     0.00000000     0.00000000     1.00000000

     9.00000000     3.00000000     1.00000000

     1.00000000     0.00000000     0.00000000
     0.00000000     1.00000000     0.00000000
     0.00000000     0.00000000     1.00000000

     6.00000000     1.50000000

     1.00000000     0.00000000
     0.00000000     1.00000000

'''
        self.assertEqual(self.bs.to_cfour(), bsstr)

    def test_printing_to_dalton(self):
        bsstr='''! test
! s functions
H   4   4
        4.0000000000        1.0000000000        0.0000000000        0.0000000000
                            0.0000000000
        2.0000000000        0.0000000000        1.0000000000        0.0000000000
                            0.0000000000
        1.0000000000        0.0000000000        0.0000000000        1.0000000000
                            0.0000000000
        0.5000000000        0.0000000000        0.0000000000        0.0000000000
                            1.0000000000
! p functions
H   3   3
        9.0000000000        1.0000000000        0.0000000000        0.0000000000
        3.0000000000        0.0000000000        1.0000000000        0.0000000000
        1.0000000000        0.0000000000        0.0000000000        1.0000000000
! d functions
H   2   2
        6.0000000000        1.0000000000        0.0000000000
        1.5000000000        0.0000000000        1.0000000000
'''
        self.assertEqual(self.bs.to_dalton(), bsstr)

    def test_printing_to_nwchem(self):

        bsstr='''BASIS "ao basis" PRINT
H s
        4.0000000000     1.00000000
H s
        2.0000000000     1.00000000
H s
        1.0000000000     1.00000000
H s
        0.5000000000     1.00000000
H p
        9.0000000000     1.00000000
H p
        3.0000000000     1.00000000
H p
        1.0000000000     1.00000000
H d
        6.0000000000     1.00000000
H d
        1.5000000000     1.00000000
END
'''
        self.assertEqual(self.bs.to_nwchem(), bsstr)
