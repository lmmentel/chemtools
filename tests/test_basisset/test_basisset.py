import unittest
from chemtools import basisset
from chemtools.basisset import BasisSet

class TestSequenceFunctions(unittest.TestCase):

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
        self.assertListEqual(basisset.welltemp(4, (0.5, 2.0, 3.0, 4.0)), [0.505859375, 1.1875, 3.8984375, 16.0])

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

class TestBasisSetInitialization(unittest.TestCase):

    def test_init_eventemp_s(self):

        f = [('s', 5, (0.5, 2.0)), ('p', 4, (1.0, 3.0))]
        bs = BasisSet.from_sequence(formula='eventemp', functs=f, element='H', name='small et')
        self.assertIsInstance(bs, BasisSet)
        self.assertEqual(bs.name, 'small et')
        self.assertEqual(bs.element, 'H')
        self.assertAlmostEqual(bs.functions['s']['e'], [0.5, 1.0, 2.0, 4.0, 8.0])

    def test_init_eventemp_sp(self):

        f = [('s', 5, (0.5, 2.0)), ('p', 4, (1.0, 3.0))]
        bs = BasisSet.from_sequence(formula='eventemp', functs=f, element='H', name='small et')
        self.assertIsInstance(bs, BasisSet)
        self.assertAlmostEqual(bs.functions['s']['e'], [0.5, 1.0, 2.0, 4.0, 8.0])
        self.assertAlmostEqual(bs.functions['p']['e'], [1.0, 3.0, 9.0, 27.0])

    def test_init_eventemp_df(self):

        f = [('d', 4, (0.5, 2.0)), ('f', 3, (1.0, 3.0))]
        bs = BasisSet.from_sequence(formula='eventemp', functs=f, element='H', name='small et')
        self.assertIsInstance(bs, BasisSet)
        self.assertAlmostEqual(bs.functions['d']['e'], [0.5, 1.0, 2.0, 4.0])
        self.assertAlmostEqual(bs.functions['f']['e'], [1.0, 3.0, 9.0])

    def test_init_welltemp_s(self):

        f = [('s', 4, (0.5, 2.0, 3.0, 4.0))]
        bs = BasisSet.from_sequence(formula='welltemp', functs=f, element='H', name='small wt')
        self.assertIsInstance(bs, BasisSet)
        self.assertAlmostEqual(bs.functions['s']['e'], [0.505859375, 1.1875, 3.8984375, 16.0])

    def test_init_welltemp_sp(self):

        wt_s = [0.55241203390786842, 1.2408224685280691, 2.7834955069665117,
                6.2130589875561473, 13.785155024015765, 30.399999999999999]
        wt_p = [0.54936270901760076, 1.8646357406075393, 6.3531011529592449,
                21.580546589187453, 72.900000000000006]
        f = [('s', 6, (0.5, 2.0, 0.9, 1.2)), ('p', 5, (0.5, 3.0, 0.8, 1.3))]
        bs = BasisSet.from_sequence(formula='welltemp', functs=f, element='H', name='small wt')
        self.assertIsInstance(bs, BasisSet)
        self.assertAlmostEqual(bs.functions['s']['e'], wt_s)
        self.assertAlmostEqual(bs.functions['p']['e'], wt_p)

    def test_init_legendre_s(self):

        le_s = [0.36787944117144233, 1.3956124250860895, 5.2944900504700287,
                20.085536923187668]

        f = [('s', 4, (1.0, 2.0))]
        bs = BasisSet.from_sequence(formula='legendre', functs=f, element='H', name='small et')
        self.assertIsInstance(bs, BasisSet)
        self.assertAlmostEqual(bs.functions['s']['e'], le_s)

    def test_init_legendre_spd(self):

        leg_s = [0.12245642825298195, 0.65376978512984729, 1.2214027581601699,
                 7.5761109446368486, 1480.2999275845475]
        leg_p = [1.8221188003905087, 4.7027533114762283, 5.2166313162210942,
                 330.29955990964862]
        leg_d = [1.8221188003905091, 1.5683121854901687, 36.59823444367801]
        f = [('s', 5, (1.0, 3.5, 1.6, 1.2)), ('p', 4, (2.0, 1.5, 1.2, 1.1)), ('d', 3, (1.0, 1.5, 1.1))]
        bs = BasisSet.from_sequence(formula='legendre', functs=f, element='H', name='small et')
        self.assertIsInstance(bs, BasisSet)
        self.assertAlmostEqual(bs.functions['s']['e'], leg_s)
        self.assertAlmostEqual(bs.functions['p']['e'], leg_p)
        self.assertAlmostEqual(bs.functions['d']['e'], leg_d)

    def test_init_mixed_formula_spd(self):

        et_s = [0.5, 1.0, 2.0, 4.0, 8.0]
        wt_p = [0.505859375, 1.1875, 3.8984375, 16.0]
        leg_d = [1.8221188003905091, 1.5683121854901687, 36.59823444367801]
        formulas = ['even', 'well', 'legendre']
        f = [('s', 5, (0.5, 2.0)), ('p', 4, (0.5, 2.0, 3.0, 4.0)), ('d', 3, (1.0, 1.5, 1.1))]
        bs = BasisSet.from_sequence(formula=formulas, functs=f, element='H', name='small et')
        self.assertIsInstance(bs, BasisSet)
        self.assertAlmostEqual(bs.functions['s']['e'], et_s)
        self.assertAlmostEqual(bs.functions['p']['e'], wt_p)
        self.assertAlmostEqual(bs.functions['d']['e'], leg_d)

if __name__ == "__main__":
    unittest.main()
