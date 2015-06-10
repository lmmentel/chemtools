import unittest
from scipy.optimize import curve_fit

from chemtools import cbs

class TestPoly(unittest.TestCase):

    def setUp(self):
        # VO Be_2 VXZ corr energies at R_e
        self.x = [2, 3, 4, 5, 6]
        self.ens = [-0.1006245329, -0.1047828536, -0.1062511383, -0.1066997691, -0.1068704347]

    def tearDown(self):
        del self.x
        del self.ens

    def test_poly_2point_default(self):
        f = cbs.poly()
        eopt, ecov = curve_fit(f, self.x[-2:], self.ens[-2:])
        self.assertAlmostEqual(-0.1071048655, eopt[0])

    def test_poly_3point(self):
        f = cbs.poly(twopoint=False)
        eopt, ecov = curve_fit(f, self.x[-3:], self.ens[-3:])

class TestExpo(unittest.TestCase):

    def setUp(self):
        self.x = [2, 3, 4, 5, 6]
        self.ens = [-29.1322987785, -29.1337749163, -29.1340309055, -29.1341443009, -29.1341692835]

    def tearDown(self):
        del self.x
        del self.ens

    def test_expo_3point(self):
        eopt, ecov = curve_fit(cbs.expo, self.x[2:], self.ens[2:])
        self.assertAlmostEqual(-29.1341763428, eopt[0])

    def test_expo_4point(self):
        eopt, ecov = curve_fit(cbs.expo, self.x[1:], self.ens[1:])
        self.assertAlmostEqual(-29.1341981701, eopt[0])

    def test_expo_5point(self):
        eopt, ecov = curve_fit(cbs.expo, self.x, self.ens)
        self.assertAlmostEqual(-29.1341546294, eopt[0])

class TestExposqrt(unittest.TestCase):

    def setUp(self):
        self.x = [2, 3, 4, 5, 6]
        self.ens = [-29.1322987785, -29.1337749163, -29.1340309055, -29.1341443009, -29.1341692835]

    def tearDown(self):
        del self.x
        del self.ens

    def test_exposqrt_2point(self):
        f = cbs.exposqrt(twopoint=True)
        eopt, ecov = curve_fit(f, self.x[-2:], self.ens[-2:])
        self.assertAlmostEqual(-29.1341567932, eopt[0])

    def test_exposqrt_3point(self):
        f = cbs.exposqrt(twopoint=False)
        eopt, ecov = curve_fit(f, self.x[-3:], self.ens[-3:])
        self.assertAlmostEqual(-29.1341783198, eopt[0])

    def test_exposqrt_4point(self):
        f = cbs.exposqrt(twopoint=False)
        eopt, ecov = curve_fit(f, self.x[-4:], self.ens[-4:])
        self.assertAlmostEqual(-29.1342137839, eopt[0])

    def test_exposqrt_5point(self):
        f = cbs.exposqrt(twopoint=False)
        eopt, ecov = curve_fit(f, self.x[-5:], self.ens[-5:])
        self.assertAlmostEqual(-29.1341729017, eopt[0])


class TestExposum(unittest.TestCase):

    def setUp(self):
        self.x = [2, 3, 4, 5, 6]
        self.ens = [-29.1322987785, -29.1337749163, -29.1340309055, -29.1341443009, -29.1341692835]

    def tearDown(self):
        del self.x
        del self.ens


    def test_exposum_default(self):
        eopt, ecov = curve_fit(cbs.exposum, self.x[-3:], self.ens[-3:])
        self.assertAlmostEqual(-29.1341837985, eopt[0])

class TestExporapolate(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_poly_2point(self):

        pass
        #cbs.extrapolate(
        #self.assertAlmostEqual(-0.1068704347,)

if __name__ == "__main__":
    unittest.main()
