import os
import unittest

from chemtools.calculators.gamessus import GamessInput 

class TestGIPonHeH2(unittest.TestCase):

    def setUp(self):
        inpfile = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data/heh2/he-h2_avdz_ormas.inp')
        self.gip = GamessInput(fname=inpfile)
        self.gip.parse()

    def test_inpdict(self):
        self.assertIsInstance(self.gip.parsed, dict)

    def test_contrl(self):
        d ={"scftyp" : "rhf",
            "cityp" : "ormas",
            "runtyp" : "energy",
            "maxit" : "30",
            "mult" : "1",
            "icut" : "30",
            "itol" : "30",
            "ispher" : "1",
            "units" : "bohr"}
        self.assertDictEqual(d, self.gip.parsed["$contrl"])

    def test_trans(self):
        d = {"cuttrf" : "1.0d-10"}
        self.assertDictEqual(d, self.gip.parsed["$trans"])

    def test_system(self):
        d = {"timlim" : "525600",
             "mwords" : "100"}
        self.assertDictEqual(d, self.gip.parsed["$system"])

    def test_title(self):
        self.assertEqual(self.gip.parsed["$data"]["title"], "He-H2 FCI")

    def test_group(self):
        self.assertEqual(self.gip.parsed["$data"]["group"], "cnv 2")

    def test_atoms(self):
        self.assertEqual(len(self.gip.parsed["$data"]["atoms"]), 2)

    def  test_atom1_xyz(self):
        t = tuple([0.0, 0.0, -3.0])
        self.assertTupleEqual(t, self.gip.parsed["$data"]["atoms"][0]["xyz"])

    def  test_atom1_atomic(self):
        self.assertEqual(self.gip.parsed["$data"]["atoms"][0]["atomic"], 2.0)

    def  test_atom1_label(self):
        self.assertEqual(self.gip.parsed["$data"]["atoms"][0]["label"], "He")

    def  test_atom2_xyz(self):
        t = tuple([0.0, 0.724368, 3.0])
        self.assertTupleEqual(t, self.gip.parsed["$data"]["atoms"][1]["xyz"])

    def  test_atom2_atomic(self):
        self.assertEqual(self.gip.parsed["$data"]["atoms"][1]["atomic"], 1.0)

    def  test_atom2_label(self):
        self.assertEqual(self.gip.parsed["$data"]["atoms"][1]["label"], "H")

class TestGIPonNe(unittest.TestCase):
    ne_no_inp = ''' $CONTRL
    scftyp=none
    RUNTYP=ENERGY
    CITYP=GUGA
    MAXIT=200
    MULT=1
    ISPHER=1
    ICUT=20
    ITOL=20
    QMTTOL=1.E-8
    ICHARG=0
 $END
 $SYSTEM TIMLIM=14400 MWORDS=100 MEMDDI=0 $END
 $SCF
    DIRSCF=.F. FDIFF=.F. VVOS=.F.
 $END
 $CIDRT
    iexcit=2
    nfzc=1
    ndoc=4
    nval=9
    group=d2h
    STSYM=AG
 $END
 $GUGDIA
    CVGTOL=1E-10
 $END
 $BASIS GBASIS=CCD $END
 $DATA
Title
DNH 2

NE 10.0       0.00000000       0.00000000      0.0000000
 $END
 $GUESS guess=moread norb=14 $END

 $VEC
 1  1 1.00040943E+00 2.94355142E-07-2.09957395E-03 0.00000000E+00 0.00000000E+00
 1  2 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00
 1  3 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00
 2  1 0.00000000E+00 0.00000000E+00 0.00000000E+00 6.96472997E-01 0.00000000E+00
 2  2 0.00000000E+00 4.55660763E-01 0.00000000E+00 0.00000000E+00 0.00000000E+00
 2  3 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00
 3  1-5.82488255E-01-1.26490410E+00 1.62236917E+00 0.00000000E+00 0.00000000E+00
 3  2 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00
 3  3 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00
 4  1 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00
 4  2 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00
 4  3 0.00000000E+00 0.00000000E+00 0.00000000E+00 1.00000000E+00 0.00000000E+00
 5  1 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 5.80423950E-01
 5  2 0.00000000E+00 0.00000000E+00 5.80423950E-01 0.00000000E+00 0.00000000E+00
 5  3 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00
 6  1-0.00000000E+00-0.00000000E+00-0.00000000E+00-0.00000000E+00-9.84522233E-01
 6  2-0.00000000E+00-0.00000000E+00 9.84522233E-01-0.00000000E+00-0.00000000E+00
 6  3-0.00000000E+00-0.00000000E+00-0.00000000E+00-0.00000000E+00-0.00000000E+00
 7  1 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00
 7  2 5.80423950E-01 0.00000000E+00 0.00000000E+00 5.80423950E-01 0.00000000E+00
 7  3 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00
 8  1-0.00000000E+00-0.00000000E+00-0.00000000E+00-0.00000000E+00-0.00000000E+00
 8  2-9.84522233E-01-0.00000000E+00-0.00000000E+00 9.84522233E-01-0.00000000E+00
 8  3-0.00000000E+00-0.00000000E+00-0.00000000E+00-0.00000000E+00-0.00000000E+00
 9  1 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00
 9  2 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00
 9  3 0.00000000E+00 0.00000000E+00 1.00000000E+00 0.00000000E+00 0.00000000E+00
10  1 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00
10  2 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00
10  3 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 1.00000000E+00
11  1-0.00000000E+00-0.00000000E+00-0.00000000E+00-0.00000000E+00-0.00000000E+00
11  2-0.00000000E+00-0.00000000E+00-0.00000000E+00-0.00000000E+00 8.66025404E-01
11  3-8.66025404E-01-0.00000000E+00-0.00000000E+00-0.00000000E+00-0.00000000E+00
12  1 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00
12  2 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00-5.00000000E-01
12  3-5.00000000E-01 1.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00
13  1-0.00000000E+00-0.00000000E+00-0.00000000E+00-9.06146430E-01-0.00000000E+00
13  2-0.00000000E+00 1.04811701E+00-0.00000000E+00-0.00000000E+00-0.00000000E+00
13  3-0.00000000E+00-0.00000000E+00-0.00000000E+00-0.00000000E+00-0.00000000E+00
14  1 2.16942447E-01 1.02369108E+00-1.09192955E-03 0.00000000E+00 0.00000000E+00
14  2 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00
14  3 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00
 $END'''

    def setUp(self):
        self.gip = GamessInput()
        self.gip.parse_from_string(TestGIPonNe.ne_no_inp)

    def test_inpdict(self):
        self.assertIsInstance(self.gip.parsed, dict)

    def test_vec(self):
        self.assertIsInstance(self.gip.parsed["$vec"], str)

if __name__ == "__main__":
    unittest.main()
