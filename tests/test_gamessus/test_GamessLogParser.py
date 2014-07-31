# -*- coding: utf-8 -*-

import unittest
import os
from chemtools.gamessus import GamessLogParser

class TestGLPonHeMini(unittest.TestCase):

    def setUp(self):
        log = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data/he/he_mini_hf.log')
        self.glp = GamessLogParser(log)

    def test_logexists(self):
        self.assertTrue(self.glp.logexists())

    def test_accomplished(self):
        self.assertTrue(self.glp.accomplished())

    def test_get_version(self):
        self.assertEqual(self.glp.get_version(), "1 MAY 2013 (R1)")

    def test_get_charge(self):
        self.assertEqual(self.glp.get_charge(), 0)

    def test_get_electrons(self):
        self.assertEqual(self.glp.get_electrons(), 2)

    def test_get_homo(self):
        self.assertEqual(self.glp.get_homo(), 0)

    def test_get_number_of_atoms(self):
        self.assertEqual(self.glp.get_number_of_atoms(), 1)

    def test_get_number_of_aos(self):
        self.assertEqual(self.glp.get_number_of_aos(), 1)

    def test_get_number_of_mos(self):
        self.assertEqual(self.glp.get_number_of_mos(), 1)

    def test_get_number_of_core_mos(self):
        self.assertEqual(self.glp.get_number_of_core_mos(), None)

    def test_get_linear_deps(self):
        self.assertIs(self.glp.get_linear_deps(), None)

    def test_get_scf_type(self):
        self.assertEqual(self.glp.get_scf_type(), "RHF")

    def test_get_cc_type(self):
        self.assertEqual(self.glp.get_cc_type(), "NONE")

    def test_get_ci_type(self):
        self.assertEqual(self.glp.get_ci_type(), "NONE")

    def test_get_mplevel(self):
        self.assertEqual(self.glp.get_mplevel(), 0)

    def test_get_hf_total_energy(self):
        self.assertAlmostEqual(self.glp.get_hf_total_energy(), -2.8356798736)

    def test_get_energy_components_hf(self):
        d = {"WAVEFUNCTION NORMALIZATION" : 1.0000000000,
             "ONE ELECTRON ENERGY" : -3.8642167967,
             "TWO ELECTRON ENERGY" : 1.0285369231,
             "NUCLEAR REPULSION ENERGY" : 0.0000000000,
             "TOTAL ENERGY" : -2.8356798736,
             "ELECTRON-ELECTRON POTENTIAL ENERGY" : 1.0285369231,
             "NUCLEUS-ELECTRON POTENTIAL ENERGY" : -6.6999422562,
             "NUCLEUS-NUCLEUS POTENTIAL ENERGY" : 0.0000000000,
             "TOTAL POTENTIAL ENERGY" : -5.6714053331,
             "TOTAL KINETIC ENERGY" : 2.8357254595,
             "VIRIAL RATIO (V/T)" : 1.9999839244}
        self.assertDictEqual(d, self.glp.get_energy_components('hf'))

    def test_get_lz_values(self):
        lz = [{"index" : 0, "shell" : 1, "lz" : 0}]
        self.assertListEqual(self.glp.get_lz_values(), lz)

    def test_get_ao_labels(self):
        ao = [{"index" : 0, "symbol" : "HE", "center" : 1, "component" : "S"}]
        for dtest, dref in zip(self.glp.get_ao_labels(), ao):
            self.assertDictEqual(dtest, dref)

class TestGLPonHeH2(unittest.TestCase):

    def setUp(self):
        log = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data/heh2/he-h2_avdz_ormas.log')
        self.glp = GamessLogParser(log)

    def test_logexists(self):
        self.assertTrue(self.glp.logexists())

    def test_accomplished(self):
        self.assertTrue(self.glp.accomplished())

    def test_get_version(self):
        self.assertEqual(self.glp.get_version(), "1 MAY 2013 (R1)")

    def test_get_charge(self):
        self.assertEqual(self.glp.get_charge(), 0)

    def test_get_electrons(self):
        self.assertEqual(self.glp.get_electrons(), 4)

    def test_get_homo(self):
        self.assertEqual(self.glp.get_homo(), 1)

    def test_get_number_of_atoms(self):
        self.assertEqual(self.glp.get_number_of_atoms(), 3)

    def test_get_number_of_aos(self):
        self.assertEqual(self.glp.get_number_of_aos(), 27)

    def test_get_number_of_mos(self):
        self.assertEqual(self.glp.get_number_of_mos(), 27)

    def test_get_number_of_core_mos(self):
        self.assertEqual(self.glp.get_number_of_core_mos(), 0)

    def test_get_linear_deps(self):
        self.assertIs(self.glp.get_linear_deps(), None)

    def test_get_scf_type(self):
        self.assertEqual(self.glp.get_scf_type(), "RHF")

    def test_get_cc_type(self):
        self.assertEqual(self.glp.get_cc_type(), "NONE")

    def test_get_ci_type(self):
        self.assertEqual(self.glp.get_ci_type(), "ORMAS")

    def test_get_mplevel(self):
        self.assertEqual(self.glp.get_mplevel(), 0)

    def test_get_hf_total_energy(self):
        self.assertAlmostEqual(self.glp.get_hf_total_energy(), -3.9842266066)

    def test_get_energy_components_hf(self):
        d = {"WAVEFUNCTION NORMALIZATION" : 1.0000000000,
             "ONE ELECTRON ENERGY" : -7.6711052217,
             "TWO ELECTRON ENERGY" : 2.3347610066,
             "NUCLEAR REPULSION ENERGY" : 1.3521176085,
             "TOTAL ENERGY" : -3.9842266066,
             "ELECTRON-ELECTRON POTENTIAL ENERGY" : 2.3347610066,
             "NUCLEUS-ELECTRON POTENTIAL ENERGY" : -11.5902333599,
             "NUCLEUS-NUCLEUS POTENTIAL ENERGY" : 1.3521176085,
             "TOTAL POTENTIAL ENERGY" : -7.9033547447,
             "TOTAL KINETIC ENERGY" :3.9191281381,
             "VIRIAL RATIO (V/T)" : 2.0166104465}
        self.assertDictEqual(d, self.glp.get_energy_components('hf'))

    def test_get_energy_components_ormas(self):
        d = {"WAVEFUNCTION NORMALIZATION" : 1.0000000000,
             "ONE ELECTRON ENERGY" : -7.6031513994,
              "TWO ELECTRON ENERGY" : 2.1965614362,
              "NUCLEAR REPULSION ENERGY" : 1.3521176085,
              "TOTAL ENERGY" : -4.0544723546,
              "ELECTRON-ELECTRON POTENTIAL ENERGY" : 2.1965614362,
              "NUCLEUS-ELECTRON POTENTIAL ENERGY" : -11.5724441222,
              "NUCLEUS-NUCLEUS POTENTIAL ENERGY" : 1.3521176085,
              "TOTAL POTENTIAL ENERGY" : -8.0237650775,
              "TOTAL KINETIC ENERGY" : 3.9692927229,
              "VIRIAL RATIO (V/T)" : 2.0214596498}
        self.assertDictEqual(d, self.glp.get_energy_components('ci'))

class TestGLPonNeDZ(unittest.TestCase):

    def setUp(self):
        log = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data/ne/ne_dz_guga.log')
        self.glp = GamessLogParser(log)

    def test_logexists(self):
        self.assertTrue(self.glp.logexists())

    def test_accomplished(self):
        self.assertTrue(self.glp.accomplished())

    def test_get_version(self):
        self.assertEqual(self.glp.get_version(), "1 MAY 2012 (R1)")

    def test_get_charge(self):
        self.assertEqual(self.glp.get_charge(), 0)

    def test_get_electrons(self):
        self.assertEqual(self.glp.get_electrons(), 10)

    def test_get_homo(self):
        self.assertEqual(self.glp.get_homo(), 4)

    def test_get_number_of_atoms(self):
        self.assertEqual(self.glp.get_number_of_atoms(), 1)

    def test_get_number_of_aos(self):
        self.assertEqual(self.glp.get_number_of_aos(), 15)

    def test_get_number_of_mos(self):
        self.assertEqual(self.glp.get_number_of_mos(), 14)

    def test_get_number_of_core_mos(self):
        self.assertEqual(self.glp.get_number_of_core_mos(), 0)

    def test_get_linear_deps(self):
        self.assertIs(self.glp.get_linear_deps(), 0)

    def test_get_scf_type(self):
        self.assertEqual(self.glp.get_scf_type(), "RHF")

    def test_get_cc_type(self):
        self.assertEqual(self.glp.get_cc_type(), "NONE")

    def test_get_ci_type(self):
        self.assertEqual(self.glp.get_ci_type(), "GUGA")

    def test_get_mplevel(self):
        self.assertEqual(self.glp.get_mplevel(), 0)

    def test_get_hf_total_energy(self):
        self.assertAlmostEqual(self.glp.get_hf_total_energy(), -128.4887755517)

    def test_get_energy_components_hf(self):
        d = {"WAVEFUNCTION NORMALIZATION" : 1.0000000000,
             "ONE ELECTRON ENERGY" : -182.6161002898,
             "TWO ELECTRON ENERGY" : 54.1273247381,
             "NUCLEAR REPULSION ENERGY" : 0.0000000000,
             "TOTAL ENERGY" : -128.4887755517,
             "ELECTRON-ELECTRON POTENTIAL ENERGY" : 54.1273247381,
             "NUCLEUS-ELECTRON POTENTIAL ENERGY" : -311.1049627743,
             "NUCLEUS-NUCLEUS POTENTIAL ENERGY" : 0.0000000000,
             "TOTAL POTENTIAL ENERGY" : -256.9776380362,
             "TOTAL KINETIC ENERGY" : 128.4888624845,
             "VIRIAL RATIO (V/T)" : 1.9999993234}

        self.assertDictEqual(d, self.glp.get_energy_components('hf'))

    def test_get_energy_components_ormas(self):
        d = {"WAVEFUNCTION NORMALIZATION" : 1.0000000000,
             "ONE ELECTRON ENERGY" : -182.3536925117,
              "TWO ELECTRON ENERGY" : 53.6772930949,
              "NUCLEAR REPULSION ENERGY" : 0.0000000000,
              "TOTAL ENERGY" : -128.6763994168,
              "ELECTRON-ELECTRON POTENTIAL ENERGY" : 53.6772930949,
              "NUCLEUS-ELECTRON POTENTIAL ENERGY" : -310.8610615681,
              "NUCLEUS-NUCLEUS POTENTIAL ENERGY" : 0.0000000000,
              "TOTAL POTENTIAL ENERGY" : -257.1837684732,
              "TOTAL KINETIC ENERGY" : 128.5073690564,
              "VIRIAL RATIO (V/T)" : 2.0013153359}
        self.maxDiff = None
        self.assertDictEqual(d, self.glp.get_energy_components('ci'))

    def test_get_lz_values_ne(self):
        lz = [{'index': 0, 'lz': 0, 'shell': 1},
            {'index': 1, 'lz': 0, 'shell': 2},
            {'index': 2, 'lz': 0, 'shell': 3},
            {'index': 3, 'lz': 1, 'shell': 3},
            {'index': 4, 'lz': 1, 'shell': 3},
            {'index': 5, 'lz': 1, 'shell': 4},
            {'index': 6, 'lz': 1, 'shell': 4},
            {'index': 7, 'lz': 0, 'shell': 4},
            {'index': 8, 'lz': 0, 'shell': 5},
            {'index': 9, 'lz': 2, 'shell': 6},
            {'index': 10, 'lz': 1, 'shell': 6},
            {'index': 11, 'lz': 0, 'shell': 6},
            {'index': 12, 'lz': 1, 'shell': 6},
            {'index': 13, 'lz': 2, 'shell': 6}]
        self.assertListEqual(self.glp.get_lz_values(), lz)

if __name__ == "__main__":
    unittest.main()
