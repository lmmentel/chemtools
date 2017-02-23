
import unittest
from chemtools import molecule


class TestAtomAN(unittest.TestCase):

    def setUp(self):
        self.h = molecule.Atom(1)

    def tearDown(self):
        del self.h

    def test_Atom_default(self):

        self.assertEqual('H', self.h.symbol)
        self.assertEqual('Hydrogen', self.h.name)
        self.assertEqual(1, self.h.atomic_number)
        self.assertEqual(1.008, self.h.mass)


class TestAtomSymbol(unittest.TestCase):

    def setUp(self):
        self.o = molecule.Atom('O')

    def tearDown(self):
        del self.o

    def test_Atom_default(self):

        self.assertEqual('O', self.o.symbol)
        self.assertEqual('Oxygen', self.o.name)
        self.assertEqual(8, self.o.atomic_number)
        self.assertEqual(15.999, self.o.mass)


if __name__ == "__main__":
    unittest.main()
