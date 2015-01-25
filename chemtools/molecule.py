'''Module for handling atoms and molecules.'''

from collections import namedtuple
from math import sqrt

class Atom(object):

    '''Basic atom class representing an atom.'''

    coords = namedtuple("XYZ", ['x', 'y', 'z'])

    _element = ['X',  'H',  'He',
                'Li', 'Be', 'B',  'C',  'N',  'O',  'F',  'Ne',
                'Na', 'Mg', 'Al', 'Si', 'P',  'S',  'Cl', 'Ar',
                'K',  'Ca', 'Sc', 'Ti', 'V',  'Cr', 'Mn', 'Fe',
                'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se',
                'Br', 'Kr', 'Rb', 'Sr', 'Y',  'Zr', 'Nb', 'Mo',
                'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
                'Sb', 'Te', 'I',  'Xe', 'Cs', 'Ba', ]

    def __init__(self, at, xyz=(0.0, 0.0, 0.0), id=None):
        self.atomic = at
        self.symbol = Atom._element[abs(self.atomic)]
        self.xyz    = Atom.coords(xyz[0], xyz[1], xyz[2])

    def move(self, x=0.0, y=0.0, z=0.0):

        '''Move atom to a set of new coordinates given in xyz'''

        self.xyz = Atom.coords(x, y, z)

    def gamess_rep(self):

        out = "{0:<10s} {1:5.1f}\t{2:15.5f}{3:15.5f}{4:15.5f}\n".format(
              self.symbol, float(self.atomic), self.xyz.x, self.xyz.y, self.xyz.z)
        return out

    def __repr__(self):
        outs = "{0:<10s} {1:5.1f}\t{2:15.5f}{3:15.5f}{4:15.5f}".format(
               self.symbol, float(self.atomic), self.xyz.x, self.xyz.y, self.xyz.z)
        return outs

class Molecule(object):

    '''Molecule class handling all the operations on single molecules'''

    def __init__(self, name="", atoms=[], unique=[], sym="", charge=0, multiplicity=1):
        self.name         = name
        self.charge       = charge
        self.multiplicity = multiplicity
        self.atoms        = atoms
        if sym == "":
            self.symmetry     = "dnh 2"
        else:
            self.symmetry     = sym
        self.electrons    = self.nele()
        #if atomList:
        #    self.set_atoms(atomList)
        if unique:
            self.unique  = unique
        else:
            self.unique  = range(len(atoms))

    def nele(self):

        '''Get the total number of electrons in a molecule.'''

        nelectrons = 0
        for atom in self.atoms:
            if atom.atomic > 0:
                nelectrons += atom.atomic
        return nelectrons - self.charge

    def set_atoms(self, atom_list):

        '''Restructure raw molecule data into molecule made of atoms.'''

        for i, (atomic, xyz) in enumerate(atom_list):
            id = i + 1
            self.atoms.append(Atom(atomic, xyz, id))

    def print_molecule(self):

        '''Print formatted molecule data.'''

        print "Molecule: {n:<10s} Charge: {c:<10d} Multiplicty: {m:<10d} Electrons: {e:<10d}\n".format(
              n=self.name, c=self.charge, m=self.multiplicity, e=self.nele())
        print "{0:<10s} {1:14s}\t{2:^10s}{3:^10s}{4:^10s}\n".format("Element", "Nuclear Charge", "x", "y", "z")
        for atom in self.atoms:
            print atom


    def get_distance(self, atom1, atom2):

        '''Calcualte the distance between two atoms.'''

        dist = 0.0
        for i in range(3):
            dist += (self.atoms[atom1].xyz[i]-self.atoms[atom2].xyz[i])**2
        return sqrt(dist)


    def gamess_rep(self):

        out = ""
        for atom in self.atoms:
            out = out + "{0:<10s} {1:5.1f}\t{2:15.5f}{3:15.5f}{4:15.5f}\n".format(
                atom.symbol, float(atom.atomic), atom.xyz.x, atom.xyz.y, atom.xyz.z)
        return out

    def molpro_rep(self):

        out = "geomtyp=xyz\ngeometry={{\n{0:4d}\n{1:<80s}\n".format(len(self.atoms), self.name)
        for atom in self.atoms:
            out = out + "{0:<10s}, {1:15.5f}, {2:15.5f}, {3:15.5f}\n".format(
                atom.symbol, atom.xyz.x, atom.xyz.y, atom.xyz.z)
        out = out + "}\n"
        return out

