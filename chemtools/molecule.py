# -*- coding: utf-8 -*-

#The MIT License (MIT)
#
#Copyright (c) 2014 Lukasz Mentel
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.

'''
Module for handling atoms and molecules.
'''

from math import sqrt
import numpy as np
from mendeleev import element


class Atom(object):

    '''Basic atom class representing an atom.'''

    def __init__(self, identifier, xyz=(0.0, 0.0, 0.0), dummy=False, id=None):

        self.xyz = np.asarray(xyz)
        self.is_dummy = dummy
        self._set_attributes(identifier)

    @property
    def xyz(self):
        return self._xyz

    @xyz.setter
    def xyz(self, values):

        if len(values) != 3:
            raise ValueError("Expecting 3 coordinates (x, y, z), got: {0:d}".format(len(values)))
        else:
            self._xyz = np.asarray(values)

    def _set_attributes(self, identifier):
        '''
        Set the attributes of an atom based on the unique "indentifier" provided.
        The attributes are read from the elements.json file.
        '''

        attrs = ["name", "symbol", "atomic_number", "mass"]
        atom = element(identifier)

        for attr in attrs:
            setattr(self, attr, getattr(atom, attr))

        if self.is_dummy:
            self.set_atomic_number(0.0)

    def set_atomic_number(self, value):
        self.atomic_number = value

    def move(self, x=0.0, y=0.0, z=0.0):

        '''Move atom to a set of new coordinates given in xyz'''

        self.xyz = np.asarray([x, y, z], dtype=self._dtxyz)

    def gamess_rep(self):

        out = "{0:<10s} {1:5.1f}\t{2:15.5f}{3:15.5f}{4:15.5f}\n".format(
                self.symbol, float(self.atomic_number), self.xyz[0], self.xyz[1], self.xyz[2])
        return out

    def __repr__(self):
        out = "{0:<10s} {1:5.1f}\t{2:15.5f}{3:15.5f}{4:15.5f}".format(
                self.symbol, float(self.atomic_number), self.xyz[0], self.xyz[1], self.xyz[2])
        return out

    def __str__(self):
        out = "{0:<10s} {1:14.2f}\t{2:15.5f}{3:15.5f}{4:15.5f}\n".format(
                self.symbol, float(self.atomic_number), self.xyz[0], self.xyz[1], self.xyz[2])
        return out


class Molecule(object):

    '''Molecule class handling all the operations on single molecules'''

    def __init__(self, name="", atoms=None, unique=None, sym="", charge=0, multiplicity=1):
        self.name = name
        self.charge = charge
        self.multiplicity = multiplicity
        self.atoms = atoms
        if sym == "":
            self.symmetry = "c1"
        else:
            self.symmetry = sym
        self.electrons = self.nele()
        if unique is not None:
            self.unique_labels = unique
        else:
            self.unique_labels = range(len(atoms))

    @property
    def atoms(self):
        return self._atoms

    @atoms.setter
    def atoms(self, values):
        self._atoms = []
        if hasattr(values, "__iter__"):
            for v in values:
                if len(v) == 1:
                    self._atoms.append(Atom(identifier=v[0]))
                elif len(v) == 2:
                    if isinstance(v[1], bool):
                        self._atoms.append(Atom(identifier=v[0], dummy=v[1]))
                    elif isinstance(v[1], (list, tuple)):
                        self._atoms.append(Atom(identifier=v[0], xyz=v[1]))
                    else:
                        raise ValueError("Second argument should be <bool> or <tuple>, not {}".format(type(v[1])))
                elif len(v) == 3:
                    self._atoms.append(Atom(identifier=v[0], xyz=v[1], dummy=v[2]))
                else:
                    raise ValueError("wrong number of Atom arguments: {}, expecting: 1, 2, 3.".format(len(v)))
        else:
            raise TypeError("{0:s} object is not iterable".format(type(values)))

    def unique(self):
        '''Get a list of unique atom specified by unique keyword'''
        return [self.atoms[i] for i in self.unique_labels]

    def nele(self):

        '''Get the total number of electrons in a molecule.'''

        nelectrons = 0
        for atom in self.atoms:
            if atom.atomic_number > 0:
                nelectrons += atom.atomic_number
        return nelectrons - self.charge


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
                atom.symbol, float(atom.atomic_number), atom.xyz[0], atom.xyz[1], atom.xyz[2])
        return out

    def molpro_rep(self):

        out = "geomtyp=xyz\ngeometry={{\n{0:4d}\n{1:<80s}\n".format(len(self.atoms), self.name)
        for atom in self.atoms:
            out = out + "{0:<10s}, {1:15.5f}, {2:15.5f}, {3:15.5f}\n".format(
                atom.symbol, atom.xyz[0], atom.xyz[1], atom.xyz[2])
        out = out + "}\n"
        return out

    def __repr__(self):

        out = "<Molecule(name={}, charge={}, multiplicity={},\n".format(self.name, self.charge, self.multiplicity)
        for atom in self.atoms:
            out += atom.__repr__() + "\n"
        out += ")>"
        return out

    def __str__(self):

        '''Print formatted molecule data.'''

        out = 'Name: {n:<10s} Charge: {c:<10d} Multiplicty: {m:<10d} Electrons: {e:<10d}\n'.format(
              n=self.name, c=self.charge, m=self.multiplicity, e=self.nele())
        out += 'Atoms:\n'
        out += '{0:<10s} {1:^14s}\t{2:^15s}{3:^15s}{4:^15s}\n'.format("Element", "Nuclear Charge", "x", "y", "z")
        for atom in self.atoms:
            out += str(atom)
        return out
