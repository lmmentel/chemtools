'''Module for handling atoms and molecules.

Atomic properties are read from the elements.json file and set as attributes
Units
* atomic_radius : pm
* atomic_volume : cm3/mol
* boiling point K
* covalent radius : pm
* density : g/cm3
* evaporation heat : kJ/mol
* first ionization : kJ/mol
* fusion heat : kJ/mol
* ionic radius : pm
* lattice constant : ang
* melting point : K
* specific heat : J/ g mol @ 20 C
* thermal conductivity : W/m K @25 C

ionization energies
Kramida, A., Ralchenko, Yu., Reader, J., and NIST ASD Team (2014). NIST Atomic Spectra Database (ver. 5.2), [Online]. Available: http://physics.nist.gov/asd [2015, April 13]. National Institute of Standards and Technology, Gaithersburg, MD. 

'''

from math import sqrt
import numpy as np
import os
import json

_project_root = os.path.abspath(os.path.dirname(__file__))
_elements_json = os.path.join(_project_root, "elements.json")

ELEMENTS = json.loads(open(_elements_json, 'rb').read().decode('ascii', errors="ignore"), encoding='utf=8')

def element_from_name(name):
    '''Get the data on the element from the ELEMENTS dictionary specified by its name'''
    for k, v in ELEMENTS.items():
        if v["name"].lower() == name.lower():
            return k, v
    else:
        raise ValueError("name {0:s} not found".format(name))

def element_from_atomic_number(atomic_number):
    '''Get the data on the element from the ELEMENTS dictionary specified by its atomic number'''
    for k, v in ELEMENTS.items():
        if v["atomic_number"] == atomic_number:
            return k, v
    else:
        raise ValueError("element with atomic number {0:d} not found".format(atomic_number))

class Atom(object):

    '''Basic atom class representing an atom.'''


    def __init__(self, identifier, xyz=(0.0, 0.0, 0.0), dummy=False, id=None):

        self._dtxyz = np.dtype([('x', 'f8'), ('y', 'f8'), ('z', 'f8')])
        self.xyz = np.asarray(xyz, dtype=self._dtxyz)
        self.is_dummy = dummy
        self._set_attributes(identifier)

    def _set_attributes(self, identifier):
        '''
        Set the attributes of an atom based on the unique "indentifier" provided.
        The attributes are read from the elements.json file.
        '''

        if isinstance(identifier, str):
            if len(identifier) > 2:
                symbol, data = element_from_name(identifier)
                setattr(self, "symbol", symbol)
                for k, v in data.items():
                    setattr(self, k, v)
            else:
                setattr(self, "symbol", identifier)
                for k, v in ELEMENTS[identifier].items():
                    setattr(self, k, v)
        elif isinstance(identifier, int):
            symbol, data = element_from_atomic_number(identifier)
            setattr(self, "symbol", symbol)
            for k, v in data.items():
                setattr(self, k, v)
        else:
            raise ValueError("wrong element identifier: {0:s}, use symbol, name or atomic number".format(identifier))

        if self.is_dummy:
            self.set_atomic_number(0.0)

    def set_atomic_number(self, value):
        self.atomic_number = value

    def move(self, x=0.0, y=0.0, z=0.0):

        '''Move atom to a set of new coordinates given in xyz'''

        self.xyz = np.asarray([x, y, z], dtype=self._dtxyz)

    def gamess_rep(self):

        out = "{0:<10s} {1:5.1f}\t{2:15.5f}{3:15.5f}{4:15.5f}\n".format(
                self.symbol, float(self.atomic_number), self.xyz['x'], self.xyz['y'], self.xyz['z'])
        return out

    def __repr__(self):
        outs = "{0:<10s} {1:5.1f}\t{2:15.5f}{3:15.5f}{4:15.5f}".format(
                self.symbol, float(self.atomic_number), self.xyz['x'], self.xyz['y'], self.xyz['z'])
        return outs

class Molecule(object):

    '''Molecule class handling all the operations on single molecules'''

    def __init__(self, name="", atoms=None, unique=None, sym="", charge=0, multiplicity=1):
        self.name = name
        self.charge = charge
        self.multiplicity = multiplicity
        self.atoms = atoms
        if sym == "":
            self.symmetry = "dnh 2"
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
                    self._atoms.append(Atom(identifier=v[0], dummy=v[1]))
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

    def print_molecule(self):

        '''Print formatted molecule data.'''

        print("Molecule: {n:<10s} Charge: {c:<10d} Multiplicty: {m:<10d} Electrons: {e:<10d}\n".format(
              n=self.name, c=self.charge, m=self.multiplicity, e=self.nele()))
        print("{0:<10s} {1:14s}\t{2:^10s}{3:^10s}{4:^10s}\n".format("Element", "Nuclear Charge", "x", "y", "z"))
        for atom in self.atoms:
            print(atom)


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
                atom.symbol, float(atom.atomic_number), atom.xyz['x'], atom.xyz['y'], atom.xyz['z'])
        return out

    def molpro_rep(self):

        out = "geomtyp=xyz\ngeometry={{\n{0:4d}\n{1:<80s}\n".format(len(self.atoms), self.name)
        for atom in self.atoms:
            out = out + "{0:<10s}, {1:15.5f}, {2:15.5f}, {3:15.5f}\n".format(
                atom.symbol, atom.xyz['x'], atom.xyz['y'], atom.xyz['z'])
        out = out + "}\n"
        return out

