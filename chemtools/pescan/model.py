import os
import math
from collections import namedtuple

from sqlalchemy import Column, Integer, String, Float, create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.ext.hybrid import hybrid_property

xyz = namedtuple("XYZ", ['x', 'y', 'z'])

Base = declarative_base()

class Atom(Base):
    __tablename__ = 'atoms'

    id = Column(Integer, primary_key=True)
    name = Column(String)
    abspath = Column(String)
    output_name = Column(String)
    basisset = Column(String)
    hf_energy = Column(Float)
    ci_energy = Column(Float)
    cc_energy = Column(Float)
    status = Column(String)

    @hybrid_property
    def output(self):
        '''
        Return the absolute path of the output.
        '''
        return os.path.join(self.abspath, self.output_name)

    @hybrid_property
    def input_name(self):
        '''
        Compose the input name from the atom name and the basis set.
        '''
        return "{m:s}_{b:s}.inp".format(m=self.name, b=self.basisset)

    def get_xyz(self, x0=0.0, y0=0.0, z0=0.0):
        '''
        Calculate (x,y,z) coordiantes from internal (r_atom1_atom2).
        '''

        return [xyz(0.0, 0.0, 0.0)]

    def __repr__(self):
        return "<Atom(name={n}, basis={b})>".format(n=self.name, b=self.basisset)

class Dimer(Base):
    __tablename__ = 'dimers'

    id = Column(Integer, primary_key=True)
    name = Column(String)
    abspath = Column(String)
    output_name = Column(String)
    basisset = Column(String)
    hf_energy = Column(Float)
    ci_energy = Column(Float)
    cc_energy = Column(Float)
    status = Column(String)
    r_atom1_atom2 = Column(Float)

    @hybrid_property
    def output(self):
        '''
        Return the absolute path of the output.
        '''
        return os.path.join(self.abspath, self.output_name)

    @hybrid_property
    def input_name(self):
        '''
        Compose the input name from the molecule name, basis set and distance.
        '''
        return "{m:s}_{b:s}_r{r:.2f}.inp".format(m=self.name,
                                                 b=self.basisset,
                                                 r=self.r_atom1_atom2)

    def get_xyz(self, x0=0.0, y0=0.0, z0=0.0):
        '''
        Calculate (x,y,z) coordiantes from internal (r_atom1_atom2).
        '''

        return [xyz(0.0, 0.0, -self.r_atom1_atom2/2.0),
                xyz(0.0, 0.0, self.r_atom1_atom2/2.0)]

    def __repr__(self):
        return "<Dimer(name={m}, basis={b}, r_atom1_atom2={raa:5.2f)}>".format(
                m=self.name, b=self.basisset, raa=self.r_atom1_atom2)

class Trimer(Base):
    __tablename__ = 'trimers'

    id = Column(Integer, primary_key=True)
    name = Column(String)
    abspath = Column(String, nullable=True)
    output_name = Column(String)
    basisset = Column(String)
    hf_energy = Column(Float, nullable=True)
    ci_energy = Column(Float, nullable=True)
    cc_energy = Column(Float, nullable=True)
    status = Column(String)
    r_atom_mol = Column(Float)
    gamma = Column(Float)
    r_atom1_atom2 = Column(Float)

    @hybrid_property
    def output(self):
        '''
        Return the absolute path of the output.
        '''
        return os.path.join(self.abspath, self.output_name)

    @hybrid_property
    def input_name(self):
        '''
        Compose the input name from the molecule name, basis set and distance.
        '''
        return "{m:s}_{b:s}_R{R:.2f}_r{r:.2f}_g{g:.2f}.inp".format(m=self.name,
                b=self.basisset, R=self.r_atom_mol, r=self.r_atom1_atom2,
                g=self.gamma)

    def get_xyz(self, x0=0.0, y0=0.0, z0=0.0):
        '''
        Calculate (x,y,z) coordiantes from internal (R_heh2, r_hh, gamma).
        '''

        atom1 = xyz(x0, y0, z0)
        cm = xyz(atom1.x, atom1.y, atom1.z + self.r_atom_mol)
        atom2 = xyz(atom1.x + 0.5*self.r_atom1_atom2*math.sin(math.radians(self.gamma)),
                                 atom1.y,
                                 cm.z + 0.5*self.r_atom1_atom2*math.cos(math.radians(self.gamma)))
        atom3 = xyz(atom1.x - 0.5*self.r_atom1_atom2*math.sin(math.radians(self.gamma)),
                                 atom1.y,
                                 cm.z - 0.5*self.r_atom1_atom2*math.cos(math.radians(self.gamma)))
        if self.gamma not in [0.0, 90.0]:
            atom1 = xyz(atom1.z, atom1.x, atom1.y)
            atom2 = xyz(atom2.z, atom2.x, atom2.y)
            atom3 = xyz(atom3.z, atom3.x, atom3.y)
        return [atom1, atom2, atom3]

    def __repr__(self):
        return "<Trimer(name={m}, basis={b}, r_atom_mol={r:5.2f}, gamma={g:5.2f}, r_atom1_atom2={raa:5.2f})>".format(
                m=self.name, b=self.basisset, r=self.r_atom_mol, g=self.gamma, raa=self.r_atom1_atom2)

class Tetramer(Base):
    __tablename__ = 'tetramers'

    id = Column(Integer, primary_key=True)
    name = Column(String)
    abspath = Column(String)
    input_name = Column(String)
    output_name = Column(String)
    basisset = Column(String)
    hf_energy = Column(Float)
    ci_energy = Column(Float)
    cc_energy = Column(Float)
    status = Column(String)
    r_mol1_mol2 = Column(Float)
    gamma = Column(Float)
    phi_1 = Column(Float)
    phi_2 = Column(Float)
    r_mol1 = Column(Float)
    r_mol2 = Column(Float)

    def get_xyz(self, r_atom1_cm23, r_atom2_atom3, gamma, x0=0.0, y0=0.0, z0=0.0):
        '''
        Calculate (x,y,z) coordiantes from internal (R_heh2, r_hh, gamma).
        '''

        pass

    def __repr__(self):
        return "r_mol1_mol2={r:5.2f}, gamma={g:5.2f}, phi_1={p1:5.2f}, phi_2={p2:5.2f}>".format(
                r=self.r_mol1_mol2, g=self.gamma, p1=self.phi_1, p2=self.phi_2)

