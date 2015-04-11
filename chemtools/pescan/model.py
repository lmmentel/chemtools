import os
import math
from collections import namedtuple
import numpy as np
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
        atom2 = xyz(cm.x + 0.5*self.r_atom1_atom2*math.sin(math.radians(self.gamma)),
                    cm.y,
                    cm.z + 0.5*self.r_atom1_atom2*math.cos(math.radians(self.gamma)))
        atom3 = xyz(cm.x - 0.5*self.r_atom1_atom2*math.sin(math.radians(self.gamma)),
                    cm.y,
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

    def get_xyz(self, x0=0.0, y0=0.0, z0=0.0):
        '''
        Calculate (x,y,z) coordiantes from internal (r_mol1_mol2, r_mol1, r_mol2, phi_1, phi_2, gamma).
        '''

        o = xyz(x0, y0, z0)
        cm1 = xyz(o.x, o.y, o.z - self.r_mol1_mol2/2.0)
        cm2 = xyz(o.x, o.y, o.z + self.r_mol1_mol2/2.0)

        atom1 = xyz(cm1.x + 0.5*self.r_mol1*math.sin(math.radians(self.phi_1)),
                    cm1.y,
                    cm1.z + 0.5*self.r_mol1*math.cos(math.radians(self.phi_1)))
        atom2 = xyz(cm1.x - 0.5*self.r_mol1*math.sin(math.radians(self.phi_1)),
                    cm1.y,
                    cm1.z - 0.5*self.r_mol1*math.cos(math.radians(self.phi_1)))
        atom3 = xyz(cm2.x + 0.5*self.r_mol2*math.sin(math.radians(self.phi_2)),
                    cm2.y,
                    cm2.z + 0.5*self.r_mol2*math.cos(math.radians(self.phi_2)))
        atom4 = xyz(cm2.x - 0.5*self.r_mol2*math.sin(math.radians(self.phi_2)),
                    cm2.y,
                    cm2.z - 0.5*self.r_mol2*math.cos(math.radians(self.phi_2)))

        # get the rotation matrix around the z axis
        Rmatrix = rotation_matrix(math.radians(self.gamma), [0, 0, 1])

        newxyz = np.dot(Rmatrix, np.array([atom3.x, atom3.y, atom3.z]))
        atom3 = xyz(newxyz[0], newxyz[1], newxyz[2])
        newxyz = np.dot(Rmatrix, np.array([atom4.x, atom4.y, atom4.z]))
        atom4 = xyz(newxyz[0], newxyz[1], newxyz[2])

        return [atom1, atom2, atom3, atom4]

    def __repr__(self):
        return "<Tetramer(name={n:15s}, basis={b:15s}, r_mol1_mol2={r:5.2f}, gamma={g:5.2f}, phi_1={p1:5.2f}, phi_2={p2:5.2f}>".format(
               n=self.name, b=self.basisset, r=self.r_mol1_mol2, g=self.gamma, p1=self.phi_1, p2=self.phi_2)

def unit_vector(data, axis=None, out=None):
    """Return ndarray normalized by length, i.e. Euclidean norm, along axis.

    >>> v0 = numpy.random.random(3)
    >>> v1 = unit_vector(v0)
    >>> numpy.allclose(v1, v0 / numpy.linalg.norm(v0))
    True
    >>> v0 = numpy.random.rand(5, 4, 3)
    >>> v1 = unit_vector(v0, axis=-1)
    >>> v2 = v0 / numpy.expand_dims(numpy.sqrt(numpy.sum(v0*v0, axis=2)), 2)
    >>> numpy.allclose(v1, v2)
    True
    >>> v1 = unit_vector(v0, axis=1)
    >>> v2 = v0 / numpy.expand_dims(numpy.sqrt(numpy.sum(v0*v0, axis=1)), 1)
    >>> numpy.allclose(v1, v2)
    True
    >>> v1 = numpy.empty((5, 4, 3))
    >>> unit_vector(v0, axis=1, out=v1)
    >>> numpy.allclose(v1, v2)
    True
    >>> list(unit_vector([]))
    []
    >>> list(unit_vector([1]))
    [1.0]

    """
    if out is None:
        data = np.array(data, dtype=np.float64, copy=True)
        if data.ndim == 1:
            data /= math.sqrt(np.dot(data, data))
            return data
    else:
        if out is not data:
            out[:] = np.array(data, copy=False)
        data = out
    length = np.atleast_1d(np.sum(data*data, axis))
    np.sqrt(length, length)
    if axis is not None:
        length = np.expand_dims(length, axis)
    data /= length
    if out is None:
        return data

def rotation_matrix(angle, direction, point=None):
    """Return matrix to rotate about axis defined by point and direction.

    >>> R = rotation_matrix(math.pi/2, [0, 0, 1], [1, 0, 0])
    >>> numpy.allclose(numpy.dot(R, [0, 0, 0, 1]), [1, -1, 0, 1])
    True
    >>> angle = (random.random() - 0.5) * (2*math.pi)
    >>> direc = numpy.random.random(3) - 0.5
    >>> point = numpy.random.random(3) - 0.5
    >>> R0 = rotation_matrix(angle, direc, point)
    >>> R1 = rotation_matrix(angle-2*math.pi, direc, point)
    >>> is_same_transform(R0, R1)
    True
    >>> R0 = rotation_matrix(angle, direc, point)
    >>> R1 = rotation_matrix(-angle, -direc, point)
    >>> is_same_transform(R0, R1)
    True
    >>> I = numpy.identity(4, numpy.float64)
    >>> numpy.allclose(I, rotation_matrix(math.pi*2, direc))
    True
    >>> numpy.allclose(2, numpy.trace(rotation_matrix(math.pi/2,
    ...                                               direc, point)))
    True

    """
    sina = math.sin(angle)
    cosa = math.cos(angle)
    direction = unit_vector(direction[:3])
    # rotation matrix around unit vector
    R = np.diag([cosa, cosa, cosa])
    R += np.outer(direction, direction) * (1.0 - cosa)
    direction *= sina
    R += np.array([[ 0.0,         -direction[2],  direction[1]],
                      [ direction[2], 0.0,          -direction[0]],
                      [-direction[1], direction[0],  0.0]])
    return R

