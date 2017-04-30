# -*- coding: utf-8 -*-

import numpy as np
import os

from itertools import product

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from chemtools.pescan.model import Base, Atom, Dimer, Trimer, Tetramer


def grid_product(dicts):
    '''
    Given a dict whose item are lists return a generator over dict with
    all possible selections from the list (cartesian product)

    >>> options = {"number": [1,2,3], "color": ["orange","blue"] }
    >>> print list( my_product(options) )
    [ {"number": 1, "color": "orange"},
      {"number": 1, "color": "blue"},
      {"number": 2, "color": "orange"},
      {"number": 2, "color": "blue"},
      {"number": 3, "color": "orange"},
      {"number": 3, "color": "blue"}
    ]
    '''
    return (dict(zip(dicts, x)) for x in product(*dicts.itervalues()))


def expand_grids(grids=[], debug=False):

    '''Expand the grid definition and calculate all the grid points.'''

    res = list()
    for grid in grids:
        array = np.empty(0)
        for start, stop, num in grid:
            array = np.concatenate((array, np.linspace(start, stop, num)), axis=0)
        res.append(array)
    return res

# controller methods


def add_atom(session, name, basisset):

    atom = Atom(name=name,
                abspath=None,
                output_name=None,
                basisset=basisset,
                hf_energy=None,
                ci_energy=None,
                cc_energy=None,
                status=None)
    session.add(atom)
    session.commit()
    return True


def add_trimer(session, name, basisset, ram, raa, gamma):
    '''
    Add a Trimer to the database.
    '''

    trimer = Trimer(
        name=name,
        abspath=None,
        output_name=None,
        basisset=basisset,
        hf_energy=None,
        ci_energy=None,
        cc_energy=None,
        status=None,
        r_atom_mol=ram,
        gamma=gamma,
        r_atom1_atom2=raa)
    session.add(trimer)
    session.commit()
    return True


def add_tetramer(session, name, basisset, r_mol1_mol2, r_mol1, r_mol2, phi_1,
                 phi_2, gamma):
    '''
    Add a Tetramer to the database.
    '''

    tetramer = Tetramer(
        name=name,
        basisset=basisset,
        r_mol1_mol2=r_mol1_mol2,
        r_mol1=r_mol1,
        r_mol2=r_mol2,
        gamma=gamma,
        phi_1=phi_1,
        phi_2=phi_2)
    session.add(tetramer)
    session.commit()
    return True


def add_dimer(session, name, basisset, raa):
    '''
    Loop over all jobs and return a list of dicts with job info

    Args:
      session:
        database session object (connection)
      name : str
        name of the system
      basissets : list
        list of strings with basis set names
      grids : numpy.array
        1-D numpy array with the values of the internuclear distances

    '''

    dimer = Dimer(name=name,
                  abspath=None,
                  output_name=None,
                  basisset=basisset,
                  hf_energy=None,
                  ci_energy=None,
                  cc_energy=None,
                  status=None,
                  r_atom1_atom2=raa)
    session.add(dimer)
    session.commit()
    return True


def create_dirs(session, modelobj=None, workdir=os.getcwd()):
    '''
    Create the directory tree associated with a given table.
    '''

    # if workdir doesn't exist create it
    # raise exception when modelobj is None

    os.chdir(workdir)

    molecules = session.query(modelobj.name).distinct().all()
    for mol in molecules:
        os.mkdir(mol.name)
        os.chdir(mol.name)
        basis = session.query(modelobj.basisset).filter(modelobj.name == mol.name).distinct().all()
        for bas in basis:
            os.mkdir(bas.basisset)
            jobpath = os.path.join(os.getcwd(), bas.basisset)
            for entry in session.query(modelobj).filter(modelobj.name == mol.name).filter(modelobj.basisset == bas.basisset).all():
                entry.abspath = jobpath
                session.add(entry)
        os.chdir("..")
    session.commit()

# database connection methods

def new_db(path):
    '''
    Create a new database under path and return the session.
    '''
    engine = create_engine("sqlite:///{path:s}".format(path=path), echo=False)
    Base.metadata.create_all(engine)
    db_session = sessionmaker(bind=engine)
    return db_session()

def get_session(dbpath):
    '''
    Connect to a database under the name stored in "path" and return the
    session.
    '''

    if not os.path.exists(dbpath):
        raise OSError("database {0} does not exist".format(path))

    engine = create_engine("sqlite:///{path:s}".format(path=dbpath), echo=False)
    db_session = sessionmaker(bind=engine)
    return db_session()
