import numpy as np
import os

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from pescan.model import Base, Dimer, Trimer, Tetramer


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

def collect(session):
    '''
    Lopp over elements form the joblist, parse the output files and insert
    the parsed value into the database.
    '''

    for job in joblist:
        parser = GamessLogParser(job.output)
        job.hf_energy = parser.get_hf_total_energy()
        job.ci_energy = parser.get_energy_components('ci')['TOTAL ENERGY']
        session.add(job)
    session.commit()

def add_trimer(session, mol_name, basisset, ram, raa, gamma):
    '''
    Loop over all jobs and return a list of dicts with job info
    '''

    trimer = Trimer(molecule=mol_name,
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

def add_dimer(session, mol_name, basisset, raa):
    '''
    Loop over all jobs and return a list of dicts with job info

    Args:
    =====
    session:
        database session object (connection)
    molecule (str):
        molecule (system) name
    basissets (list):
        list of strings with basis set names
    grids (np.array):
        1-D numpy array with the values of the internuclear distances

    Returns:
    ========
    '''

    dimer = Dimer(molecule=mol_name,
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

    molecules = session.query(modelobj.molecule).distinct().all()
    for mol in molecules:
        os.mkdir(mol.molecule)
        os.chdir(mol.molecule)
        basis = session.query(modelobj.basisset).filter(modelobj.molecule == mol.molecule).distinct().all()
        for bas in basis:
            os.mkdir(bas.basisset)
            jobpath = os.path.join(os.getcwd(), bas.basisset)
            for entry in session.query(modelobj).filter(modelobj.molecule == mol.molecule).filter(modelobj.basisset == bas.basisset).all():
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

def connect_db(path):
    '''
    Connect to a database under the name stored in "path" and return the
    session.
    '''

    if not os.path.exists(path):
        raise OSError("database {0} does not exist".format(path))

    engine = create_engine("sqlite:///{path:s}".format(path=path), echo=False)
    db_session = sessionmaker(bind=engine)
    return db_session()
