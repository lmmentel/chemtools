
from chemtools.molecule import Molecule
from chemtools.calculators import molpro
from chemtools.basisopt import opt_shell_by_nf, opt_multishell
import pytest


@pytest.fixture
def be2():
    return Molecule(name="Be2", atoms=[("Be", (0.0, 0.0, 1.2268016705), False),
                                       ("Be", (0.0, 0.0, -1.2268016705), False)],
                    charge=0, multiplicity=1)


@pytest.fixture
def be_atom():
    return Molecule(name="Be", atoms=[("Be", (0.0, 0.0, 0.0), False)],
                    charge=0, multiplicity=1)


@pytest.fixture
def optimization():
    return {
        "method": "Nelder-Mead",
        "lambda": 10.0,
        "tol": 1.0e-3,
        "options": {"maxiter": 500,
                    "disp": True}
    }


@pytest.fixture
def calculator():
    return molpro.Molpro(exevar="MOLPRO_EXE",
                         runopts=["-s", "-n", "4", "-d", "."])


@pytest.mark.skip(reason="no way of currently testing this")
def test_optimize_s(calculator, optimization, be):

    job = {"method": "hf",
           "objective": "total energy",
           "core": [0, 0, 0, 0, 0, 0, 0, 0],
           "inpname": "be_hf.inp",
           "inpdata": '',
           "verbose": True,
           }

    nfs = range(2, 20)
    bsopt = {"typ": "legendre",
             "name": "nowa baza",
             "kind": "hf",
             "params": [(1.37, 2.83, 0.4)],
             "element": "Be"}

    opt_shell_by_nf(shell='s',
                    nfs=nfs,
                    max_params=5,
                    opt_tol=1.0e-2,
                    save=True,
                    code=calculator,
                    job=job,
                    mol=be,
                    bsnoopt=None,
                    bsopt=bsopt,
                    opt=optimization)


@pytest.mark.skip(reason="no way of currently testing this")
def test_optimize_multishell(calculator, optimization, be2):

    shells = ['s', 'p', 'd']
    nfps = [range(4, 15), range(1, 10), range(1, 5)]
    guesses = [[1.37, 2.83, -0.4, 0.3], [2.5], [0.5]]

    bsopt = {"typ": "legendre",
             "name": "nowa baza",
             "kind": "hf",
             "element": "Be"}

    job = {"method": "hf",
           "objective": "total energy",
           "core": [0, 0, 0, 0, 0, 0, 0, 0],
           "inpname": "be_hf.inp",
           "inpdata": '',
           "verbose": True,
           }

    # "params"   : [(1.37, 2.83, -0.4, 0.3),],
    opt_multishell(shells=shells, nfps=nfps, guesses=guesses, max_params=5,
                   opt_tol=1.0e-3,
                   save=True, bsopt=bsopt,
                   code=calculator,
                   job=job,
                   mol=be2,
                   bsnoopt=None,
                   opt=optimization)
