
from chemtools.molecule import Atom, Molecule
from chemtools import molpro
from chemtools.basisset import BasisSet
from chemtools.basisopt import opt_shell_by_nf, opt_multishell


be2 = Molecule(name="Be2", atoms=[("Be", (0.0, 0.0,  1.2268016705), False),
                                  ("Be", (0.0, 0.0, -1.2268016705), False)], charge=0, multiplicity=1)

be = Molecule(name="Be", atoms=[("Be", (0.0, 0.0, 0.0), False)], charge=0, multiplicity=1)

optimization = {"method"  : "Nelder-Mead",
                "lambda"  : 10.0,
                "tol"     : 1.0e-3,
                "options" : {"maxiter" : 500,
                             "disp"    : True,
                            }
               }

job = {"method"    : "hf",
       "objective" : "total energy",
       "core"      : [0,0,0,0,0,0,0,0],
       "inpname"   : "be_hf.inp",
       "inpdata"   : open("molpro.inp", 'r').read(),
       "verbose"   : True,
       }

mp = molpro.Molpro(
        name="MOLPRO",
        execpath="/share/apps/molprop_2012_1_Linux_x86_64_i8/bin/molpro",
        runopts=["-s", "-n", "1", "-d", "."],
            )

def test_optimize_s():

    nfs = range(2, 20)
    bsopt = {"typ"      : "legendre",
             "name"     : "nowa baza",
             "kind"     : "hf",
             "params"   : [(1.37, 2.83, 0.4),],
             "element"  : "Be",
    }
    opt_shell_by_nf(shell='s',
                    nfs=nfs,
                    max_params=5,
                    opt_tol=1.0e-2,
                    save=True,
                    code=mp,
                    job=job,
                    mol=be,
                    bsnoopt=None,
                    bsopt=bsopt,
                    opt=optimization)

def test_optimize_multishell():

    shells = ['s', 'p', 'd']
    nfps = [range(4, 15), range(1, 10), range(1, 5)]
    guesses = [[1.37, 2.83, -0.4, 0.3], [2.5], [0.5]]

    bsopt = {"typ"      : "legendre",
             "name"     : "nowa baza",
             "kind"     : "hf",
             "element"  : "Be",
    }

    #"params"   : [(1.37, 2.83, -0.4, 0.3),],
    opt_multishell(shells=shells, nfps=nfps, guesses=guesses, max_params=5, opt_tol=1.0e-3,
                   save=True, bsopt=bsopt,
                   code=mp,
                   job=job,
                   mol=be2,
                   bsnoopt=None,
                   opt=optimization)
