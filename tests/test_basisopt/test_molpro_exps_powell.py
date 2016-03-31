
from chemtools.basisopt import BSOptimizer
from chemtools.molecule import Molecule
from chemtools.calculators.molpro import Molpro
from chemtools.basisset import BasisSet
import pytest
import os

@pytest.mark.skipif(os.getenv('MOLPRO_EXE') is None, reason="<MOLPRO_EXE> undefined")
def test_molpro_exps_opt(tmpdir):

    tmpdir.chdir()

    template = '''***,he
memory,100,m
GTHRESH,THROVL=1.0e-9

%geometry

%basis

%core

{rhf; wf,2,1,0}

'''

    he = Molecule(name="He", atoms=[('He',)], charge=0, multiplicity=1)

    optimization = {"method"  : "Powell",
        "tol"     : 1.0e-5,
        "lambda"  : 10.0,
        "jacob"   : False,
        "options" : {"maxiter" : 100, "disp" : True},
    }

    mp = Molpro(exevar="MOLPRO_EXE", runopts=["-s", "-n", "1", "-d", "."])

    exps = (1.00473157e+03,   3.51102141e+02,   1.35742439e+02,
            5.69386800e+01,   2.51416924e+01,   1.15836670e+01,
            5.23629402e+00,   2.38685192e+00,   7.60836609e-01)

    sfuncts = {'He' : [('s', 'exp', 9, exps)]}

    bso = BSOptimizer(objective='hf total energy', template=template, code=mp, mol=he,
                  fsopt=sfuncts, staticbs=None, core=[0,0,0,0,0,0,0,0], fname='molpro_powell.inp',
                  verbose=True, uselogs=True, optalg=optimization)


    bso.run()

    energy = -2.8614232306
    assert abs(bso.result.fun - energy) < 1.0e-8, 'wrong objective'
