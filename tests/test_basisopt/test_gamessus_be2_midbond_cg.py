
import os
import pytest
from chemtools.basisopt import BSOptimizer
from chemtools.molecule import Molecule
from chemtools.calculators.gamessus import GamessUS
from chemtools.basisset import BasisSet


@pytest.mark.skipif(os.getenv('GAMESS_EXE') is None,
                    reason="<GAMESS_EXE> undefined")
def test_optimize(tmpdir):

    tmpdir.chdir()

    gamesstemp = ''' $CONTRL
    scftyp=rhf
    cityp=guga
    runtyp=energy
    maxit=30
    mult=1
    ispher=1
    units=bohr
    exetyp=run
 $END
 $SYSTEM
    timlim=525600
    mwords=100
 $END
 $CIDRT
   iexcit=2
   nfzc=2
   ndoc=2
   nval=40
   group=d2h
 $END

%basis
'''

    be2X = Molecule(name="Be2", atoms=[('Be', (0.0, 0.0, -1.5)),
                                       ('H', (0.0, 0.0, 0.0), True),
                                       ('Be', (0.0, 0.0, 1.5))],
                    sym="dnh 2", charge=0, multiplicity=1, unique=[0, 1])

    bsstr = '''basis={
s,Be,2.940000E+03,4.412000E+02,1.005000E+02,2.843000E+01,9.169000E+00,3.196000E+00,1.159000E+00,1.811000E-01,5.890000E-02;
c,1.9,6.800000E-04,5.236000E-03,2.660600E-02,9.999300E-02,2.697020E-01,4.514690E-01,2.950740E-01,1.258700E-02,-3.756000E-03;
c,1.9,-1.230000E-04,-9.660000E-04,-4.831000E-03,-1.931400E-02,-5.328000E-02,-1.207230E-01,-1.334350E-01,5.307670E-01,5.801170E-01;
c,9.9,1.000000E+00;

p,Be,3.619000E+00,7.110000E-01,1.951000E-01,6.018000E-02;
c,1.4,2.911100E-02,1.693650E-01,5.134580E-01,4.793380E-01;
c,4.4,1.000000E+00;

d,Be,2.354000E-01;
c,1.1,1.000000E+00;
}
'''

    vdz = BasisSet.from_str(bsstr, fmt='molpro')

    optimization = {"method": "CG",
                    "tol": 1.0e-4,
                    "jacob": False,
                    "options": {"maxiter": 100, "disp": True, 'eps': 1.0e-3}}

    gams = GamessUS(exevar="GAMESS_EXE", runopts=['1'])

    midbond = {'H': [('s', 'et', 4, (0.05, 2.0)),
                     ('p', 'et', 4, (0.04, 2.0))]}

    bso = BSOptimizer(objective='cisd total energy', template=gamesstemp,
                      code=gams, mol=be2X, fsopt=midbond, staticbs={'Be': vdz},
                      core=2, fname='midbond.inp', verbose=True,
                      optalg=optimization)

    bso.run()
    energy = -29.1326831928
    assert abs(bso.result.fun - energy) < 1.0e-8, 'wrong objective'
