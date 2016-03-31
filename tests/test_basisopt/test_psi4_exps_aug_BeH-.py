from chemtools.basisopt import BSOptimizer
from chemtools.molecule import Molecule
from chemtools.calculators.psi4 import Psi4
from chemtools.basisset import BasisSet
import pytest
import os

@pytest.mark.skipif(os.getenv('PSI4_EXE') is None, reason="<PSI4_EXE> undefined")
def test_optimize(tmpdir):

    tmpdir.chdir()

    psi = Psi4(exevar="PSI4_EXE", runopts=["-n", "4"])

    beh = Molecule(name="BeH-", atoms=[('Be',), ('H', (0.0, 0.0, 2.724985))],
                   sym="cnv 2", charge=-1, multiplicity=1)


    psitemp = '''#! test BeH-

memory 100 mb

molecule beh {
    Be  0.0000  0.0000  0.0000
    H   0.0000  0.0000  1.4420
    units angstrom
}

beh.set_molecular_charge(-1)
beh.set_multiplicity(1)


basis = {

assign mybasis

[ mybasis ]
spherical

%basisset

}

set {
freeze_core True
e_convergence 13
d_convergence 10
}

energy('cisd')
'''

    optimization = {"method"  : "Nelder-Mead",
                    "tol"     : 1.0e-4,
                    "options" : {"maxiter" : 100,
                                 "disp"    : True,
                                }
                    }

    ccpvdz_str = '''basis={
s,Be,2.940000E+03,4.412000E+02,1.005000E+02,2.843000E+01,9.169000E+00,3.196000E+00,1.159000E+00,1.811000E-01,5.890000E-02;
c,1.9,6.800000E-04,5.236000E-03,2.660600E-02,9.999300E-02,2.697020E-01,4.514690E-01,2.950740E-01,1.258700E-02,-3.756000E-03;
c,1.9,-1.230000E-04,-9.660000E-04,-4.831000E-03,-1.931400E-02,-5.328000E-02,-1.207230E-01,-1.334350E-01,5.307670E-01,5.801170E-01;
c,9.9,1.000000E+00;

p,Be,3.619000E+00,7.110000E-01,1.951000E-01,6.018000E-02;
c,1.4,2.911100E-02,1.693650E-01,5.134580E-01,4.793380E-01;
c,4.4,1.000000E+00;

d,Be,2.354000E-01;
c,1.1,1.000000E+00;

s, H , 13.0100000, 1.9620000, 0.4446000, 0.1220000, 0.0297400
c, 1.3, 0.0196850, 0.1379770, 0.4781480
c, 4.4, 1
c, 5.5, 1
p, H , 0.7270000, 0.1410000
c, 1.1, 1
c, 2.2, 1
}
'''

    ccpvdz = BasisSet.from_str(ccpvdz_str, fmt='molpro')

    augfs = {'Be' : [('s', 'exp', 1, (0.02,)),
                     ('p', 'exp', 1, (0.01,)),
                     ('d', 'exp', 1, (0.07,))]}

    bso = BSOptimizer(objective='cisd total energy', template=psitemp, code=psi, mol=beh,
                      fsopt=augfs, staticbs=ccpvdz, core=[1,0,0,0,0,0,0,0], verbose=True,
                      optalg=optimization, fname='beh_aug.dat')

    bso.run()

    energy = -15.204232302468
    assert abs(bso.result.fun - energy) < 1.0e-8, 'wrong objective'
