
from chemtools.basisopt import BSOptimizer
from chemtools.molecule import Molecule
from chemtools.calculators.molpro import Molpro
from chemtools.basisset import BasisSet

def test_optimize(tmpdir):

    tmpdir.chdir()

    be = Molecule(name="Be", atoms=[('Be',)], charge=0, multiplicity=1)

    template = '''***,be ccpvdz tight fs test
memory,100,m

%geometry

%basis

%core

{rhf; wf,4,1,0}
cisd
'''

    mp = Molpro(executable="/home/lmentel/Programs/molprop_2012_1_Linux_x86_64_i8/bin/molpro",
                runopts=["-s", "-n", "1", "-d", "/home/lmentel/scratch"])


    corefs = {'Be' :  [('s', 'exp', 1, (1.1,)), ('p', 'exp', 1, (4.2,))]}

    ccpvdz_str = '''basis={
!BAZA Petersona
s,Be,2.940000E+03,4.412000E+02,1.005000E+02,2.843000E+01,9.169000E+00,3.196000E+00,1.159000E+00,1.811000E-01,5.890000E-02;
c,1.9,6.800000E-04,5.236000E-03,2.660600E-02,9.999300E-02,2.697020E-01,4.514690E-01,2.950740E-01,1.258700E-02,-3.756000E-03;
c,1.9,-1.230000E-04,-9.660000E-04,-4.831000E-03,-1.931400E-02,-5.328000E-02,-1.207230E-01,-1.334350E-01,5.307670E-01,5.801170E-01;
c,9.9,1.000000E+00;

p,Be,3.619000E+00,7.110000E-01,1.951000E-01,6.018000E-02;
c,1.4,2.911100E-02,1.693650E-01,5.134580E-01,4.793380E-01;
c,4.4,1.000000E+00;
!
!d,Be,2.354000E-01;
!c,1.1,1.000000E+00;
}
'''
    ccpvdz = BasisSet.from_str(ccpvdz_str, fmt='molpro')

    bso = BSOptimizer(objective='cisd total energy', core=[[1,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0]], template=template,
                    verbose=True, code=mp, mol=be, fsopt=corefs, staticbs=ccpvdz, uselogs=True, runcore=True,
                    fname='be_tight.inp')

    bso.run()
    coreenergy = -0.031510764610001019
    assert abs(bso.result.fun - coreenergy) < 1.0e-8, 'wrong objective'
