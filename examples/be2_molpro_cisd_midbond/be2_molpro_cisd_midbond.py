from chemtools.basisopt import driver
from chemtools.molecule import Atom, Molecule
from chemtools.molpro import Molpro
from chemtools import molpro
from chemtools.basisset import BasisSet, get_x0

# the dummy center will be on element 'H'

inpstr = '''***,h2o test
memory,100,m                            !allocate 500 MW dynamic memory

geometry

basis

dummy, H

core

{rhf; wf,8,1,0}
cisd

'''

be2 = Molecule(name="Be2", atoms=[Atom(at=4, xyz=(0.0, 0.0, -1.5)),
                                  Atom(at=1, xyz=(0.0, 0.0,  0.0)),
                                  Atom(at=4, xyz=(0.0, 0.0,  1.5))], sym="dnh 2", charge=0, multiplicity=1)

optimization = {"method"  : "Nelder-Mead",
                "tol"     : 1.0e-4,
                "lambda"  : 10.0,
                "options" : {"maxiter" : 100,
                             "disp"    : True,
                            }
               }

job = {"method"    : "cisd",
       "objective" : "total energy",
       "core"      : [2,0,0,0,0,0,0,0],
       "inpname"   : "be2_cisd_midbond.inp",
       "inpdata"   : inpstr,
       "verbose"   : True,
       }

mp = Molpro(
        name="MOLPRO",
        executable="/home/lmentel/Programs/MOLPRO/molprop_2012_1_Linux_x86_64_i8/bin/molpro",
        runopts=["-s", "-n", "1", "-d", "/home/lmentel/scratch"],
            )

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

auxbs = molpro.parse_basis(bsstr)

be_ccpvdz = BasisSet.from_dict({
            "name"      : "cc-pvdz",
            "element"   : "Be",
            "functions" : auxbs["Be"]})

mb = {"typ"      : "eventemp",
      "name"     : "mb44",
      "params"   : [(0.05, 2.0), (0.04, 2.0)],
      "nfpshell" : [4, 4],
      "element"  : "H",
      "type"     : "midbond",
     }

def main():

    res = driver(code=mp, job=job, mol=be2, bsnoopt=[be_ccpvdz], bsopt=mb, opt=optimization)
    print res.success
    print res.x
    print res.fun

if __name__ == "__main__":
    main()
