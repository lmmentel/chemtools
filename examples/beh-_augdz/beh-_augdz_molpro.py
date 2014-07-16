from chemtools.basisopt import driver
from chemtools.molecule import Atom, Molecule
from chemtools.molpro import Molpro
from chemtools import molpro
from chemtools.basisset import BasisSet, get_x0

inpstr = '''***,h2o test
memory,100,m                            !allocate 500 MW dynamic memory

geometry

basis

core

{rhf; wf,6,1,0}

'''

beh = Molecule(name="Be", atoms=[Atom(at=4),
                                 Atom(at=1, xyz=(0.0, 0.0, 2.724985))], sym="cnv 2", charge=-1, multiplicity=1)

optimization = {"method"  : "CG",
                "tol"     : 1.0e-4,
                "lambda"  : 10.0,
                "options" : {"maxiter" : 100,
                             "disp"    : True,
                            }
               }

job = {"method"    : "hf",
       "objective" : "total energy",
       "core"      : [1,0,0,0,0,0,0,0],
       "inpname"   : "be_cisd_augdz.inp",
       "inpdata"   : inpstr,
       "verbose"   : True,
       }

mp = Molpro(
        name="MOLPRO",
        execpath="/home/lmentel/Programs/MOLPRO/molprop_2012_1_Linux_x86_64_i8/bin/molpro",
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

s, H , 13.0100000, 1.9620000, 0.4446000, 0.1220000, 0.0297400
c, 1.3, 0.0196850, 0.1379770, 0.4781480
c, 4.4, 1
c, 5.5, 1
p, H , 0.7270000, 0.1410000
c, 1.1, 1
c, 2.2, 1
}
'''

auxbs = molpro.parse_basis(bsstr)

be_ccpvdz = BasisSet.from_dict({"element" : "Be",
            "functions" : auxbs["Be"]})
h_ccpvdz = BasisSet.from_dict({"element" : "H",
            "functions" : auxbs["H"]})

augfs = {"typ"      : "direxp",
         "name"     : "aug",
         "params"   : [(0.02,), (0.01,), (0.07,)],
         "nfpshell" : [1, 1, 1],
         "element"  : "Be",
         "type"     : "diffuse",
        }

def main():

    bsnopt = [be_ccpvdz, h_ccpvdz]

    res = driver(code=mp, job=job, mol=beh, bsnoopt=bsnopt, bsopt=augfs, opt=optimization)
    print res.success
    print res.x
    print res.fun

if __name__ == "__main__":
    main()
