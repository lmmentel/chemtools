
from chemtools.basisopt import driver
from chemtools.molecule import Atom, Molecule
from chemtools.molpro import Molpro
from chemtools import molpro
from chemtools.basisset import BasisSet

be = Molecule(name="Be", atoms=[Atom(at=4)], charge=0, multiplicity=1)

optimization = {"method"  : "Nelder-Mead",
                "lambda"  : 10.0,
                "tol"     : 1.0e-4,
                "options" : {"maxiter" : 500,
                             "disp"    : True,
                            }
               }

job = {"method"    : "cisd",
       "objective" : "regexp",
       "regexp"    : r"\s+MP2 total energy:\s+(\-?\d+\.\d+)",
       "core"      : [1,0,0,0,0,0,0,0],
       "inpname"   : "be_cisd.inp",
       "inpdata"   : open("molpro.inp", 'r').read(),
       "verbose"   : True,
       }

ccpvdz = BasisSet.from_dict({"element"   : "Be", "name" : "cc-pvdz",
          "functions" : molpro.parse_basis(open("ccpvdz_molpro.bas", "r").read())["Be"],
         })

mp = Molpro(
        name="MOLPRO",
        executable="/home/lmentel/Programs/MOLPRO/molprop_2012_1_Linux_x86_64_i8/bin/molpro",
        runopts=["-s", "-n", "1", "-d", "/home/lmentel/scratch"],
            )

welltemp = {"typ"      : "welltemp",
            "name"     : "pol",
            "params"   : [(1.45, 4.717, 0.41, 0.5)],
            "nfpshell" : [0,5],
            "element"  : "Be",
           }


def main():

    res = driver(code=mp, job=job, mol=be, bsnoopt=ccpvdz, bsopt=welltemp, opt=optimization)

    print res.success
    print res.x
    print res.fun

if __name__ == "__main__":
    main()
