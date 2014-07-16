
from chemtools.basisopt import driver
from chemtools.molecule import Atom, Molecule
from chemtools.molpro import Molpro
from chemtools import molpro
from chemtools.basisset import BasisSet, get_x0

be = Molecule(name="Be", atoms=[Atom(at=4)], charge=0, multiplicity=1)

optimization = {"method"  : "Nelder-Mead",
                "lambda"  : 10.0,
                "tol"     : 1.0e-4,
                "options" : {"maxiter" : 500,
                             "disp"    : True,
                            }
               }

job = {"method"    : "cisd",
       "objective" : "core energy",
       "core"      : [[1,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0]],
       "inpname"   : "be_cisd_coredz.inp",
       "inpdata"   : open("molpro.inp", 'r').read(),
       "verbose"   : True,
       }

mp = Molpro(
        name="MOLPRO",
        execpath="/home/lmentel/Programs/MOLPRO/molprop_2012_1_Linux_x86_64_i8/bin/molpro",
        runopts=["-s", "-n", "1", "-d", "/home/lmentel/scratch"],
            )

# a single element tuple must be followed by comma since otherwise its not
# treate as a tuple for example (1.0) should be (1.0,)

corefs = {"typ"      : "direct",
          "name"     : "aug",
          "params"   : [(1.1,), (4.2,)],
          "nfpshell" : [1, 1],
          "element"  : "Be",
          "type"     : "tight",
         }

ccpvdz = BasisSet.from_dict({"element"   : "Be",
          "functions" : molpro.parse_basis(open("ccpvdz_molpro.bas", "r").read())["Be"],
         })


def main():

    res = driver(code=mp, job=job, mol=be, bsnoopt=ccpvdz, bsopt=corefs, opt=optimization)
    print res.success
    print res.x
    print res.fun

if __name__ == "__main__":
    main()
