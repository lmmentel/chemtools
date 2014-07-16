from chemtools.basisopt import driver
from chemtools.molecule import Atom, Molecule
from chemtools.molpro import Molpro
from chemtools import molpro
from chemtools.basisset import BasisSet, get_x0

inpstr = '''***,h2o test
memory,100,m                            !allocate 500 MW dynamic memory
GTHRESH,THROVL=1.0e-9

geometry

basis

core

{rhf; wf,2,1,0}

'''

he = Molecule(name="He", atoms=[Atom(at=2)], charge=0, multiplicity=1)

optimization = {"method"  : "Nelder-Mead",
                "tol"     : 1.0e-4,
                "lambda"  : 10.0,
                "options" : {"maxiter" : 100,
                             "disp"    : True,
                            }
               }

job = {"method"    : "hf",
       "objective" : "total energy",
       "core"      : [0,0,0,0,0,0,0,0],
       "inpname"   : "he_hf.inp",
       "inpdata"   : inpstr,
       "verbose"   : True,
       }

mp = Molpro(
        name="MOLPRO",
        execpath="/home/lmentel/Programs/MOLPRO/molprop_2012_1_Linux_x86_64_i8/bin/molpro",
        runopts=["-s", "-n", "1", "-d", "/home/lmentel/scratch"],
            )

legendre = {"typ"      : "legendre",
         "name"     : "orbital",
         "params"   : [(3.55090649,  6.54677461,  1.02692142,  0.31719821)],
         "nfpshell" : [14,],
         "element"  : "He",
         "type"     : "hf",
        }

def main():

    # Energy from numberical HF solution:
    # total energy:   -2.861680026576101

    res = driver(code=mp, job=job, mol=he, bsopt=legendre, opt=optimization)
    print res.success
    print res.x
    print res.fun

if __name__ == "__main__":
    main()
