
from basisopt import driver
from molecule import Atom, Molecule
from molpro import Molpro

be = Molecule(name="Be", atoms=[Atom(at=4)], charge=0, multiplicity=1)

optimization = {"method"  : "Nelder-Mead",
                "lambda"  : 10.0,
                "tol"     : 1.0e-4,
                "options" : {"maxiter" : 500,
                             "disp"    : True,
                            }
               }

job = {"method"    : "hf",
       "objective" : "total energy",
       "core"      : [0,0,0,0,0,0,0,0],
       "inpname"   : "be_hf_10s5p.inp",
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

event  = {"typ"      : "eventemp",
          "name"     : "aug",
          "params"   : [(0.9, 2.03),(1.0, 3.0)],
          "nfpshell" : [10,4],
          "atomic"   : 4
         }

def main():


    res = driver(code=mp, job=job, mol=be, bsopt=event, opt=optimization)

    print res.success
    print res.x
    print res.fun

if __name__ == "__main__":
    main()
