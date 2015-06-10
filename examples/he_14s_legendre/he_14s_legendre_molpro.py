from chemtools.basisopt import Job
from chemtools.molecule import Atom, Molecule
from chemtools.molpro import Molpro
from chemtools import molpro
from chemtools.basisset import BasisSet, get_x0

template = '''***,h2o test
memory,100,m                            !allocate 500 MW dynamic memory
GTHRESH,THROVL=1.0e-9

%geometry

%basis

%core

{rhf; wf,2,1,0}

'''

he = Molecule(name="He", atoms=[('He',)], charge=0, multiplicity=1)

nm = {"method"  : "Nelder-Mead",
      "tol"     : 1.0e-4,
      "lambda"  : 15.0,
      "options" : {"maxiter" : 100, "disp" : True},
     }

molpro = Molpro(
            name="MOLPRO",
            executable="/home/lmentel/Programs/molprop_2012_1_Linux_x86_64_i8/bin/molpro",
            runopts=["-s", "-n", "1", "-d", "/home/lmentel/scratch"],
            )

legendre = {"typ"      : "legendre",
         "name"     : "orbital",
         "params"   : [(3.55090649,  6.54677461,  1.02692142,  0.31719821)],
         "nfpshell" : [14,],
         "element"  : "He",
         "type"     : "hf",
        }

lege_sp = {"typ"     : "legendre",
          "name"     : "cc-pV6Z tight sp",
          "params"   : [(2.9623643, 1.39661345, 0.5), (3.04754371, 1.55118609, 0.2),],
          "nfpshell" : [5, 5,],
          "element"  : "He",
          "type"     : "tight",
         }

job = Job(calcmethod='hf', objective='total energy', core=[0,0,0,0,0,0,0,0],
          template=template, code=molpro, optalg=nm, mol=he, bsopt=legendre)

def main():

    # Energy from numerical HF solution:
    # total energy:   -2.861680026576101

    job.run()
    print(job.results)

if __name__ == "__main__":
    main()
