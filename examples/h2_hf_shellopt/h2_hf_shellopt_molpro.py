
from chemtools.basisopt import opt_shell_by_nf
from chemtools.molecule import Atom, Molecule
from chemtools import molpro
from chemtools.basisset import BasisSet

from scipy.constants import angstrom, physical_constants

def angs_to_bohr(value):
    return value*angstrom/physical_constants['Bohr radius'][0]

def bohr_to_angs(value):
    return value*physical_constants['Bohr radius'][0]/angstrom


r_a0 = 1.448736
r_angs = bohr_to_angs(r_a0)

h2 = Molecule(name="H2", atoms=[('H', (0.0, 0.0,  r_angs/2.0)),
                                ('H', (0.0, 0.0, -r_angs/2.0))],
              charge=0, multiplicity=1, unique=[0])

optimization = {"method"  : "Nelder-Mead",
                "lambda"  : 10.0,
                "tol"     : 1.0e-4,
                "options" : {"maxiter" : 500,
                             "disp"    : True,
                            }
               }

mpinpstr = '''***,h2
memory,100,m

%geometry

%basis

%core

gthresh,energy=1.0e-9
{rhf; wf,2,1,0}

'''

job = {"method"    : "hf",
       "objective" : "total energy",
       "core"      : [0,0,0,0,0,0,0,0],
       "inpname"   : "H2_hf_molpro.inp",
       "template"  : mpinpstr,
       "verbose"   : True,
       }

mp = molpro.Molpro(
        executable="/home/lmentel/Programs/molprop_2012_1_Linux_x86_64_i8/bin/molpro",
        runopts=["-s", "-n", "1", "-d", "."],
            )

def main():

    OPT_TOL = 1.0e-6
    MAX_PARAMS = 2

    bsopt = {"typ"      : "eventemp",
             "name"     : "hf",
             "kind"     : "cisd",
             "params"   : [( 0.33, )],
             "element"  : "H",
    }

    bsnoopt = BasisSet.from_dict(eval(open("optimized__legendre_2_11s.bas", "r").read()))
    bsp = BasisSet.from_dict(eval(open("optimized__legendre_2_5p.bas", "r").read()))
    bsd = BasisSet.from_dict(eval(open("optimized__legendre_2_2d.bas", "r").read()))
    bsf = BasisSet.from_dict(eval(open("optimized__legendre_2_1f.bas", "r").read()))

    bsnoopt.add(bsp)
    bsnoopt.add(bsd)
    bsnoopt.add(bsf)

    opt_shell_by_nf(shell='g',
                    nfs=range(1, 20),
                    max_params=MAX_PARAMS,
                    opt_tol=OPT_TOL,
                    save=True,
                    code=mp,
                    job=job,
                    mol=h2,
                    bsnoopt=bsnoopt,
                    bsopt=bsopt,
                    opt=optimization)

if __name__ == "__main__":
    main()
