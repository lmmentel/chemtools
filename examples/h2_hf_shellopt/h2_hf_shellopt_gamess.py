
from chemtools.basisopt import opt_shell_by_nf
from chemtools.molecule import Atom, Molecule
from chemtools import gamessus
from chemtools.basisset import BasisSet

h2 = Molecule(name="H2", atoms=[('H', (0.0, 0.0,  0.724368)),
                                ('H', (0.0, 0.0, -0.724368))],
              charge=0, multiplicity=1, unique=[0])

optimization = {"method"  : "Nelder-Mead",
                "lambda"  : 10.0,
                "tol"     : 1.0e-4,
                "options" : {"maxiter" : 500,
                             "disp"    : True,
                            }
               }

template = {'$contrl': {'icut': '30', 'ispher': '1', 'itol': '30',
                        'maxit': '30', 'mult': '1', 'qmttol': '1.0d-8',
                        'runtyp': 'energy', 'scftyp': 'rhf', 'units': 'bohr'},
            '$scf': {'conv': '1.0d-9'},
            '$system': {'mwords': '100', 'timlim': '525600'},
            '$trans': {'cuttrf': '1.0d-14'},
            '$data' : {"title" : "h2", "group" : "dnh 2" }
}


job = {"method"    : "hf",
       "objective" : "total energy",
       "core"      : "",
       "inpname"   : "H2_hf_gamess.inp",
       "template"  : template,
       "verbose"   : True,
       }

gamess = gamessus.Gamess(
            executable="/home/lmentel/Programs/gamess-us-dec2014/rungms",
            scratch="/home/lmentel/scratch",
            runopts = {"remove_dat" : True, "nproc" : 2})

def main():

    OPT_TOL = 1.0e-6
    MAX_PARAMS = 2

    bsopt = {"typ"      : "eventemp",
             "name"     : "hf",
             "kind"     : "cisd",
             "params"   : [( 0.33, 2.13),],
             "element"  : "H",
    }

    opt_shell_by_nf(shell='s',
                    nfs=range(4, 20),
                    max_params=MAX_PARAMS,
                    opt_tol=OPT_TOL,
                    save=True,
                    code=gamess,
                    job=job,
                    mol=h2,
                    bsnoopt=None,
                    bsopt=bsopt,
                    opt=optimization)

if __name__ == "__main__":
    main()
