#!/usr/bin/env python

# lmentel's packages
from gamessus import Gamess
from molecule import Atom, Molecule
from molpro import Molpro
import basisset as bas

# standard python packages
from scipy.optimize import minimize
from subprocess import Popen, PIPE
import datetime
import numpy as np
import os
import re
import sys
import time

_shells = ["S", "P", "D", "F", "G", "H", "I"]

def write_gamess_input(fname=None, core=None, bs=None, code=None, mol=None):

    '''write gamess input stored in the dictionary <dinput> to the file
    <fname> '''

    dinput = parse_gamess_input(code["input"])
    dinput = set_gamess_input(dinput, mol, bs, code, core)

    end = "$end"

    if not isinstance(dinput, dict):
        raise TypeError("expected a dictionary but got {0:s}".format(type(dinput)))
    inp = open(fname, 'w')

    for key, value in sorted(dinput.items()):
        inp.write(" {g:<s}\n".format(g=key))
        for fkey, fvalue in sorted(value.items()):
            inp.write("    {k:s}={v:s}\n".format(k=fkey, v=str(fvalue)))
        inp.write(" {g:<s}\n".format(g=end))
    # write basis set and atom/molecule
    inp.write(" $data\n")
    inp.write("title\n")
    inp.write("{0:s}\n\n".format(mol.symmetry))
    for atom in [a for a in mol.atoms if mol.atoms.index(a) in mol.unique]:
        inp.write(atom.gamess_rep())
        for function in [f for f in bs if f["atomic"] == atom.atomic]:
            inp.write(" {s:<1s} {n:>3d}\n".format(s=_shells[function["shell"]], n=len(function["exps"])))
            for i, prim in enumerate(zip(function["exps"], function["coeffs"])):
                inp.write("{0:>3d}{z:>20.10e}{c:>20.10f}\n".format(i+1, z=prim[0], c=prim[1]))
        inp.write("\n")
    inp.write(" $end")
    inp.close()

def set_gamess_input(dinp, mol, bs, code, core):

    if "$contrl" in dinp.keys():
        dinp["$contrl"]["icharg"] = mol.charge
        dinp["$contrl"]["mult"]   = mol.multiplicity
        if mol.multiplicity == 1:
            dinp["$contrl"]["scftyp"] = "rhf"
        elif mol.multiplicity > 1:
            dinp["$contrl"]["scftyp"] = "rohf"
        if code["method"].lower() == "hf":
            dinp["$contrl"] = {key:value for key, value in dinp["$contrl"].items() if key != "cityp"}
    else:
        sys.exit("no $contrl group in the gamess input string")
    if "$cidet" in dinp.keys():
        dinp["$cidet"]["nact"]  = bas.get_no_functions(bs) - core
        dinp["$cidet"]["ncore"] = core
        dinp["$cidet"]["nels"]  = mol.electrons - core*2
    if "$ormas" not in dinp.keys():
        dinp["$ormas"] = {}
    if code["method"] == "cisd":
        # sum of singly and doubly occupied  orbitals
        noccupied = mol.multiplicity - 1 + (mol.electrons-(mol.multiplicity - 1))/2
        dinp["$ormas"]["nspace"] = 2
        dinp["$ormas"]["mstart(1)"] = ','.join([str(core+1),str(noccupied+1)])
        dinp["$ormas"]["mine(1)"] = ",".join([str((mol.electrons-core*2)-2), str(0)])
        dinp["$ormas"]["maxe(1)"] = ",".join([str(mol.electrons-core*2), str(2)])
    elif code["method"] == "mrcisd":
        dinp["$ormas"]["nspace"] = 3
        dinp["$ormas"]["mstart(1)"] = ','.join([str(core+1),str(mol.electrons-core+1)])
        dinp["$ormas"]["mine(1)"] = ",".join([str(mol.electrons-2), str(0)])
        dinp["$ormas"]["maxe(1)"] = ",".join([str(mol.electrons), str(2)])

    return dinp

def parse_gamess_input(contents):
    '''
    Parse gamess input file into a dictionary of dictionaries, where the
    highest level entries are gamess namelist fileds and that contain
    dictionaries of options. All key are converted to lowercase. For example if
    the following input was parsed:

    " $CONTRL scftyp=rhf units=bohr
          runtyp=energy   cityp=ormas $END
      $SYSTEM TIMLIM=1000 mwords=500 $END
      ..."

    the follwoing dictionary will be produced:

    {"$contrl" : {"scftyp" : "rhf",
                  "units"  : "bohr",
                  "runtyp" : "energy",
                  "cityp"  : "ormas"},
     "$system" : {"timlim" : "1000",
                  "mwords" : "500"},
    ...
    }
    '''

    dontparse = ["$data", "$vec", ]

    pat = re.compile(r'(?P<group>\$[a-zA-Z]{3,6})\s+(?P<entries>.*?)\$END', flags=re.S)

    dinput = {}

    iterator = pat.finditer(contents)

    for match in iterator:
        if match.group("group").lower() not in dontparse:
            dinput[match.group("group").lower()] = {}
            fields = [s.strip() for s in match.group("entries").split("\n")]
            for field in fields:
                if not field.startswith("!"):
                    for line in field.split():
                        key, value = line.split("=")
                        dinput[match.group("group").lower()][key.lower()] = value
        elif match.group("group").lower() == "$data":
            bs = parse_basis(match.group("entries").split("\n"))
    if "$data" in [match.group("group").lower() for match in iterator]:
        return dinput, bs
    else:
        return dinput


def run_gamess(fname, code=None):

    '''Run a single gamess job interactively - without submitting to the queue.'''

    datfile = os.path.splitext(fname)[0]+".dat"
    if os.path.exists(os.path.join(code["opts"]["scratch"],datfile)):
        os.remove(os.path.join(code["opts"]["scratch"], datfile))

    process = Popen([code["opts"]["exec"], fname, code["opts"]["version"], code["opts"]["nproc"]], stdout=PIPE, stderr=PIPE)
    out, err = process.communicate()
    output_file = os.path.splitext(fname)[0]+"_{0:d}.log".format(process.pid)
    fout = open(output_file, 'w')
    fout.write(out)
    fout.write("{0:s}\n{1:^80s}\n{0:s}\n".format("="*80, "Error messages:"))
    fout.write(err)
    fout.close()
    return output_file


def driver(code=None, job=None, mol=None, bsnoopt=None, bsopt=None, opt=None):
    '''
    Driver for the basis set optimization
    '''

    print "Script started at {0}".format(datetime.datetime.today())

    starttime = time.time()

    # set default optimization options if opts not given
    if not opt:
        opt = {"method"  : "BFGS",
                "lambda"  : 10.0,
                "tol"     : 1.0e-4,
                "options" : {"maxiter" : 50,
                             "disp"    : True,
                             "eps"     : 0.01
                            }
               }
    else:
        if opt["method"].lower() == "nelder-mead":
            jacob = None
        else:
            jacob = False
        if "lambda" not in opt.keys():
            opt["lambda"] = 10.0

    if not bsnoopt:
        bsnooptobj = []

    if job:
        if job["objective"].lower() == "core energy":
            if len(job["core"]) == 2:
                if job["core"][0] > job["core"][1]:
                    job["core"].reverse()
            else:
                sys.exit("<driver>: core should have 2 entries, but it has {0:4d}".format(len(job["core"])))
            function = run_core_energy
        elif job["objective"] in ["total energy", "correlation energy", "regexp"]:
            function = run_total_energy
        else:
            sys.exit("<driver>: wrong objective in code dictionary")
        if isinstance(code, Gamess):
            if bsnoopt and len(bsnoopt) > 0:
                bsnooptobj = bas.parse_gamess_basis(bsnoopt)
        elif isinstance(code, Molpro):
            if bsnoopt and len(bsnoopt) > 0:
                bsnooptobj = bas.parse_molpro_basis(bsnoopt)
        else:
            sys.exit("<driver>: unknown code object given")
    else:
        sys.exit("<driver>: no dictionary describing code given")

    if not bsopt:
        sys.exit("<driver>: no dictionary describing basis set to be optimized given")

    if not mol:
        sys.exit("<driver>: no molecule object specified")

    x0 = bas.get_x0(bsopt)

    res = minimize(function, x0,
            args=(bsopt, bsnooptobj, code, job, mol, opt,),
            method=opt["method"],
            jac=jacob,
            tol=opt["tol"],
            options=opt["options"])
    print res
    print "Elapsed time : {0:>20.3f} sec".format(time.time()-starttime)
    # write a nice printer of the optimized exponents if converged
    return res

def run_total_energy(x0, *args):
    '''
    Funtion for running a single point calculation and parsing the
    resulting energy (or property) as specified by the objective
    function.

    Input
        x0: list/numpy 1d array
            contains a list of parameters to be optimized, may be
            explicit exponents or parametrized exponents in terms
            of some polynomial
        args: tuple of dictionaries
            bsopt, bsnoopt, code, job, mol, opt, needed for writing
            input and parsing output
    Output
        parsed result of the single point calculation as speficied by the
        objective function in the "job" dictionary
    '''

    # unpack the args tuple for code readability
    bsopt, bsnoopt, code, job, mol, opt = args

    if bsopt["typ"] in ["direxp", "direct", "exps", "exponents", "event", "eventemp"]:
        #penalty = sum(min(0, x)**2 for x in x0)
        penalty = 0.0
        if any(x < 0 for x in x0):
            x0 = [abs(x) for x in x0]
    else:
        penalty = 0.0

    bs2opt = bas.get_basis(x0, bsopt)

    code.write_input(job["inpname"], job["core"], bs=bsnoopt+bs2opt, mol=mol, inpdata=job["inpdata"])
    code.run_single(job["inpname"])
    if code.isok():
        objective = code.parse(job["method"], job["objective"], job.get("regexp", None))
        if job["verbose"]:
            print "{0:<s}".format("Job Terminated without errors")
            print "Current exponents"
            bas.print_functions(bs2opt)
            print "x0 : ", ", ".join([str(x) for x in x0])
            print "\n{0:<20s} : {1:>30s}".format("Output", code.outfile)
            print "{0:<20s} : {1:>30.10f}".format("Obejctive", objective + opt["lambda"]*penalty)
            print "="*80
        return objective +opt["lambda"]*penalty
    else:
        sys.exit("something went wrong, check output {0:s}".format(code.outfile))

def run_core_energy(x0, *args):

    # unpack the args tuple for code readability
    bsopt, bsnoopt, code, job, mol, opt = args

    if bsopt["typ"] in ["direxp", "direct", "exps", "exponents", "event", "eventemp"]:
        #penalty = sum(min(0, x)**2 for x in x0)
        penalty = 0.0
        if any(x < 0 for x in x0):
            x0 = [abs(x) for x in x0]
    else:
        penalty = 0.0

    bs2opt = bas.get_basis(x0, bsopt)

    citote = []
    stats  = []
    base   = os.path.splitext(job["inpname"])[0]
    inputs = [base+"_core"+str(sum(x))+".inp" for x in job["core"]]

    for inpname, core in zip(inputs, job["core"]):
        code.write_input(inpname, core=core, bs=bsnoopt+bs2opt, mol=mol, inpdata=job["inpdata"])

    outs = code.run_multiple(inputs)
    for out in outs:
        citote.append(code.parse_tote(out))
        stats.append(code.isok(out))

    if stats[0] and stats[1]:
        if job["verbose"]:
            print "Current exponents"
            bas.print_functions(bs2opt)
            print "x0 : ", ", ".join([str(x) for x in x0])
            print "{0:<20s} : {1:>30s} {2:>30s}".format("Terminated OK", str(stats[0]), str(stats[1]))
            print "{0:<20s} : {1:>30.10f} {2:>30.10f}".format("CI total energy", citote[0], citote[1])
            print "-"*84
        coreenergy = citote[0] - citote[1]
        if job["verbose"]:
            print "{0:<20s} : {1:>30.10f}".format("Core energy", coreenergy)
            print "{0:<20s} : {1:>30.10f}".format("Objective", coreenergy)
            print "="*84
        return coreenergy
    else:
        print "{0:<20s} : {1:>30s}".format("Job terminated with ERRORS, check output.")
        sys.exit("something went wrong, check outputs {0:s}".format(", ".join(outs)))
