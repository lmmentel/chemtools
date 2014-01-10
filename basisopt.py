#!/usr/bin/env python

# lmentel's packages 
from gamessus import GamessParser
from molecule import Atom, Molecule
from molpro import OutputParser
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

def write_molpro_input(fname=None, core=None, bs=None, code=None, mol=None):

    inpt = code["input"]

    inpt = re.sub('geometry', mol.molpro_rep(), inpt, flags=re.I)
    inpt = re.sub('basis', bas.write_molpro_basis(bs), inpt, flags=re.I)
    inpt = re.sub("core","core,{},0,0,0,0,0,0,0\n".format(core), inpt, flags=re.I)

    inp = open(fname, 'w')

    inp.write(inpt)

    inp.close()


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

def run_molpro(fname, code=None):

    '''Run a single molpro job interactively - without submitting to the queue.'''

    opts = []
    opts.append(code["exec"])
    opts.append(fname)
    opts.extend(code["opts"])
    process = Popen(opts, stdout=PIPE, stderr=PIPE)
    out, err = process.communicate()
    error_file = os.path.splitext(fname)[0]+".err"
    ferr = open(error_file, 'w')
    ferr.write(out)
    ferr.write("{0:s}\n{1:^80s}\n{0:s}\n".format("="*80, "Error messages:"))
    ferr.write(err)
    ferr.close()

def driver(code=None, mol=None, bsnoopt=None, bsopt=None, opts=None):
    '''
    Driver for the basis set optimization
    '''

    print "Script started at {0}".format(datetime.datetime.today())

    starttime = time.time()

    # set default optimization options if opts not given
    if not opts:
        opts = {"method"  : "BFGS",
                "lambda"  : 10.0,
                "tol"     : 1.0e-4,
                "options" : {"maxiter" : 50,
                             "disp"    : True,
                             "eps"     : 0.01
                            }
               }
    else:
        if opts["method"].lower() == "nelder-mead":
            jacob = None
        else:
            jacob = False
        if "lambda" not in opts.keys():
            opts["lambda"] = 10.0

    if not bsnoopt:
        bsnooptobj = []

    if code:
        if code["objective"].lower() == "coreenergy":
            if len(code["core"]) == 2:
                if code["core"][0] > code["core"][1]:
                    code["core"].reverse()
            else:
                sys.exit("<driver>: core should have 2 entries, but it has {0:4d}".format(len(code["core"])))
            function = run_core_energy
        elif code["objective"].lower() == "totalenergy":
            function = run_total_energy
        else:
            sys.exit("<driver>: wrong objective in code dictionary")
        if code["name"].lower() in ["gamess", "gamessus", "gamess-us"]:
            if bsnoopt and len(bsnoopt) > 0:
                bsnooptobj = bas.parse_gamess_basis(bsnoopt)
        elif code["name"].lower() in ["molpro"]:
            if bsnoopt and len(bsnoopt) > 0:
                bsnooptobj = bas.parse_molpro_basis(bsnoopt)
        else:
            sys.exit("<driver>: wrong code name: {}".format(code["name"]))
    else:
        sys.exit("<driver>: no dictionary describing code given")

    if not bsopt:
        sys.exit("<driver>: no dictionary describing basis set to be optimized given")

    if not mol:
        sys.exit("<driver>: no molecule object specified")

    x0 = bas.get_x0(bsopt)

    res = minimize(function, x0,
            args=(bsopt, bsnooptobj, code, mol, opts,),
            method=opts["method"],
            jac=jacob,
            tol=opts["tol"],
            options=opts["options"])
    print res
    print "Elapsed time : {0:>20.3f} sec".format(time.time()-starttime)

def run_total_energy(x0, *args):

    if args[0]["typ"] in ["direxp", "direct", "exps", "exponents", "event", "eventemp"]:
        #penalty = sum(min(0, x)**2 for x in x0)
        penalty = 0.0
        if any(x < 0 for x in x0):
            x0 = [abs(x) for x in x0]
    else:
        penalty = 0.0

    bs2opt = bas.get_basis(x0, args[0])

    fname = args[2]["inputname"]
    core  = args[2]["core"][0]

    if args[2]["name"].lower() in ["molpro"]:
        write_molpro_input(fname, core, code=args[2], bs=args[1]+bs2opt, mol=args[3])
        run_molpro(fname, args[2])
        output = os.path.splitext(fname)[0] + ".out"
        parser = OutputParser(out=output)
        hfenergy = parser.get_hf_total_energy()
        cienergy = parser.get_cisd_total_energy()
    elif args[2]["name"].lower() in ["gamess", "gamess-us", "gamessus"]:
        write_gamess_input(fname, core, code=args[2], bs=args[1]+bs2opt, mol=args[3])
        output = run_gamess(fname, args[2])
        parser = GamessParser(log=output)
        hfenergy = parser.get_hf_total_energy()
        cienergy = parser.get_ormas_total_energy()
    else:
        sys.exit("<run_total_energy>: unknown code")

    print "Current exponents"
    bas.print_functions(bs2opt)
    print "x0 : ", ", ".join([str(x) for x in x0])
    print "\n{0:<20s} : {1:>30s}".format("Output", output)

    if parser.terminatedOK():
        termok = "YES"
        print "{0:<20s} : {1:>30s}".format("Terminated OK", termok)
        if args[2]["method"] == "hf":
            print "{0:<20s} : {1:>30.10f}".format("HF total energy", hfenergy)
            print "{0:<20s} : {1:>30.10f}".format("Obejctive", hfenergy + args[4]["lambda"]*penalty)
            print "="*80
            return hfenergy + args[4]["lambda"]*penalty
        if args[2]["method"] in ["cisd", "mrcisd"]:
            print "{0:<20s} : {1:>30.10f}".format("HF total energy", hfenergy)
            print "{0:<20s} : {1:>30.10f}".format("CI total energy", cienergy)
            print "{0:<20s} : {1:>30.10f}".format("Objective", cienergy + args[4]["lambda"]*penalty)
            print "="*80
            return cienergy + args[4]["lambda"]*penalty
    else:
        termok = "NO"
        print "{0:<20s} : {1:>30s}".format("Terminated OK", termok)
        if args[2]["name"] == "gamess" and not isinstance(parser.get_linear_deps(), type(None)):
            return parser.get_hf_total_energy() + args[4]["lambda"]*parser.get_linear_deps()
        else:
            sys.exit("something went wrong, check output {0:s}".format(output))

def run_core_energy(x0, *args):

    if args[0]["typ"] in ["direxp", "direct", "exps", "exponents"]:
        penalty = sum(min(0, x)**2 for x in x0)
        if any(x < 0 for x in x0):
            x0 = [abs(x) for x in x0]
    else:
        penalty = 0.0

    bs2opt = bas.get_basis(x0, args[0])

    nb = os.path.splitext(args[2]["inputname"])[0]
    outs = []
    pars = []
    hfe  = []
    cie  = []
    for fname, core in zip([nb+"_core-"+str(x)+".inp" for x in args[2]["core"]], args[2]["core"]):
        if args[2]["name"].lower() in ["molpro"]:
            write_molpro_input(fname, core, code=args[2], bs=args[1]+bs2opt, mol=args[3])
            run_molpro(fname, args[2])
            output = os.path.splitext(fname)[0] + ".out"
            outs.append(output)
            parser = OutputParser(out=output)
            pars.append(parser)
            hfe.append(parser.get_hf_total_energy())
            cie.append(parser.get_cisd_total_energy())
        elif args[2]["name"].lower() in ["gamess", "gamess-us", "gamessus"]:
            write_gamess_input(fname, core, code=args[2], bs=args[1]+bs2opt, mol=args[3])
            output = run_gamess(fname, args[2])
            outs.append(output)
            parser = GamessParser(log=output)
            pars.append(parser)
            hfe.append(parser.get_hf_total_energy())
            cie.append(parser.get_ormas_total_energy())
        else:
            sys.exit("<run_total_energy>: unknown code")

    print "Current exponents"
    bas.print_functions(bs2opt)
    print "x0 : ", ", ".join([str(x) for x in x0])
    print "\n{0:<20s} : {1:>30s} {2:>30s}".format("Outputs", outs[0], outs[1])

    termok = []
    for parser in pars:
        if parser.terminatedOK():
            termok.append("YES")
        else:
            termok.append("NO")

    print "{0:<20s} : {1:>30s} {2:>30s}".format("Terminated OK", termok[0], termok[1])
    if pars[0].terminatedOK() and pars[1].terminatedOK():
        print "{0:<20s} : {1:>30.10f} {2:>30.10f}".format("HF total energy", hfe[0], hfe[1])
        print "{0:<20s} : {1:>30.10f} {2:>30.10f}".format("CI total energy", cie[0], cie[1])
        coreenergy = cie[0] - cie[1]
        print "-"*84
        print "{0:<20s} : {1:>30.10f}".format("Core energy", coreenergy)
        print "{0:<20s} : {1:>30.10f}".format("Objective", coreenergy)
        print "="*84
        return coreenergy
    else:
        if args[2]["name"] == "gamess" and not isinstance(gpsc.get_linear_deps(), type(None)):
            return gpsc.get_hf_total_energy() + args[4]["lambda"]*gp.get_linear_deps()
        else:
            sys.exit("something went wrong, check output {0:s}".format(output))
