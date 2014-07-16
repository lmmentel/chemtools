# standard python packages
from scipy.optimize import minimize
import datetime
import numpy as np
import os
import re
import sys
import time
# chemtools packages
from basisset import BasisSet, get_x0

def driver(code=None, job=None, mol=None, bsnoopt=None, bsopt=None, opt=None):
    '''
    Driver for the basis set optimization
    '''

    print("Script started at {0}".format(datetime.datetime.today()))

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
    else:
        sys.exit("<driver>: no dictionary describing code given")

    if not bsopt:
        sys.exit("<driver>: no dictionary describing basis set to be optimized given")

    if not mol:
        sys.exit("<driver>: no molecule object specified")

    x0 = get_x0(bsopt)
    res = minimize(function, x0,
            args=(bsopt, bsnoopt, code, job, mol, opt,),
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
    Funtion for running a single point calculation and parsing the resulting
    energy (or property) as specified by the objective function.

    Args:
        x0: list/numpy 1d array
            contains a list of parameters to be optimized, may be
            explicit exponents or parametrized exponents in terms
            of some polynomial
        args: tuple of dictionaries
            bsopt, bsnoopt, code, job, mol, opt, needed for writing
            input and parsing output
    Returns:
        parsed result of the single point calculation as speficied by the
        objective function in the "job" dictionary

    '''

    # unpack the args tuple for code readability
    bsopt, bsnoopt, code, job, mol, opt = args

    if bsopt["typ"] in ["direxp", "direct", "exps", "exponents", "event", "eventemp", "well", "welltemp"]:
        #penalty = sum(min(0, x)**2 for x in x0)
        penalty = 0.0
        if any(x < 0 for x in x0):
            x0 = [abs(x) for x in x0]
    else:
        penalty = 0.0

    bslist = []
    bs2opt = BasisSet.from_optdict(x0, bsopt)
    if isinstance(bsnoopt, list):
        for bs in bsnoopt:
            if bs2opt.element == bs.element:
                bs2opt.add(bs)
            else:
                bslist.append(bs)
        bslist.append(bs2opt)
    else:
        bs2opt.add(bsnoopt)
        bslist.append(bs2opt)
    code.write_input(job["inpname"], job["core"], bs=bslist, mol=mol, inpdata=job["inpdata"])
    code.run(job["inpname"])
    if code.accomplished():
        objective = code.parse(job["method"], job["objective"], job.get("regexp", None))
        if job["verbose"]:
            print("{0:<s}".format("Job Terminated without errors"))
            print("Current exponents")
            bs2opt.print_exponents()
            print("x0 : ", ", ".join([str(x) for x in x0]))
            print("\n{0:<20s} : {1:>30s}".format("Output", code.outfile))
            print("{0:<20s} : {1:>30.10f}".format("Objective", objective + opt["lambda"]*penalty))
            print("="*80)
        return objective +opt["lambda"]*penalty
    else:
        sys.exit("something went wrong, check output {0:s}".format(code.outfile))

def run_core_energy(x0, *args):
    '''
    Funtion for running two single point calculations and parsing the resulting
    energy (or property) as specified by the objective function, primarily
    designed to extract core energy.

    Args:
        x0: list/numpy 1d array
            contains a list of parameters to be optimized, may be
            explicit exponents or parametrized exponents in terms
            of some polynomial
        args: tuple of dictionaries
            bsopt, bsnoopt, code, job, mol, opt, needed for writing
            input and parsing output
    Returns:
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

    bs2opt = BasisSet.from_optdict(x0, bsopt)

    citote = []
    stats  = []
    base   = os.path.splitext(job["inpname"])[0]
    inputs = [base+"_core"+str(sum(x))+".inp" for x in job["core"]]
    if isinstance(bsnoopt, list):
        for bs in bsnoopt:
            bs2opt.add(bs)
    else:
        bs2opt.add(bsnoopt)

    for inpname, core in zip(inputs, job["core"]):
        code.write_input(inpname, core=core, bs=bs2opt, mol=mol, inpdata=job["inpdata"])

    outs = code.run_multiple(inputs)
    for out in outs:
        citote.append(code.parse_tote(out))
        stats.append(code.accomplished(out))

    if stats[0] and stats[1]:
        if job["verbose"]:
            print("Current exponents")
            bs2opt.print_exponents()
            print("x0 : ", ", ".join([str(x) for x in x0]))
            print("{0:<20s} : {1:>30s} {2:>30s}".format("Terminated OK", str(stats[0]), str(stats[1])))
            print("{0:<20s} : {1:>30.10f} {2:>30.10f}".format("CI total energy", citote[0], citote[1]))
            print("-"*84)
        coreenergy = citote[0] - citote[1]
        if job["verbose"]:
            print("{0:<20s} : {1:>30.10f}".format("Core energy", coreenergy))
            print("{0:<20s} : {1:>30.10f}".format("Objective", coreenergy))
            print("="*84)
        return coreenergy
    else:
        print("Job terminated with ERRORS, check output.")
        sys.exit("something went wrong, check outputs {0:s}".format(", ".join(outs)))
