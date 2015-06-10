
""" basisopt is a module for optimizing primitive exponents of basis sets"""

from __future__ import print_function
from scipy.optimize import minimize
import datetime
import os
import sys
import time
import pprint
import string
import random

# chemtools packages
from chemtools.basisset import BasisSet, get_x0

class Job(object):

    def __init__(self, calcmethod=None, objective=None, core=None, template=None,
                 regexp=None, verbose=False, code=None, optalg=None, mol=None,
                 bsopt=None, bsnoopt=None, fname=None):
        '''

        Args:
          calcmethod : str
            Electornic structure method to use in optimzation
          objective : str
            Obejctive on which the optimization will be evaluated
          code : chemtools.code.Code
            Instance of the `Code` subclass specifying which program to use
          mol : chemtools.molecule.Molecule
            Instance of the `Molecule` class specifying the system
          bsnoopt : chemtools.basisset.BasisSet
            Instance of the `BasisSet` class or list of those intances with a
            basis set whose exponents are not going to be optimized
          bsopt : dict
            A dictionary specifying the basis set to be optimized
         optalg : dict
            A dictionary specyfing the optimization algorithm and its options

        '''
        self.bsopt = bsopt
        self.bsnoopt = bsnoopt
        self.calcmethod = calcmethod
        self.objective = objective
        self.regexp = regexp
        self.template = template
        self.verbose = verbose
        self.core = core
        self.code = code
        self.optalg = optalg
        self.mol = mol
        self.results = []
        self.fname = fname
        #self.parameters = parameters
        self.function = self.get_function()

    @property
    def fname(self):
        return self._fname

    @fname.setter
    def fname(self, value):
        if value is None:
            self._fname = ''.join(random.choice(string.ascii_letters) for _ in range(10)) + '.inp'
        else:
            self._fname = value

    @property
    def optalg(self):
        return self._optalg

    @optalg.setter
    def optalg(self, value):
        '''set default optimization options if opts not given'''
        if value is None:
            self._optalg = {"method"  : "BFGS",
                            "lambda"  : 10.0,
                            "tol"     : 1.0e-4,
                            "options" : {"maxiter" : 50,
                                            "disp"    : True,
                                            "eps"     : 0.01
                                        }
                            }
        elif isinstance(value, dict):
            print('value is dict')
            if value.get('method').lower() == "nelder-mead":
                value['jacob'] = None
            else:
                value['jacob'] = False
            if not value.has_key('lambda'):
                value["lambda"] = 10.0
            self._optalg = value
        else:
            raise ValueError("optalg should be a <dict>, got: {}".format(type(value)))

    @property
    def objective(self):
        return self._objective

    @objective.setter
    def objective(self, value):

        if value in ['total energy', 'core energy', 'correlation energy', 'regexp']:
            self._objective = value
        else:
            raise ValueError('unknown objective: {}'.format(value))

    @property
    def bsopt(self):
        return self._bsopt

    @bsopt.setter
    def bsopt(self, value):

        if value is None:
            raise ValueError("no dictionary describing basis set to be optimized given")
        else:
            self._bsopt = value

    def get_function(self):
        if self.objective == 'core energy':
            return run_core_energy
        elif self.objective in ['total energy', 'correlation energy', 'regexp']:
            return run_total_energy

    def header(self):
        '''
        Return the basic information about the execution environment and
        optimization settings.
        '''

        script = sys.argv[0]
        path = os.path.abspath(os.path.dirname(script))

        out = "\n".join(["Script name : {0}".format(script),
                         "Workdir     : {0}".format(path),
                         "Start time  : {0}".format(datetime.datetime.today()),
                        ])

        return out

    def jobinfo(self):
        '''Return the information on the optimization objects and parameters'''

        #for name, obj in [("code", code), ("job", job), ("mol", mol),
        #                ("bsnoopt", bsnoopt), ("bsopt", bsopt), ("opt", opt)]:
        attrs = ['code', 'mol', 'optalg']
        out = ''
        for attr in attrs:
            out += attr.upper().center(80, "=") + '\n'
            obj = getattr(self, attr)
            if isinstance(obj, dict):
                out += pprint.pformat(obj)
            else:
                out += str(obj)
        return out

    def run(self):
        '''
        Driver for the basis set optimization

        Returns:
            res : OptimizeResult
            An instance of the ``scipy.optimize.OptimizeResult`` class
        '''

        print(self.header())
        print("\n"+"="*80, "STARTING OPTIMIZATION".center(80, "="), "="*80, end="\n\n", sep='\n')
        print(self.jobinfo())

        starttime = time.time()

        x0 = get_x0(self.bsopt)
        res = minimize(self.function, x0,
                    args=(self,),
                    method=self.optalg["method"],
                    jac=self.optalg['jacob'],
                    tol=self.optalg["tol"],
                    options=self.optalg["options"])
        print(res)
        print("Elapsed time : {0:>20.3f} sec".format(time.time()-starttime))
        self.results.append(res)

def run_total_energy(x0, *args):
    '''
    Funtion for running a single point calculation and parsing the resulting
    energy (or property) as specified by the objective function.

    Args:
        x0 : list or numpy.array
            contains a list of parameters to be optimized, may be
            explicit exponents or parametrized exponents in terms
            of some polynomial
        args : tuple of dicts
            bsopt, bsnoopt, code, job, mol, opt, needed for writing
            input and parsing output

    Returns:
        parsed result of the single point calculation as speficied by the
        objective function in the "job" dictionary

    '''

    # unpack the args tuple for code readability
    job = args[0]

    if job.bsopt["typ"] in ["direxp", "direct", "exps", "exponents", "event", "eventemp", "well", "welltemp"]:
        #penalty = sum(min(0, x)**2 for x in x0)
        penalty = 0.0
        if any(x < 0 for x in x0):
            x0 = [abs(x) for x in x0]
    else:
        penalty = 0.0

    bslist = []
    bs2opt = BasisSet.from_optdict(x0, job.bsopt)
    if job.verbose:
        print("Current exponents")
        bs2opt.print_exponents()
    if isinstance(job.bsnoopt, list):
        for bs in job.bsnoopt:
            if bs2opt.element == bs.element:
                bs2opt.add(bs)
            else:
                bslist.append(bs)
        bslist.append(bs2opt)
    else:
        bs2opt.add(job.bsnoopt)
        bslist.append(bs2opt)
    job.code.write_input(fname=job.fname, template=job.template, bs=bslist, mol=job.mol, core=job.core)
    output = job.code.run(job.fname)
    if job.code.accomplished(output):
        objective = job.code.parse(output, job.calcmethod, job.objective, job.regexp)
        if objective is None:
            raise ValueError("Unable to parse the objective, check output")
        if job.verbose:
            print("{0:<s}".format("Job Terminated without errors"))
            print("x0 : ", ", ".join([str(x) for x in x0]))
            print("\n{0:<20s} : {1:>30s}".format("Output", output))
            print("{0:<20s} : {1:>30.10f}".format("Objective", objective + job.optalg["lambda"]*penalty))
            print("="*80)
        return objective + job.optalg["lambda"]*penalty
    else:
        sys.exit("something went wrong, check output {0:s}".format(output))

def run_core_energy(x0, *args):
    '''
    Funtion for running two single point calculations and parsing the resulting
    energy (or property) as specified by the objective function, primarily
    designed to extract core energy.

    Args:
        x0: list or numpy.array
            contains a list of parameters to be optimized, may be
            explicit exponents or parametrized exponents in terms
            of some polynomial
        args: tuple of dicts
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
    if job["verbose"]:
        print("Current exponents")
        bs2opt.print_exponents()

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

    outputs = code.run_multiple(inputs)
    for output in outputs:
        citote.append(code.parse(output, job["method"], job["objective"], job.get("regexp", None)))
        stats.append(code.accomplished(output))

    if stats[0] and stats[1]:
        if job["verbose"]:
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
        sys.exit("something went wrong, check outputs {0:s}".format(", ".join(outputs)))

def opt_shell_by_nf(shell=None, nfs=None, max_params=5, opt_tol=1.0e-4, save=False, bsopt=None, **kwargs):
    '''
    For a given shell optimize the functions until the convergence criterion is reached
    the energy difference for two consecutive function numbers is less than the threshold

    Kwargs:
        shell : string
            string label for the shell to be optimized
        nfs : list of ints
            list of integers representing the number of basis functions to be
            inceremented in the optimization,
        max_params : int
            maximal number of parameters to be used in the legendre expansion,
            (length of the expansion)
        opt_tol : float
            threshold controlling the termination of the shell optimization,
            if energy difference between two basis sets with subsequent number
            of functionsis larger than this threshold, another function is
            added to this shell and parameters are reoptimized,
        save : bool
            a flag to trigger saving all the optimized basis functions for each
            shell,
        kwargs : dict
            options for the basis set optimization driver, see driver function
            from the basisopt module

    Returns:
        BasisSet object instance with optimized functions for the specified shell

    Raises:
        ValueError:
            when `shell` is different from `s`, `p`, `d`, `f`, `g`, `h`, `i`
            when number of parameters equals 1 and there are more functions
            when there are more parameters than functions to optimize
    '''

    _shells = ['s', 'p', 'd', 'f', 'g', 'h', 'i']

    if shell not in _shells:
        raise ValueError('shell must be one of the following: {}'.format(", ".join(_shells)))

    if len(bsopt['params'][0]) == 1 and min(nfs) != 1:
        raise ValueError("1 parameter and {0} functions doesn't make sense".format(min(nfs)))

    if min(nfs) < len(bsopt['params'][0]):
        raise ValueError("more parameters ({0}) than functions to optimize ({1})".format(len(bsopt['params'][0]), min(nfs)))

    bsopt["typ"] = "legendre"
    e_last = 0.0
    x_last = bsopt['params'][0]

    for nf in nfs:
        bsopt['nfpshell'] = [0]*_shells.index(shell) + [nf]
        res = driver(bsopt=bsopt, **kwargs)

        print("Completed optimization for {0:d} {1}-type functions".format(nf, shell))

        print("{s:<25s} : {v:>20.10f}".format(s="Current Function value", v=res.fun))
        print("{s:<25s} : {v:>20.10f}".format(s="Previous Function value", v=e_last))
        print("{s:<25s} : {v:>20.10f}".format(s="Difference", v=abs(res.fun-e_last)))

        if abs(res.fun - e_last) < opt_tol:
            print("Basis saturated with respect to threshold")
            bsopt['nfpshell'] = nfps_last
            if save:
                save_basis(x_last, bsopt)
            return BasisSet.from_optdict(x_last, bsopt)
        else:
            print("Threshold exceeded, continuing optimization.")
            x_last = res.x.tolist()
            nfps_last = bsopt['nfpshell']
            if len(bsopt['params'][0]) < max_params:
                print("adding more parameters")
                restup = tuple(res.x)
                # this assumption should be revised, improved or justified
                # suggestion: maybe change to average of existing parameters
                restup += tuple([abs(min(restup)*0.25)])
                bsopt['params'] = [restup,]
            else:
                print("not adding more parameters")
                bsopt['params'] = [tuple(res.x),]
            e_last = res.fun
    else:
        print("Supplied number of functions exhausted but the required accuracy was not reached")
        if save:
            save_basis(res.x.tolist(), bsopt)
        return BasisSet.from_optdict(res.x.tolist(), bsopt)

def opt_multishell(shells=None, nfps=None, guesses=None, max_params=5, opt_tol=1.0e-4, save=False, bsopt=None, **kwargs):
    '''
    Optimize a basis set by saturating the function space shell by shell

    Kwargs:
        shells (list of strings):
            list of shells to be optimized, in the order the optimization should
            be performed,

        nfps (list of lists of integers):
            list specifying a set of function numbers to be scanned per each
            shell,

        guesses (list of lists of floats):
            list specifying a set of starting parameters per each shell,

        max_params (int)
            maximal number of parameters to be used in the legendre expansion,
            (length of the expansion)

        opt_tol (float):
            threshold controlling the termination of the shell optimization,
            if energy difference between two basis sets with subsequent number
            of functionsis larger than this threshold, another function is
            added to this shell and parameters are reoptimized

        kwargs:
            options for the basis set optimization driver, see driver function
            from the basisopt module
    '''

    # bsnoopt needs to exists since it will be appended with optimized
    # functions after each shell is optimized
    if kwargs['bsnoopt'] is None:
        kwargs['bsnoopt'] = BasisSet.from_dict({"element"   : "Be",
                                      "functions" : {}})

    # begin the main loop over shells, nr f per shell and guesses
    for shell, nfs, guess in zip(shells, nfps, guesses):

        header = " Beginning optimization for {s:s} shell ".format(s=shell)
        print("="*100)
        print(format(header, '-^100'))
        print("="*100)

        bsopt['params'] = [tuple(guess)]
        # begin the optimization for a given shell and store the optimized
        # functions under optimized_shell
        optimized_shell = opt_shell_by_nf(shell, nfs, max_params=max_params,
                                          opt_tol=opt_tol, save=save,
                                          bsopt=bsopt, **kwargs)
        # add the optimized shell to the total basis set
        kwargs['bsnoopt'].add(optimized_shell)

    # save the final complete basis set
    basis_dict = vars(kwargs['bsnoopt'])
    with open("final.bas", 'wb') as ff:
        ff.write(str(basis_dict))

