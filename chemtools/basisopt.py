
""" basisopt is a module for optimizing primitive exponents of basis sets"""

from __future__ import division, print_function
from scipy.optimize import minimize
import datetime
import os
import sys
import time
import pprint
import string
import random
import numpy as np

# chemtools packages
from chemtools.basisset import BasisSet

class BSOptimizer(object):

    def __init__(self, method=None, objective=None, core=None, template=None,
                 regexp=None, verbose=False, code=None, optalg=None, mol=None,
                 fsopt=None, staticbs=None, fname=None):

        '''
        Args:
          method : str
            Electornic structure method to use in optimzation
          objective : str
            Obejctive on which the optimization will be evaluated
          code : chemtools.code.Code
            Instance of the `Code` subclass specifying which program to use
          mol : chemtools.molecule.Molecule
            Instance of the `Molecule` class specifying the system
          staticbs : chemtools.basisset.BasisSet
            Instance of the `BasisSet` class or list of those intances with a
            basis set whose exponents are not going to be optimized
          fsopt : dict
            A dictionary specifying the basis set to be optimized
          optalg : dict
            A dictionary specyfing the optimization algorithm and its options
        '''

        self.fsopt = fsopt
        self.staticbs = staticbs
        self.method = method
        self.objective = objective
        self.regexp = regexp
        self.template = template
        self.verbose = verbose
        self.core = core
        self.code = code
        self.optalg = optalg
        self.mol = mol
        self.result = None
        self.fname = fname
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

        #defaults
        tol = 1.0e-4
        lbd = 10.0
        jac = False
        mxi = 100

        if value is None:
            self._optalg = {"method"  : "Nelder-Mead",
                            "lambda"  : lbd,
                            "tol"     : tol,
                            "jacob"   : None,
                            "options" : {"maxiter" : 100,
                                         "disp"    : True,
                                        }
                            }
        elif isinstance(value, dict):
            if value.get('method').lower() == "nelder-mead":
                value['jacob'] = None
            else:
                value['jacob'] = jac
            if not value.has_key('lambda'):
                value["lambda"] = lbd
            if not value.has_key('tol'):
                value["lambda"] = tol
            if not value['options'].has_key('maxiter'):
                value['options']['maxiter'] = mxi
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
    def fsopt(self):
        return self._fsopt

    @fsopt.setter
    def fsopt(self, value):

        if value is None:
            raise ValueError("no dictionary describing basis set to be optimized given")
        else:
            self._fsopt = value

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

    def get_x0(self):
        '''
        Collect all the parameters in a consecutive list of elements.
        '''

        x0 = list()
        for fperatom in self.fsopt.values():
            for fs in fperatom:
                x0.extend(fs[-1])
        return x0

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

        x0 = self.get_x0()
        self.result = minimize(self.function, x0,
                        args=(self,),
                        method=self.optalg["method"],
                        jac=self.optalg['jacob'],
                        tol=self.optalg["tol"],
                        options=self.optalg["options"])
        print(self.result)
        print("Elapsed time : {0:>20.3f} sec".format(time.time()-starttime))

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
    bso = args[0]

    for atom, functs in bso.fsopt.items():
        ni = 0; nt = 0
        for shell, seq, nf, params in functs:
            nt += nf
            if seq not in ['le', 'legendre']:
                x0[ni:nt] = np.abs(x0[ni:nt])
            ni += nf

    penalty = 0.0

    bsdict = dict()
    for atom, functs in bso.fsopt.items():
        bsdict[atom] = BasisSet.from_optpars(x0, functs=functs, name='opt', element=atom)

    if bso.verbose:
        print("Current exponents being optimized:")
        for atom, basis in bsdict.items():
            print(atom, basis.print_functions())

    if bso.staticbs is not None:
        if isinstance(bso.staticbs, dict):
            common_atoms = set(bso.staticbs.keys()) & set(bsdict.keys())
            if common_atoms:
                for atom in common_atoms:
                    bsdict[atom].append(bso.staticbs[atom])
            diff_atoms = set(bso.staticbs.keys()) - set(bsdict.keys())
            if diff_atoms:
                for atom in diff_atoms:
                    bsdict[atom] = bso.staticbs[atom]

    bso.code.write_input(fname=bso.fname, template=bso.template, bs=bsdict.values(), mol=bso.mol, core=bso.core)
    output = bso.code.run(bso.fname)
    if bso.code.accomplished(output):
        objective = bso.code.parse(output, bso.method, bso.objective, bso.regexp)
        if objective is None:
            raise ValueError("Unable to parse the objective, check output")
        if bso.verbose:
            print("{0:<s}".format("Job Terminated without errors"))
            print("x0 : ", ", ".join([str(x) for x in x0]))
            print("\n{0:<20s} : {1:>30s}".format("Output", output))
            print("{0:<20s} : {1:>30.10f}".format("Objective", objective + bso.optalg["lambda"]*penalty))
            print("="*80)
        return objective + bso.optalg["lambda"]*penalty
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
    bso = args[0]

    for atom, functs in bso.fsopt.items():
        ni = 0; nt = 0
        for shell, seq, nf, params in functs:
            nt += nf
            if seq not in ['le', 'legendre']:
                x0[ni:nt] = np.abs(x0[ni:nt])
            ni += nf

    penalty = 0.0

    bsdict = dict()
    for atom, functs in bso.fsopt.items():
        bsdict[atom] = BasisSet.from_optpars(x0, functs=functs, name='opt', element=atom)

    if bso.verbose:
        print("Current exponents being optimized:")
        for atom, basis in bsdict.items():
            print(atom, basis.print_functions())

    if bso.staticbs is not None:
        if isinstance(bso.staticbs, dict):
            common_atoms = set(bso.staticbs.keys()) & set(bsdict.keys())
            if common_atoms:
                for atom in common_atoms:
                    bsdict[atom].append(bso.staticbs[atom])
            diff_atoms = set(bso.staticbs.keys()) - set(bsdict.keys())
            if diff_atoms:
                for atom in diff_atoms:
                    bsdict[atom] = bso.staticbs[atom]

    citote = []
    stats  = []
    base   = os.path.splitext(bso.fname)[0]
    inputs = [base+"_core"+str(sum(x))+".inp" for x in bso.core]

    for inpname, core in zip(inputs, bso.core):
        bso.code.write_input(fname=inpname, core=core, bs=bsdict.values(), mol=bso.mol, template=bso.template)

    outputs = bso.code.run_multiple(inputs)
    for output in outputs:
        citote.append(bso.code.parse(output, bso.method, bso.objective, bso.regexp))
        stats.append(bso.code.accomplished(output))

    if stats[0] and stats[1]:
        if bso.verbose:
            print("x0 : ", ", ".join([str(x) for x in x0]))
            print("{0:<20s} : {1:>30s} {2:>30s}".format("Terminated OK", str(stats[0]), str(stats[1])))
            print("{0:<20s} : {1:>30.10f} {2:>30.10f}".format("CI total energy", citote[0], citote[1]))
            print("-"*84)
        coreenergy = citote[0] - citote[1]
        if coreenergy > 0.0:
            coreenergy = -1.0*coreenergy
        if bso.verbose:
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

