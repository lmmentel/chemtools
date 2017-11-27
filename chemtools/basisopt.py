# -*- coding: utf-8 -*-

#The MIT License (MIT)
#
#Copyright (c) 2014 Lukasz Mentel
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.

"""\
``basisopt`` is a module containing flexible methods for optimizing primitive
exponents of basis sets"""

from __future__ import division, print_function

import atexit
import datetime
import os
import pprint
import random
import string
import sys
import time
from collections import OrderedDict
import numpy as np
from scipy.optimize import minimize

from chemtools.basisparse import ORBITALS
from chemtools.basisset import BasisSet, get_num_params, sliceinto


class BSOptimizer(object):
    '''
    Basis Set Optimizer class is a convenient wrapper for optimizing primitive
    exponents of Gaussian basis sets using different code interfaces for
    performing the actual electronic structure calculations.

    Args:
        objective : str or callable
            Name of the objective using which the optimization will be
            evaluated, if it is a callable/function it should take the output
            name as argument

        code : :py:class:`Calculator <chemtools.calculators.calculator.Calculator>`
            Subclass of the :py:class:`Calculator <chemtools.calculators.calculator.Calculator>`
            specifying the electronic structure program to use

        mol : :py:class:` Molecule <chemtools.molecule.Molecule>`
            Molecule specifying the system

        staticbs : dict or :py:class:`BasisSet <chemtools.basisset.BasisSet>`
            Basis set or ``dict`` of basis set with basis sets whose exponents
            are not going to be optimized

        fsopt : dict
            A dictionary specifying the basis set to be optimized, the keys
            should be element/atom symbol and the values should contain lists
            of 4-tuples composed of: shell, type, number of functions and
            parameters/exponents, for example,

            >>> fs = {'H' : [('s', 'et', 10, (0.5, 2.0))]}

            describes 10 ``s``-type exponents for H generated from even
            tempered formula with parameters 0.5 and 2.0

        optalg : dict
            A dictionary specifying the optimization algorithm and its options

        fname : str
            Name of the job/input file for the single point calculator

        uselogs : bool
            Use natural logarithms of exponents in optimization rather than
            their values, this option is only relevant if functions asre given
            as ``exp`` or ``exponents``

        regexp : str
            Regular expression to use in search for the objective if
            ``objective`` is regexp

        runcore : bool
            Flag to mark wheather to run separate single point jobs with
            different numbers of frozen core orbitals to calculate core energy

        penalize : bool
            Flag enabling the use of penalty function, when two exponent in
            any shell are too close the objective is multiplied by a penalty
            factor calculated in :py:func:`get_penalty
            <chemtools.basisopt.get_penalty>`

        penaltykwargs : dict
            Keyword arguments for the penalty function, default
            ``{'alpha' : 25.0, 'smallestonly' : True}``
    '''

    def __init__(self, objective=None, core=None, template=None,
                 regexp=None, verbose=False, code=None, optalg=None, mol=None,
                 fsopt=None, staticbs=None, fname=None, uselogs=True,
                 runcore=False, penalize=None, penaltykwargs=None,
                 logfile=None):

        self.fsopt = fsopt
        self.staticbs = staticbs
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
        self.uselogs = uselogs
        self.penalize = penalize
        self.penaltykwargs = penaltykwargs
        self.logfile = logfile

        if runcore:
            self.function = run_core_energy
        else:
            self.function = run_total_energy

    @property
    def fname(self):
        return self._fname

    @fname.setter
    def fname(self, value):
        '''
        Set the name of the job/input file for the code used for single point
        calculations.
        '''
        if value is None:
            self._fname = 'bso_' + ''.join(random.choice(string.ascii_letters)
                                           for _ in range(10)) + '.inp'
        else:
            self._fname = value

    @property
    def optalg(self):
        return self._optalg

    @optalg.setter
    def optalg(self, value):
        '''
        set default optimization options if opts not given
        '''

        # defaults
        tol = 1.0e-4
        jac = False
        mxi = 100

        if value is None:
            self._optalg = {"method": "Nelder-Mead",
                            "tol": tol,
                            "jacob": None,
                            "options": {"maxiter": 100,
                                        "disp": True,
                                        }
                            }

        elif isinstance(value, dict):
            if value.get('method').lower() == "nelder-mead":
                value['jacob'] = None
            else:
                value['jacob'] = jac
            if 'tol' not in value.keys():
                value['tol'] = tol
            if 'maxiter' not in value['options'].keys():
                value['options']['maxiter'] = mxi
            self._optalg = value
        else:
            raise ValueError("optalg should be a <dict>, got: {}".format(type(value)))

    @property
    def fsopt(self):
        return self._fsopt

    @fsopt.setter
    def fsopt(self, value):

        seqs = ["et", "even", "eventemp", "even tempered",
                "wt", "well", "welltemp", "well tempered",
                "le", "legendre", "exp", "exponents"]

        if value is None:
            raise ValueError("<fsopt> cannot be None")
        else:
            if isinstance(value, dict):
                for symbol, funlist in value.items():
                    for funtuple in funlist:
                        if funtuple[0].lower() not in ORBITALS:
                            msg = '<{}>: <shell> '.format(symbol) +\
                                  ' should be one of {}'.format(', '.join(ORBITALS)) +\
                                  ', got: {}'.format(funtuple[0])
                            raise ValueError(msg)
                        if funtuple[1] not in seqs:
                            msg = '<{}>: <seq> '.format(symbol) +\
                                  ' should be one of {}'.format(', '.join(seqs)) +\
                                  ', got: {}'.format(funtuple[1])
                            raise ValueError(msg)
            else:
                ValueError('<fsopt> should be a dict, got: {}'.format(type(value)))

            self._fsopt = OrderedDict(sorted(value.items(),
                                             key=lambda t: t[0]))

    @property
    def penaltykwargs(self):
        return self._penaltykwargs

    @penaltykwargs.setter
    def penaltykwargs(self, value):

        if value is None:
            self._penaltykwargs = {'alpha': 25.0, 'smallestonly': True}
        else:
            self._penaltykwargs = value

    @property
    def logfile(self):
        return self._logfile

    @logfile.setter
    def logfile(self, value):
        '''
        if a filename is specified open a text file for writing with
        line buffering
        '''

        if value is None:
            self._logfile = None
            self.log = sys.stdout
        else:
            self._logfile = value
            self.log = open(value, 'wt', buffering=1)
            atexit.register(self.finalize)

    def finalize(self):

        if self.log is not sys.stdout:
            self.log.close()

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

        # for name, obj in [("code", code), ("job", job), ("mol", mol),
        #                ("bsnoopt", bsnoopt), ("bsopt", bsopt), ("opt", opt)]:
        # TODO add printing of static basis set
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
        Collect all the parameters in an array of consecutive elements.
        '''

        x0 = np.empty(0)
        for fperatom in self.fsopt.values():
            for fs in fperatom:
                if fs[1].lower() in ['exp', 'exponents'] and self.uselogs:
                    x0 = np.concatenate((x0, np.log(fs[-1])), axis=0)
                else:
                    x0 = np.concatenate((x0, np.array(fs[-1])), axis=0)
        return x0

    def run(self):
        '''
        Start the basis set optimization

        Returns:
            res : OptimizeResult
                An instance of the ``scipy.optimize.OptimizeResult`` class
        '''

        self.log.write(self.header())
        self.log.write("\n" + "\n".join(["=" * 80,
            "STARTING OPTIMIZATION".center(80, "="), "=" * 80]) + "\n\n")
        self.log.write(self.jobinfo())
        self.log.flush()

        starttime = time.time()

        x0 = self.get_x0()
        self.result = minimize(self.function, x0,
                               args=(self,),
                               method=self.optalg["method"],
                               jac=self.optalg['jacob'],
                               tol=self.optalg["tol"],
                               options=self.optalg["options"])
        self.log.write(str(self.result))
        self.log.write("Elapsed time : {0:>20.3f} sec".format(time.time() -
                                                              starttime))

    def get_basis(self, name=None, element=None):
        '''
        Construct the BasisSet object from the result of the
        optimized exponents and function definition.

        Args:
            name : str
                Name to be assigned to the basis set
            element : str
                Element symbol for the basis set

        Returns:
            basis : chemtools.basisset.BasisSet
                :py:class:`BasisSet <chemtools.basisset.BasisSet>` object with
                the optimized functions
        '''

        bsdict = {}
        # get number of parameters per atom
        npars = [sum(get_num_params(t) for t in funs) for funs in self.fsopt.values()]
        x0peratom = sliceinto(self.result.x, npars)
        for (atom, funs), xpars in zip(self.fsopt.items(), x0peratom):

            bsdict[atom] = BasisSet.from_optpars(xpars, funs, name=name,
                                                 element=element,
                                                 explogs=self.uselogs)

        if self.staticbs is not None:
            if isinstance(self.staticbs, dict):
                common_atoms = set(self.staticbs.keys()) & set(bsdict.keys())
                if common_atoms:
                    for atom in common_atoms:
                        bsdict[atom].append(self.staticbs[atom])
                diff_atoms = set(self.staticbs.keys()) - set(bsdict.keys())
                if diff_atoms:
                    for atom in diff_atoms:
                        bsdict[atom] = self.staticbs[atom]
            elif isinstance(self.staticbs, BasisSet):
                if self.staticbs.element in set(bsdict.keys()):
                    bsdict[atom].append(self.staticbs)
                else:
                    bsdict[self.staticbs.element] = self.staticbs

        return bsdict


def get_basis_dict(bso, x0):
    '''
    Return a dictionary with :py:class:`BasisSet` objects as values and
    element symbols as keys. The dictionary is composed based on the current
    parameters ``x0`` and attributes of the :py:class:`BSOptimizer` including
    the ``staticbs``
    '''

    bsdict = dict()
    for atom, functs in bso.fsopt.items():
        bsdict[atom] = BasisSet.from_optpars(x0, funs=functs, name='opt',
                                             element=atom, explogs=bso.uselogs)

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
        elif isinstance(bso.staticbs, BasisSet):
            if bso.staticbs.element in set(bsdict.keys()):
                bsdict[atom].append(bso.staticbs)
            else:
                bsdict[bso.staticbs.element] = bso.staticbs

    return bsdict


def get_penalty(bsdict, alpha=25.0, smallestonly=True):
    '''
    For a given dict of basis sets calculate the penalty for pairs of
    exponents being too close together.

    Args:
        bsdict : dict
            Dictionary of :py:class:`BasisSet` objects

        alpha : float
            Parameter controlling the magnitude and range of the penalty

        smallestonly : bool
            A flag to mark whether to use only the smallest ratio to calculate
            the penalty or all smallest ratios from each shell and calculate
            the penalty as a product of individual penalties

    For each basis and shell within the basis ratios between pairs of sorted
    exponents are calculated. The minimal ratio (closest to 1.0) is taken to
    calculate the penalty according to the formula

    .. math::

       P = 1 + \\exp(-\\alpha (r - 1))

    where :math:`r` is the ratio between the two closest exponents (>1).
    '''

    minratios = []
    for basis in bsdict.values():
        for functions in basis.functions.values():
            exps = np.sort(functions['e'])
            ratios = exps[1:] / exps[:-1]
            minratios.append(np.min(ratios))

    if smallestonly:
        return 1 + np.exp(-alpha * (np.min(np.array(minratios)) - 1))
    else:
        return np.prod(1 + np.exp(-alpha * (np.array(minratios) - 1)))


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
        ni = 0
        nt = 0
        for shell, seq, nf, params in functs:
            nt += nf
            if seq not in ['le', 'legendre'] and not bso.uselogs:
                x0[ni:nt] = np.abs(x0[ni:nt])
            ni += nf

    bsdict = get_basis_dict(bso, x0)

    # set the penalty value
    if bso.penalize:
        penalty = get_penalty(bsdict, **bso.penaltykwargs)
    else:
        penalty = 1.0

    if bso.verbose:
        bso.log.write("Current exponents being optimized:\n")
        for atom, functs in bso.fsopt.items():
            basis = BasisSet.from_optpars(x0, funs=functs, name='opt',
                                          element=atom, explogs=bso.uselogs)
            bso.log.write(atom + basis.print_functions())
            bso.log.flush()

    bso.code.write_input(fname=bso.fname, template=bso.template, basis=bsdict,
                         mol=bso.mol, core=bso.core)
    output = bso.code.run(bso.fname)
    if bso.code.accomplished(output):

        if callable(bso.objective):
            objective = bso.objective(output)
        else:
            objective = bso.code.parse(output, bso.objective, bso.regexp)

        if objective is None:
            raise ValueError("Unable to parse the objective, check output")
        if bso.verbose:
            bso.log.write("{0:<s}".format("Job Terminated without errors\n"))
            bso.log.write("x0 : " + ", ".join([str(x) for x in x0]) + "\n")
            bso.log.write("\n{0:<20s} : {1:>30s}\n".format("Output", output))
            bso.log.write("{0:<20s} : {1:>30.10f}\n".format(str(bso.objective), objective))
            bso.log.write("{0:<20s} : {1:>30.10f}\n".format("Objective", objective * penalty))
            bso.log.write("=" * 80 + "\n")
            bso.log.flush()
        return objective * penalty
    else:
        raise ValueError("something went wrong, check output {0:s}".format(output))


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
        ni = 0
        nt = 0
        for shell, seq, nf, params in functs:
            nt += nf
            if seq not in ['le', 'legendre'] and not bso.uselogs:
                x0[ni:nt] = np.abs(x0[ni:nt])
            ni += nf

    bsdict = get_basis_dict(bso, x0)

    # set the penalty value
    if bso.penalize:
        penalty = get_penalty(bsdict, bso.penaltykwargs)
    else:
        penalty = 1.0

    if bso.verbose:
        bso.log.write("Current exponents being optimized:\n")
        for atom, functs in bso.fsopt.items():
            basis = BasisSet.from_optpars(x0, funs=functs, name='opt',
                                          element=atom, explogs=bso.uselogs)
            bso.log.write(atom + "\n" + basis.print_functions())
            bso.log.flush()

    citote = []
    stats = []
    base = os.path.splitext(bso.fname)[0]
    inputs = [base + "_core" + str(sum(x)) + ".inp" for x in bso.core]

    for inpname, core in zip(inputs, bso.core):
        bso.code.write_input(fname=inpname, core=core, basis=bsdict,
                             mol=bso.mol, template=bso.template)

    outputs = bso.code.run_multiple(inputs)
    for output in outputs:
        citote.append(bso.code.parse(output, bso.objective, bso.regexp))
        stats.append(bso.code.accomplished(output))

    if stats[0] and stats[1]:
        if bso.verbose:
            bso.log.write("x0 : " + ", ".join([str(x) for x in x0]) + "\n")
            bso.log.write("{0:<20s} : {1:>30s} {2:>30s}\n".format("Terminated OK",
                                                                str(stats[0]),
                                                                str(stats[1])))
            bso.log.write("{0:<20s} : {1:>30.10f} {2:>30.10f}\n".format(str(bso.objective),
                                                              citote[0],
                                                              citote[1]))
            bso.log.write("-" * 84)
        coreenergy = citote[0] - citote[1]
        if coreenergy > 0.0:
            coreenergy = -1.0 * coreenergy
        if bso.verbose:
            bso.log.write("{0:<20s} : {1:>30.10f}\n".format("Core energy",
                                                            coreenergy))
            bso.log.write("{0:<20s} : {1:>30.10f}\n".format("Objective",
                                                            coreenergy * penalty))
            bso.log.write("=" * 84)
        bso.log.flush()
        return coreenergy * penalty
    else:
        raise ValueError("something went wrong, check outputs {0:s}".format(", ".join(outputs)))


def opt_shell_by_nf(shell=None, nfs=None, max_params=5, opt_tol=1.0e-4,
                    save=False, bsopt=None, **kwargs):
    '''
    For a given shell optimize the functions until the convergence criterion
    is reached the energy difference for two consecutive function numbers is
    less than the threshold

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
        BasisSet object instance with optimized functions for the specified
        shell

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
