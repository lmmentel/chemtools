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

from __future__ import division, print_function

from copy import deepcopy
import numpy as np
from itertools import chain
from collections import OrderedDict
from scipy.linalg import sqrtm, inv
from scipy.special import factorial, factorial2, binom
import pickle
import re
import os
from mendeleev import element
from chemtools.basisparse import parse_basis, CFDTYPE, SHELLS

def read_pickle(fname):
    '''Read a pickled BasisSet object from file

    Args:
      fname : str
        File name containing the BasisSet
    '''

    with open(fname, 'rb') as fil:
        return pickle.load(fil)

class BasisSet(object):
    '''
    Basis set module supporting basic operation on basis sets and can be used
    as a API for mongoDB basis set repository.
    '''

    def __init__(self, name, element, family=None, kind=None,
                 functions=None):
        '''
        Args:
          name : str
            Name of the basis set
          element : str
            Symbol of the element

        Kwargs:
          kind : str
            Classification of the basis set functions, *diffuse*, *tight*
          family : str
            basis set family
          functions : dict
            Dict of functions with *s*, *p*, *d*, *f*, ... as keys
          params : list of dicts
            Parameters for generating the functions according to the model
        '''

        self.name = name
        self.element = element
        self.family = family
        self.kind = kind
        self.functions = functions

    @classmethod
    def from_optpars(cls, x0, functs=None, name=None, element=None):
        '''
        Return a basis set object generated from a sequence based on the specified
        arguments.

        Args:
          x0 : list or numpy.array
            Parameters to generate the basis set as a contonuous list or array
          functs : list of tuples
            A list of tuple specifying the shell type, number of functions and
            parameters, e.g. `[('s', 'et', 4, (0.5, 2.0)), ('p', 'et', 3, (1.0, 3.0))]`
          name : str
            Name of the basis set
          element : str
            Chemical symbol of the element

        Returns:
          out : BasisSet
        '''

        funs = dict()
        ni = 0; nt = 0
        for shell, seq, nf, params in functs:
            funs[shell] = dict()
            if seq in ["exp", "exponents"]:
                nt += nf
                funs[shell]['e'] = generate_exponents(seq, nf, x0[ni:nt])
                ni += nf
            else:
                nt += len(params)
                funs[shell]['e'] = generate_exponents(seq, nf, x0[ni:nt])
                ni += len(params)
        functions = OrderedDict(sorted(funs.items(), key=lambda x: SHELLS.index(x[0])))
        bs = cls(name=name, element=element, functions=functions)
        bs.uncontract()
        return bs

    @classmethod
    def from_sequence(cls, functs=None, name=None, element=None):
        '''
        Return a basis set object generated from a sequence based on the specified
        arguments.

        Args:
          functs : list of tuples
            A list of tuple specifying the shell type, number of functions and
            parameters, e.g. `[('s', 'et', 4, (0.5, 2.0)), ('p', 'et', 3, (1.0, 3.0))]`
          name : str
            Name of the basis set
          element : str
            Chemical symbol of the element

        Returns:
          out : BasisSet
        '''

        funs = dict()

        for shell, seq, nf, params in functs:
            funs[shell] = dict()
            funs[shell]['e'] = generate_exponents(seq, nf, params)
        functions = OrderedDict(sorted(funs.items(), key=lambda x: SHELLS.index(x[0])))
        bs = cls(name=name, element=element, functions=functions)
        bs.uncontract()
        bs.sort()
        return bs

    @classmethod
    def from_file(cls, fname=None, fmt=None, name=None):
        '''Read and parse a basis set from file and return a BasisSet object

        Args:
          fname : str
            File name
          fmt : str
            Format of the basis set in the file (*molpro*, *gamessus*)
          name : str
            Name of the basis set

        Returns:
          out : BasisSet or dict
            Basisset object parsed from file or dictionary of BasisSet objects
        '''

        if name is None:
            name = os.path.splitext(fname)[0]

        with open(fname, 'r') as fobj:
            basstr = fobj.read()

        return BasisSet.from_str(basstr, fmt=fmt, name=name)

    @classmethod
    def from_str(cls, string, fmt=None, name=None):
        '''
        Parse a basis set from string

        Args:
          string : str
            A string with the basis set
          fmt : str
            Format of the basis set in the file: *molpro*, *gamessus*
          name : str
            Name of the basis set

        Returns:
          out : BasisSet or dict
            Basisset object parsed from string or dictionary of BasisSet
            objects
        '''

        res = parse_basis(string, fmt=fmt)
        # return a BasisSet object if only one basis parsed or dict of BasisSet
        # objects with atomic symbol as key

        out = dict()
        if len(res) == 1:
            atom, fs = res.items()[0]
            return cls(name=name, element=res.keys()[0],
                    functions=OrderedDict(sorted(fs.items(), key=lambda x: SHELLS.index(x[0]))))
        else:
            for atom, fs in res.items():
                out[atom] = cls(name=name, element=atom,
                        functions=OrderedDict(sorted(fs.items(), key=lambda x: SHELLS.index(x[0]))))
            return out

    def __repr__(self):
        keys = ['name', 'element', 'family', 'kind']

        res = "<BasisSet(\n"
        for key in keys:
            res += "\t{k:<20s} = {v}\n".format(k=key, v=getattr(self, key))
        res += self.print_functions()
        res += ")>"
        return res

    def __str__(self):
        keys = ['name', 'element', 'family', 'kind']

        res = ''
        for key in keys:
            res += "{k:<20s} = {v}\n".format(k=key.capitalize(), v=getattr(self, key))
        res += 'Functions:\n'
        res += self.print_functions()
        return res

    def __add__(self, other):
        '''Add functions from another BasisSet object

        Args:
          other : BasisSet
            BasisSet object whose functions will be added to the existing ones

        Returns:
          BasisSet instance with functions from `self` and `other` merged
        '''

        selffuncs = deepcopy(self.functions)
        otherfuncs = deepcopy(other.functions)
        for oshell, ofs in otherfuncs.items():
            if oshell.lower() in selffuncs.keys():
                i = self.functions[oshell]['e'].size
                selffuncs[oshell]['e'] = np.concatenate((selffuncs[oshell]['e'], ofs['e']))
                for cf in ofs['cf']:
                    cf['idx'] = cf['idx'] + i
                    selffuncs[oshell]['cf'].append(cf)
            else:
                self.functions[oshell] = ofs
        return BasisSet(name=self.name, element=self.element, functions=selffuncs)

    def append(self, other):
        '''Append functions from another BasisSet object

        Args:
          other : BasisSet
            BasisSet object whose functions will be added to the existing ones
        '''

        otherfuncs = deepcopy(other.functions)
        for oshell, ofs in otherfuncs.items():
            if oshell.lower() in self.functions.keys():
                i = self.functions[oshell]['e'].size
                self.functions[oshell]['e'] = np.concatenate((self.functions[oshell]['e'], ofs['e']))
                for cf in ofs['cf']:
                    cf['idx'] = cf['idx'] + i
                    self.functions[oshell]['cf'].append(cf)
            else:
                self.functions[oshell] = ofs

    def get_exponents(self, shell=None):
        '''
        Return the exponents of a given shell or if the shell isn't specified
        return all of the available exponents
        '''

        if shell is None:
            return chain.from_iterable([self.functions[k]['e'] for k in self.functions.keys()])
        else:
            return self.functions[shell]['e']

    def uncontract(self, copy=False):
        '''
        Uncontract the basis set. This replaces the contraction coefficients in
        the current object.

        Args:
          copy : bool
            If `True` return an uncontracted copy of the basis set rather than
            uncontracting in place, default is `False`.
        '''

        if copy:
            res = deepcopy(self)
            for shell, fs in res.functions.items():
                fs['cf'] = [np.array([tuple([i, 1.0])], dtype=CFDTYPE) for i, _ in enumerate(fs['e'])]
            return res
        else:
            for shell, fs in self.functions.items():
                fs['cf'] = [np.array([tuple([i, 1.0])], dtype=CFDTYPE) for i, _ in enumerate(fs['e'])]

    def to_cfour(self, comment="", efmt="15.8f", cfmt="15.8f"):
        '''
        Return a string with the basis set in (new) CFOUR format.

        Args:
          comment : str
            comment string
          efmt : str
            string describing output format for the exponents, default: "20.10f"
          cfmt : str
            string describing output format for the contraction coefficients,
            default: "15.8f"

        Returns:
          res : str
            basis set string in Cfour format
        '''

        am, ne, cf = [], [], []
        for shell, shellfs in sorted(self.functions.items(), key=lambda x: SHELLS.index(x[0])):
            am.append(SHELLS.index(shell))
            ne.append(len(shellfs["e"]))
            cf.append(len(shellfs["cf"]))

        res = "\n{e}:{s}\n{c}\n\n".format(e=self.element, s=self.name, c=comment)
        res += "{0:3d}\n".format(len(self.functions.keys()))
        res += "".join(["{0:5d}".format(x) for x in am]) + "\n"
        res += "".join(["{0:5d}".format(x) for x in cf]) + "\n"
        res += "".join(["{0:5d}".format(x) for x in ne]) + "\n"
        res += "\n"

        for shell, fs in self.functions.items():
            for lst in splitlist(fs['e'], 5):
                res += "".join(["{0:>{efmt}}".format(e, efmt=efmt) for e in lst])  + "\n"
            res += "\n"
            # create an array with all the contraction coefficients for a given shell
            cc = self.contraction_matrix(shell)
            for row in cc:
                for splitrow in splitlist(row, 5):
                    res += "{c}".format(c="".join(["{0:{cfmt}}".format(c, cfmt=cfmt) for c in splitrow])) + "\n"
            res += "\n"
        return res

    def to_dalton(self, efmt="20.10f", cfmt="15.8f"):
        '''
        Return a string with the basis set in DALTON format.

        Args:
          efmt : str
            string describing output format for the exponents, default: "20.10f"
          cfmt : str
            string describing output format for the contraction coefficients,
            default: "15.8f"

        Returns:
          res : str
            basis set string in Dalton format
        '''

        res = "! {s}\n".format(s=self.name)
        for shell, fs in self.functions.items():
            res += "! {s} functions\n".format(s=shell)
            res += "{f:1s}{p:>4d}{c:>4d}\n".format(f="F", p=len(fs['e']), c=len(fs['cf']))
            # create an array with all the contraction coefficients for a given shell
            cc = self.contraction_matrix(shell)
            for expt, cfc in zip(fs['e'], cc):
                res += "{e:>{efmt}}{c}".format(e=expt, efmt=efmt, c="".join(["{0:{cfmt}}".format(c, cfmt=cfmt) for c in cfc])) + "\n"
        return res

    def normalze(self):
        '''
        Normalize contraction coefficients for each contracted functions based
        on the primitive overlaps so that the norm is equal to 1.
        '''

        for shell, fs in self.functions.items():
            cc = self.contraction_matrix(shell)
            po = primitive_overlap(SHELLS.index(shell), fs['e'], fs['e'])
            for col in range(cc.shape[1]):
                norm2 = np.dot(cc[:, col], np.dot(po, cc[:, col]))
                fs['cf'][col]['cc'] = cc[fs['cf'][col]['idx'], col]/np.sqrt(norm2)

    def normalzation(self):
        '''
        For each function (contracted) calculate the norm and return a list of
        tuples containing the shell, function index and the norm respectively.
        '''

        out = list()
        for shell, fs in self.functions.items():
            cc = self.contraction_matrix(shell)
            po = primitive_overlap(SHELLS.index(shell), fs['e'], fs['e'])
            for col in range(cc.shape[1]):
                norm2 = np.dot(cc[:, col], np.dot(po, cc[:, col]))
                out.append(tuple([shell, col, norm2]))
        return out

    def contraction_matrix(self, shell):
        '''
        Return the contraction coefficients for a given shell in a matrix form
        with size ne*nc, where ne is the number of exponents and nc is the number
        of contracted functions

        Args:
          shell : str
            shell label, *s*, *p*, *d*, ...

        Returns:
          out : 2D numpy.array
            2D array with contraction coefficients
        '''

        if shell in self.functions.keys():
            fs = self.functions[shell]
            out = np.zeros((fs['e'].size, len(fs['cf'])))
            for i, cf in enumerate(fs['cf']):
                out[cf['idx'], i] = cf['cc']
            return out
        else:
            raise ValueError("shell '{}' is not present in the BasisSet".format(shell))

    def to_gamessus(self, efmt="20.10f", cfmt="15.8f"):
        '''
        Return a string with the basis set in GAMESS(US) format.

        Args:
          efmt : str
            string describing output format for the exponents, default: "20.10f"
          cfmt : str
            string describing output format for the contraction coefficients,
            default: "15.8f"

        Returns:
          res : str
            basis set string in Gamess(US) format
        '''

        res = "{0:s}\n".format(self.element)
        for shell, fs in self.functions.items():
            for contraction in fs["cf"]:
                res += "{s:<1s}{n:>3d}\n".format(s=shell.upper(), n=len(contraction))
                for i, (idx, coeff) in enumerate(contraction, start=1):
                    res += "{i:3d}{e:>{efmt}}{c:>{cfmt}}".format(i=i, e=fs['e'][idx], efmt=efmt, c=coeff, cfmt=cfmt)+ "\n"
        return res + "\n"

    def to_gaussian(self, efmt="20.10f", cfmt="15.8f"):
        '''
        Return a string with the basis set in Gaussian format.

        Args:
          efmt : str
            string describing output format for the exponents, default: "20.10f"
          cfmt : str
            string describing output format for the contraction coefficients,
            default: "15.8f"

        Returns:
          res : str
            basis set string in Gaussian format
        '''

        res = "****\n{e:<7s}0\n".format(e=self.element)
        for shell, fs in self.functions.items():
            for contraction in fs["cf"]:
                res += "{s:<1s}{n:>4d}{i:>7.2f}\n".format(s=shell.upper(), n=len(contraction),i=1.0)
                for idx, coeff in contraction:
                    res += "{e:>{efmt}}{c:>{cfmt}}".format(e=fs['e'][idx], efmt=efmt, c=coeff, cfmt=cfmt)+ "\n"
        return res + "****\n"

    def to_molpro(self, pars=False, efmt="20.10f", cfmt="15.8f"):
        '''
        Return a string with the basis set in MOLPRO format.

        Args:
          pars : bool
            a flag to indicate whether to wrap the basis with "basis={ }" string
          efmt : str
            string describing output format for the exponents, default: "20.10f"
          cfmt : str
            string describing output format for the contraction coefficients,
            default: "15.8f"

        Returns:
          res : str
            basis set string
        '''

        res = ""
        for shell, fs in self.functions.items():
            exps = ", ".join(["{0:>{efmt}}".format(e, efmt=efmt).lstrip() for e in fs['e']])
            res += "{s:>s}, {e:>s}, ".format(s=shell, e=self.element) + exps + '\n'
            for cf in fs['cf']:
                coeffs = ", ".join(["{0:>{cfmt}}".format(cc, cfmt=cfmt).lstrip() for cc in cf['cc']])
                res += "c, {0:d}.{1:d}, ".format(cf['idx'].min()+1, cf['idx'].max()+1) + coeffs + '\n'
        if pars:
            res = 'basis={\n' + res + '}'
        return res

    def to_nwchem(self, efmt="20.10f", cfmt="15.8f"):
        '''
        Return a string with the basis set in NWChem format.

        Args:
          efmt : str
            string describing output format for the exponents, default: "20.10f"
          cfmt : str
            string describing output format for the contraction coefficients,
            default: "15.8f"

        Returns:
          res : str
            basis set string in NwChem format
        '''

        res = 'BASIS "ao basis" PRINT\n'
        for shell, fs in self.functions.items():
            # create an array with all the contraction coefficients for a given shell
            cc = self.contraction_matrix(shell)

            # select columns with more than 1 non zero coefficients
            nonzerocolmask = np.array([np.count_nonzero(col) > 1 for col in cc.T])
            nonzerorowmask = np.array([np.count_nonzero(row) > 0 for row in cc[:, nonzerocolmask]])
            if np.any(nonzerocolmask):
                res += "{e} {s}\n".format(e=self.element, s=shell)
                for expt, cfc in zip(fs['e'][nonzerorowmask], cc[np.ix_(nonzerorowmask, nonzerocolmask)]):
                    res += "{e:>{efmt}}{c}".format(e=expt, efmt=efmt, c="".join(["{0:{cfmt}}".format(c, cfmt=cfmt) for c in cfc])) + "\n"

            if np.any(np.logical_not(nonzerocolmask)):
                for colidx in np.where(np.logical_not(nonzerocolmask))[0]:
                    res += "{e} {s}\n".format(e=self.element, s=shell)
                    nonzerorowmask = np.array([np.count_nonzero(row) > 0 for row in cc[:, colidx]])
                    e = fs['e'][nonzerorowmask][0]
                    c = cc[nonzerorowmask, :][0][colidx]
                    res += "{e:>{efmt}}{c:>{cfmt}}".format(e=e, efmt=efmt, c=c, cfmt=cfmt) + "\n"
        return res + "END\n"

    def print_functions(self, efmt="20.10f", cfmt="15.8f"):
        '''
        Return a string with the basis set.

        Args:
          efmt : str
            string describing output format for the exponents, default: "20.10f"
          cfmt : str
            string describing output format for the contraction coefficients,
            default: "15.8f"

        Returns:
          res : str
            basis set string
        '''

        res = ''
        for shell, fs in self.functions.items():
            # create an array with all the contraction coefficients for a given shell
            res += "\n" + "{s} shell".format(s=shell).center(40, '=') + "\n"
            cc = self.contraction_matrix(shell)
            count = 0
            # select columns with more than 1 non zero coefficients
            nonzerocolmask = np.array([np.count_nonzero(col) > 1 for col in cc.T])
            nonzerorowmask = np.array([np.count_nonzero(row) > 0 for row in cc[:, nonzerocolmask]])
            if np.any(nonzerocolmask):
                res += 'Contracted:\n'
                for expt, cfc in zip(fs['e'][nonzerorowmask], cc[np.ix_(nonzerorowmask, nonzerocolmask)]):
                    count += 1
                    res += "{i:5d}{e:>{efmt}}{c}".format(i=count, e=expt, efmt=efmt, c="".join(["{0:{cfmt}}".format(c, cfmt=cfmt) for c in cfc])) + "\n"

            if np.any(np.logical_not(nonzerocolmask)):
                res += 'Uncontracted:\n'
                for colidx in np.where(np.logical_not(nonzerocolmask))[0]:
                    count += 1
                    nonzerorowmask = np.array([np.count_nonzero(row) > 0 for row in cc[:, colidx]])
                    e = fs['e'][nonzerorowmask][0]
                    c = cc[nonzerorowmask, :][0][colidx]
                    res += "{i:5d}{e:>{efmt}}{c:>{cfmt}}".format(i=count, e=e, efmt=efmt, c=c, cfmt=cfmt) + "\n"
        return res

    def to_pickle(self, fname=None):
        '''Save the basis set in pickle format under the filename `fname`

        Args:
          fname : str
            File name
        '''

        if fname is None:
            fname = self.name.strip().replace(' ', '_') + '.bas.pkl'

        with open(fname, 'wb') as fbas:
            pickle.dump(self, fbas)

    def nf(self, spherical=True):
        '''
        Calculate the number of basis functions

        Args:
          spherical : bool
            flag indicating if spherical or cartesian functions should be used,
            default: True

        Returns:
          res : int
            number of basis functions
        '''

        if spherical:
            return sum(nspherical(get_l(shell))*len(fs['cf']) for shell, fs in self.functions.items())
        else:
            return sum(ncartesian(get_l(shell))*len(fs['cf']) for shell, fs in self.functions.items())

    def nprimitive(self, spherical=True):
        '''
        Return the number of primitive functions assuming sphrical or cartesian
        gaussians.

        Args:
          spherical : bool
            A flag to select either spherical or cartesian gaussians

        Returns:
          out : int
            Number of primitive function in the basis set
        '''

        if spherical:
            # calculate the number of spherical components per shell
            ncomp = [nspherical(get_l(shell)) for shell in self.functions.keys()]
        else:
            # calculate the number of cartesian components per shell
            ncomp = [ncartesian(get_l(shell)) for shell in self.functions.keys()]

        return sum([prim*nc for prim, nc in zip(self.primitives_per_shell(), ncomp)])


    def contraction_scheme(self):
        '''
        Return a string describing the contraction scheme.
        '''

        cs, ec = [], []
        for shell, fs in self.functions.items():
            cs.append((shell, len(fs['e']), len(fs['cf'])))
            ec.append([len(cfs) for cfs in fs['cf']])
        return "({p:s}) -> [{c:s}] : {{{ec}}}".format(
                p="".join(["{0:d}{1:s}".format(c[1], c[0]) for c in cs]),
                c="".join(["{0:d}{1:s}".format(c[2], c[0]) for c in cs]),
                ec="/".join([" ".join(["{0:d}".format(c) for c in x]) for x in ec]))

    def primitives_per_shell(self):
        '''
        Calculate how many primitive functions are in each shell.

        Returns:
          out : list of ints
        '''

        return [len(f['e']) for s, f in self.functions.items()]

    def contractions_per_shell(self):
        '''
        Calculate how many contracted functions are in each shell.

        Returns:
          out : list of ints
        '''

        return [len(f['cf']) for s, f in self.functions.items()]

    def primitives_per_contraction(self):
        '''
        Calculate how many primities are used in each contracted function.

        Returns:
          out : list of ints
        '''
        return [[len(cc) for cc in f['cf']] for s, f in self.functions.items()]

    def contraction_type(self):
        '''
        Try to determine the contraction type: segmented, general, uncontracted, unknown.
        '''

        pps = self.primitives_per_shell()
        ppc = self.primitives_per_contraction()

        if any(x > 1 for x in pps):
            if all(all(x == 1 for x in shell) for shell in ppc):
                return "uncontracted"
            elif all(all(pinc == np for pinc in shell) for np, shell in zip(pps, ppc)):
                return "general"
            else:
                return "unknown"

        # one function per shell case
        if all(all(x == 1 for x in shell) for shell in ppc):
            return "uncontracted 1fps"

    def sort(self, reverse=False):
        '''
        Sort shells in the order of increasing angular momentum and for each
        shell sort the exponents.

        Args:
          reverse : bool
            If `False` sort the exponents in each shell in the descending order
            (default), else sort exponents in ascending order
        '''

        self.functions = OrderedDict(sorted(self.functions.items(), key=lambda x: SHELLS.index(x[0])))

        for shell, fs in self.functions.items():
            if reverse:
                idx = np.argsort(fs['e'])
            else:
                idx = np.argsort(fs['e'])[::-1]

            # actually sort the exponents and coefficients
            fs['e'] = fs['e'][idx]
            for cf in fs['cf']:
                cf['idx'] = idx[cf['idx']]
                cf.sort(order='idx')

    def partial_wave_expand(self):
        """
        From a given basis set with shells spdf... return a list of basis sets that
        are subsets of the entered basis set with increasing angular momentum functions
        included [s, sp, spd, spdf, ...]
        """
        res = list()
        shells = self.functions.keys()
        for i in range(1, len(shells)+1):
            bscopy = copy.copy(self)
            bscopy.functions = {k:v for k, v in self.functions.items() if k in shells[:i]}
            res.append(bscopy)
        return res

    def completeness_profile(self, zetas):
        '''
        Calculate the completeness profile of each shell of the basis set

        Args:
          zetas : numpy.array
            Scaning exponents

        Returns:
          out : numpy.array
            numpy array with values of the profile (shells, zetas)
        '''

        out = np.zeros((zetas.size, len(self.functions.keys())))

        for i, (shell, fs) in enumerate(self.functions.items()):
            cc = self.contraction_matrix(shell)
            # calculate the overlap matrix
            S = self.shell_overlap(shell)
            # calculate the inverse square root of S
            X = inv(sqrtm(S))
            SO = primitive_overlap(SHELLS.index(shell), fs['e'], zetas)
            J = np.dot(np.dot(cc, X).T, SO)
            Y = np.sum(np.square(J), axis=0)
            out[:, i] = Y
        return out

    def shell_overlap(self, shell):
        '''
        Calculate the overlap integrals for a given shell

        Args:
          shell : str
            Shell

        Returns:
          out : numpy.array
            Overlap integral matrix
        '''

        exps = self.functions[shell]['e']
        S = primitive_overlap(SHELLS.index(shell), exps, exps)
        cc = self.contraction_matrix(shell)
        return np.dot(cc.T, np.dot(S, cc))

def get_l(shell):
    '''Return the angular momentum value of a given shell'''

    return SHELLS.index(shell)

def merge(first, other):
    '''
    Merge functions from two BasisSet objects

    Args:
      first : BasisSet
        First `BasisSet` object to merge
      other : BasisSet
        Second `BasisSet` object to merge

    Returns:
      our : BasisSet
        BasisSet instance with functions from `first` and `other` merged
    '''

    selffuncs = deepcopy(first.functions)
    otherfuncs = deepcopy(other.functions)
    for oshell, ofs in otherfuncs.items():
        if oshell.lower() in selffuncs.keys():
            exps, idxs, idxo = merge_exponents(selffuncs[oshell]['e'], ofs['e'])
            selffuncs[oshell]['e'] = exps
            for cf in selffuncs[oshell]['cf']:
                cf['idx'] = idxs[cf['idx']]
            for cf in ofs['cf']:
                cf['idx'] = idxo[cf['idx']]
            selffuncs[oshell]['cf'].extend(ofs['cf'])
        else:
            selffuncs[oshell] = ofs
    return BasisSet(name=first.name, element=first.element, functions=selffuncs)

def primitive_overlap(l, a, b):
    '''
    Calculate the overlap integrals for a given shell `l` and two sets of exponents

    .. math::

       S(\zeta_i, \zeta_j) = \\frac{2^{l + \\frac{3}{2}} \left(\zeta_{1} \zeta_{2}\\right)^{\\frac{l}{2} + \\frac{3}{4}}}{\left(\zeta_{1} + \zeta_{2}\\right)^{l + \\frac{3}{2}}}

    Args:
      l : int
        Angular momentum quantum number of the shell
      a (M,) : numpy.array
        First vector of exponents
      b (N,) : numpy.array
        Second vector of exponents

    Returns:
      out (M, N) : numpy.array
        Overlap integrals
    '''

    return np.power(2.0*np.sqrt(np.outer(a, b))/np.add.outer(a, b), l + 1.5)

def nspherical(l):
    '''
    Calculate the number of spherical components of a function with a given angular
    momentum value *l*.
    '''
    return 2*l+1

def ncartesian(l):
    '''
    Calculate the number of cartesian components of a function with a given angular
    momentum value *l*.
    '''
    return int((l+1)*(l+2)/2)

def zetas2legendre(zetas, kmax):
    '''
    From a set of exponents (`zetas`), using least square fit calculate the
    expansion coefficients into the legendre polynomials of the order `kmax`.

    Args:
      kmax : int
        length of the legendre expansion
      zetas : list
        list of exponents (floats) to be fitted

    Returns:
      coeff : numpy.array
        numpy array of legendre expansion coeffcients of length kmax
    '''

    # exponents should be sorted in the ascending order
    zetas = sorted(zetas)
    c = np.asarray([1.0]*kmax, dtype=float)

    leg = np.polynomial.Legendre(c)
    a = np.zeros((len(zetas), kmax))

    for j in range(len(zetas)):
        for k in range(kmax):
            arg = (2.0*(j+1.0)-2.0)/(len(zetas)-1.0)-1.0
            a[j, k] = leg.basis(k)(arg)

    return np.linalg.lstsq(a, np.log(zetas))[0]

def zetas2eventemp(zetas):
    '''
    From a set of exponents (`zetas`), using least square fit calculate the
    approximate $\\alpha$ and $\\beta$ parameters from the even tempered
    expansion.

    Args:
      zetas : list
        list of exponents (floats) to be fitted

    Returns:
      coeff : numpy.array
        numpy array of legendre expansion coeffcients of length kmax
    '''

    return min(zetas), np.power(max(zetas)/min(zetas), 1.0/(len(zetas)-1))

def generate_exponents(formula, nf, params):
    '''
    Generate a sequence of exponents from a specified formula

    Args:
      formula : str
        name of the sequence from which the exponents are generated, one of:
          - *even*, *eventemp*, *even tempered*
          - *well*, *welltemp*, *well tempered*
          - *legendre*
          - *exp*, *exponents*
      nf : int
        number of exponents
      params : tuple
        parameters for the formula

    Returns:
      numpy.array with exponents
    '''

    if formula.lower() in ["et", "even", "eventemp", "even tempered"]:
        return eventemp(nf, params)
    elif formula.lower() in ["wt", "well", "welltemp", "well tempered"]:
        return welltemp(nf, params)
    elif formula.lower() in ["le", "legendre"]:
        return legendre(nf, params)
    elif formula.lower() in ["exp", "exponents"]:
        return np.sort(np.array(params))[::-1]
    else:
        raise ValueError('unknown sequence name: {}'.format(formula))

def eventemp(nf, params):
    '''
    Generate a sequence of nf even tempered exponents accodring to
    the even tempered formula

    .. math::
       \\zeta_i = \\alpha \cdot \\beta^{i-1}

    Args:
      nf : int
        number fo exponents to generate
      params : tuple of floats
        alpha and beta parameters
    Returns:
      res : numpy array
        array of generated exponents (floats)
    '''

    if not isinstance(nf, int):
        raise TypeError('"nf" variable should be of "int" type, got: {}'.format(type(nf)))
    if len(params) !=  2:
        raise ValueError('"params" tuple should have exactly 2 entries, got {}'.format(len(params)))

    alpha, beta = params
    zetas = alpha * np.power(beta, np.arange(nf))
    return zetas[::-1]

def welltemp(nf, params):
    '''
    Generate a sequence of nf well tempered exponents accodring to
    the well tempered fromula

    .. math::

       \\zeta_i = \\alpha \cdot \\beta^{i-1} \cdot \\left[1 + \\gamma \cdot \\left(\\frac{i}{N}\\right)^{\delta}\\right]

    Args:
      nf : int
        number fo exponents to generate
      params : tuple of floats
        alpha, beta, gamma and delta parameters

    Returns:
      res : numpy.array
        array of generated exponents (floats)
    '''
    if not isinstance(nf, int):
        raise TypeError('"nf" variable should be of "int" type, got: {}'.format(type(nf)))
    if len(params) !=  4:
        raise ValueError('"params" tuple should have exactly 4 entries, got {}'.format(len(params)))

    alpha, beta, gamma, delta = params
    zetas = alpha*np.power(beta, np.arange(nf))*(1+gamma*np.power(np.arange(1, nf+1)/nf, delta))
    zetas.sort()
    return zetas[::-1]

def legendre(nf, coeffs):
    '''
    Generate a sequence of nf exponents from expansion in the orthonormal
    legendre polynomials as described in: Peterson, G. A. et.al J. Chem. Phys.,
    Vol. 118, No. 3 (2003), eq. (7).

    .. math::
       \ln \\zeta_i = \\sum^{k_{\max}}_{k=0} A_k P_k \\left(\\frac{2j-2}{N-1}-1\\right)

    Args:
      nf : int
        number fo exponents to generate
      params : tuple of floats
        polynomial coefficients (expansion parameters)

    Returns:
      res : numpy array
        array of generated exponents (floats)
    '''

    if not isinstance(nf, int):
        raise TypeError('"nf" variable should be of "int" type, got: {}'.format(type(nf)))
    if len(coeffs) <  1:
        raise ValueError('"coeffs" tuple should have at least 1 entry, got {}'.format(len(coeffs)))

    # special case for one function
    if len(coeffs) == 1:
        return [np.exp(coeffs[0])]

    poly = np.polynomial.legendre.Legendre(coeffs)
    zetas = [poly(((2.0*(i+1.0)-2.0)/(nf-1.0))-1.0) for i in range(nf)]
    return np.exp(zetas[::-1])

def splitlist(l, n):
    if len(l) % n == 0:
        splits = len(l)//n
    elif len(l) % n != 0 and len(l) > n:
        splits = len(l)//n+1
    else:
        splits = 1

    for i in range(splits):
        yield l[n*i:n*i+n]

def sliceinto(l, s):
    '''slice a list into chunks with sizes defined in s'''

    if len(l) != sum(s):
        raise ValueError('cannot slice list, size mismatch {}!={}'.format(len(l), sum(s)))

    si = [sum(s[:i]) for i in range(0, len(s))]
    return [l[i:i+s] for i, s in zip(si, sizes)]

def have_equal_floats(a, b):
    '''
    Check if the arrays a and b have equal items

    Pairwise compare all the item from a with all the items form b

    Args:
      a : numpy.array
      b : numpy.array

    Returns:
      res : bool
        - True if there are common items in `a` and `b`
        - False if there are no common items
    '''

    x = np.array([x for x in product(a, b)])
    return np.any(np.isclose(x[:, 0], x[:, 1]))


def xyzlist(l):
    '''
    Generate an array of :math:`l_x`, :math:`l_y`, :math:`l_z` components of
    cartesian gaussian with a given angular momentum value in canonical order.

    For exampe:
      - ``l = 0`` generates the result ``array([[0, 0, 0]])``
      - ``l = 1`` generates ``array([[1, 0, 0], [0, 1, 0], [0, 0, 1])``
      - etc.

    The functions are coded by triples of powers of *x*, *y*, *z*, namely
    ``[1, 2, 3]`` corresponds to :math:`xy^{2}z^{3}`.

    Args:
      l : int
        Angular momentum value

    Returns:
      out ((l+1)*(l+2)/2, 3) : numpy.array
        Array of monomial powers, where columns correspond to *x*, *y*, and *z*
        respectively and rows correspond to functions.
    '''

    ncart = (l + 1)*(l + 2) // 2
    out = np.zeros((ncart, 3), dtype=np.int32)
    index = 0
    for i in range(l + 1):
        for j in range(i + 1):
            out[index,:] = [l - i, i - j, j]
            index += 1
    return out

def zlmtoxyz(l):
    '''
    Generates the expansion coefficients of the real spherical harmonics in
    terms of products of cartesian components. Method based on [1]_

    .. [1] Schlegel, H. B., & Frisch, M. J. (1995). "Transformation between
       Cartesian and pure spherical harmonic Gaussians". International Journal
       of Quantum Chemistry, 54(2), 83–87. `doi:10.1002/qua.560540202 <http:www.dx.doi.org/10.1002/qua.560540202>`_

    Args:
      l : int
        Angular momentum value

    Returns:
      out :math:`((l+1)*(l+2)/2, 2*l + 1)` : numpy.array
        Expansion coefficients of real spherical harmonics in terms of cartesian
        gaussians
    '''

    ncart = (l + 1)*(l + 2)//2
    nspher = 2*l + 1

    out = np.zeros((ncart, nspher))

    cartc = xyzlist(l)

    m = 0
    for icart in range(ncart):
        kx = cartc[icart, 0]
        ky = cartc[icart, 1]
        kz = cartc[icart, 2]
        if kx + ky % 2 == 1:
            continue
        j = (kx + ky) / 2
        tmpc = factorial2(2*kx - 1)*factorial2(2*ky - 1)*factorial2(2*kz - 1)/factorial2(2*l - 1)
        tmpc = np.sqrt(tmpc)/factorial(j)
        for i in range(l//2 + 1):
            if j > i:
                continue
            if kx % 2 == 1:
                continue
            k = kx//2
            if k > j:
                continue
            tmpi = tmpc*factorial2(2*l - 2*i - 1)/factorial(l - 2*i)/factorial(i - j)/2.0**i
            if i % 2 == 1:
                tmpi = -tmpi
            out[icart, l] += binom(j, k)*tmpi

    for m in range(1, l + 1):
        for icart in range(ncart):
            kx = cartc[icart, 0]
            ky = cartc[icart, 1]
            kz = cartc[icart, 2]
            jj = kx + ky - m
            if jj % 2 == 1:
                continue
            if  jj < 0:
                continue
            j = jj//2
            tmpc = factorial2(2*kx - 1)*factorial2(2*ky - 1)*factorial2(2*kz - 1)/factorial2(2*l - 1)
            tmpc = np.sqrt(2.0*tmpc*factorial(l-m)/factorial(l+m))/factorial(j)
            for i in range((l-m)//2 + 1):
                if j > i:
                    continue
                tmpi = tmpc*factorial2(2*l - 2*i-1)/factorial(l-m-2*i)/factorial(i - j)/2.0**i
                if i % 2 == 1:
                    tmpi = -tmpi
                for k in range(j + 1):
                    kk = kx - 2*k
                    if kk < 0 or kk > m:
                        continue
                    tmpk = tmpi*binom(j, k)*binom(m, kk)
                    kkk = m - kk
                    if kkk % 4 == 0:
                        out[icart, m + l] += tmpk
                    elif kkk % 4 == 1:
                        out[icart, -m + l] += tmpk
                    elif kkk % 4 == 2:
                        out[icart, m + l] -= tmpk
                    else:
                        out[icart, -m + l] -= tmpk
    return out

