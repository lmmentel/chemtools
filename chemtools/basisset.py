#import pymongo

from __future__ import division, print_function

from copy import deepcopy
import numpy as np
from itertools import chain
from collections import OrderedDict
from scipy.linalg import sqrtm, inv
import pickle
import re
import os
from elements import element

SHELLS = ["s", "p", "d", "f", "g", "h", "i", "j", "k"]
CFDTYPE = [('idx', np.int32), ('cc', np.float64)]

def read_pickle(fname):
    '''Read a pickled BasisSet object from file

    Args:
      fname : str
        File name containing the BasisSet
    '''

    with open(fname, 'rb') as fil:
        return pickle.load(fil)

def read_dict(dct):
    '''
    Initialize the BasisSet object from a general dictionary with an
    assumed structure of "functions" entry. This method can be used to
    initialize the BasisSet object directly from mongoDB database.

    Args:
      dct : dict
        Dictionary with the basis set data, required keys are:
          - name
          - element
          - functions
    '''
    required = set(['name', 'element', 'functions'])
    if isinstance(dct, dict):
        if not required.issubset(set(d.keys())):
            missing = list(required.difference(set(dct.keys())))
            raise ValueError('Following keys are required: {0}'.format(", ".join(missing)))
        bs = BasisSet()
        for key, val in dct.items():
            if key == 'functions':
                setattr(bs, key, OrderedDict(sorted(val.items(), key=lambda x: SHELLS.index(x[0]))))
            else:
                setattr(bs, key, val)
    return bs

class BasisSet:
    '''
    Basis set module supporting basic operation on basis sets and can be used
    as a API for mongoDB basis set repository.
    '''


    def __init__(self, name, element, family=None, kind=None,
                 functions=None, params=None):
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
            Parameters for generating the functins according to the model
        '''

        self.name = name
        self.element = element
        self.family = family
        self.kind = kind
        self.params = params
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
        return bs

    @classmethod
    def from_file(cls, fname=None, fmt=None, name=None):
        '''Read and parse a basis set from file and return a BasisSet object

        Args:
          fname : str
            File name
          fmt : str
            Format of the basis set in the file: *molpro*, *gamessus*
         name : str
           Name of the basis set
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
        '''

        formats = ['molpro', 'gamessus']

        if fmt.lower() not in formats:
            raise ValueError("<fmt> should be one of: {}".format(", ".join(formats)))

        if fmt == 'molpro':
            bsd = parse_molpro_basis(string)
        elif fmt == 'gamessus':
            bsd = parse_gamessus_basis(string)

        # return a BasisSet object if only one basis parsed or dict of BasisSet
        # objects with atomic symbol as key

        res = dict()
        if len(bsd) == 1:
            atom, fs = bsd.items()[0]
            return cls(name=name, element=bsd.keys()[0],
                    functions=OrderedDict(sorted(fs.items(), key=lambda x: SHELLS.index(x[0]))))
        else:
            for atom, fs in bsd.items():
                res[atom] = cls(name=name, element=atom,
                        functions=OrderedDict(sorted(fs.items(), key=lambda x: SHELLS.index(x[0]))))
            return res

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

        res = ""
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
            return sum(nspherical(SHELLS.index(shell))*len(fs['cf']) for shell, fs in self.functions.items())
        else:
            return sum(ncartesian(SHELLS.index(shell))*len(fs['cf']) for shell, fs in self.functions.items())

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

def merge(first, other):
    '''Merge functions from two BasisSet objects

    Args:
      first : BasisSet
      other : BasisSet
        BasisSet object whose functions will be added to the existing ones

    Returns:
        BasisSet instance with functions from `self` and `other` merged
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

def merge_exponents(a, b):
    '''
    Concatenate the arrays `a` and `b` using only the unique items from both arrays

    Args:
      a : numpy.array
      b : numpy.array

    Returns:
      res : 3-tuple of numpy arrays
        - res[0] sorted union of `a` and `b`
        - res[1] indices of `a` items in res[0]
        - res[2] indices of `b` items in res[0]
    '''

    u = np.union1d(a, b)[::-1]
    idxa = np.where(np.in1d(u, a))[0]
    idxb = np.where(np.in1d(u, b))[0]
    return u, idxa, idxb

def parse_molpro_basis(string):
    '''
    Parse basis set from a string in Molpro format.
    '''

    bas_re = re.compile(r'basis\s*=\s*\{(.*?)\}', flags=re.DOTALL|re.IGNORECASE)

    m = bas_re.search(string)
    if m:
        lines = m.group(1).split("\n")
    else:
        raise ValueError("basis string not found")

    start = []
    for i, line in enumerate(lines):
        if line.split(",")[0].lower() in SHELLS:
            start.append(i)
    if len(start) == 0:
        return None

    startstop = []
    for i in range(len(start)-1):
        startstop.append((start[i], start[i+1]))
    startstop.append((start[-1], len(lines)))

    bs = {}
    for i in startstop:
        at_symbol, shell = parse_molpro_shell(lines[i[0]], lines[i[0]+1:i[1]])
        if at_symbol in bs.keys():
            bs[at_symbol] = dict(list(bs[at_symbol].items()) + list(shell.items()))
        else:
            bs[at_symbol] = shell
    return bs

def parse_molpro_shell(expsline, coeffs):
    '''
    Parse functions of one shell in molpro format.
    '''

    real = re.compile(r'[dD]')
    # remove empty strings and whitespace line breaks and tabs
    coeffs = [x.strip() for x in coeffs if x.strip() not in ['', '\n', '\t']]

    fs = {}

    shell  = expsline.split(",")[0].lower()
    at_symbol = expsline.split(",")[1].strip().capitalize()
    exps   = np.array([float(real.sub('E', x)) for x in expsline.rstrip(";").split(",")[2:]])

    fs[shell] = {'e' : exps, 'cf' : []}
    if len(coeffs) != 0:
        for line in coeffs:
            lsp = line.rstrip(";").split(",")
            if lsp[0] == "c":
                i, j = [int(x) for x in lsp[1].split(".")]
                coeffs = [float(real.sub('E', x)) for x in lsp[2:]]
                fs[shell]['cf'].append(np.array(list(zip(list(range(i-1, j)), coeffs)), dtype=CFDTYPE))
    else:
        for i in range(len(exps)):
            fs[shell]['cc'].append(np.array([tuple([i, 1.0])], dtype=CFDTYPE))
    return at_symbol, fs

def parse_ecp(ecpstring):

    ecp_re = re.compile(r'\!\s*Effective core Potentials.*-{25}\s*\n(.*?)\n\s*\n', flags=re.DOTALL)

    lines = ecpstring.split("\n")

    start = []
    for i, line in enumerate(lines):
        if line.split(",")[0].lower() == 'ecp':
            start.append(i)

    if len(start) == 0:
        return None

    startstop = []
    for i in range(len(start)-1):
        startstop.append((start[i], start[i+1]))
    startstop.append((start[-1], len(lines)))

    ecp = {}
    for i in startstop:
        ecp = dict(list(ecp.items()) + list(parse_coeffs(lines[i[0] : i[1]]).items()))
    return ecp

def parse_coeffs(lines):

    firstl = lines[0].replace(';', '').split(',')
    element = firstl[1].strip().capitalize()
    nele = firstl[2]
    lmax = firstl[3]

    liter = iter(x for x in lines[1:] if x != '')

    ecp = {element : {"nele" : nele, "lmax" : lmax, "shells" : []}}

    while True:
        try:
            temp = next(liter)
        except StopIteration as err:
            break
        nlmax = int(temp.split(";")[0])
        comment = temp.split(";")[1].replace("!", "")
        tt = {'comment' : comment, 'parameters' : []}
        for i in range(nlmax):
            param = next(liter).replace(";", "").split(",")
            tt['parameters'].append({'m' : float(param[0]), 'gamma' : float(param[1]),'c' : float(param[2])})
        ecp[element]['shells'].append(tt)
    return ecp

def parse_gamessus_basis(string):
    '''
    Parse the basis set into a list of dictionaries from a string in
    gamess format.
    '''

    bas_re = re.compile(r'\$DATA\n(.*?)\$END', flags=re.DOTALL|re.IGNORECASE)

    m = bas_re.search(string)
    if m:
        basisstr = m.group(1)
    else:
        raise ValueError("basis string not found")

    pat = re.compile(r'^\s*(?P<shell>[SPDFGHIspdfghi])\s*(?P<nf>[1-9]+)')
    res = dict()

    for item in basisstr.split('\n\n'):
        if len(item) > 0:
            atom, basis = item.split('\n', 1)
            if atom != "":
                functions = dict()
                elem = element(atom.capitalize())
                bslines = basis.split("\n")
                for i, line in enumerate(bslines):
                    match = pat.search(line)
                    if match:
                        shell, nf = match.group("shell").lower(), match.group("nf")
                        exps, indxs, coeffs = parse_gamessus_function(bslines[i+1:i+int(nf)+1])
                        if shell in functions.keys():
                            sexp, idxs, idxo = merge_exponents(functions[shell]['e'], exps)
                            functions[shell]['e'] = sexp
                            for cf in functions[shell]['cf']:
                                cf['idx'] = idxs[cf['idx']]
                            newcf = np.array(zip(idxo[indxs-1], coeffs), dtype=CFDTYPE)
                            functions[shell]['cf'].append(newcf)
                        else:
                            functions[shell] = dict()
                            functions[shell]['cf'] = list()
                            functions[shell]['e'] = exps
                            functions[shell]['cf'].append(np.array(zip(indxs-1, coeffs), dtype=CFDTYPE))
                res[elem.symbol] = functions
    return res

def parse_gamessus_function(los):
    '''
    Parse a basis set function information from list of strings into
    three lists containg: exponents, indices, coefficients.

    Remeber that python doesn't recognise the `1.0d-3` format where `d` or
    `D` is used to the regex subsitution has to take care of that.
    '''

    real = re.compile(r'[dD]')

    indxs = np.array([int(item.split()[0]) for item in los], dtype=np.int32)
    exps = np.array([float(real.sub('E', item.split()[1])) for item in los], dtype=np.float64)
    coeffs = np.array([float(real.sub('E', item.split()[2])) for item in los], dtype=np.float64)

    return (exps, indxs, coeffs)


