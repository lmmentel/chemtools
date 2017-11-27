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

import json
import os
import pickle
import re
from collections import OrderedDict
from copy import copy, deepcopy
from itertools import chain
import numpy as np
from scipy.linalg import sqrtm, inv
from scipy.special import factorial, factorial2, binom
from chemtools.basisparse import (parse_basis, merge_exponents, CFDTYPE, get_l,
                                  NumpyEncoder)


class BasisSet(object):
    '''
    Basis set module supporting basic operation on basis sets.

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

    def __init__(self, name, element, family=None, kind=None,
                 functions=None, info=None):

        self.name = name
        self.element = element
        self.family = family
        self.kind = kind
        self.functions = functions
        self.info = info

    def __add__(self, other):
        '''Add functions from another BasisSet object

        Args:
          other : BasisSet
            BasisSet object whose functions will be added to the existing ones

        Returns:
          BasisSet instance with functions from `self` and `other` merged
        '''

        newf = deepcopy(self.functions)
        for oshell, ofs in other.functions.items():
            if oshell.lower() in self.functions.keys():
                i = self.functions[oshell]['e'].size
                newf[oshell]['e'] = np.concatenate((self.functions[oshell]['e'],
                                                    ofs['e']))
                newcf = [np.copy(cf) for cf in ofs['cf']]
                for cf in newcf:
                    cf['idx'] = cf['idx'] + i
                newf[oshell]['cf'].extend(newcf)
            else:
                newf[oshell] = deepcopy(ofs)

        return BasisSet(name=self.name, element=self.element, functions=newf)

    @staticmethod
    def from_pickle(fname, **kwargs):
        '''Read a pickled BasisSet object from a pickle file

        Args:
            fname : str
                File name containing the BasisSet
            kwargs : dict
                Extra arguments for the `pickle.load` method

        Raises:
            UnicodeDecodeError
                When you try to read a python2 pickle using python3, to fix
                that use `encoding='latin1'` option
        '''

        with open(fname, 'rb') as fil:
            return pickle.load(fil, **kwargs)

    @classmethod
    def from_optpars(cls, x0, funs=None, name=None, element=None,
                     explogs=False):
        '''
        Return a basis set object generated from a sequence based on the
        specified arguments.

        Args:
          x0 : list or numpy.array
            Parameters to generate the basis set as a continuous list or array

          funs : list of tuples
            A list of tuple specifying the shell type, number of functions and
            parameters, e.g. `[('s', 'et', 4, (0.5, 2.0)), ('p', 'et', 3,
            (1.0, 3.0))]`

          name : str
            Name of the basis set

          element : str
            Chemical symbol of the element

        Returns:
          out : BasisSet
        '''

        functions = dict()
        ni = 0
        nt = 0
        for shell, seq, nf, params in funs:
            functions[shell] = dict()
            if seq in ["exp", "exponents"]:
                nt += nf
                if explogs:
                    functions[shell]['e'] = generate_exponents(seq, nf, np.exp(x0[ni: nt]))
                else:
                    functions[shell]['e'] = generate_exponents(seq, nf, x0[ni: nt])
                ni += nf
            else:
                # TODO:
                # params shouldn't be here since it can have a dummy,
                # real values are taken from x0
                nt += len(params)
                functions[shell]['e'] = generate_exponents(seq, nf, x0[ni:nt])
                ni += len(params)
        functions = OrderedDict(sorted(functions.items(),
                                key=lambda x: get_l(x[0])))
        bs = cls(name=name, element=element, functions=functions)
        bs.uncontract()
        return bs

    @classmethod
    def from_sequence(cls, funs=None, name=None, element=None):
        '''
        Return a basis set object generated from a sequence based on the
        specified arguments.

        Args:
          funs : list of tuples
            A list of tuple specifying the shell type, number of functions and
            parameters, e.g. `[('s', 'et', 4, (0.5, 2.0)), ('p', 'et', 3,
            (1.0, 3.0))]`

          name : str
            Name of the basis set

          element : str
            Chemical symbol of the element

        Returns:
          out : BasisSet
        '''

        functions = dict()

        for shell, seq, nf, params in funs:
            functions[shell] = dict()
            functions[shell]['e'] = generate_exponents(seq, nf, params)
        functions = OrderedDict(sorted(functions.items(),
                                key=lambda x: get_l(x[0])))
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
            atom, fs = list(res.items())[0]
            return cls(name=name, element=list(res.keys())[0],
                       functions=OrderedDict(sorted(fs.items(),
                                             key=lambda x: get_l(x[0]))))
        else:
            for atom, fs in res.items():
                out[atom] = cls(name=name, element=atom,
                                functions=OrderedDict(sorted(fs.items(),
                                                      key=lambda x: get_l(x[0]))))
            return out

    @classmethod
    def from_json(cls, jsonstring, **kwargs):
        '''
        Instantiate the `BasisSet` object from a JSON string

        Args:
            jsonstring: str
                A JSON serialized string with the basis set
        '''

        def json_numpy_obj_hook(dct):
            '''
            Decodes a previously encoded numpy ndarray with proper shape
            and dtype.

            Args:
                dct: (dict) json encoded ndarray

            Returns:
                (ndarray) if input was an encoded ndarray
            '''

            if isinstance(dct, dict) and 'dtype' in dct:
                if dct['dtype'] == 'CFDTYPE':
                    return np.array([tuple(r) for r in dct['data']],
                                    dtype=CFDTYPE)
                else:
                    return np.array(dct['data'], dtype=dct['dtype'])
            return dct

        loaded = json.loads(jsonstring, object_hook=json_numpy_obj_hook,
                            **kwargs)

        return cls(**loaded)

    def append(self, other):
        '''Append functions from another BasisSet object

        Args:
          other : BasisSet
            BasisSet object whose functions will be added to the existing ones
        '''

        for oshell, ofs in other.functions.items():
            if oshell.lower() in self.functions.keys():
                i = self.functions[oshell]['e'].size
                self.functions[oshell]['e'] = np.concatenate((self.functions[oshell]['e'], ofs['e']))
                newcf = [np.copy(cf) for cf in ofs['cf']]
                for cf in newcf:
                    cf['idx'] = cf['idx'] + i
                self.functions[oshell]['cf'].extend(newcf)
            else:
                self.functions[oshell] = deepcopy(ofs)

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
            SO = primitive_overlap(get_l(shell), fs['e'], zetas)
            J = np.dot(np.dot(cc, X).T, SO)
            Y = np.sum(np.square(J), axis=0)
            out[:, i] = Y
        return out

    def contraction_matrix(self, shell):
        '''
        Return the contraction coefficients for a given shell in a
        matrix form with size `ne * nc`, where `ne` is the number of
        exponents and `nc` is the number of contracted functions

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

    def contractions_per_shell(self):
        '''
        Calculate how many contracted functions are in each shell.

        Returns:
          out : list of ints
        '''

        return [len(f['cf']) for s, f in self.functions.items()]

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

    def contraction_type(self):
        '''
        Try to determine the contraction type: segmented, general,
        uncontracted, unknown.
        '''

        pps = self.primitives_per_shell()
        ppc = self.primitives_per_contraction()

        if any(x > 1 for x in pps):
            if all(all(x == 1 for x in shell) for shell in ppc):
                return "uncontracted"
            elif all(all(pinc == np for pinc in shell)
                     for np, shell in zip(pps, ppc)):
                return "general"
            else:
                return "unknown"

        # one function per shell case
        if all(all(x == 1 for x in shell) for shell in ppc):
            return "uncontracted 1fps"

    def get_exponents(self, asdict=True):
        '''
        Return the exponents of a given shell or if the shell isn't
        specified return all of the available exponents

        Args:
            asdict (bool): if `True` a dict with exps per shell is
                           returned
        '''

        if asdict:
            return {shell: f['e'] for shell, f in self.functions.items()}
        else:
            exps = [self.functions[k]['e'] for k in self.functions.keys()]
            array = np.fromiter(chain.from_iterable(exps), dtype=np.float)
            return np.sort(array)[::-1]

    def get_filename(self, ext=None):

        name = self.name.strip().replace(' ', '_') + '-' + self.element
        if ext:
            return name + '.' + ext
        return name

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
            return sum(nspherical(get_l(shell)) * len(fs['cf'])
                       for shell, fs in self.functions.items())
        else:
            return sum(ncartesian(get_l(shell)) * len(fs['cf'])
                       for shell, fs in self.functions.items())

    def normalize(self):
        '''
        Normalize contraction coefficients for each contracted functions based
        on the primitive overlaps so that the norm is equal to 1.
        '''

        for shell, fs in self.functions.items():
            cc = self.contraction_matrix(shell)
            po = primitive_overlap(get_l(shell), fs['e'], fs['e'])
            for col in range(cc.shape[1]):
                norm2 = np.dot(cc[:, col], np.dot(po, cc[:, col]))
                fs['cf'][col]['cc'] = cc[fs['cf'][col]['idx'], col] / np.sqrt(norm2)

    def normalization(self):
        '''
        For each function (contracted) calculate the norm and return a list of
        tuples containing the shell, function index and the norm respectively.
        '''

        out = list()
        for shell, fs in self.functions.items():
            cc = self.contraction_matrix(shell)
            po = primitive_overlap(get_l(shell), fs['e'], fs['e'])
            for col in range(cc.shape[1]):
                norm2 = np.dot(cc[:, col], np.dot(po, cc[:, col]))
                out.append(tuple([shell, col, norm2]))
        return out

    def nprimitive(self, spherical=True):
        '''
        Return the number of primitive functions assuming sphrical or
        cartesian gaussians.

        Args:
          spherical : bool
            A flag to select either spherical or cartesian gaussians

        Returns:
          out : int
            Number of primitive function in the basis set
        '''

        if spherical:
            # calculate the number of spherical components per shell
            ncomp = [nspherical(get_l(shell))
                     for shell in self.functions.keys()]
        else:
            # calculate the number of cartesian components per shell
            ncomp = [ncartesian(get_l(shell))
                     for shell in self.functions.keys()]

        return sum([prim * nc for prim, nc in zip(self.primitives_per_shell(),
                                                  ncomp)])

    def partial_wave_expand(self):
        '''
        From a given basis set with shells spdf... return a list of basis sets
        that are subsets of the entered basis set with increasing angular
        momentum functions included [s, sp, spd, spdf, ...]
        '''

        res = list()
        shells = self.functions.keys()
        for i in range(1, len(shells) + 1):
            bscopy = copy(self)
            bscopy.functions = {k: v for k, v in self.functions.items()
                                if k in shells[:i]}
            res.append(bscopy)
        return res

    def primitives_per_shell(self):
        '''
        Calculate how many primitive functions are in each shell.

        Returns:
          out : list of ints
        '''

        return [len(f['e']) for s, f in self.functions.items()]

    def primitives_per_contraction(self):
        '''
        Calculate how many primities are used in each contracted function.

        Returns:
          out : list of ints
        '''
        return [[len(cc) for cc in f['cf']] for s, f in self.functions.items()]

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
                fs['cf'] = [np.array([tuple([i, 1.0])], dtype=CFDTYPE)
                            for i, _ in enumerate(fs['e'])]
            return res
        else:
            for shell, fs in self.functions.items():
                fs['cf'] = [np.array([tuple([i, 1.0])], dtype=CFDTYPE)
                            for i, _ in enumerate(fs['e'])]

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
        S = primitive_overlap(get_l(shell), exps, exps)
        cc = self.contraction_matrix(shell)
        return np.dot(cc.T, np.dot(S, cc))

    def sort(self, reverse=False):
        '''
        Sort shells in the order of increasing angular momentum and for each
        shell sort the exponents.

        Args:
          reverse : bool
            If `False` sort the exponents in each shell in the descending order
            (default), else sort exponents in ascending order
        '''

        self.functions = OrderedDict(sorted(self.functions.items(),
                                            key=lambda x: get_l(x[0])))

        for shell, fs in self.functions.items():
            if reverse:
                idx = np.argsort(fs['e'])
            else:
                idx = np.argsort(fs['e'])[::-1]

            # actually sort the exponents and coefficients
            fs['e'] = fs['e'][idx]
            for cf in fs['cf']:
                cf['idx'] = np.asarray([np.nonzero(idx == y)[0][0]
                                        for y in cf['idx']])
                cf.sort(order='idx')

    def to_cfour(self, comment="", efmt="15.8f", cfmt="15.8f"):
        '''
        Return a string with the basis set in (new) CFOUR format.

        Args:
          comment : str
            comment string
          efmt : str
            string describing output format for the exponents,
            default: "20.10f"
          cfmt : str
            string describing output format for the contraction
            coefficients, default: "15.8f"

        Returns:
          res : str
            basis set string in Cfour format
        '''

        am, ne, cf = [], [], []
        for shell, shellfs in sorted(self.functions.items(),
                                     key=lambda x: get_l(x[0])):
            am.append(get_l(shell))
            ne.append(len(shellfs["e"]))
            cf.append(len(shellfs["cf"]))

        res = "\n{e}:{s}\n{c}\n\n".format(e=self.element, s=self.name,
                                          c=comment)
        res += "{0:3d}\n".format(len(self.functions.keys()))
        res += "".join(["{0:5d}".format(x) for x in am]) + "\n"
        res += "".join(["{0:5d}".format(x) for x in cf]) + "\n"
        res += "".join(["{0:5d}".format(x) for x in ne]) + "\n"
        res += "\n"

        for shell, fs in self.functions.items():
            for lst in splitlist(fs['e'], 5):
                res += "".join(["{0:>{efmt}}".format(e, efmt=efmt)
                                for e in lst]) + "\n"
            res += "\n"
            # create an array with all the contraction coefficients for a given shell
            cc = self.contraction_matrix(shell)
            for row in cc:
                for splitrow in splitlist(row, 5):
                    res += "{c}".format(c="".join(["{0:{cfmt}}".format(c, cfmt=cfmt) for c in splitrow])) + "\n"
            res += "\n"
        return res

    def to_dalton(self, fmt='prec'):
        '''
        Return a string with the basis set in DALTON format.

        Args:
          fmt : str
            string describing output format for the exponents and coefficents
                - `prec` "20.10f"
                - `default` "10.4f"
                - or python format string e.g. "15.8f"

        Returns:
          res : str
            basis set string in Dalton format
        '''

        formats = {'prec': '20.10f', 'default': '10.4f'}

        ffmt = formats.get(fmt, fmt)

        if fmt == 'prec':
            ffmt = '20.10f'
            nitems = 3
            fmtlabel = 'H'
        elif fmt == 'default':
            ffmt = '10.4f'
            nitems = 7
            fmtlabel = ' '
        else:
            ffmt = fmt
            cwidth = int(ffmt.split('.')[0])
            nitems = (80 - cwidth) // cwidth
            fmtlabel = "{0:d}F{1}.{2}".format(nitems, ffmt.split('.')[0], re.sub('[A-Za-z]+', '', ffmt.split('.')[1]))

        res = "! {s}\n".format(s=self.name)
        for shell, fs in self.functions.items():
            res += "! {s} functions\n".format(s=shell)
            res += "{f:1s}{p:>4d}{c:>4d}\n".format(f=fmtlabel, p=len(fs['e']),
                                                   c=len(fs['cf']))
            # create an array with all the contraction coefficients for a given shell
            cc = self.contraction_matrix(shell)
            for expt, row in zip(fs['e'], cc):
                it = splitlist(row, nitems)
                res += "{e:>{efmt}}{c}".format(e=expt, efmt=ffmt, c="".join(["{0:{cfmt}}".format(c, cfmt=ffmt) for c in next(it)])) + "\n"

                for row in it:
                    res += ' ' * int(ffmt.split('.')[0]) + \
                        "".join(["{0:{cfmt}}".format(c, cfmt=ffmt) for c in row]) + "\n"
        return res

    def to_gamessus(self, efmt="20.10f", cfmt="15.8f"):
        '''
        Return a string with the basis set in GAMESS(US) format.

        Args:
          efmt : str
            string describing output format for the exponents,
            default: "20.10f"
          cfmt : str
            string describing output format for the contraction
            coefficients, default: "15.8f"

        Returns:
          res : str
            basis set string in Gamess(US) format
        '''

        res = ''
        for shell, fs in self.functions.items():
            for contraction in fs["cf"]:
                res += "{s:<1s}{n:>3d}\n".format(s=shell.upper(),
                                                 n=len(contraction))
                for i, (idx, coeff) in enumerate(contraction, start=1):
                    res += "{i:3d}{e:>{efmt}}{c:>{cfmt}}".format(i=i,
                                                                 e=fs['e'][idx],
                                                                 efmt=efmt,
                                                                 c=coeff,
                                                                 cfmt=cfmt) + "\n"
        return res + "\n"

    def to_gaussian(self, efmt="20.10f", cfmt="15.8f"):
        '''
        Return a string with the basis set in Gaussian format.

        Args:
          efmt : str
            string describing output format for the exponents,
            default: "20.10f"
          cfmt : str
            string describing output format for the contraction
            coefficients, default: "15.8f"

        Returns:
          res : str
            basis set string in Gaussian format
        '''

        res = "****\n{e:<7s}0\n".format(e=self.element)
        for shell, fs in self.functions.items():
            for contraction in fs["cf"]:
                res += "{s:<1s}{n:>4d}{i:>7.2f}\n".format(s=shell.upper(),
                                                          n=len(contraction),
                                                          i=1.0)
                for idx, coeff in contraction:
                    res += "{e:>{efmt}}{c:>{cfmt}}".format(e=fs['e'][idx],
                                                           efmt=efmt, c=coeff,
                                                           cfmt=cfmt) + "\n"
        return res + "****\n"

    def to_json(self, fname=None, **kwargs):
        '''
        Serizalize the basis set object to JSON format
        '''

        if fname:
            with open(fname, 'w') as fjson:
                json.dump(self.__dict__, fjson, **kwargs)

        return json.dumps(self.__dict__, cls=NumpyEncoder, **kwargs)

    def to_latex(self, efmt="20.10f", cfmt="15.8f"):
        '''
        Return a string with the basis set as LaTeX table/

        Args:
            efmt : str
                Output format for the exponents, default: "20.10f"
            cfmt : str
                Output format for the contraction coefficients,
                default: "15.8f"

        Returns:
            res : str
                basis set string in LaTeX format
        '''

        # get the number of contracted functions per shell
        ncf = [sum([len(cfs) > 1 for cfs in fs['cf']])
               for sh, fs in self.functions.items()]
        idx = np.argmax(ncf)
        nccols = ncf[idx]

        out = '\\begin{{tabular}}{{{}}}\n'.format('r' * (nccols + 2))
        out += 'No. & \multicolumn{1}{c}{Exponent} & ' +\
               '\multicolumn{{{0:d}}}{{c}}{{Coefficients }} \\\ \n'.format(nccols)

        for shell, fs in self.functions.items():
            out += '\hline \n'
            out += '\multicolumn{{{0:d}}}{{c}}{{ {1:s} shell }} \\\ \hline \n'.format(nccols + 2, shell)

            cc = self.contraction_matrix(shell)
            count = 0
            # select columns with more than 1 non zero coefficients
            nonzerocolmask = np.array([np.count_nonzero(col) > 1
                                       for col in cc.T])
            nonzerorowmask = np.array([np.count_nonzero(row) > 0
                                       for row in cc[:, nonzerocolmask]])
            if np.any(nonzerocolmask):
                for expt, cfc in zip(fs['e'][nonzerorowmask],
                                     cc[np.ix_(nonzerorowmask,
                                               nonzerocolmask)]):
                    count += 1
                    coeffrow = " & ".join(["{0:{cfmt}}".format(c, cfmt=cfmt)
                                           for c in cfc])
                    out += "{i:5d} & {e:>{efmt}} & {c}".format(i=count, e=expt,
                                                               efmt=efmt,
                                                               c=coeffrow) + " \\\ \n"

            if np.any(np.logical_not(nonzerocolmask)):
                for colidx in np.where(np.logical_not(nonzerocolmask))[0]:
                    count += 1
                    nonzerorowmask = np.array([np.count_nonzero(row) > 0
                                               for row in cc[:, colidx]])
                    e = fs['e'][nonzerorowmask][0]
                    c = cc[nonzerorowmask, :][0][colidx]
                    out += "{i:5d} & {e:>{efmt}} & \\\ \n".format(i=count, e=e,
                                                                  efmt=efmt)
        out += '\end{tabular}'
        return out

    def to_molpro(self, withpars=False, efmt="20.10f", cfmt="15.8f"):
        '''
        Return a string with the basis set in MOLPRO format.

        Args:
            withpars : bool
                A flag to indicate whether to wrap the basis with
                `basis={ }` string
            efmt : str
                Output format for the exponents, default: "20.10f"
            cfmt : str
                Output format for the contraction coefficients,
                default: "15.8f"

        Returns:
            res : str
                basis set string
        '''

        # reorder indices of cf if necessary without modiyfing the instance
        funs = deepcopy(self.functions)
        for shell, fs in funs.items():
            if not has_consecutive_indices(fs):
                funs[shell] = reorder_shell_to_consecutive(fs)

        res = ""
        for shell, fs in funs.items():
            exps = ", ".join(["{0:>{efmt}}".format(e, efmt=efmt).lstrip()
                              for e in fs['e']]) + '\n'
            res += "{s:>s}, {e:>s}, ".format(s=shell, e=self.element) + exps

            for icol, cf in enumerate(fs['cf']):
                if cf.size == 1:
                    coeffs = ", ".join(["{0:>{cfmt}}".format(cc,
                        cfmt=cfmt).lstrip() for cc in cf['cc']])
                    res += "c, {0:d}.{0:d}, ".format(cf['idx'][0] + 1) + \
                        coeffs + '\n'
                else:
                    # indices should be consecutive
                    coeffs = ", ".join(["{0:>{cfmt}}".format(cc,
                        cfmt=cfmt).lstrip() for cc in cf['cc']])
                    res += "c, {0:d}.{1:d}, ".format(cf['idx'].min() + 1,
                                                     cf['idx'].max() + 1) +\
                        coeffs + '\n'
        if withpars:
            res = 'basis={\n' + res + '}'
        return res

    def to_nwchem(self, efmt="20.10f", cfmt="15.8f"):
        '''
        Return a string with the basis set in NWChem format.

        Args:
          efmt : str
            string describing output format for the exponents,
            default: "20.10f"
          cfmt : str
            string describing output format for the contraction
            coefficients, default: "15.8f"

        Returns:
          res : str
            basis set string in NwChem format
        '''

        res = 'BASIS "ao basis" PRINT\n'
        for shell, fs in self.functions.items():
            # create an array of contraction coefficients for a given shell
            cc = self.contraction_matrix(shell)

            # select columns with more than 1 non zero coefficients
            nonzerocolmask = np.array([np.count_nonzero(col) > 1
                                       for col in cc.T])
            nonzerorowmask = np.array([np.count_nonzero(row) > 0
                                       for row in cc[:, nonzerocolmask]])
            if np.any(nonzerocolmask):
                res += "{e} {s}\n".format(e=self.element, s=shell)
                for expt, cfc in zip(fs['e'][nonzerorowmask],
                                     cc[np.ix_(nonzerorowmask, nonzerocolmask)]):
                    res += "{e:>{efmt}}{c}".format(e=expt, efmt=efmt, c="".join(["{0:{cfmt}}".format(c, cfmt=cfmt) for c in cfc])) + "\n"

            if np.any(np.logical_not(nonzerocolmask)):
                for colidx in np.where(np.logical_not(nonzerocolmask))[0]:
                    res += "{e} {s}\n".format(e=self.element, s=shell)
                    nonzerorowmask = np.array([np.count_nonzero(row) > 0
                                               for row in cc[:, colidx]])
                    e = fs['e'][nonzerorowmask][0]
                    c = cc[nonzerorowmask, :][0][colidx]
                    res += "{e:>{efmt}}{c:>{cfmt}}".format(e=e, efmt=efmt, c=c,
                                                           cfmt=cfmt) + "\n"
        return res + "END\n"

    def to_pickle(self, fname=None):
        '''Save the basis set in pickle format under the filename `fname`

        Args:
          fname : str
            File name
        '''

        if fname is None:
            fname = self.get_filename('pkl')

        with open(fname, 'wb') as fbas:
            pickle.dump(self, fbas)

    def print_functions(self, efmt="20.10f", cfmt="15.8f"):
        '''
        Return a string with the basis set.

        Args:
          efmt : str
            string describing output format for the exponents,
            default: "20.10f"
          cfmt : str
            string describing output format for the contraction
            coefficients, default: "15.8f"

        Returns:
          res : str
            basis set string
        '''

        res = ''
        for shell, fs in self.functions.items():
            # create an array with all the contraction coefficients
            # for a given shell
            res += "\n" + "{s} shell".format(s=shell).center(40, '=') + "\n"
            cc = self.contraction_matrix(shell)
            count = 0
            # select columns with more than 1 non zero coefficients
            nonzerocolmask = np.array([np.count_nonzero(col) > 1
                                       for col in cc.T])
            nonzerorowmask = np.array([np.count_nonzero(row) > 0
                                       for row in cc[:, nonzerocolmask]])
            if np.any(nonzerocolmask):
                res += 'Contracted:\n'
                for expt, cfc in zip(fs['e'][nonzerorowmask],
                                     cc[np.ix_(nonzerorowmask,
                                               nonzerocolmask)]):
                    count += 1
                    coeffrow = "".join(["{0:{cfmt}}".format(c, cfmt=cfmt)
                                        for c in cfc])
                    res += "{i:5d}{e:>{efmt}}{c}".format(i=count, e=expt,
                                                         efmt=efmt,
                                                         c=coeffrow) + "\n"

            if np.any(np.logical_not(nonzerocolmask)):
                res += 'Uncontracted:\n'
                for colidx in np.where(np.logical_not(nonzerocolmask))[0]:
                    count += 1
                    nonzerorowmask = np.array([np.count_nonzero(row) > 0
                                               for row in cc[:, colidx]])
                    e = fs['e'][nonzerorowmask][0]
                    c = cc[nonzerorowmask, :][0][colidx]
                    res += "{i:5d}{e:>{efmt}}{c:>{cfmt}}".format(i=count, e=e,
                                                                 efmt=efmt,
                                                                 c=c,
                                                                 cfmt=cfmt) + "\n"
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
            res += "{k:<20s} = {v}\n".format(k=key.capitalize(),
                                             v=getattr(self, key))
        res += 'Functions:\n'
        res += self.print_functions()
        return res


def has_consecutive_indices(shell):
    '''
    Check if all the contracted functions have consecutive indices

    Args:
        shell : dict
            Basis functions for a given shell asa dict with structure
            ``{'e' : np.array(), 'cf': [np.array(), np.array(), ...]}``
    '''

    for i, cf in enumerate(shell['cf']):
        if len(cf) > 1:
            if not np.all(np.sort(cf['idx'])[1:] -
                          np.sort(cf['idx'])[:-1] == 1):
                return False
    else:
        return True


def reorder_shell_to_consecutive(shell):
    '''
    Reorder the exponents so that the indices of the contracted functions
    have consecutive inidices.

    Args:
        shell : dict
            Basis functions for a given shell asa dict with structure
            ``{'e' : np.array(), 'cf': [np.array(), np.array(), ...]}``

    Returns:
        shell : dict
            Same shell as on input but with reordered exponents and relabelled
            contracted functions
    '''

    # get the index of the largest contracted function
    cfidx = np.argmax([len(x) for x in shell['cf']])

    # get the new indices where the ones corresponding to the contracted
    # function are consecutive
    oldidx = shell['cf'][cfidx]['idx']
    conidx = np.arange(shell['e'].size, dtype=int)
    restidx = np.setdiff1d(conidx, oldidx)
    newidx = np.concatenate((oldidx, restidx))

    funs = {'e': shell['e'][newidx], 'cf': []}

    newidxlist = newidx.tolist()

    for cf in shell['cf']:
        funs['cf'].append(np.array([(newidxlist.index(i), c) for i, c in cf],
                                   dtype=CFDTYPE))

    # check if all the contracted functions have consecutive indices
    if not has_consecutive_indices(funs):
        raise ValueError('cannot reorder functions, check #cf: {}'.format(i))

    return funs


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

    newf = deepcopy(first.functions)
    for oshell, ofs in other.functions.items():
        if oshell.lower() in first.functions.keys():
            exps, idxs, idxo = merge_exponents(first.functions[oshell]['e'], ofs['e'])
            newf[oshell]['e'] = exps
            for cf in newf[oshell]['cf']:
                cf['idx'] = idxs[cf['idx']]
            newcf = [np.copy(cf) for cf in ofs['cf']]
            for cf in newcf:
                cf['idx'] = idxo[cf['idx']]
            newf[oshell]['cf'].extend(newcf)
        else:
            newf[oshell] = deepcopy(ofs)

    return BasisSet(name=first.name, element=first.element, functions=newf)


def primitive_overlap(l, a, b):
    '''
    Calculate the overlap integrals for a given shell `l` and two sets of
    exponents

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

    return np.power(2.0 * np.sqrt(np.outer(a, b)) / np.add.outer(a, b), l + 1.5)


def nspherical(l):
    '''
    Calculate the number of spherical components of a function with a given
    angular momentum value `l`.
    '''
    return 2 * l + 1


def ncartesian(l):
    '''
    Calculate the number of cartesian components of a function with a given
    angular momentum value *l*.
    '''
    return int((l + 1) * (l + 2) / 2)


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
    c = np.asarray([1.0] * kmax, dtype=float)

    leg = np.polynomial.Legendre(c)
    a = np.zeros((len(zetas), kmax))

    for j in range(len(zetas)):
        for k in range(kmax):
            arg = (2.0 * (j + 1.0) - 2.0) / (len(zetas) - 1.0) - 1.0
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

    return min(zetas), np.power(max(zetas) / min(zetas), 1.0 / (len(zetas) - 1))


def generate_exponents(formula, nf, params):
    '''
    Generate a sequence of exponents from a specified formula

    Args:
      formula : str
        name of the sequence from which the exponents are generated, one of:
          - `et`, `even`, `eventemp`, `even tempered`
          - `wt`, `well`, `welltemp`, `well tempered`
          - `le`, `legendre`
          - `exp`, `exponents`
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


def get_num_params(ftuple):
    '''
    Get the number of parameter for a given 4-tuple desciribing functions

    >>> get_num_params(('s', 'et', 20, (0.4, 2.0)))
    2
    '''

    formula = ftuple[1]
    if formula.lower() in ["et", "even", "eventemp", "even tempered"]:
        return 2
    elif formula.lower() in ["wt", "well", "welltemp", "well tempered"]:
        return 4
    elif formula.lower() in ["le", "legendre"]:
        return ftuple[2]
    elif formula.lower() in ["exp", "exponents"]:
        return ftuple[2]
    else:
        raise ValueError('unknown sequence name: {}'.format(formula))


def eventemp(numexp, params):
    '''
    Generate a sequence of nf even tempered exponents accodring to
    the even tempered formula

    .. math::
       \\zeta_i = \\alpha \cdot \\beta^{i-1}

    Args:
        numexp : int
            Number fo exponents to generate
        params : tuple of floats
            Alpha and beta parameters
    Returns:
        res : numpy array
            Array of generated exponents (floats)
    '''

    if not isinstance(numexp, int):
        raise TypeError('"nf" variable should be of "int" type, got: {}'.format(type(numexp)))
    if len(params) != 2:
        raise ValueError('"params" tuple should have exactly 2 entries, got {}'.format(len(params)))

    alpha, beta = params
    zetas = alpha * np.power(beta, np.arange(numexp))
    return zetas[::-1]


def welltemp(numexp, params):
    '''
    Generate a sequence of nf well tempered exponents accodring to
    the well tempered fromula

    .. math::

       \\zeta_i = \\alpha \cdot \\beta^{i-1} \cdot \\left[1 + \\gamma \cdot \\left(\\frac{i}{N}\\right)^{\delta}\\right]

    Args:
        numexp : int
            Number fo exponents to generate
        params : tuple of floats
            Alpha, beta, gamma and delta parameters

    Returns:
        res : numpy.array
            Array of generated exponents (floats)
    '''

    if not isinstance(numexp, int):
        raise TypeError('"nf" variable should be of "int" type, got: {}'.format(type(numexp)))
    if len(params) != 4:
        raise ValueError('"params" tuple should have exactly 4 entries, got {}'.format(len(params)))

    alpha, beta, gamma, delta = params
    zetas = alpha * np.power(beta, np.arange(numexp)) * \
            (1 + gamma * np.power(np.arange(1, numexp + 1) / numexp, delta))
    zetas.sort()
    return zetas[::-1]


def legendre(numexp, coeffs):
    '''
    Generate a sequence of nf exponents from expansion in the orthonormal
    legendre polynomials as described in: Peterson, G. A. et.al J. Chem. Phys.,
    Vol. 118, No. 3 (2003), eq. (7).

    .. math::
       \ln \\zeta_i = \\sum^{k_{\max}}_{k=0} A_k P_k \\left(\\frac{2j-2}{N-1}-1\\right)

    Args:
        numexp : int
            Number fo exponents to generate
        params : tuple of floats
            Polynomial coefficients (expansion parameters)

    Returns:
        res : numpy array
            Array of generated exponents (floats)
    '''

    if not isinstance(numexp, int):
        raise TypeError('"nf" variable should be of "int" type, got: {}'.format(type(numexp)))
    if len(coeffs) < 1:
        raise ValueError('"coeffs" tuple should have at least 1 entry, got {}'.format(len(coeffs)))

    # special case for one function
    if len(coeffs) == 1:
        return [np.exp(coeffs[0])]

    poly = np.polynomial.legendre.Legendre(coeffs)
    zetas = [poly(((2.0 * (i + 1.0) - 2.0) / (numexp - 1.0)) - 1.0) for i in range(numexp)]
    return np.exp(zetas[::-1])


def splitlist(l, n):
    '''
    Split a list into sublists of size `n`
    '''

    if len(l) % n == 0:
        splits = len(l) // n
    elif len(l) % n != 0 and len(l) > n:
        splits = len(l) // n + 1
    else:
        splits = 1

    for i in range(splits):
        yield l[n * i:n * i + n]


def sliceinto(toslice, chunk_sizes):
    '''
    Slice a list ``toslice`` into chunks of sizes defined in ``chunk_sizes``

    >>> sliceinto(list(range(10)), (6, 4))
    [[0, 1, 2, 3, 4, 5], [6, 7, 8, 9]]
    '''

    if len(toslice) != sum(chunk_sizes):
        raise ValueError('cannot slice list, size mismatch {} != {}'.format(len(toslice),
                                                                            sum(chunk_sizes)))

    sidx = [sum(chunk_sizes[:i]) for i in range(0, len(chunk_sizes))]
    return [toslice[i:i + s] for i, s in zip(sidx, chunk_sizes)]


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

    ncart = (l + 1) * (l + 2) // 2
    out = np.zeros((ncart, 3), dtype=np.int32)
    index = 0
    for i in range(l + 1):
        for j in range(i + 1):
            out[index, :] = [l - i, i - j, j]
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
        Expansion coefficients of real spherical harmonics in terms of
        cartesian gaussians
    '''

    ncart = (l + 1) * (l + 2) // 2
    nspher = 2 * l + 1

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
        tmpc = factorial2(2 * kx - 1) * factorial2(2 * ky - 1) * factorial2(2 * kz - 1) / factorial2(2 * l - 1)
        tmpc = np.sqrt(tmpc) / factorial(j)
        for i in range(l // 2 + 1):
            if j > i:
                continue
            if kx % 2 == 1:
                continue
            k = kx // 2
            if k > j:
                continue
            tmpi = tmpc * factorial2(2 * l - 2 * i - 1) / factorial(l - 2 * i) / factorial(i - j) / 2.0**i
            if i % 2 == 1:
                tmpi = -tmpi
            out[icart, l] += binom(j, k) * tmpi

    for m in range(1, l + 1):
        for icart in range(ncart):
            kx = cartc[icart, 0]
            ky = cartc[icart, 1]
            kz = cartc[icart, 2]
            jj = kx + ky - m
            if jj % 2 == 1:
                continue
            if jj < 0:
                continue
            j = jj // 2
            tmpc = factorial2(2 * kx - 1) * factorial2(2 * ky - 1) * factorial2(2 * kz - 1) / factorial2(2 * l - 1)
            tmpc = np.sqrt(2.0 * tmpc * factorial(l - m) / factorial(l + m)) / factorial(j)
            for i in range((l - m) // 2 + 1):
                if j > i:
                    continue
                tmpi = tmpc * factorial2(2 * l - 2 * i - 1) / factorial(l - m - 2 * i) / factorial(i - j) / 2.0**i
                if i % 2 == 1:
                    tmpi = -tmpi
                for k in range(j + 1):
                    kk = kx - 2 * k
                    if kk < 0 or kk > m:
                        continue
                    tmpk = tmpi * binom(j, k) * binom(m, kk)
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
