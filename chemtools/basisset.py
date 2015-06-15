#import pymongo

from __future__ import division

import copy
import decimal
import numpy as np
from itertools import chain
from collections import OrderedDict
import pickle
import re
import os
from elements import element

SHELLS = ["s", "p", "d", "f", "g", "h", "i", "j", "k"]

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
    def from_optdict(cls, x0, functs, name=None, element=None):

        total_nf = sum([x[1] for x in functs])
        if len(x0) != total_nf:
            raise ValueError('size mismatch {0} != {1}'.format(len(x0), total_nf))

        funs = dict()
        ni = 0; nt = 0
        for shell, nf, _ in functs:
            nt += nf
            funs[shell] = dict()
            funs[shell]['e'] = x0[ni:nt]
            ni += nf
        functions = OrderedDict(sorted(funs.items(), key=lambda x: SHELLS.index(x[0])))
        bs = cls(name=name, element=element, functions=functions)
        return bs

    @classmethod
    def from_sequence(cls, formula=None, functs=None, name=None, element=None):
        '''
        Return a basis set object generated from a sequence based on the specified
        arguments.

        Args:
          formula : str
            Generative formula for the exponents, allowed values are:
                - *eventemp*, *even*, *even tempered*
                - *welltemp*, *well*, *well tempered*
                - *legendre*
          functs : list of tuples
            A list of tuple specifying the shell type, number of functions and
            parameters, e.g. `[('s', 4, (0.5, 2.0)), ('p', 3, (1.0, 3.0))]`
          name : str
            Name of the basis set
          element : str
            Chemical symbol of the element
        '''

        funs = dict()

        if isinstance(formula, str):
            formulas = [formula]*len(functs)
        elif isinstance(formula, (list, tuple)):
            formulas = formula
        else:
            raise ValueError('formula should be <str>, <list> or <tuple>, got: {}'.format(type(formula)))

        for seqf, (shell, nf, params) in zip(formulas, functs):
            funs[shell] = dict()
            funs[shell]['e'] = generate_exponents(seqf, nf, params)
        functions = OrderedDict(sorted(funs.items(), key=lambda x: SHELLS.index(x[0])))
        bs = cls(name=name, element=element, functions=functions)
        return bs

    @classmethod
    def from_file(cls, fname=None, fmt=None, name=None):
        '''Read and parse a basis set from file and return a BasisSet object

        Args:
          fname : str
            File name
          fmt : str
            Format of the basis set in the file: *molpro*, *gamess*
        '''

        if name is None:
            name = os.path.splitext(fname)[0]

        formats = ['molpro', 'gamessus']

        if fmt.lower() not in formats:
            raise ValueError("<fmt> should be one of: {}".format(", ".join(formats)))

        with open(fname, 'r') as fobj:
            basstr = fobj.read()

        res = dict()
        if fmt == 'molpro':
            bsd = parse_molpro_basis(basstr)
        elif fmt == 'gamessus':
            bsd = parse_gamess_basis(basstr)
        # return a BasisSet object if only one basis parsed or dict of BasisSet
        # objects with atomic symbol as key

        if len(bsd.keys()) == 1:
            return BasisSet(name=name, element=atom,
                    functions=OrderedDict(sorted(fs.items(), key=lambda x: SHELLS.index(x[0]))))
        else:
            for atom, fs in bsd.items():
                res[atom] = BasisSet(name=name, element=atom,
                        functions=OrderedDict(sorted(fs.items(), key=lambda x: SHELLS.index(x[0]))))
            return res

    def get_exponents(self, shell=None):
        '''
        Return the exponents of a given shell or if the shell isn't specified
        return all of the available exponents
        '''

        if shell is None:
            return chain.from_iterable([self.functions[k]['e'] for k in self.functions.keys()])
        else:
            return self.functions[shell]['e']

    def uncontract(self):
        '''
        Uncontract the basis set. This replaces the contraction coefficients in
        the current object.
        '''

        for shell, fs in self.functions.items():
            fs['cc'] = [[tuple([i, 1.0])] for i, _ in enumerate(fs['e'])]

    def __repr__(self):
        res = "<BasisSet(\n"
        for key, val in self.__dict__.items():
            res += "\t{k:<20s} = {v}\n".format(k=key, v=val)
        res += ")>"
        return res

    def __str__(self):
        res = "<BasisSet(\n"
        for key, val in self.__dict__.items():
            res += "\t{k:<20s} = {v}\n".format(k=key, v=val)
        res += ")>"
        return res

    def __add__(self, bs2add):

        if not isinstance(bs2add, BasisSet):
            raise TypeError('only instances of "BasisSet" class can be added')

        for oshell, ofs in bs2add.functions.items():
            if oshell.lower() in self.functions.keys():
                if self.common_floats(self.functions[oshell]['exponents'], ofs['exponents']):
                    exps, indx = self.merge_list_of_floats(self.functions[oshell.lower()]['exponents'], ofs['exponents'])
                    self.functions[oshell.lower()]['exponents'] = exps
                    for cf in ofs['contractedfs']:
                        newc = copy.copy(cf)
                        newc['indices'] = [indx[i] for i in newc['indices']]
                        self.functions[oshell.lower()]['contractedfs'].append(newc)
                else:
                    nexp = len(self.functions[oshell.lower()]['exponents'])
                    self.functions[oshell.lower()]['exponents'].extend(ofs['exponents'])
                    for cf in ofs['contractedfs']:
                        newc = copy.copy(cf)
                        newc['indices'] = [i + nexp for i in newc['indices']]
                        self.functions[oshell.lower()]['contractedfs'].append(newc)
            else:
                self.functions[oshell] = ofs

    def add(self, bs2add):
        '''
        Merge functions to the BasisSet object "functions" attribute from
        another BasisSet object if the "element" atrribute has the same value.
        '''

        if bs2add is not None and self.element == bs2add.element:
            for oshell, ofs in bs2add.functions.items():
                if oshell.lower() in self.functions.keys():
                    if self.common_floats(self.functions[oshell]['exponents'], ofs['exponents']):
                        exps, indx = self.merge_list_of_floats(self.functions[oshell.lower()]['exponents'], ofs['exponents'])
                        self.functions[oshell.lower()]['exponents'] = exps
                        for cf in ofs['contractedfs']:
                            newc = copy.copy(cf)
                            newc['indices'] = [indx[i] for i in newc['indices']]
                            self.functions[oshell.lower()]['contractedfs'].append(newc)
                    else:
                        nexp = len(self.functions[oshell.lower()]['exponents'])
                        self.functions[oshell.lower()]['exponents'].extend(ofs['exponents'])
                        for cf in ofs['contractedfs']:
                            newc = copy.copy(cf)
                            newc['indices'] = [i + nexp for i in newc['indices']]
                            self.functions[oshell.lower()]['contractedfs'].append(newc)
                else:
                    self.functions[oshell] = ofs
        else:
            return

    def print_exponents(self):

        for shell, fs in sorted(self.functions.items(), key=lambda x: SHELLS.index(x[0])):
            for exp in fs['e']:
                print("{s:10s}  {e:>25.10f}".format(s=shell, e=exp))


    @staticmethod
    def merge_list_of_floats(lof, lof2add, prec=6):
        '''
        merge a list of floats "lof2add" into a reference list
        of floats "lof" ommitting duplicate floats and return
        merged list and indices of newly added floats in the reference list.

        Args:
            lof (list of floats)
                reference list of floats to which new floats will be appended
            lof2add (list of floats)
                list of floats to be added to the reference list
            prec (int)
                precision (number of significant digits used in float comparison
                default=6
        Returns:
            lof (list of floats)
                appended list of floats
            newidx (list of integers)
                list of integers corresponding to the indices of newly added floats
        '''

        c = decimal.Context(prec=prec)
        decimal.setcontext(c)
        D = decimal.Decimal

        expref  = [D(x)*1 for x in lof]
        exp2add = [D(x)*1 for x in lof2add]
        nexp = len(lof)
        newidx = list()
        for addexp in exp2add:
            if addexp in expref:
                newidx.append(expref.index(addexp))
            else:
                expref.append(addexp)
                newidx.append(expref.index(addexp))
        for addidx, addexp in zip(newidx, lof2add):
            if addidx > nexp:
                lof.append(addexp)
        return lof, newidx

    def to_cfour(self, comment="", efmt="15.8f", cfmt="15.8f"):
        '''
        Return a string with the basis set in (new) CFOUR format.

        Args:
            comment (str)
                comment string
            efmt (str)
                string describing output format for the exponents, default:
                "20.10f"
            cfmt (str)
                string describing output format for the contraction
                coefficients, default: "15.8f"

        Returns:
            res (str)
                basis set string
        '''

        am, ne, cf = [], [], []
        for shell, shellfs in sorted(self.functions.items(), key=lambda x: SHELLS.index(x[0])):
            am.append(SHELLS.index(shell))
            ne.append(len(shellfs["e"]))
            cf.append(len(shellfs["cc"]))

        res = "\n{e}:{s}\n{c}\n\n".format(e=self.element, s=self.name, c=comment)
        res += "{0:3d}\n".format(len(self.functions.keys()))
        res += "".join(["{0:5d}".format(x) for x in am]) + "\n"
        res += "".join(["{0:5d}".format(x) for x in ne]) + "\n"
        res += "".join(["{0:5d}".format(x) for x in cf]) + "\n"
        res += "\n"

        for shell, fs in self.functions.items():
            for lst in splitlist(fs['e'], 5):
                res += "".join(["{0:>{efmt}}".format(e, efmt=efmt) for e in lst])  + "\n"
            res += "\n"
            for i, expt in enumerate(fs['e']):
                cc = [t[1] if t[0] == i else 0.0 for csf in fs['cc'] for t in csf]
                for lst in splitlist(cc, 5):
                    res += "{c}".format(c="".join(["{0:{cfmt}}".format(c, cfmt=cfmt) for c in lst])) + "\n"
            res += "\n"
        return res

    def to_dalton(self, efmt="20.10f", cfmt="15.8f"):
        '''
        Return a string with the basis set in DALTON format.

        Args:
            efmt (str)
                string describing output format for the exponents, default:
                "20.10f"
            cfmt (str)
                string describing output format for the contraction
                coefficients, default: "15.8f"

        Returns:
            res (str)
                basis set string
        '''

        res = "! {s}\n".format(s=self.name)
        for shell, fs in self.functions.items():
            res += "! {s} functions\n".format(s=shell)
            res += "{f:1s}{p:>4d}{c:>4d}\n".format(f="F", p=len(fs['e']), c=len(fs['cc']))
            for i, expt in enumerate(fs['e']):
                cc = [t[1] if t[0] == i else 0.0 for csf in fs['cc'] for t in csf]
                res += "{e:>{efmt}}{c}".format(e=expt, efmt=efmt, c="".join(["{0:{cfmt}}".format(c, cfmt=cfmt) for c in cc])) + "\n"
        return res

    def to_gamess(self, efmt="20.10f", cfmt="15.8f"):
        '''
        Return a string with the basis set in GAMESS(US) format.

        Args:
            efmt (str)
                string describing output format for the exponents, default:
                "20.10f"
            cfmt (str)
                string describing output format for the contraction
                coefficients, default: "15.8f"

        Returns:
            res (str)
                basis set string
        '''

        res = ""
        for shell, fs in self.functions.items():
            for contraction in fs["cc"]:
                res += "{s:<1s}{n:>3d}\n".format(s=shell.upper(), n=len(contraction))
                for i, (idx, coeff) in enumerate(contraction, start=1):
                    res += "{i:3d}{e:>{efmt}}{c:>{cfmt}}".format(i=i, e=fs['e'][idx], efmt=efmt, c=coeff, cfmt=cfmt)+ "\n"
        return res + "\n"

    def to_molpro(self, efmt="20.10f", cfmt="15.8f"):
        '''
        Return a string with the basis set in MOLPRO format.

        Args:
            efmt (str)
                string describing output format for the exponents, default:
                "20.10f"
            cfmt (str)
                string describing output format for the contraction
                coefficients, default: "15.8f"

        Returns:
            res (str)
                basis set string
        '''

        res = ""
        for shell, fs in self.functions.items():
            exps = ", ".join(["{0:>{efmt}}".format(e, efmt=efmt).lstrip() for e in fs['e']])
            res += "{s:>s}, {e:>s}, ".format(s=shell, e=self.element) + exps + '\n'
            for contraction in fs['cc']:
                indices = [c[0] for c in contraction]
                coeffs = ", ".join(["{0:>{cfmt}}".format(c[1], cfmt=cfmt).lstrip() for c in contraction])
                res += "c, {0:d}.{1:d}, ".format(min(indices)+1, max(indices)+1) + coeffs + '\n'
        return res


    def to_nwchem(self, efmt="20.10f", cfmt="15.8f"):
        '''
        Return a string with the basis set in NWChem format.

        Args:
            efmt (str)
                string describing output format for the exponents, default:
                "20.10f"
            cfmt (str)
                string describing output format for the contraction
                coefficients, default: "15.8f"

        Returns:
            res (str)
                basis set string
        '''

        res = 'BASIS "ao basis" PRINT\n'
        for shell, fs in self.functions.items():
            res += "{e} {s}\n".format(e=self.element, s=shell)
            for i, expt in enumerate(fs['e']):
                cc = [t[1] if t[0] == i else 0.0 for csf in fs['cc'] for t in csf]
                res += "{e:>{efmt}}{c}".format(e=expt, efmt=efmt, c="".join(["{0:{cfmt}}".format(c, cfmt=cfmt) for c in cc])) + '\n'
        return res + "END\n"

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

    def contraction_scheme(self):
        '''
        Return a string describing the contraction scheme.
        '''

        cs, ec = [], []
        for shell, fs in self.functions.items():
            cs.append((shell, len(fs['e']), len(fs['cc'])))
            ec.append([len(cfs) for cfs in fs['cc']])
        return "({p:s}) -> [{c:s}] : {{{ec}}}".format(
                p="".join(["{0:d}{1:s}".format(c[1], c[0]) for c in cs]),
                c="".join(["{0:d}{1:s}".format(c[2], c[0]) for c in cs]),
                ec="/".join([" ".join(["{0:d}".format(c) for c in x]) for x in ec]))

    def nf(self, spherical=True):
        '''
        Calculate the number of basis functions

        Args:
            spherical : bool
                flag indicating if spherical or cartesian functions should be
                used, default: True

        Returns:
            res : int
                number of basis functions
        '''

        if spherical:
            return sum(nspherical(SHELLS.index(shell))*len(fs['cc']) for shell, fs in self.functions.items())
        else:
            return sum(ncartesian(SHELLS.index(shell))*len(fs['cc']) for shell, fs in self.functions.items())

    @classmethod
    def uncontracted(cls, bs):
        '''
        Return a new BasisSet object with uncotracted version of the basis.
        '''
        bsnew = copy.deepcopy(bs)
        for shell, shellfs in sorted(bsnew.functions.items(), key=lambda x: SHELLS.index(x[0])):
            shellfs["contractedfs"] = [{"indices" : [i], "coefficients" : [1.0]} for i, e in enumerate(shellfs["exponents"])]
        return bsnew

    def primitives_per_shell(self):
        return [len(f['e']) for s, f in self.functions.items()]

    def contractions_per_shell(self):
        return [len(f['cc']) for s, f in self.functions.items()]

    def primitives_per_contraction(self):
        return [[len(cc) for cc in f['cc']] for s, f in self.functions.items()]

    def contraction_type(self):
        '''
        Determine the contraction type: segmented, general, uncontracted, unknown
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
    From a set of exponents "zetas", using least square fit calculate the
    expansion coefficients into the legendre polynomials of the order "kmax".

    Args:
        kmax (int)
            length of the legendre expansion
        zetas (list)
            list of exponents (floats) to be fitted

    Returns:
        coeff (np.array)
            numpy array of legendre expansion coeffcients of length kmax
    '''

    # exponents should be sorted in the acsending order
    zetas = sorted(zetas)
    c = np.asarray([1.0]*kmax, dtype=float)

    leg = np.polynomial.Legendre(c)
    a = np.zeros((len(zetas), kmax))

    for j in range(len(zetas)):
        for k in range(kmax):
            arg = (2.0*(j+1.0)-2.0)/(len(zetas)-1.0)-1.0
            a[j, k] = leg.basis(k)(arg)

    return np.linalg.lstsq(a, np.log(zetas))[0]

def generate_exponents(formula, nf, params):

    if formula.lower() in ["even", "eventemp", "even tempered"]:
        return eventemp(nf, params)
    elif formula.lower() in ["well", "welltemp", "well tempered"]:
        return welltemp(nf, params)
    elif formula.lower() in ["legendre"]:
        return legendre(nf, params)
    else:
        raise ValueError('unknown sequence name: {}'.format(formula))

def eventemp(nf, params):
    '''
    Generate a sequence of nf even tempered exponents accodring to
    the even tempered formula

    :math:`\\zeta_i = \\alpha \cdot \\beta^{i-1}`.

    Args:
        nf : int
            number fo exponents to generate
        params : tuple of floats
            alpha and beta parameters
    Returns:
        res : list
            list of generated exponents (floats)
    '''
    if not isinstance(nf, int):
        raise TypeError('"nf" variable should be of "int" type, got: {}'.format(type(nf)))
    if len(params) !=  2:
        raise ValueError('"params" tuple should have exactly 2 entries, got {}'.format(len(params)))

    alpha, beta = params
    return [alpha * np.power(beta, i) for i in range(nf)]

def welltemp(nf, params):
    '''
    Generate a sequence of nf well tempered exponents accodring to
    the well tempered fromula

    :math:`\\zeta_i = \\alpha \cdot \\beta^{i-1} \cdot \\left[1 + \\gamma \cdot \\left(\\frac{i}{N}\\right)^{\delta}\\right]`.

    Args:
        nf : int
            number fo exponents to generate
        params : tuple of floats
            alpha, beta, gamma and delta parameters
    Returns:
        res : list of floats
            list of generated exponents (floats)
    '''
    if not isinstance(nf, int):
        raise TypeError('"nf" variable should be of "int" type, got: {}'.format(type(nf)))
    if len(params) !=  4:
        raise ValueError('"params" tuple should have exactly 4 entries, got {}'.format(len(params)))

    alpha, beta, gamma, delta = params
    return [alpha*np.power(beta, i)*(1+gamma*np.power((i+1)/nf, delta)) for i in range(nf)]

def legendre(nf, coeffs):
    '''
    Generate a sequence of nf exponents from expansion in the orthonormal
    legendre polynomials as described in: Peterson, G. A. et.al J. Chem. Phys.,
    Vol. 118, No. 3 (2003), eq. (7).

    :math:`\ln \\zeta_i = \\sum^{k_{\max}}_{k=0} A_k P_k \\left(\\frac{2j-2}{N-1}-1\\right)`

    Args:
        nf : int
            number fo exponents to generate
        params : tuple of floats
            polynomial coefficients (expansion parameters)
    Returns:
        res : list of floats
            list of generated exponents (floats)
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
    return [np.exp(x) for x in zetas]

def group(lst, n):
    """group([0,3,4,10,2,3], 2) => [(0,3), (4,10), (2,3)]

    Group a list into consecutive n-tuples. Incomplete tuples are
    discarded e.g.

    >>> group(range(10), 3)
    [(0, 1, 2), (3, 4, 5), (6, 7, 8)]
    """
    return list(zip(*[lst[i::n] for i in range(n)]))

def get_x0(basisopt):
    '''
    Collect all the parameters in a consecutive list of elements.
    '''

    return list(chain.from_iterable(basisopt["params"]))

def splitlist(l, n):
    if len(l) % n == 0:
        splits = len(l)/n
    elif len(l) % n != 0 and len(l) > n:
        splits = len(l)/n+1
    else:
        splits = 1

    for i in range(splits):
        yield l[n*i:n*i+n]

def common_floats(lof, reflof):

    c = decimal.Context(prec=6)
    decimal.setcontext(c)
    D = decimal.Decimal

    l    = [D(x)*1 for x in lof]
    refl = [D(x)*1 for x in reflof]
    return any(D('0') == x.compare(y) for x in l for y in refl)

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
        at_symbol, shell = parse_shell(lines[i[0]], lines[i[0]+1:i[1]])
        if at_symbol in bs.keys():
            bs[at_symbol] = dict(list(bs[at_symbol].items()) + list(shell.items()))
        else:
            bs[at_symbol] = shell
    return bs

def parse_shell(expsline, coeffs):
    '''
    Parse functions of one shell in molpro format.
    '''

    # remove empty strings and whitespace line breaks and tabs
    coeffs = [x.strip() for x in coeffs if x.strip() not in ['', '\n', '\t']]

    fs = {}

    shell  = expsline.split(",")[0]
    at_symbol = expsline.split(",")[1].strip().capitalize()
    exps   = [float(x) for x in expsline.rstrip(";").split(",")[2:]]

    fs[shell.lower()] = {'e' : exps, 'cc' : []}
    if len(coeffs) != 0:
        for line in coeffs:
            lsp = line.rstrip(";").split(",")
            if lsp[0] == "c":
                i, j = [int(x) for x in lsp[1].split(".")]
                coeffs = [float(x) for x in lsp[2:]]
                fs[shell.lower()]['cc'].append(list(zip(list(range(i-1, j)), coeffs)))
    else:
        for i in range(len(exps)):
            fs[shell.lower()]['cc'].append([tuple([i, 1.0])])
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

def parse_gamess_basis(string):
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
                print(elem.symbol)
                bslines = basis.split("\n")
                for i, line in enumerate(bslines):
                    match = pat.search(line)
                    if match:
                        shell, nf = match.group("shell").lower(), match.group("nf")
                        exps, indxs, coeffs = parse_function(bslines[i+1:i+int(nf)+1])
                        if shell in functions.keys():
                            lasti = len(functions[shell]['e'])
                            functions[shell]['e'].extend(exps)
                            functions[shell]['cc'].append(list(zip([lasti + i - 1 for i in indxs], coeffs)))
                        else:
                            functions[shell] = dict()
                            functions[shell]['e'] = exps
                            functions[shell]['cc'] = list()
                            functions[shell]['cc'].append(list(zip([i - 1 for i in indxs], coeffs)))
                res[elem.symbol] = functions
    return res

def parse_function(los):
    '''
    Parse a basis set function information from list of strings into
    three lists containg: exponents, indices, coefficients.

    Remeber that python doesn't recognise the `1.0d-3` format where `d` or
    `D` is used to the regex subsitution has to take care of that.
    '''

    real = re.compile(r'[dD]')

    indxs = [int(item.split()[0]) for item in los]
    exps = [float(real.sub('E', item.split()[1])) for item in los]
    coeffs = [float(real.sub('E', item.split()[2])) for item in los]

    return (exps, indxs, coeffs)


