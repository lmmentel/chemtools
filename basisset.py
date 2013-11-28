import os
import re
import sys
import periodic
import numpy as np
from operator import itemgetter
from itertools import groupby, chain


_shells = ["S", "P", "D", "F", "G", "H", "I"]

class Basis(object):
    '''A class for handling basis sets and parising basis sets from files.'''


    _angular = {"S" : 0, "P" : 1, "D" : 2, "F" : 3, "G" :  4, "H" : 5, "I" : 6}

    def __init__(self, path, name):
        '''To initialize just pass the path under which there is the basis set repository.'''
        self.path = path
        self.name = name
        self.exists()

    def __repr__(self):
        return "Basis set repository is located in: {}".format(self.path)

    def exists(self):
        '''Check if the basis set repository exists.'''
        if os.path.exists(self.path):
            return True
        else:
            sys.exit("Basis set repository doesn't exist under '{}'".format(self.path))

    def get_basis(self, atomic):
        '''Retrieve the basis set from the basis set repository file. The location of
           the repository is give when the basis object is initialized.

            The repository file can be simply downloaded from EMSL and saved
            as text file for further processing.
        '''

        with open(os.path.join(self.path, Basis._basisfiles[self.name]), 'r') as f:
            lines  = f.readlines()
        start = 0
        for i, line in enumerate(lines):
            if Basis._elements[abs(atomic)] in line:
                start = i + 1
                break

        if start == 0:
            sys.exit("No {0:s} basis for {1:s} found in the reposository, exiting...".format(self.name, Basis._elements[abs(atomic)]))

        for i, line in enumerate(lines[start::]):
            if line == '\n':
                stop = i
                break
        basisout = lines[start:start + stop]
        return "".join(basisout)

def uncontract(bs):
    '''
    Uncontract the contracted basis set
    '''

    # check for duplicate exps in the shell
    newfs = []
    for key, group in groupby(bs, lambda x: x["shell"]):
        exps = []
        for f in group:
            for zeta in f["exps"]:
                if zeta not in exps:
                    exps.append(zeta)
                    newdict = f.copy()
                    newdict["exps"] = [zeta]
                    newdict["coeffs"] = [1.0]
                    newfs.append(newdict)
    return newfs

def get_spherical(l):
    ''' Calculate number of spherical gaussian functions for a
        given angular momentum value "l".
    '''
    return 2*l+1

def get_cartesian(l):
    ''' Calculate number of cartesian gaussian functions for a
        given angular momentum value "l".
    '''
    return (l+1)*(l+2)/2

def get_no_functions(bs, spherical=True):
    '''
    Calculate total number of functions
    '''

    if spherical:
        nof = sum([get_spherical(f["shell"]) for f in bs])
    else:
        nof = sum([get_cartesian(f["shell"]) for f in bs])
    return nof

def print_functions(fs):

    print "{0:<10s}{1:<10s}{2:<10s}{3:^20s}{4:^20s}".format("Atom",
                                                            "No.",
                                                            "Shell",
                                                            "Exponent",
                                                            "Coefficient")
    for i, row in enumerate(sorted(fs, key=itemgetter("atomic", 'shell'))):
        for f in xrange(len(row["exps"])):
            print "{0:<10s}{1:<10d}{2:<10s}{3:>20.10f}{4:>20.10f}".format(
                    periodic.element(row["atomic"]).symbol,
                    i+1,
                    _shells[row["shell"]],
                    float(row["exps"][f]),
                    float(row["coeffs"][f]))


def group(lst, n):
    """group([0,3,4,10,2,3], 2) => [(0,3), (4,10), (2,3)]

    Group a list into consecutive n-tuples. Incomplete tuples are
    discarded e.g.

    >>> group(range(10), 3)
    [(0, 1, 2), (3, 4, 5), (6, 7, 8)]
    """
    return zip(*[lst[i::n] for i in range(n)])

def get_x0(basisopt):
    '''Collect all the parameters in a consecutive list of elements.'''

    x0 = list(chain.from_iterable(basisopt["params"]))
    return x0

def get_basis(x0, bsoptdict):

    icount = 0
    basis = []
    for lqn, nf in enumerate(bsoptdict["nfpshell"]):
        if bsoptdict["typ"].lower() in ["direxp", "direct", "exps", "exponents"]:
            exps = x0[icount:icount+nf]
            icount += nf
        elif bsoptdict["typ"].lower() in ["even", "eventemp", "eventempered"]:
            groups = group(x0, 2)
            exps = eventemp(nf, groups[lqn])
        elif bsoptdict["typ"].lower() in ["well", "welltemp", "welltempered"]:
            groups = group(x0, 4)
            exps = welltemp(nf, groups[lqn])
        elif bsoptdict["typ"].lower() in ["legendre"]:
            npar = len(bsoptdict["params"][lqn])
            pars = x0[icount:icount+npar]
            exps = legendre(nf, pars)

        for function in xrange(nf):
            basis.append({"shell"  : lqn,
                          "typ"    : bsoptdict["name"],
                          "exps"   : [exps[function]],
                          "coeffs" : [1.0],
                          "atomic" : bsoptdict["atomic"],
                         })
    return basis

def get_params_and_basis(bsopt):

    x0    = get_x0(bsopt)
    basis = []
    for lqn, nf in enumerate(bsopt["nfpshell"]):
        if bsopt["typ"].lower() in ["direxp", "direct", "exps", "exponents"]:
            exps = bsopt["params"][lqn]
        elif bsopt["typ"].lower() in ["even", "eventemp", "eventempered"]:
            exps = eventemp(nf, bsopt["params"][lqn])
        elif bsopt["typ"].lower() in ["well", "welltemp", "welltempered"]:
            exps = welltemp(nf, bsopt["params"][lqn])
        elif bsopt["typ"].lower() in ["legendre"]:
            exps = legendre(nf, bsopt["params"][lqn])
        else:
            sys.exit("wrong specification of basis optimzation parameters type : {0}".format(bsopt["typ"]))
        for function in xrange(nf):
            basis.append({"shell"  : lqn,
                            "typ"    : bsopt["name"],
                            "exps"   : [exps[function]],
                            "coeffs" : [1.0],
                            "atomic" : bsopt["atomic"],
                            })
    return x0, basis

def eventemp(nf, params):
    '''
    Generate a sequence of nf even tempered exponents accodring to

    zeta_i = alpha * beta**(i-1)
    '''
    alpha, beta = params
    return [alpha * np.power(beta, i) for i in xrange(nf)]

def welltemp(nf, params):
    '''
    Generate a sequence of nf well tempered exponents accodring to

    zeta_i = alpha * beta**(i-1) * [1 + gamma * (i/nf)**delta]
    '''

    alpha, beta, gamma, delta = params
    return [alpha*np.power(beta, i)*(1+gamma*np.power((i+1)/nf, delta)) for i in xrange(nf)]

def legendre(nf, coeffs):
    '''
    Generate a sequence of nf exponents from expansion in the orthonormal
    legendre polynomials as described in:
    Peterson, G. A. et.al J. Chem. Phys., Vol. 118, No. 3 (2003).
    '''

    poly = np.polynomial.legendre.Legendre(coeffs)
    zetas = [poly(((2.0*(i+1.0)-2.0)/(nf-1.0))-1.0) for i in xrange(nf)]
    return [np.exp(x) for x in zetas]

# MolPro parser

def parse_shell(expsline, coeffs):
    '''
    Parse functions of one shell in molpro format.
    '''

    fs = []

    shell  = expsline.split(",")[0]
    element = periodic.element(expsline.split(",")[1].strip())
    symbol = expsline.split(",")[1]
    exps   = [float(x) for x in expsline.rstrip(";").split(",")[2:]]

    for line in coeffs:
        lsp = line.rstrip(";").split(",")
        if lsp[0] == "c":
            i, j = [int(x) for x in lsp[1].split(".")]
            coeffs = [float(x) for x in lsp[2:]]
            fs.append({"shell"  : _shells.index(shell.upper()),
                       "typ"    : "",
                       "exps"   : exps[i-1:j],
                       "coeffs" : coeffs,
                       "atomic" : element.atomic,
                      })
    return fs


def parse_molpro_basis(basisstring):
    '''
    Parse basis set from a string in MolPro format.
    '''
    # belowe is the regular expression that can be used for parsing
    #basis set from  molpro input files
    #>>> basre = re.compile(r'basis\s*=\s*\{(.*?)\}', flags=re.S)

    lines = basisstring.split("\n")

    start = []
    for i, line in enumerate(lines):
        if line.split(",")[0] in ["s","p", "d", "f", "g", "h", "i"]:
            start.append(i)

    startstop = []
    for i in range(len(start)-1):
        startstop.append((start[i], start[i+1]))
    startstop.append((start[-1], len(lines)))

    bs = []
    for i in startstop:
        bs.extend(parse_shell(lines[i[0]], lines[i[0]+1:i[1]]))
    return bs

def write_molpro_basis(basisset):
    '''
    Write basis set in molpro format

    This little function is quite dirty and would benefit from  rewriting!
    '''

    newd = []

    for atom, atomgroup in groupby(basisset, lambda x: x["atomic"]):
        for shell, shellgroup in groupby(atomgroup, lambda x: x["shell"]):
            exps    = []
            indices = []
            coeffs  = []
            for function in shellgroup:
                for zeta in function["exps"]:
                    if zeta not in exps:
                        exps.append(zeta)
                istart = exps.index(function["exps"][0]) + 1
                istop  = exps.index(function["exps"][-1]) + 1
                indices.append((istart, istop))
                coeffs.append(function["coeffs"])
            newd.append({"atomic" : atom,
                            "shell"  : shell,
                            "exps"   : exps,
                            "indices": indices,
                            "coeffs" : coeffs})
    outstring = "basis={\n"
    for f in sorted(newd, key=itemgetter("atomic", "shell")):
        elem = periodic.element(f["atomic"])
        outstring = outstring + "{0}, {1}, {2}\n".format(
                _shells[f["shell"]].lower(),
                elem.symbol,
                ", ".join([str(x) for x in f["exps"]]))
        for i, item in enumerate(f["coeffs"]):
            outstring = outstring + "{0}, {1}, {2}\n".format(
                "c",
                ".".join([str(x) for x in f["indices"][i]]),
                ", ".join([str(x) for x in item]))
    outstring = outstring + "}\n"
    return outstring

# Gamess-US parser

def parse_gamess_basis(basisset):
    '''
    Parse the basis set into a list of dictionaries from a string in
    gamess format.
    '''

    pat = re.compile(r'^\s*(?P<shell>[SPDFGHIspdfghi])\s*(?P<nf>[1-9])')

    bs = []
    for atomic, rawbs in basisset.items():
        bslines = rawbs.split("\n")

        for i, line in enumerate(bslines):
            match = pat.search(line)
            if match:
                bs.extend(parse_function(match.group("shell"),
                                         atomic,
                                         bslines[i+1:i+int(match.group("nf"))+1]))
    return bs

def parse_function(shell, atomic, los):

    '''Parse basis info from list of strings'''

    exps   = [float(item.split()[1]) for item in los]
    coeffs = [float(item.split()[2]) for item in los]

    bs =[{"shell"  : _shells.index(shell.upper()),
          "typ"    : "",
          "exps"   : exps,
          "coeffs" : coeffs,
          "atomic" : atomic,
         }]
    return bs
