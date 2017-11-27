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
import re

import numpy as np
from mendeleev import element

CFDTYPE = [('idx', np.int32), ('cc', np.float64)]
ORBITALS = ('s', 'p', 'd', 'f', 'g', 'h', 'i', 'k')


def get_l(shell):
    'Return the orbital angular momentum quantum number for a given subshell'

    if shell in ORBITALS:
        return ORBITALS.index(shell.lower())
    else:
        raise ValueError('"{}" is not a proper shell label'.format(shell))


def parse_basis(string, fmt=None):
    '''
    A wrapper for parsing the basis sets in different formats.

    Args:
      string : str
        A string with the basis set

      fmt : str
        Format in which the basis set is specified

    Returns:
      out : dict
        A dictionary of parsed basis sets with element symbols as keys and
        basis set functions as values
    '''

    formats = ['molpro', 'gamessus', 'gaussian']

    if fmt == 'molpro':
        return parse_molpro_basis(string)
    elif fmt == 'gamessus':
        return parse_gamessus_basis(string)
    elif fmt == 'gaussian':
        return parse_gaussian_basis(string)
    else:
        raise ValueError("Can only parse <fmt>: {}".format(", ".join(formats)))


def parse_molpro_basis(string):
    '''
    Parse basis set from a string in Molpro format.
    '''

    bas_re = re.compile(r'basis\s*=\s*\{(.*?)\}', flags=re.DOTALL | re.I)

    m = bas_re.search(string)
    if m:
        lines = m.group(1).split("\n")
    else:
        raise ValueError('basis string: "basis={*}" not found')

    start = []
    for i, line in enumerate(lines):
        if line.split(",")[0].lower() in ORBITALS:
            start.append(i)
    if len(start) == 0:
        return None

    startstop = []
    for i in range(len(start) - 1):
        startstop.append((start[i], start[i + 1]))
    startstop.append((start[-1], len(lines)))

    bs = {}
    for i in startstop:
        at_symbol, shell = parse_molpro_shell(lines[i[0]],
                                              lines[i[0] + 1: i[1]])
        if at_symbol in bs.keys():
            bs[at_symbol] = dict(list(bs[at_symbol].items()) +
                                 list(shell.items()))
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

    shell = expsline.split(",")[0].lower()
    at_symbol = expsline.split(",")[1].strip().capitalize()
    exps = np.array([float(real.sub('E', x))
                     for x in expsline.rstrip(";").split(",")[2:]])

    fs[shell] = {'e': exps, 'cf': []}
    if len(coeffs) != 0:
        for line in coeffs:
            lsp = line.rstrip(";").split(",")
            if lsp[0] == "c":
                i, j = [int(x) for x in lsp[1].split(".")]
                coeffs = [float(real.sub('E', x)) for x in lsp[2:]]
                fs[shell]['cf'].append(np.array(list(zip(list(range(i - 1, j)),
                                                coeffs)), dtype=CFDTYPE))
    else:
        for i in range(len(exps)):
            fs[shell]['cf'].append(np.array([tuple([i, 1.0])], dtype=CFDTYPE))
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
    for i in range(len(start) - 1):
        startstop.append((start[i], start[i + 1]))
    startstop.append((start[-1], len(lines)))

    ecp = {}
    for i in startstop:
        ecp = dict(list(ecp.items()) +
                   list(parse_coeffs(lines[i[0]: i[1]]).items()))
    return ecp


def parse_coeffs(lines):

    firstl = lines[0].replace(';', '').split(',')
    element = firstl[1].strip().capitalize()
    nele = firstl[2]
    lmax = firstl[3]

    liter = iter(x for x in lines[1:] if x != '')

    ecp = {element: {"nele": nele, "lmax": lmax, "shells": []}}

    while True:
        try:
            temp = next(liter)
        except StopIteration:
            break
        nlmax = int(temp.split(";")[0])
        comment = temp.split(";")[1].replace("!", "")
        tt = {'comment': comment, 'parameters': []}
        for i in range(nlmax):
            param = next(liter).replace(";", "").split(",")
            tt['parameters'].append({'m': float(param[0]),
                                     'gamma': float(param[1]),
                                     'c': float(param[2])})
        ecp[element]['shells'].append(tt)
    return ecp


def parse_gaussian_basis(string):
    '''
    Parse the basis set into a list of dictionaries from a string in
    gaussian format.
    '''

    shellre = re.compile(r'^\s*(?P<shells>[SPDFGHILspdfghil]+)\s*(?P<nf>[1-9]+)\s*(?P<scale>\-?\d+\.\d+)')
    out = dict()
    for i, item in enumerate(string.split('****\n')):
        if len(item) > 5:
            atomline, basis = item.split('\n', 1)
            atom = atomline.split()[0]
            bslines = basis.split('\n')
            functions = dict()
            for i, line in enumerate(bslines):
                match = shellre.search(line)
                if match:
                    shells, nf, scale = match.group('shells').lower(), match.group('nf'), match.group('scale')
                    exps, indxs, coeffs = parse_gaussian_function(bslines[i + 1: i + int(nf) + 1])
                    for shell, cc in zip(shells, coeffs.T):
                        if shell in functions.keys():
                            sexp, idxs, idxo = merge_exponents(functions[shell]['e'], exps)
                            functions[shell]['e'] = sexp * float(scale)**2
                            for cf in functions[shell]['cf']:
                                cf['idx'] = idxs[cf['idx']]
                            newcf = np.array(list(zip(idxo[indxs], cc)),
                                             dtype=CFDTYPE)
                            functions[shell]['cf'].append(newcf)
                        else:
                            functions[shell] = dict()
                            functions[shell]['cf'] = list()
                            functions[shell]['e'] = exps
                            functions[shell]['cf'].append(np.array(list(zip(indxs, cc)), dtype=CFDTYPE))
                out[atom] = functions
    return out


def parse_gaussian_function(lines):
    '''
    Parse a basis set function information from list of strings into
    three lists containg: exponents, indices, coefficients.

    Remeber that python doesn't recognise the `1.0d-3` format where `d` or
    `D` is used to the regex subsitution has to take care of that.
    '''

    real = re.compile(r'[dD]')

    indxs = np.arange(len(lines), dtype=np.int32)
    exps = np.array([float(real.sub('E', line.split()[0])) for line in lines],
                    dtype=np.float64)
    coeffs = np.array([[float(real.sub('E', x))
                        for x in line.split()[1:]] for line in lines],
                      dtype=np.float64)
    return (exps, indxs, coeffs)


def parse_gamessus_basis(string):
    '''
    Parse the basis set into a list of dictionaries from a string in
    gamess format.
    '''

    bas_re = re.compile(r'\$DATA\n(.*?)\$END', flags=re.DOTALL | re.IGNORECASE)

    m = bas_re.search(string)
    if m:
        basisstr = m.group(1)
    else:
        raise ValueError("basis not found, should be inside '$DATA' and '$END' group")

    pat = re.compile(r'^\s*(?P<shells>[SPDFGHILspdfghil]+)\s*(?P<nf>[1-9]+)')
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
                        shells, nf = match.group("shells").lower(), match.group("nf")
                        exps, indxs, coeffs = parse_gamessus_function(bslines[i + 1: i + int(nf) + 1])
                        if shells in ['L', 'l']:
                            shells = ['s', 'p']
                        for shell, cc in zip(shells, coeffs.T):
                            if shell in functions.keys():
                                sexp, idxs, idxo = merge_exponents(functions[shell]['e'], exps)
                                functions[shell]['e'] = sexp
                                for cf in functions[shell]['cf']:
                                    cf['idx'] = idxs[cf['idx']]
                                newcf = np.array(list(zip(idxo[indxs-1], cc)),
                                                 dtype=CFDTYPE)
                                functions[shell]['cf'].append(newcf)
                            else:
                                functions[shell] = dict()
                                functions[shell]['cf'] = list()
                                functions[shell]['e'] = exps
                                functions[shell]['cf'].append(np.array(list(zip(indxs - 1, cc)), dtype=CFDTYPE))
                res[elem.symbol] = functions
    return res


def parse_gamessus_function(lines):
    '''
    Parse a basis set function information from list of strings into
    three lists containg: exponents, indices, coefficients.

    Remeber that python doesn't recognise the `1.0d-3` format where `d` or
    `D` is used to the regex subsitution has to take care of that.
    '''

    real = re.compile(r'[dD]')

    indxs = np.array([int(line.split()[0]) for line in lines], dtype=np.int32)
    exps = np.array([float(real.sub('E', line.split()[1])) for line in lines],
                    dtype=np.float64)
    coeffs = np.array([[float(real.sub('E', x)) for x in line.split()[2:]]
                        for line in lines], dtype=np.float64)
    return (exps, indxs, coeffs)


def merge_exponents(a, b):
    '''
    Concatenate the arrays `a` and `b` using only the unique items from both
    arrays

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


class NumpyEncoder(json.JSONEncoder):

    def default(self, obj):
        '''
        If input object is an ndarray it will be converted into a dict
        holding the data and dtype.
        '''

        if isinstance(obj, np.ndarray):

            if obj.dtype == CFDTYPE:
                dtype = 'CFDTYPE'
            else:
                dtype = str(obj.dtype)

            return dict(data=obj.tolist(), dtype=dtype)
        # Let the base class default method raise the TypeError
        return json.JSONEncoder.default(self, obj)
