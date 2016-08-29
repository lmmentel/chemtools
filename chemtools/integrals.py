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

from __future__ import print_function

import numpy as np
from scipy.special import factorial2
from chemtools.basisset import ncartesian, get_l
from numba import jit, int32, float64
from scipy.linalg import sqrtm, inv


@jit(float64(int32, float64))
def norm(n, a):
    '''
    Calculate the normalization factor for a cartesian gaussian of the form

    .. math::

       G_{n}(x, a, A_{x}) = x^{n}_{A}\exp(-ax^{2}_{A})

    Args:
      n : int
        Power of the preexponential factor
      a : float
        Exponent

    Returns:
      out : float
        Value of the normalization factor
    '''

    # first option
    #return np.sqrt(np.power(2.0*a, n + 0.5, dtype=np.float64)/gamma(n + 0.5), dtype=np.float64)
    # second option
    return np.sqrt(np.power(4.0*a, n)*np.sqrt(2.0*a/np.pi, dtype=np.float64)/factorial2(2*n - 1), dtype=np.float64)
    # third option
    #return np.sqrt(np.power(8.0*a, n)*factorial(n)*np.sqrt(2.0*a/np.pi)/factorial(2*n), dtype=np.float64)

@jit(float64(int32, int32, float64, float64, float64, float64))
def obara_saika(i, j, a, b, Ax, Bx):
    '''
    Calculates the overlap integral between two  cartesian gaussian components using
    the Obara-Saika recurrence relations.

    Args:
      i : int
        Power of the preexponential monomial in the first gaussian
      j : int
        Power of the preexponential monomial in the second gaussian
      a : float
        Exponent of the first gaussian
      b : float
        Exponent of the second gaussian
      Ax : float
        Coordinate of the first gaussian
      Bx : float
        Coordiante of the second gaussian

    Returns:
      out : float
        The value of the :math:`S_{i,j}` integral
    '''

    p = a + b
    mu = a * b / p
    Px = (a * Ax + b * Bx) / p
    X_AB = Ax - Bx
    X_PA = Px - Ax
    X_PB = Px - Bx

    if i == 0 and j == 0:
        return np.sqrt(np.pi / p) * np.exp(-mu * X_AB**2)
    elif i == 0 and j == 1:
        return X_PB * obara_saika(0, 0, a, b, Ax, Bx)
    elif i == 1 and j == 0:
        return X_PA * obara_saika(0, 0, a, b, Ax, Bx)
    elif i == 1 and j == 1:
        return X_PA * obara_saika(0, 1, a, b, Ax, Bx) + \
                obara_saika(0, 0, a, b, Ax, Bx)/(2.0*p)
    elif i == 0 and j > 1:
        return X_PB * obara_saika(0, j - 1, a, b, Ax, Bx) + \
                (j-1)*obara_saika(0, j - 2, a, b, Ax, Bx)/(2.0*p)
    elif i == 1 and j > 1:
        return X_PB * obara_saika(1, j - 1, a, b, Ax, Bx) + \
                (obara_saika(0, j-1, a, b, Ax, Bx) + \
                (j - 1) * obara_saika(1, j-2, a, b, Ax, Bx))/(2.0*p)
    elif i > 1 and j == 0:
        return X_PA * obara_saika(i-1, 0, a, b, Ax, Bx) + \
                (i - 1)*obara_saika(i-2, 0, a, b, Ax, Bx)/(2.0*p)
    elif i > 1 and j == 1:
        return X_PA*obara_saika(i-1, 1, a, b, Ax, Bx) + \
                ((i-1)*obara_saika(i-2, 1, a, b, Ax, Bx) + \
                obara_saika(i-1, 0, a, b, Ax, Bx))/(2.0*p)
    else:
        return X_PA*obara_saika(i-1, j, a, b, Ax, Bx) + \
                ((i-1)*obara_saika(i-2, j, a, b, Ax, Bx) + \
                j*obara_saika(i-1, j-1, a, b, Ax, Bx))/(2.0*p)


def get_basinfo(bases, positions, order='canonical'):
    '''
    Compose a record array with all of the exponents/functions, The record has
    the follwing structure

      * ``exp``,  ``float``, the exponent
      * ``shell``, ``str``, the shell
      * ``l``, ``int``, the angular momentum quantum number
      * ``ix``, ``int``, the power fo the ``x`` monomial
      * ``iy``, ``int``, the power fo the ``y`` monomial
      * ``iz``, ``int``, the power fo the ``z`` monomial
      * ``atom``, ``str``, symbol of the element
      * ``x``, ``float``, the x coordinate
      * ``y``, ``float``, the y coordinate
      * ``z``, ``float``, the z coordinate

    Args:
      bases : list
        A list of ``BasisSet`` objects
      positions : list
        List of 3-tuples with the (x, y, z) coordinates
      order : str
        Order in which the components are sorted, default is ``canonical``, meaning
        [*xx*, *xy*, *xz*, *yy*, *yz*, *zz*] for the *d* shell, other choice is ``gamessus``
        where the highest powers are listed first [*xx*, *yy*, *zz*, *xy*, *xz*, *yz*]

    Returns:
      out : numpy.recarray
        Numpy record array with all available functions
    '''

    dtype = [('exp', np.float64), ('shell', 'S1'), ('l', np.int32),
             ('ix', np.int32), ('iy', np.int32), ('iz', np.int32),
             ('atom', 'S3'), ('x', np.float64), ('y', np.float64), ('z', np.float64)]

    nf = sum([bs.nprimitive(spherical=False) for bs in bases])
    out = np.recarray((nf,), dtype=dtype)

    index = -1
    for basis, atom in zip(bases, positions):
        for shell, fs in basis.functions.items():
            l = get_l(shell)
            for ex in fs['e']:
                for ixyz in xyzlist(l, order=order):
                    index += 1
                    out[index] = (ex, shell, l, ixyz[0], ixyz[1], ixyz[2],
                                  basis.element, atom[0], atom[1], atom[2])
    return out


def xyzlist(l, order='canonical'):
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
      order : str
        Order in which the components are sorted, default is ``canonical``, meaning
        [*xx*, *xy*, *xz*, *yy*, *yz*, *zz*] for the *d* shell, other choice is ``gamessus``
        where the highest powers are listed first [*xx*, *yy*, *zz*, *xy*, *xz*, *yz*]

    Returns:
      out ((l+1)*(l+2)/2, 3) : numpy.array
        Array of monomial powers, where columns correspond to *x*, *y*, and *z*
        respectively and rows correspond to functions.
    '''

    ncart = ncartesian(l)
    out = np.zeros((ncart, 3), dtype=np.int32)

    index = 0
    for i in range(l + 1):
        for j in range(i + 1):
            out[index,:] = [l - i, i - j, j]
            index += 1
    if order == 'canonical' or l == 0:
        return out
    elif order == 'gamessus':
        # indices of the second largest value in each row
        imax2 = np.argsort(out, axis=1)[:, 1]
        # indices of largest elements in each row
        imax = np.argmax(out, axis=1)
        # negative sum of two largest values
        nmax2 = -np.sum(np.sort(out, axis=1)[:, 1:], axis=1)
        # negative largest values in each row
        nmax = -np.max(out, axis=1)
        return out[np.lexsort((imax2, imax, nmax2, nmax))]
    else:
        raise ValueError('Unknown <order> value: {}'.format(order))


@jit
def primitive_overlap(exps):
    '''
    Calculate the overlaps between the primitive cartesian gaussians

    Args:
      exps (N,) : numpy.recarray
        Numpy record array with available functions, for details see ``get_exps``

    Returns:
      S (N, N) : numpy.array
        Overlap integral matrix between the primitives
    '''

    # initialize the overlap matrix
    S = np.zeros((exps.size, exps.size), dtype=np.float64)

    # this should be done is some smarter way avoiding the explicit loop indexing
    # maybe np.nditer?

    for i in range(exps.size):

        irow = exps[i]
        iix, iiy, iiz, iexp, xi, yi, zi = irow['ix'], irow['iy'], irow['iz'], irow['exp'], irow['x'], irow['y'], irow['z']

        normi = norm(iix, iexp)*norm(iiy, iexp)*norm(iiz, iexp)

        for j in range(i + 1):
            normj = norm(exps[j]['ix'], exps[j]['exp'])*\
                    norm(exps[j]['iy'], exps[j]['exp'])*\
                    norm(exps[j]['iz'], exps[j]['exp'])

            x = obara_saika(iix, exps[j]['ix'], iexp, exps[j]['exp'], xi, exps[j]['x'])
            y = obara_saika(iiy, exps[j]['iy'], iexp, exps[j]['exp'], yi, exps[j]['y'])
            z = obara_saika(iiz, exps[j]['iz'], iexp, exps[j]['exp'], zi, exps[j]['z'])

            S[i, j] = normi*normj*x*y*z

    S = S + S.T
    np.fill_diagonal(S, S.diagonal()/2.0)
    return S


def primitive_overlap_between(bas1, bas2):

    # initialize the overlap matrix
    S = np.zeros((bas1.size, bas2.size), dtype=np.float64)

    for i in range(bas1.size):
        normi = norm(bas1[i]['ix'], bas1[i]['exp'])*norm(bas1[i]['iy'], bas1[i]['exp'])*norm(bas1[i]['iz'], bas1[i]['exp'])
        for j in range(bas2.size):
            normj = norm(bas2[j]['ix'], bas2[j]['exp'])*norm(bas2[j]['iy'], bas2[j]['exp'])*norm(bas2[j]['iz'], bas2[j]['exp'])

            x = obara_saika(bas1[i]['ix'], bas2[j]['ix'], bas1[i]['exp'], bas2[j]['exp'], bas1[i]['x'], bas2[j]['x'])
            y = obara_saika(bas1[i]['iy'], bas2[j]['iy'], bas1[i]['exp'], bas2[j]['exp'], bas1[i]['y'], bas2[j]['y'])
            z = obara_saika(bas1[i]['iz'], bas2[j]['iz'], bas1[i]['exp'], bas2[j]['exp'], bas1[i]['z'], bas2[j]['z'])

            S[i, j] = normi*normj*x*y*z

    return S


def contraction_matrix(bases):
    '''
    Construct and return the matrix of contraction coefficents

    Args:
      bases : list
        List of BasisSet object

    Returns:
      cc : numpy.array
        Matrix with contraction coefficients
    '''

    # calcualte the number of primitives
    nprim = sum([bs.nprimitive(spherical=False) for bs in bases])
    # calculate the number of contracted functions
    nc = sum([bs.nf(spherical=False) for bs in bases])

    cc = np.zeros((nprim, nc), dtype=np.float64)
    icol = 0

    irow = 0
    for basis in bases:
        for shell, fs in basis.functions.items():
            l = get_l(shell)
            ncart = ncartesian(l)
            nprs = len(fs['e']) * ncart
            for cf in fs['cf']:
                # calculate the indices of rows for all the components
                rowidx = np.asarray([[i * ncart + n  for i in cf['idx']] for n in range(ncart)]).T
                # reshape the contraction coefficients to assign them to elements of the contraction matrix
                cc[rowidx + irow, np.arange(ncart) + icol] = np.repeat(cf['cc'], ncart).reshape(cf.size, ncart)
                icol += ncart
            irow += nprs
    return cc


def completeness_profile(basinfo, cc, zetas):
    '''
    Calculate the completeness profile of each shell of the basis set

    Args:
      zetas : numpy.array
        Scaning exponents

    Returns:
      out : numpy.array
        numpy array with values of the profile (shells, zetas)
    '''

    # calculate the overlap matrix
    Sprim = primitive_overlap(basinfo)
    S = np.dot(cc.T, np.dot(Sprim, cc))
    # calculate the inverse square root of S
    X = inv(sqrtm(S))
    SO = primitive_overlap_between(basinfo, zetas)
    J = np.dot(np.dot(cc, X).T, SO)
    out = np.sum(np.square(J), axis=0)

    return out
