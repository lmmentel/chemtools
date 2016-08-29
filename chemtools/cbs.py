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

'''Module for Complete Basis Set (CBS) Extrapolations.'''

from scipy.optimize import curve_fit
import numpy as np


def extrapolate(x, energy, method, **kwargs):
    '''
    An interface for performing CBS extrapolations using various methods.

    Args:
      x : numpy.array
        A vector of basis set cardinal numbers
      energy : numpy.array
        A vector of corresponding energies
      method : str
        Method/formula to use to perform the extrapolation

    Kwargs:
      Keyword arguments to be passed to the requested extrapolation function
      using the `method` argument
    '''

    methods = {
        'poly': poly,
        'expo': expo,
        'exposqrt': exposqrt,
        'exposum': exposum,
        'uste': uste,
    }

    if len(x) != len(energy):
        raise ValueError("x and energy should have the same size")

    if method in methods.keys():
        eopt, ecov = curve_fit(methods[method](**kwargs), x, energy)
        return eopt
    else:
        raise ValueError("wrong method: {0}, accepted values are: {1}".format(method,
            ", ".join(list(methods.keys()))))


def uste(method="CI"):
    '''
    CBS extrapolation using uniform singlet and triplet pair extrapolation
    (USTE) scheme [1]_.

    .. [1] Varandas, A. J. C. (2007). "Extrapolating to the one-electron
       basis-set limit in electronic structure calculations. The Journal of
       Chemical Physics, 126(24), 244105. `doi:10.1063/1.2741259 <http://www.dx.doi.org/10.1063/1.2741259>`_

    Args:
      x : int
        Cardinal number of the basis set
      e_cbs : float
        Approximation to the energy value at the CBS limit
      a : float
        Empirical A3 parameter
      method : str
        One of: *ci*, *cc*

    Returns:
      function object
    '''
    def uste_ci(x, e_cbs, a):
        # parameters calibrated for MRCI(Q)
        params = {
            "A05": 0.003769,
            "c": -1.1784771,
            "m": 1.25,
            "alpha": -0.375}

        a5 = params["A05"] + params['c'] * np.power(a, params["m"])
        v = e_cbs + a / np.power(x + params["alpha"], 3) + a5 / np.power(x + params["alpha"], 5)
        return v

    def uste_cc(x, e_cbs, a):
        # parameters calibrated for CC
        params = {
            "A05": 0.1660699,
            "c": -1.4222512,
            "m": 1.0,
            "alpha": -0.375}

        a5 = params["A05"] + params['c'] * np.power(a, params["m"])
        v = e_cbs + a / np.power(x + params["alpha"], 3) + a5 / np.power(x + params["alpha"], 5)
        return v

    if method.lower() == "ci":
        return uste_ci
    elif method.lower() == "cc":
        return uste_cc
    else:
        ValueError("wrong method: {}, accepted values are: 'ci', 'cc'".format(method))


def exposqrt(twopoint=True):
    '''
    Three-point formula for extrapolating the HF reference energy [2]_.

    .. [2] Karton, A., & Martin, J. M. L. (2006). Comment on: “Estimating the
       Hartree-Fock limit from finite basis set calculations” [Jensen F (2005)
       Theor Chem Acc 113:267]. Theoretical Chemistry Accounts, 115, 330–333.
       `doi:10.1007/s00214-005-0028-6 <http://www.dx.doi.org/10.1007/s00214-005-0028-6>`_

    .. math::

       E^{HF}(X) = E_{CBS} + a\cdot \exp(-b\sqrt{X})

    Args:
      twpoint : bool
        A flag marking the use of two point extrapolation with `b=9.0`

    Returns:
      funtion object
    '''

    def exposqrt2(x, e_cbs, a):
        '''
        Two-point formula for extrapolating the HF reference energy, as
        proposed by A. Karton and J. M. L. Martin, Theor. Chem. Acc. 115, 330.
        (2006)

        :math:`E^{HF}(X) = E_{CBS} + A\cdot \exp(-9\sqrt{X})`

        Args:
        x : int
            Cardinal number of the basis set
        e_cbs : float
            Approximation to the energy value at the CBS limit
        a : float
            Pre-exponential empirical parameter
        '''
        return e_cbs + a * np.exp(-9.0 * np.sqrt(x))

    def exposqrt3(x, e_cbs, a, b):
        '''
        Three-point formula for extrapolating the HF reference energy, as
        proposed by A. Karton and J. M. L. Martin, Theor. Chem. Acc. 115, 330.
        (2006)

        :math:`E^{HF}(X) = E_{CBS} + A\cdot \exp(-B\sqrt{X})`

        Args:
        x : int
            Cardinal number of the basis set
        e_cbs : float
            Approximation to the energy value at the CBS limit
        a : float
            Pre-exponential empirical parameter
        b : float
            Exponential empirical parameter
        '''
        return e_cbs + a * np.exp(-b * np.sqrt(x))

    if twopoint:
        return exposqrt2
    else:
        return exposqrt3


def expo():
    '''
    CBS extrapolation formula by exponential Dunning-Feller relation.

    .. math::

       E^{HF}(X) = E_{CBS} + a\cdot\exp(-bX)

    Returns:
      function object
    '''

    def exponential(x, e_cbs, a, b):
        '''
        CBS extrapolation formula by exponential Dunning-Feller relation.

        Args:
        x : int
            Cardinal number of the basis set
        e_cbs : float
            Approximation to the energy value at the CBS limit
        a : float
            Pre-exponential empirical parameter
        b : float
            Exponential empirical parameter
        '''

        return e_cbs + b * np.exp(-a * x)

    return exponential


def exposum():
    '''
    Three point extrapolation through sum of exponentials expression

    .. math::

       E(X) = E_{CBS} + a \cdot\exp(-(X-1)) + b\cdot\exp(-(X-1)^2)

    '''

    def exponentialssum(x, e_cbs, a, b):
        '''
        Three point extrapolation through sum of exponentials expression

        :math:`E(X) = E_{CBS} + A \cdot\exp(-(X-1)) + B\cdot\exp(-(X-1)^2)`

        Args:
        x : int
            Cardinal number of the basis set
        e_cbs : float
            Approximation to the energy value at the CBS limit
        a : float
            Empirical coefficent
        b : float
            Empirical coefficent
        '''
        return e_cbs + a * np.exp(-(x - 1)) + b * np.exp(-(x - 1)**2)

    return exponentialssum


def poly(p=0.0, z=3.0, twopoint=True):
    '''
    CBS extrapolation by polynomial relation.

    .. math::

       E(X) = E_{CBS} + \sum_{i}a_{i}\cdot (X + P)^{-b_{i}}

    Kwargs:
      twpoint : bool
        A flag for choosing the two point extrapolation
      z : float or list of floats
        Order of the polynomial, *default=3.0*
      p : float
        A parameter modifying the cardinal number, *default=0.0*
    '''

    def poly3(x, e_cbs, a, z):
        '''
        Args:
          x : int
            Cardinal number of the basis set
          e_cbs : float
            Approximation to the energy value at the CBS limit
          a : float or list of floats
            Polynomial coefficient
        '''

        a = np.array(a)
        z = np.array(z)
        return e_cbs + np.dot(a, np.power((x + p), -z))

    def poly2(x, e_cbs, a):
        '''
        Args:
          x : int
            Cardinal number of the basis set
          e_cbs : float
            Approximation to the energy value at the CBS limit
          a : float or list of floats
            Polynomial coefficient
        '''

        a = np.array(a)
        zeta = np.array(z)
        return e_cbs + np.dot(a, np.power((x + p), -zeta))

    if twopoint:
        return poly2
    else:
        return poly3
