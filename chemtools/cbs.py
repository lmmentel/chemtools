
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
    '''

    methods = {
        'poly' : poly,
        'expo' : expo,
        'exposqrt' : exposqrt,
        'exposum' : exposum,
        'uste' : uste,
    }

    if len(x) != len(energy):
        raise ValueError("x and energy should have the same size")

    if method in methods.keys():
        eopt, ecov = curve_fit(methods[method](**kwargs), x, energy)
        return eopt
    else:
        raise ValueError("wrong method: {0}, accepted values are: {1}".format(method, ", ".join(list(methods.keys()))))


def uste(method="CI"):
    '''
    CBS extrapolation using USTE shceme based on
    A. J. C. Varandas, JPCA 114, 8505-8516 (2010).

    Args:
      x : int
        Cardinal number of the basis set
      e_cbs : float
        Approximation to the energy value at the CBS limit
      a : float
        Empirical A3 parameter
      method : str
        One of: *ci*, *cc*
    '''
    def uste_ci(x, e_cbs, a):
        # parameters calibrated for MRCI(Q)
        params = {
            "A05" :  0.003769,
            "c"   : -1.1784771,
            "m"   : 1.25,
            "alpha" : -0.375}

        a5 = params["A05"] + params['c']*np.power(a, params["m"])
        v = e_cbs + a/np.power(x+params["alpha"], 3) +a5/np.power(x+params["alpha"], 5)
        return v

    def uste_cc(x, e_cbs, a):
        # parameters calibrated for CC
        params = {
            "A05" :  0.1660699,
            "c"   : -1.4222512,
            "m"   : 1.0,
            "alpha" : -0.375}

        a5 = params["A05"] + params['c']*np.power(a, params["m"])
        v = e_cbs + a/np.power(x+params["alpha"], 3) +a5/np.power(x+params["alpha"], 5)
        return v

    if method.lower() == "ci":
        return uste_ci
    elif method.lower() == "cc":
        return uste_cc
    else:
        ValueError("wrong method: {}, accepted values are: 'ci', 'cc'".format(method))

def exposqrt(twopoint=True):

    def exposqrt2(x, e_cbs, a):
        '''
        Two-point formula for extrapolating the HF reference energy, as proposed
        by A. Karton and J. M. L. Martin, Theor. Chem. Acc. 115, 330.  (2006)

        :math:`E^{HF}(X) = E_{CBS} + A\cdot \exp(-9\sqrt{X})`

        Args:
        x : int
            Cardinal number of the basis set
        e_cbs : float
            Approximation to the energy value at the CBS limit
        a : float
            Pre-exponential empirical parameter
        '''
        return e_cbs + a*np.exp(-9.0*np.sqrt(x))

    def exposqrt3(x, e_cbs, a, b):
        '''
        Three-point formula for extrapolating the HF reference energy, as proposed
        by A. Karton and J. M. L. Martin, Theor. Chem. Acc. 115, 330.  (2006)

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
        return e_cbs + a*np.exp(-b*np.sqrt(x))

    if twopoint:
        return exposqrt2
    else:
        return exposqrt3

def expo(x, e_cbs, a, b):
    '''
    CBS extrapolation formula by exponential Dunning-Feller relation.

    :math:`E^{HF}(X) = E_{CBS} + A\cdot\exp(-BX)`

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

    return e_cbs + b * np.exp(-a*x)

def exposum(x, e_cbs, a, b):
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
    return e_cbs + a*np.exp(-(x-1)) + b*np.exp(-(x-1)**2)

def poly(p=0.0, z=3.0, twopoint=True):
    '''
    CBS extrapolation by polynomial relation.

    :math:`E(X) = E_{CBS} + \sum_{i}A_{i}\cdot (X + P)^{-B_{i}}`

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

