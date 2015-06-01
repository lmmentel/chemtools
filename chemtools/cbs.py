
'''Module for Complete Basis Set Extrapolations.'''

import numpy as np

def uste(x, e_cbs, A3, method="CI"):

    '''CBS extrapolation using USTE shceme based on
    A. J. C. Varandas, JPCA 114, 8505-8516 (2010).

    Args:
      x : int
        Cardinal number of the basis set
      e_cbs : float
        Approximation to the energy value at the CBS limit
      A3 : float
        Empirical parameter
      method : str
        One of: *CI*, *CC*
    '''

    # parameters calibrated for MRCI(Q)
    ci_params = {
        "A05" :  0.003769,
        "c"   : -1.1784771,
        "m"   : 1.25,
        "alpha" : -0.375}
    # parameters calibrated for CC
    cc_params = {
        "A05" :  0.1660699,
        "c"   : -1.4222512,
        "m"   : 1.0,
        "alpha" : -0.375}

    if method == "CI":
        dd = ci_params.copy()
    elif method == "CC":
        dd = cc_params.copy()

    A5 = dd["A05"] + c*math.pow(A3, dd["m"])

    v = e_cbs + A3/np.power(x+dd["alpha"], 3) +A5/np.power(x+dd["alpha"], 5)

    return v

def km2(x, e_cbs, a):
    '''
    Two-point formula for extrapolating the HF reference energy, as proposed by
    A. Karton and J. M. L. Martin, Theor. Chem. Acc. 115, 330.  (2006)

    :math:`E^{HF}(X) = E_{CBS} + a\cdot (X+1)\cdot \exp(-9\sqrt{X})`

    Args:
      x : int
        Cardinal number of the basis set
      e_cbs : float
        Approximation to the energy value at the CBS limit
      a : float
        Pre-exponential empirical parameter
    '''

    return e_cbs + a* (x + 1)*np.exp(-9.0*np.sqrt(x))

def km3(x, e_cbs, a, b):
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

def helgaker_X3(x, a, b):

    '''Inverse power extrapolation based on T. Helgaker et. al. JCP 106(23),
    9639 (1997).'''

    return a + b * np.power(x,-3)

def cbs_helgaker(x, e_cbs, a, b):

    '''Helgaker type extrapolationn to CBS with 3 parameters.

    :math:`E^{HF}(X) = E_{CBS} + A\cdot X^{-B}`

   Args:
       x : int
        Cardinal number of the basis set
      e_cbs : float
        Approximation to the energy value at the CBS limit
      b : float
        Polynomial coefficient
      a : float
        Order of the polynomial
    '''

    return e_cbs + a*np.power(x, -b)

def scf_exp(x, e_cbs, a, b):
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

def scf_pol(x, e_cbs, b, a):
    '''
    CBS extrapolation by polynomial relation.

    :math:`E^{HF}(x) = E_{CBS} + b\cdot x^{-a}`

    Args:
      x : int
        Cardinal number of the basis set
      e_cbs : float
        Approximation to the energy value at the CBS limit
      b : float
        Polynomial coefficient
      a : float
        Order of the polynomial
    '''

    return e_cbs + b * np.power(x, -a)

