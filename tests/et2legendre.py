
import numpy as np
import basisset as bs

def et2legendre(kmax, zetas):
    '''
    From a set of exponents "zetas", using least square fit calculate the
    expansion coefficients into the legendre polynomials of the order "kmax".
    '''

    c = np.asarray([1.0]*kmax, dtype=float)

    leg = np.polynomial.Legendre(c)
    a = np.zeros((len(etzetas), kmax))

    for j in xrange(len(etzetas)):
        for k in xrange(kmax):
            arg = (2.0*(j+1.0)-2.0)/(len(etzetas)-1.0)-1.0
            a[j, k] = leg.basis(k)(arg)

    return np.linalg.lstsq(a, np.log(etzetas))[0]

######### TEST


# generate exponents from even tempered series

nf, a, b = 5, 0.01, 2.0
etzetas = bs.eventemp(nf, (a, b))
print 'eventemp zetas : ', np.around(etzetas, decimals=6)

# obtain the expansion coefficients for the 4-th order legendre expansion

Ak = et2legendre(4, etzetas)
print 'fitted coeffs :  ', Ak
print 'legendre zetas : ', np.around(bs.legendre(nf, Ak), decimals=6)
