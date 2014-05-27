
import numpy as np
import chemtools.basisset as bs

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

######### TEST 1

# generate exponents from even tempered series

nf, a, b = 5, 0.01, 2.0
etzetas = bs.eventemp(nf, (a, b))
print 'eventemp zetas : ', np.around(etzetas, decimals=6)

# obtain the expansion coefficients for the 4-th order legendre expansion

Ak = et2legendre(4, etzetas)
print 'fitted coeffs :  ', Ak
print 'legendre zetas : ', np.around(bs.legendre(nf, Ak), decimals=6)


######### TEST 2

# get parameters for Yb s-type exponents

yb_s = [9934253.300000000, 10471742.50000000, 2747395.910000000, 804886.2050000000,
        258307.4000000000,  89168.3545000000,  32759.5992000000,  12669.4607000000,
          5106.3288200000,   2132.9940100000,    921.1886680000,    407.4479870000,
           186.4254560000,     88.5472069000,     43.1802249000,     21.4929947000,
            10.4678731000,      4.9894697500,      2.2641934000,      0.8266660100,
             0.3306664000,      0.1322665600,      0.0529066000,      0.0211626400,
             0.0084650500]

print 'fitted coeffs (k=2):  ', et2legendre(2, yb_s)
print 'fitted coeffs (k=3):  ', et2legendre(3, yb_s)
print 'fitted coeffs (k=4):  ', et2legendre(4, yb_s)
