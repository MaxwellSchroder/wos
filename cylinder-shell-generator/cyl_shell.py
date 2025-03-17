# Author: Rowan J. Gollan
# Date: 2024-10-17
#
# This Python script computes a semi-analytical result for the steady-state temperature
# distribution in a cylindrical shell.
#
# Solution is described in:
# Ch. 4 as Example 4-4 in Hahn and Ozisik (2012)
#
# Hahn and Ozisik (2012)
# Heat Conduction, 3rd edition
# Wiley, Hoboken N.J.
#

import sympy as sp
from math import log, cos, pi

# we make these sympy symbols global because we may want to
# interact with them inside and outside of the class
PHI = sp.Symbol('phi')

class CylindricalShell:

    _phi_0 = pi/2
    _n_terms = 10

    def __init__(self, R_i, R_o, T_i, k, flux):
        self._a = R_i
        self._b = R_o
        self._T_i = T_i
        self._k = k
        self._f = flux
        self._init_lambdas()
        self._init_constants()
        return

    def _init_lambdas(self):
        self._lambda = [ i*pi/self._phi_0 for i in range(self._n_terms) ]

    def _init_constants(self):
        self._C = list(range(self._n_terms))
        # C0 is special case
        self._C[0] = self._b/(self._k*self._phi_0)*sp.integrate(self._f, (PHI, 0.0, self._phi_0))
        # all other Cs
        for i in range(1, self._n_terms):
            l_n = self._lambda[i]
            numer = (2*self._b/(self._k*l_n*self._phi_0)) * sp.integrate(self._f * sp.cos(l_n*PHI), (PHI, 0.0, self._phi_0))
            b_on_a = self._b/self._a
            denom = pow(b_on_a, l_n) + pow(b_on_a, -l_n)
            self._C[i] = numer/denom

    def temperature(self, r, theta):
        r_on_a = r/self._a
        phi = theta - pi/2
        T = self._C[0] * log(r_on_a)
        for i in range(1, self._n_terms):
            l_n = self._lambda[i]
            T += self._C[i] * (pow(r_on_a, l_n) - pow(r_on_a, -l_n)) * cos(l_n * phi)
            
        T += self._T_i
        return T
        
        

def test():
    import numpy as np
    import matplotlib.pyplot as plt
    # INPUTS
    flux = 20_000*(10*sp.sin(PHI) - 1)
    T_i =  300.0   # K
    R_i = 2.52e-02 # m, inner nose radius, 
    R_o = 3.81e-02 # m, outer nose radius
    k = 16.24 # W/(m.K), for stainless steel 321

    cyl = CylindricalShell(R_i, R_o, T_i, k, flux)

    rs = np.linspace(R_i, R_o, 20, endpoint=True)
    thetas = np.linspace(pi/2, pi, 20, endpoint=True)
    Ts = []
    for theta in thetas:
        for r in rs:
            Ts.append(cyl.temperature(r, theta))
    Ts = np.array(Ts, dtype=float)
    Ts = Ts.reshape(len(thetas), len(rs))
    print(str(Ts))

    r, theta = np.meshgrid(rs, thetas)
    fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
    cax = ax.contourf(theta, r, Ts, 30, cmap='hot')
    ax.set_ylim(0, 5e-2)
    ax.set_xlim(pi/2, pi)
    cb = fig.colorbar(cax, location='bottom')
    cb.set_label("temperature, K")

    plt.savefig('cyl-T-dist2.png', dpi=600)
    

if __name__ == '__main__':
    test()

    






