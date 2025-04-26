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
from math import log, cos, pi, sin, tan
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.tri as tri


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
        

def setup_problem():
    flux = 20_000 * (10 * sp.sin(PHI) - 1)
    T_i = 300.0  # K
    R_i = 2.52e-02  # m
    R_o = 3.81e-02  # m
    k = 16.24  # W/(m.K)
    cyl = CylindricalShell(R_i, R_o, T_i, k, flux)
    return cyl, R_i, R_o


def compute_temperature_field(cyl, rs, thetas):
    Ts = []
    for theta in thetas:
        for r in rs:
            Ts.append(cyl.temperature(r, theta))
    Ts = np.array(Ts, dtype=float).reshape(len(thetas), len(rs))
    return Ts


def plot_polar_temperature(rs, thetas, Ts):
    r, theta = np.meshgrid(rs, thetas)
    fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
    cax = ax.contourf(theta, r, Ts, 30, cmap='hot')
    ax.set_ylim(0, 5e-2)
    ax.set_xlim(pi/2, pi)
    cb = fig.colorbar(cax, location='bottom')
    cb.set_label("temperature, K")
    plt.savefig('cyl-T-dist.png', dpi=600)

def convert_to_cartesian(rs, thetas, Ts):
    xyt_solution = []
    for i_theta, theta in enumerate(thetas):
        for i_r, r in enumerate(rs):
            x = r * cos(theta)
            y = r * sin(theta)
            t = Ts[i_theta][i_r]
            xyt_solution.append([x, y, t])
    return xyt_solution

def extract_boundary_cartesian(rs, thetas, Ts):
    xyt_boundary = []

    # 1. Bottom edge (theta = π/2, radius OUT to IN)
    theta = thetas[0]
    for r in reversed(rs):
        x = r * cos(theta)
        y = r * sin(theta)
        t = Ts[0][np.where(rs == r)[0][0]]
        xyt_boundary.append([x, y, t])

    # 2. Inner arc (r = R_i, theta π/2 to π)
    r_index = 0
    for theta in thetas:
        x = rs[r_index] * cos(theta)
        y = rs[r_index] * sin(theta)
        t = Ts[np.where(thetas == theta)[0][0]][r_index]
        xyt_boundary.append([x, y, t])

    # 3. Top edge (theta = π, radius IN to OUT)
    theta = thetas[-1]
    for r in rs:
        x = r * cos(theta)
        y = r * sin(theta)
        t = Ts[-1][np.where(rs == r)[0][0]]
        xyt_boundary.append([x, y, t])

    # 4. Outer arc (r = R_o, theta π to π/2)
    r_index = -1
    for theta in reversed(thetas):
        x = rs[r_index] * cos(theta)
        y = rs[r_index] * sin(theta)
        t = Ts[np.where(thetas == theta)[0][0]][r_index]
        xyt_boundary.append([x, y, t])

    return xyt_boundary

def write_to_csv(data, filename='analytical_solution.csv'):
    with open(filename, 'w') as f:
        for x, y, t in data:
            f.write(f"{x},{y},{t}\n")
    print(f"Writing some calculated solution to {filename}. Format: x,y,t")

def plot_cartesian_temperature(xyt_solution):
    x_vals = [x for x, y, t in xyt_solution]
    y_vals = [y for x, y, t in xyt_solution]
    t_vals = [t for x, y, t in xyt_solution]

    plt.figure()
    plt.scatter(x_vals, y_vals, c=t_vals, cmap='hot')
    plt.gca().set_aspect('equal')
    plt.colorbar(label='Temperature (K)')
    plt.title("Temperature distribution in Cartesian (x, y) space")
    plt.xlabel('x [m]')
    plt.ylabel('y [m]')
    plt.show()

def test():
    cyl, R_i, R_o = setup_problem()
    rs = np.linspace(R_i, R_o, 20)
    thetas = np.linspace(pi/2, pi, 20)

    Ts = compute_temperature_field(cyl, rs, thetas)
    # plot_polar_temperature(rs, thetas, Ts)

    interior_temperature_solution_xyt = convert_to_cartesian(rs, thetas, Ts)
    boundary_coordinates = extract_boundary_cartesian(rs, thetas, Ts)

    write_to_csv(interior_temperature_solution_xyt, filename="interior_T_solution.csv")
    write_to_csv(boundary_coordinates, filename="boundary_representation.csv")

    plot_cartesian_temperature(boundary_coordinates)

if __name__ == '__main__':
    test()

    






