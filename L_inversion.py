#! /usr/bin/env python
from sympy import symbols, pprint, exp, I, latex, Matrix
import sympy.abc
x = sympy.abc.x
z = sympy.abc.z
w = sympy.abc.omega
rho = sympy.abc.rho
alpha = sympy.abc.alpha
beta = sympy.abc.beta
xi = sympy.abc.xi
gamma = sympy.abc.gamma
mu = sympy.abc.mu

L = Matrix([
    [I*xi, I*xi, -I*beta, I*beta],
    [I*alpha, -I*alpha, I*xi, I*xi],
    [2*mu*xi*gamma, 2*mu*xi*gamma, -2*mu*xi*beta, 2*mu*xi*beta],
    [-2*mu*xi*alpha, 2*mu*xi*alpha, -2*mu*xi*gamma, -2*mu*xi*gamma]
    ])

Linv = L.inv()
Linv.simplify()

pprint(L)
pprint(Linv)
