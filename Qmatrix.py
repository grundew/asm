#! /usr/bin/env python
from sympy import symbols, diff, pprint, exp, I, latex
import sympy.abc
x = sympy.abc.x
z = sympy.abc.z
w = sympy.abc.omega
rho = sympy.abc.rho
mu = sympy.abc.mu
k = sympy.abc.k
t = sympy.abc.t

Dp = symbols("Dp", each_char=False)
Up = symbols("Up", each_char=False)
Us = symbols("Us", each_char=False)
Ds = symbols("Ds", each_char=False)

gs = symbols('gs', each_char=False)
gp = symbols('gp', each_char=False)

vs = symbols('vs', each_char=False)
vp = symbols('vp', each_char=False)

phi = (Dp*exp(-I*gp*z) + Up*exp(I*gp*z))*exp(I*(w*t-k*x))
psi = (Ds*exp(-I*gs*z) + Us*exp(I*gs*z))*exp(I*(w*t-k*x))

mu = rho*vs**2

sigmazz = rho*diff(phi, t, 2) - 2*mu*(diff(phi, x, 2) + diff(psi, x, z))
sigmaxz = rho*diff(psi, t, 2) + 2*mu*(diff(phi, x, z) + diff(psi, z, 2))
sigmaxz2 = -rho*diff(psi, t, 2) + 2*mu*(diff(phi, x, z) + diff(psi, x, 2))

ux = diff(phi, x) - diff(psi, z)
uz = diff(phi, z) - diff(psi, x)

pprint(sigmazz)
pprint(sigmaxz)
pprint(ux)
pprint(uz)

#print latex(sigmazz, mode='equation')
#print latex(sigmaxz, mode='equation')
#print latex(uz, mode='equation')
#print latex(ux, mode='equation')


