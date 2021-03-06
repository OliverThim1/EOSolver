# This program is based on a program by Christian Forssen.

import numpy as np
import scipy
import matplotlib.pyplot as plt
import radialSolver.radialSolver as rs

def heaviside(x):
    return 0.5*(np.sign(x) + 1)

def pot_step(r,V0=1,c=1):
    return -V0 * heaviside(c-r)

def pot_coulomb(r,Z=1):
    return -Z/r

def pot_gauss(r,V0=1,c=1):
    return -V0*np.exp(-r**2/(2*c**2))

R = np.logspace(-16,5.75,5000,1,np.exp(1))  # Logarithmic mesh for our integration: [10^-7,...10^2.5]
#R = np.linspace(1e-7, 250,10000)
Z=1                      # Charge of the nucleous
nmax = 5            # Maximum principal state
lmax = 1                 # Maximum angular momentum

def pot(r):
    return pot_coulomb(r,Z)+pot_gauss(r)
    #return pot_coulomb(r,Z)+pot_step(r)

# Computes eigenstates requested by lmax and nmax
#bound_states = bs.FindBoundStates(R,lmax,nmax,pot)


"""# Test av Numerov Solve.
def g(x):
    return 1
def s(x):
    return 0

y = rs._NumerovSolve(g,s,np.linspace(0,10,1000),0.01)

plt.plot(R, y)
plt.show()
"""

bound_states = rs.FindBoundStates(R,lmax,nmax,pot)

