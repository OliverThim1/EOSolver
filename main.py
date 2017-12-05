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

R = np.logspace(-7,2.5,5000) # Logarithmic mesh for our integration: [10^-7,...10^2.5]
Z=1                      # Charge of the nucleous
nmax = 5                 # Maximum principal state
lmax = 0                 # Maximum angular momentum

def pot(r):
    #return pot_coulomb(r,Z)+pot_gauss(r)
    return pot_coulomb(r,Z)+pot_step(r)

# Computes eigenstates requested by lmax and nmax
bound_states = bs.FindBoundStates(R,lmax,nmax,pot)