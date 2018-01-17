# This program is based on a program by Christian Forssen.

import numpy as np
import pylab as p
from scipy import *
from scipy import optimize
import matplotlib as plt


def _NumerovSolve(g,s,R,h,start1=0,start2=1e-4):
    """ Solves differential equation of the form y'' = -g(x)y + s(x)
        for every point in R using numerovs method.
        <Start2> only scales the solution.

        TODO: Improve array system to numpy.array.
    """
    # Solves backwards so flip R.
    R0 = R[::-1]

    y=[start1,start2] # Function value array start from the rear.

    for i in range(1,len(R)-1):
        xn = R[i]
        xnm = R[i-1]
        xnp = R[i+1]

        # Numerovs method.
        x = (  2.*y[i]*(1-((5/12)*(h**2))*g(xn)) - y[i-1]*(1+((h**2)/12)*g(xnm)) + ((h**2)/12)*(s(xnp) + 10*s(xn) + s(xnm))  ) / (1 + ((h**2)/12)*g(xnp))
        y.append(x)

    # Take the norm of the y-vector; divide with it's maximum.
    norm = max(y[:int(len(R) * 0.9)])

    y = y[::-1] # Flips back.
    u = [x / norm for x in y] # Normalizes.
    return u




def _SolveSchroedinger(E,R,l,pot):
    """ Integrates Schroedinger equation given a energy E, and angular momentum
    eigenvalue l, and potential pot(r).
    For stability, we integrate Schroedinger equation starting at infinity
    and integrating down to zero.
    We use initial condition at infinity : u(infinity)=0, and u'(infinity)=small.
    In order that solution is of the order of unity, we normalize it to (roughly) the peak
    value of the function.

    TODO: Find functions g(x), s(x), so that we can solve the diff.equation with _NumerovSolve.
    """

    """ After the variable change u(r)/r = R(r) we have a differential equation of the form

    If we consider m/ \bar(h)^2 = 1 we get

    u'' = (2*(V-E) + l(l+1)/r^2)*u
    _NumerovSolve wants a differential equation of the form y'' = -g(x)y + s(x), so we construct g(x), s(x).
    """
    h = 0.01 # Steglängd

    def g(r):
        return (2 * (E - pot(r)) - l*(l + 1) / r ** 2)

    def s(r):
        return 0

    return _NumerovSolve(g, s, R, h)


def _Shoot(E, R, l, pot):
    """ Given a energy E (momentum eigenvalue l, and potential pot(r)), and initial condition
        at infinity (u(infinity)=0, u'(infinity)=small) it finds the value of u(r=0).
        For the state to be bound state, we have the following condition: u(r=0)==0.
    """

    u = _SolveSchroedinger(E, R, l, pot)

    # Vi vill ha värdet i r=0, därför interpolerar vi sista biten i.o.m att
    # R[0] inte riktigt.
    return u[0] + (u[1] - u[0]) * (0.0 - R[0]) / (R[1] - R[0])

def FindBoundStates(R,lmax,nmax,pot,MaxSteps=1000,Etol=1e-17):
    """ For each momentum eigenvalue l=[0,...lmax], it find first nmax-l bound-states.
    It first brackets roots of function Shoot using a logarithmic mesh dense at zero.
    It later refinds the roots by optimize.brentq method.
    """
    Eb=[]
    for l in range(lmax+1): # For each momentum eigenvalue l
        nfound=0
        while nfound < nmax-l: # Find nmax bound states
            E_coul = -1/(2.*(nfound+1+l)**2) #Why? Pure coulomb energy
            E0=E_coul
            try: dE = E_coul - 0.99*Eb0
            except: dE = 3*np.abs(E_coul)
            i=0
            uo = _Shoot(E0, R, l, pot)
            un = _Shoot(E0-dE, R, l, pot)
            while un*uo > 0 and i < MaxSteps:
                i+=1
                E0-=dE
                un = _Shoot(E0-dE, R, l, pot)

            Eb0 = optimize.brentq(_Shoot, E0-dE, E0, args=(R, l, pot), xtol=Etol)
            Eb.append([l,Eb0])
            print ('Bound state # %2i (l=%2i): E=%12.9f (Pure Coulomb: E_c=%12.9f; ratio=%6.4f)' \
            % (nfound, l, Eb0, E_coul, Eb0 / E_coul))
            nfound+=1
    return Eb




        # Find a higher and lower boundry for E.

        # Find E so u(r=0) = 0.

        # Print result


    # Return the bound states
