# -*- coding: utf-8 -*-
"""
Created on Wed Jan  7 14:02:28 2015

@author: dvartany
"""

import numpy as np
from scipy.optimize import brentq
from scipy.interpolate import RectBivariateSpline

from constants import *
from utils import pickle_results

def lam(r, q, f):
    return f*q**2*G*M/a*(a/(r-a))**4

def Omega(r):
    return np.sqrt(G*M/r**3)

def ftid(r, Sigma, q, f, off=False):
    if off:
        return 0.0
    else:
        return 0.5*(OmegaIn - Omega(r))*lam(r, q, f)*Sigma

def fv(r, T, Sigma):
    return 1.125 * Omega(r)*alpha*k*T/mu * Sigma

def Tirr(r):
    return (((eta/7.0)*L/(4*np.pi*sigma))**2* k/(G*M*mu))**(1.0/7.0)*r**(-3.0/7.0)

def func(T, r, Sigma, q, f):
    return sigma*T**4 - 3*(kappa0*T**0.5*Sigma*0.125 + 2/(kappa0*Sigma*T**0.5))*(ftid(r,Sigma,q,f)\
           + fv(r,T,Sigma)) - sigma*Tirr(r)**4

def func1(T, r, Sigma, q, f):
    return sigma*T**4 - 3*(kappa1*T**(-7)*Sigma*0.125 + 2/(kappa1*Sigma*T**(-7)))*(ftid(r,Sigma,q,f)\
           + fv(r,T,Sigma)) - sigma*Tirr(r)**4

def func2(T, r, Sigma, q, f):
    return sigma*T**4 - 3*(kappa2*T**2*Sigma*0.125 + 2/(kappa2*Sigma*T**2))*(ftid(r,Sigma,q,f)\
           + fv(r,T,Sigma)) - sigma*Tirr(r)**4

def Tfin(Tcheck, r, Sigma, q, f):
    if Tcheck > 202.677:
        return Tcheck
    else:
        return brentq(func1,1,204,args=(r,Sigma,q,f))
        if 166.81 < brentq(func1,1,204,args=(r,Sigma,q,f)) < 202.677:
            return brentq(func1,1,204,args=(r,Sigma,q,f))
        else:
            return brentq(func2,1,204,args=(r,Sigma,q,f))

def buildTempTable(rGrid, q=1.0, f=0.001, Tmin=1.0, Tmax=500000, Sigmin=1.0e-3, Sigmax=1500, Sigres=1000):
    """
    Return a table of precomputed temperatures as a function of radius and surface density.

    Arguments:
    rGrid: Grid of radii to compute the temperatures for, it has to be in cgs units.

    Keywords:
    q: Ratio of the masses of the stars in the binary, between 0 and 1.
    f: Fudge factor used to tune the magnitude of the tidal torque.
    Tmax: Maximum temperature to use as a starting point for the upper bound of brentq.
    Tmin: Minimum temperature to use as a starting point for the lower bound of brentq.
    Sigmin: Minimum value for the surface density to compute temperatures for.
    Sigmax: Maximum value for the surface density to compute temperatures for.
    Sigres: Number of points to use for surface density in the table, it will be logarthmically spaced
            in surface density.
    """

    SigmaGrid = np.logspace(np.log10(Sigmin), np.log10(Sigmax), Sigres)
    temp = np.zeros((len(rGrid), Sigres)) #create m x n array for temp
    for i, r in enumerate(rGrid):
        for j, Sigma in enumerate(SigmaGrid):
            Tcheck = brentq(func, Tmin, Tmax, args=(r,Sigma,q,f))
            temp[i,j] = Tfin(Tcheck, r, Sigma, q, f)
    # Return values in logspace for interpolation
    return np.log10(rGrid), np.log10(SigmaGrid), np.log10(temp)

@pickle_results("interpolator.pkl")
def buildInterpolator(r, gamma, q, fudge, mDisk, **kargs):
    # Keep in mind that buildTemopTable() returns the log10's of the values
    rGrid, SigmaGrid, temp = buildTempTable(r*a*gamma, q=q, f=fudge, **kargs)
    # Go back to dimensionless units
    rGrid -= np.log10(a*gamma)
    SigmaGrid -= np.log10(mDisk*M/gamma**2/a**2)
    # Interpolate in the log of dimensionless units
    return RectBivariateSpline(rGrid, SigmaGrid, temp)
