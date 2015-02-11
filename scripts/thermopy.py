# -*- coding: utf-8 -*-
"""
    Created on Wed Jan  7 14:02:28 2015
    
    @author: dvartany
    """

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq
import scipy.optimize as optimize
# Constants in cgs

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

def func(T, r, Sigma, q, f, kappa):
    return sigma*T**4 - (3*op(T , r, Sigma, kappa)*T**0.5*Sigma*0.0625 + 2/(op(T, r, Sigma,kappa)*Sigma*T**0.5))*(ftid(r,Sigma,q,f) + fv(r,T,Sigma)) - sigma*Tirr(r)**4

def op(T, r, Sigma, kappa):
    if kappa == 1:
        return 0.0125987 * T**1.5
    elif kappa == 2:
        return 1.96231 * 10**8 *(Omega(r) * Sigma *(k/mu)**0.5)**0.0949916
    elif kappa == 3:
        return 0.00144313* T**1.5
    elif kappa == 4:
        return 5.01187*10**17/T**5.832
    elif kappa == 5:
        return 1.14868 *10**-6 *T**2.129
    elif kappa == 6:
        return 1.2583*10**135*(Omega(r) * Sigma*(k/mu)**0.5)**1.312/T**42.324
    elif kappa == 7:
        return 9.71628 *10**-16*T**4.0625
    elif kappa == 8:
        return 8.51138*10**58*(Omega(r) * Sigma *(k/mu)**0.5)**0.676/T**18.142
    elif kappa == 9:
        return 1.01158*10**-14*(Omega(r) * Sigma *(k/mu)**0.5)**0.498*T**3.154
    elif kappa == 10:
        return 1.15878*10**-41*(Omega(r) * Sigma * (k/mu)**0.5)**0.382*T**10.371
    elif kappa == 11:
        return 1.0617*10**12 *(Omega(r) * Sigma * (k/mu)**0.5)**0.928/T**3.824
    elif kappa == 12:
        return 10**-.48
    else:
        raise ValueError("Check your kappa input")


def Tfin(r, Sigma, q, f, idx):
    Tmin = 1e-3
    Tmax = 1e7
    try:
        T = brentq(func, Tmin, Tmax, args=(r,Sigma,q,f,idx), maxiter=200)
    except ValueError:
        Tfin(r, Sigma, q, f, idx+1)
    if rightregime(T, Sigma, r, idx):
        return T
    else:
        return Tfin(r, Sigma, q, f, idx+1)

def rightregime(T, Sigma, r, idx):
    
    if idx == 1:
        return T < 144.958* (Omega(r) * Sigma * (k/mu)**0.5)**0.019172
    
    elif idx == 2:
        return 144.958* (Omega(r) * Sigma * (k/mu)**0.5)**0.019172 <= T <= 171.54*(Omega(r) * Sigma *(k/mu)**0.5)**0.019172
    
    elif idx == 3:
        return 171.54*(Omega(r) * Sigma *(k/mu)**0.5)**0.019172 <= T <= 617.376
    
    elif idx == 4:
        return 617.376 < T < 931.773
    
    elif idx == 5:
        return 931.773 <= T <= 1584.42 *(Omega(r) * Sigma * (k/mu)**0.5)**0.027182
    
    elif idx == 6:
        return 1584.42 *(Omega(r) * Sigma * (k/mu)**0.5)**0.027182 < T < 1719.07 * (Omega(r) * Sigma * (k/mu)**0.5)**0.028398
    
    elif idx == 7:
        return 1719.07 * (Omega(r) * Sigma * (k/mu)**0.5)**0.028398 <= T <= 2137.71 * (Omega(r) * Sigma * (k/mu)**0.5)**0.030457
    
    elif idx == 8:
        return 2137.71 * (Omega(r) * Sigma * (k/mu)**0.5)**0.030457 < T < 2656.1 * (Omega(r) * Sigma*(k/mu)**0.5)**0.0083548
    
    elif idx == 9:
        return 2656.1 * (Omega(r) * Sigma*(k/mu)**0.5)**0.0083548 <= T <= 5345.15 * (Omega(r) * Sigma * (k/mu)**0.5)**0.0151134
    
    elif idx == 10:
        return 5345.15 * (Omega(r) * Sigma * (k/mu)**0.5)**0.0151134 < T < 9769.78 * (Omega(r) * Sigma * (k/mu)**0.5)**0.040816
    
    elif idx == 11:
        return 9769.78 * (Omega(r) * Sigma * (k/mu)**0.5)**0.040816 <= T <= 19529.8 *(Omega(r) * Sigma * (k/mu)**0.5)**0.32558
    
    elif idx == 12:
        return T > 19529.8 *(Omega(r) * Sigma * (k/mu)**0.5)**0.32558
    
    else:
        raise ValueError("Opacity index out of range")

def buildTempTable(rGrid, q=1.0, f=0.001, Tmin=202.6769, Tmax=1e7, Sigmin=1.0e-5, Sigmax=1.0e4, Sigres=2000, **kargs):
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
            try:
                temp[i,j] = Tfin(r, Sigma, q, f, 1)
            except ValueError:
                print "No solution found for any opacities"
    # Return values in logspace for interpolation
    return np.log10(rGrid), np.log10(SigmaGrid), np.log10(temp)

