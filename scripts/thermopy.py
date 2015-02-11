# -*- coding: utf-8 -*-
"""
    Created on Wed Jan  7 14:02:28 2015
    
    @author: dvartany
    """

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq

# Constants in cgs

alpha = 1.0e-2
eta = 3.5
c = 29979245800.0
beta = 0.5
G = 6.674e-8
kappa0 = 0.1
kappa_es = 0.4
sigma = 5.6704e-5
k = 1.3806e-16
mu = 2*1.673e-24
AU = 1.49597871e13
a = 0.2*AU # Semimajor axis of the binary system
L = 3.839e33
M = 1.9891e33
OmegaIn = (G*M/a**3)**0.5
cs = 1.0e5

Tmin = 1e-3
Tmax = 5e6

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
    elif kappa ==   2:
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
        return "Check your kappa input"

def Tfin(r, Sigma, q, f, idx):
    try:
        T = brentq(func, Tmin, Tmax, args=(r,Sigma,q,f,idx), maxiter=200)
    except ValueError:
        Tfin(r, Sigma, q, f, idx+1)
    if rightRegime(T, idx):
        return T
    else:
        Tfin(r, Sigma, q, f, idx+1)

def rightRegima(T, idx):
    # If it's in the right regime return True, False otherwise
    try:
        T1 = brentq(func, Tmin, Tmax, args=(r,Sigma,q,f,1), maxiter=200)
    except ValueError:
        return brentq(func, Tmin, Tmax, args=(r,Sigma,q,f,2), maxiter=200)
    if T1 < 144.958* (Omega(r) * Sigma * (k/mu)**0.5)**0.019172:
        return T1
    else:
        try:
            T2 =  brentq(func,1e-3,618,args=(r,Sigma,q,f,2), maxiter=200)
        except ValueError:
            return brentq(func, Tmin, Tmax, args=(r,Sigma,q,f,3), maxiter=200)
        if 144.958* (Omega(r) * Sigma * (k/mu)**0.5)**0.019172 <= T2 <= 171.54*(Omega(r) * Sigma *(k/mu)**0.5)**0.019172:
            return T2
        elif 144.958* (Omega(r) * Sigma * (k/mu)**0.5)**0.019172  > T2 or  T2 > 171.54*(Omega(r) * Sigma *(k/mu)**0.5)**0.019172:
            try:
                T3 = brentq(func, Tmin, Tmax, args=(r,Sigma,q,f,3), maxiter=200)
            except ValueError:
                return brentq(func, Tmin, Tmax, args=(r,Sigma,q,f,4), maxiter=200)
            if 171.54*(Omega(r) * Sigma *(k/mu)**0.5)**0.019172 <= T3 <= 617.376:
                return T3
        elif 171.54*(Omega(r) * Sigma *(k/mu)**0.5)**0.019172  > T3 or T3 > 617.376:
            try:
                T4 = brentq(func, Tmin, Tmax, args=(r,Sigma,q,f,4), maxiter=200)
            except ValueError:
                return brentq(func, Tmin, Tmax, args=(r,Sigma,q,f,5), maxiter=200)
            if 617.376 < T4 < 931.773:
                return T4
        elif 617.376 > T4 or T4 > 931.773:
            try:
                T5 = brentq(func, Tmin, Tmax, args=(r,Sigma,q,f,5), maxiter=200)
            except ValueError:
                return brentq(func, Tmin, Tmax, args=(r,Sigma,q,f,6), maxiter=200)
            if 931.773 <= T5 <= 1584.42 *(Omega(r) * Sigma * (k/mu)**0.5)**0.027182:
                return T5
        elif 931.773 > T5 or T5 > 1584.42 *(Omega(r) * Sigma * (k/mu)**0.5)**0.027182:
            try:
                T6 = brentq(func, Tmin, Tmax, args=(r,Sigma,q,f,6), maxiter=200)
            except ValueError:
                return brentq(func, Tmin, Tmax, args=(r,Sigma,q,f,7), maxiter=200)
            if 1584.42 *(Omega(r) * Sigma * (k/mu)**0.5)**0.027182 < T6 < 1719.07 * (Omega(r) * Sigma * (k/mu)**0.5)**0.028398:
                return T6
        elif 1584.42 *(Omega(r) * Sigma * (k/mu)**0.5)**0.027182 > T6 or T6 > 1719.07 * (Omega(r) * Sigma * (k/mu)**0.5)**0.028398:
            try:
                T7 = brentq(func, Tmin, Tmax, args=(r,Sigma,q,f,7), maxiter=200)
            except ValueError:
                return brentq(func, Tmin, Tmax, args=(r,Sigma,q,f,8), maxiter=200)
            if 1719.07 * (Omega(r) * Sigma * (k/mu)**0.5)**0.028398 <= T7 <= 2137.71 * (Omega(r) * Sigma * (k/mu)**0.5)**0.030457:
                return T7
        elif 1719.07 * (Omega(r) * Sigma * (k/mu)**0.5)**0.028398 > T7 or T7 > 2137.71 * (Omega(r) * Sigma * (k/mu)**0.5)**0.030457:
            try:
                T8 = brentq(func, Tmin, Tmax, args=(r,Sigma,q,f,8), maxiter=200)
            except ValueError:
                return brentq(func, Tmin, Tmax, args=(r,Sigma,q,f,9), maxiter=200)
            if 2137.71 * (Omega(r) * Sigma * (k/mu)**0.5)**0.030457 < T8 < 2656.1 * (Omega(r) * Sigma*(k/mu)**0.5)**0.0083548:
                return T8
        elif 2137.71 * (Omega(r) * Sigma * (k/mu)**0.5)**0.030457 > T8 or T8 > 2656.1 * (Omega(r) * Sigma*(k/mu)**0.5)**0.0083548:
            try:
                T9 = brentq(func, Tmin, Tmax, args=(r,Sigma,q,f,9), maxiter=200)
            except ValueError:
                return brentq(func, Tmin, Tmax, args=(r,Sigma,q,f,10), maxiter=200)
            if 2656.1 * (Omega(r) * Sigma*(k/mu)**0.5)**0.0083548 <= T9 <= 5345.15 * (Omega(r) * Sigma * (k/mu)**0.5)**0.0151134:
                return T9
        elif 2656.1 * (Omega(r) * Sigma*(k/mu)**0.5)**0.0083548 >T9 or T9 > 5345.15 * (Omega(r) * Sigma * (k/mu)**0.5)**0.0151134:
            try:
                T10 = brentq(func, Tmin, Tmax, args=(r,Sigma,q,f,10), maxiter=200)
            except ValueError:
                return brentq(func, Tmin, Tmax, args=(r,Sigma,q,f,11), maxiter=200)
            if 5345.15 * (Omega(r) * Sigma * (k/mu)**0.5)**0.0151134 < T10 < 9769.78 * (Omega(r) * Sigma * (k/mu)**0.5)**0.040816:
                return T10
        elif 5345.15 * (Omega(r) * Sigma * (k/mu)**0.5)**0.0151134 > T10 or T10 > 9769.78 * (Omega(r) * Sigma * (k/mu)**0.5)**0.040816:
            try:
                T11 = brentq(func, Tmin, Tmax, args=(r,Sigma,q,f,11), maxiter=200)
            except ValueError:
                return brentq(func, Tmin, Tmax, args=(r,Sigma,q,f,12), maxiter=200)
            if 9769.78 * (Omega(r) * Sigma * (k/mu)**0.5)**0.040816 <= T11 <= 19529.8 *(Omega(r) * Sigma * (k/mu)**0.5)**0.32558:
                return T11
        else:
            return brentq(func, Tmin, Tmax, args=(r,Sigma,q,f,12), maxiter=200)



def buildTempTable(rGrid, q=1.0, f=0.001, Tmin=202.6769, Tmax=5.0e6, Sigmin=1.0e-5, Sigmax=1.0e4, Sigres=2000, **kargs):
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
                T1 = brentq(func, Tmin, Tmax, args=(r,Sigma,q,f,1), maxiter=200)
                temp[i,j] = Tfin(T1, r, Sigma, q, f)
            except ValueError:
                try:
                    temp[i,j] = Tfin(1.0, r, Sigma, q, f)
                except ValueError:
                    print "args2=", r, Sigma, q, f
                    Ts = np.linspace(1.0e-3, 166.813, num=1000)
                    ys = func2(Ts, r, Sigma, q, f)
                    plt.semilogx(Ts, np.log10(ys+1.0))
                    plt.semilogx(Ts, np.log10(-ys+1.0))
                    #plt.semilogx(Ts, ys)
                    plt.show()
                    raise
    # Return values in logspace for interpolation
    return np.log10(rGrid), np.log10(SigmaGrid), np.log10(temp)