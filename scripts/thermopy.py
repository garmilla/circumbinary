# -*- coding: utf-8 -*-
"""
Created on Wed Jan  7 14:02:28 2015

@author: dvartany
"""

import numpy as np
import scipy.optimize as optimize
import matplotlib.pyplot as plt
from matplotlib import rc

alpha = 1.0e-2
eta = 3.5
G = 6.674e-8
sigma = 5.6704e-5
k = 1.3806e-16
mu = 2*1.673e-24
AU = 1.49597871e13
a = 0.2*AU # Semimajor axis of the binary system
L = 3.839e33
M = 1.9891e33
OmegaIn = (G*M/a**3)**0.5
f = 0.01
q = 1

kappa0 = 0.1
kappa1 = 2*10**16
kappa2 = 2*10**(-4)

Trange = np.linspace(1,50000,10)

rmin = 12.6021 #log scale, 4e12 cm
rmax = 13
#rmax = 17.4771 #log scale, 3e17 cm
ngrid = 5 #no. cells

Sigmin = 1000 #min density in cgs
Sigmax = 1500 #max density +1 in cgs
Sigres = 1000 #density resolution

temp = [] #create m x n array for temp
s = np.logspace(rmin,rmax,ngrid) #dummy grid to reference for indicing


for r in np.logspace(rmin,rmax,ngrid):
    temp.append([])
    for Sigma in range(Sigmin,Sigmax,Sigres):
        def lam(r):
            return f*q**2*G*M/a*(a/(r-a))**4
        
        def Omega(r):
            return np.sqrt(G*M/r**3)
        
        #option to turn off tidal heating
#        def ftid(r,Sigma):
#            return 0
        def ftid(r,Sigma):
            return 0.5*(OmegaIn - Omega(r))*lam(r)*Sigma
        
        def fv(r,T,Sigma):
            return 1.125 * Omega(r)*alpha*k*T/mu * Sigma
        
        def Tirr(r):
            return (((eta/7.0)*L/(4*np.pi*sigma))**2* k/(G*M*mu))**(1.0/7.0)*r**(-3.0/7.0)
            
        def func(T):
            return sigma*T**4 - 3*(kappa0*T**0.5*Sigma*0.125 + 2/(kappa0*Sigma*T**0.5))*(ftid(r,Sigma) + fv(r,T,Sigma)) - sigma*Tirr(r)**4 
        
        def func1(T):
            return sigma*T**4 - 3*(kappa1*T**(-7)*Sigma*0.125 + 2/(kappa1*Sigma*T**(-7)))*(ftid(r,Sigma) + fv(r,T,Sigma)) - sigma*Tirr(r)**4 
        
        def func2(T):
            return sigma*T**4 - 3*(kappa2*T**2*Sigma*0.125 + 2/(kappa2*Sigma*T**2))*(ftid(r,Sigma) + fv(r,T,Sigma)) - sigma*Tirr(r)**4 
        
        Tcheck = optimize.bisect(func,1,50000)
        
        def Tfin(Tcheck):
            if Tcheck > 202.677:
                return Tcheck
            else:
                return optimize.bisect(func1,1,204)
                if 166.81 < optimize.bisect(func1,1,204) < 202.677:
                    return optimize.bisect(func1,1,204)
                else:
                    return optimize.bisect(func2,1,204)
        
  
        z = np.where(s==r)[0][0]
   
        temp[np.where(s==r)[0][0]].append(Tfin(Tcheck))

#        print func(Trange), optimize.bisect(func,1,10000)
print temp

#fix minsig, minrad

#print Tfin(Tcheck)
plt.loglog(s/AU, temp)
plt.xlabel(r'radius [AU]')
plt.ylabel(r'T [K]')
plt.show()
#
