#This code was copied from the astroML package
#https://github.com/astroML/astroML/blob/master/astroML/decorators.py
import pickle
import numpy as np
from scipy.integrate import trapz
import matplotlib.pyplot as plt

from constants import *
import thermopy as thm

import convection as conv

_colors=['b', 'g', 'r', 'c', 'm', 'y', 'k']

def plotTmap(circ, Sigres=2000):
    SigmaGrid = np.logspace(np.log10(circ.Sigmin), np.log10(circ.Sigmax), Sigres)
    SigmaGrid = np.flipud(SigmaGrid)
    Tmap = np.zeros((Sigres, len(circ.r)))
    for i, Sigma in enumerate(SigmaGrid):
        SigmaArr = np.ones(circ.r.shape)*Sigma
        Tmap[i,:] = circ._bellLinT(SigmaArr)

    fig = plt.figure()

    plt.imshow(np.log10(Tmap),
               extent=(np.log10(circ.r[0]*a*circ.gamma/AU), np.log10(circ.r[-1]*a*circ.gamma/AU),
               np.log10(circ.dimensionalSigma(circ.Sigmin)), np.log10(circ.dimensionalSigma(circ.Sigmax))),
               interpolation='nearest', aspect='auto')

    clb = plt.colorbar()
    clb.set_label("log10(T) K")
    plt.xlabel("log10(r) AU")
    plt.ylabel("log10(Sigma) g/cm^2")

    return fig

def overPlotRuns(circ, fig=None, times=None, nTimes=4, **kargs):
    if fig is None:
        fig = plotTmap(circ, **kargs)

    ax = fig.get_axes()[0]
    xlim = ax.get_xlim(); ylim = ax.get_ylim()

    if times == None:
        times = np.logspace(np.log10(circ.times[0]), np.log10(circ.times[-1]), nTimes)
        print "You didn't specify times, I'll plot the times: {0}".format(circ.dimensionalTime(times))
    else:
        times = circ.dimensionlessTime(np.array(times))

    for i, t in enumerate(times):
        circ.loadTime(t)
        print "I'm plotting snapshot {0} yr".format(circ.dimensionalTime())
        Sigma = np.maximum(circ.dimensionalSigma(), circ.dimensionalSigma(circ.Sigmin))
        Sigma = np.log10(Sigma)
        r = np.log10(circ.r*circ.gamma*a/AU)
        ax.plot(r, Sigma, color='k')

    ax.set_xlim(xlim); ax.set_ylim(ylim)
    return fig

def plotSTF(circ, xlim=None, times=None, nTimes=4, logLog=True, sigMin=0.0001, FMin=1.0e35):
    """
    Plot panel with Sigma, temperature, and FJ
    """
    fig = plt.figure()
    
    if times == None:
        times = np.logspace(np.log10(circ.times[0]), np.log10(circ.times[-1]), nTimes)
        print "You didn't specify times, I'll plot the times: {0}".format(circ.dimensionalTime(times))
    else:
        times = circ.dimensionlessTime(np.array(times))
    
    if xlim==None:
        xlim=(circ.r[0], 1.0e5*circ.r[0])
    
    axSigma = plt.subplot(3, 1, 1)
    axT = plt.subplot(3, 1, 2)
    axFJ = plt.subplot(3, 1, 3)

    axSigma.set_ylabel("Sigma")
    axT.set_ylabel("T")
    axFJ.set_ylabel("FJ")
    axFJ.set_xlabel("r/r0")
    
    axSigma.set_xlim(xlim)
    axT.set_xlim(xlim)
    axFJ.set_xlim(xlim)
    
    for i, t in enumerate(times):
        circ.loadTime(t)
        print "I'm plotting snapshot {0} yr".format(circ.dimensionalTime())
        Sigma = circ.dimensionalSigma()
        FJ = circ.dimensionalFJ()
        if logLog:
            axSigma.loglog(circ.r, np.maximum(sigMin, Sigma), color=_colors[i%7])
            axT.loglog(circ.r, circ.T.value, color=_colors[i%7])
            axFJ.loglog(circ.r, np.maximum(FMin, FJ), color=_colors[i%7])
        else:
            axSigma.semilogx(circ.r, Sigma, color=_colors[i%7])
            axT.semilogx(circ.r, circ.T.value, color=_colors[i%7])
            axFJ.semilogx(circ.r, FJ, color=_colors[i%7])
    return fig

def plotaspect(circ, xlim=None, times=None, nTimes=4, logLog=True, sigMin=0.0001):
    """
    Plot panel with Sigma, temperature, and FJ
    """
    fig = plt.figure()
    
    if times == None:
        times = np.logspace(np.log10(circ.times[0]), np.log10(circ.times[-1]), nTimes)
        print "You didn't specify times, I'll plot the times: {0}".format(circ.dimensionalTime(times))
    else:
        times = circ.dimensionlessTime(np.array(times))
    
    if xlim==None:
        xlim=(circ.r[0], 1.0e5*circ.r[0])
    
    axaspect = plt.subplot(1, 1, 1)

    axaspect.set_ylabel("h/r")
    axaspect.set_xlabel("r/r0")
    
    axaspect.set_xlim(xlim)
    
    for i, t in enumerate(times):
        circ.loadTime(t)
        print "I'm plotting snapshot {0} yr".format(circ.dimensionalTime())
        if logLog:
            axaspect.loglog(circ.r, (k*circ.T.value*circ.r*a*circ.gamma/G/M/mu)**0.5, color=_colors[i%7])
        else:
            axaspect.semilogx(circ.r, (k*circ.T.value*circ.r*a*circ.gamma/G/M/mu)**0.5, color=_colors[i%7])
    return figc
    
def plotSTOp(circ, xlim=None, times=None, nTimes=4, logLog=True, sigMin=0.0001):
    """
        Plot panel with opacity
        """
    fig = plt.figure()
    
    if times == None:
        times = np.logspace(np.log10(circ.times[0]), np.log10(circ.times[-1]), nTimes)
        print "You didn't specify times, I'll plot the times: {0}".format(circ.dimensionalTime(times))
    else:
        times = circ.dimensionlessTime(np.array(times))
    
    if xlim==None:
        xlim=(circ.r[0], 1.0e5*circ.r[0])
    
    axSigma = plt.subplot(3, 2, 1)
    axT = plt.subplot(3, 2, 2)
    axOp = plt.subplot(3, 2, 3)
    axTau = plt.subplot(3, 2, 4)
    axidx = plt.subplot(3 ,2, 5)

    axSigma.set_ylabel("Sigma")
    axSigma.set_xlabel("r/r0")
    axT.set_ylabel("T")
    axT.set_xlabel("r/r0")
    axOp.set_ylabel("$\\kappa$")
    axOp.set_xlabel("r/r0")
    axTau.set_ylabel("$\\tau$")
    axTau.set_xlabel("r/r0")
    axidx.set_ylabel("Index")
    axidx.set_xlabel("r/r0")
    
    axSigma.set_xlim(xlim)
    axT.set_xlim(xlim)
    axOp.set_xlim(xlim)
    axTau.set_xlim(xlim)
    axidx.set_xlim(xlim)
    
    for i, t in enumerate(times):
        circ.loadTime(t)
        print "I'm plotting snapshot {0} yr".format(circ.dimensionalTime())
        Sigma = circ.dimensionalSigma()
        Sigma = np.maximum(sigMin, Sigma)
        r = circ.r*circ.gamma*a # Dimensional radius
        T = circ.T.value
        kappa = np.zeros(T.shape)
        solved = np.zeros(T.shape, dtype=bool)
        index = np.zeros(T.shape)
        for idx in range(1, 13):
            Tmin, Tmax = thm.getBracket(r, Sigma, idx)
            good = np.logical_and(True, T > Tmin)
            good = np.logical_and(good, T < Tmax)
            update = np.logical_and(good, np.logical_not(solved))
            kappa[update] = thm.op(T[update], r[update], Sigma[update], idx)
            index[update] = idx
            solved[update] = True
        
        if logLog:
            axSigma.loglog(circ.r, np.maximum(sigMin, Sigma), color=_colors[i%7])
            axT.loglog(circ.r, circ.T.value, color=_colors[i%7])
            axOp.loglog(circ.r, kappa, color=_colors[i%7])
            axTau.loglog(circ.r, np.maximum(kappa*sigMin, kappa*Sigma), color=_colors[i%7])
            axidx.loglog(circ.r,index, color=_colors[i%7])
        else:
            axSigma.semilogx(circ.r, Sigma, color=_colors[i%7])
            axT.semilogx(circ.r, circ.T.value, color=_colors[i%7])
            axOp.semilogx(circ.r, kappa, color=_colors[i%7])
            axTau.semilogx(circ.r, kappa*Sigma, color=_colors[i%7])
            axidx.loglog(circ.r,index, color=_colors(i%7))

    return fig

def plotTVI(circ, xlim=None, times=None, nTimes=4, logLog=True, sigMin=0.0001):
    """
        Plot panel with various heating terms
        """
    fig = plt.figure()
    
    if times == None:
        times = np.logspace(np.log10(circ.times[0]), np.log10(circ.times[-1]), nTimes)
        print "You didn't specify times, I'll plot the times: {0}".format(circ.dimensionalTime(times))
    else:
        times = circ.dimensionlessTime(np.array(times))
    
    if xlim==None:
        xlim=(circ.r[0], 1.0e5*circ.r[0])
    
    axheat= plt.subplot(1,1,1)


    axheat.set_ylabel("Heating Terms (erg/cm^2/s)")
    axheat.set_xlabel("r/r0")
    
    
    axheat.set_xlim(xlim)
    
    axheat.set_ylim(1.0e-0, 1.0e8)
    
    
    for i, t in enumerate(times):
        circ.loadTime(t)
        print "I'm plotting snapshot {0} yr".format(circ.dimensionalTime())
        Sigma = circ.dimensionalSigma()
        T = circ.T.value
        r = circ.r*circ.gamma*a # Dimensional radius
        
        if logLog:
            axheat.loglog(circ.r,thm.ftid(r,Sigma,circ.q,circ.fudge),color= 'r')
            axheat.loglog(circ.r,thm.fv(r,T,Sigma),color='b')
            axheat.loglog(circ.r,sigma*thm.Tirr(r, circ.q)**4,color='g')
        else:
            axheat.semilogx(circ.r,sigma*thm.Tirr(r, circ.q)**4,color ='r')
            axheat.semilogx(circ.r,thm.fv(r,T,Sigma),color='b')
            axheat.semilogx(circ.r,sigma*thm.Tirr(r, circ.q)**4,color='g')

    return fig

def plotdz(circ, xlim=None, Sigdz = None, times=None, nTimes=4, logLog=True, sigMin=0.0001):
    """
    Plot iceline and deadzone
    """
    fig = plt.figure()
    
    if times == None:
        times = np.logspace(np.log10(circ.times[0]), np.log10(circ.times[-1]), nTimes)
        print "You didn't specify times, I'll plot the times: {0}".format(circ.dimensionalTime(times))
    else:
        times = circ.dimensionlessTime(np.array(times))
    
    if xlim==None:
        xlim=(circ.r[0], 1.0e5*circ.r[0])
    
    if Sigdz == None:
        Sigdz = 20
    else:
        Sigdz = Sigdz
        
    axdz = plt.subplot(2, 1, 1)
    axT = plt.subplot(2, 1, 2)


    axT.set_ylabel("T (K)")
    axT.set_xlabel("r/r0")
    axdz.set_ylabel("Sigma")
    axdz.set_xlabel("r/r0")
    
    axdz.set_xlim(xlim)
    axT.set_xlim(xlim)
    
    for i, t in enumerate(times):
        circ.loadTime(t)
        print "I'm plotting snapshot {0} yr".format(circ.dimensionalTime())
        Sigma = circ.dimensionalSigma()
        Sigma = np.maximum(sigMin, Sigma)
        r = circ.r*circ.gamma*a # Dimensional radius
        T = circ.T.value
        kappa = np.zeros(T.shape)
        idxtab = np.zeros(T.shape)
        solved = np.zeros(T.shape, dtype=bool)
        for idx in range(1, 13):
            Tmin, Tmax = thm.getBracket(r, Sigma, idx)
            good = np.logical_and(True, T > Tmin)
            good = np.logical_and(good, T < Tmax)
            update = np.logical_and(good, np.logical_not(solved))
            kappa[update] = thm.op(T[update], r[update], Sigma[update], idx)
            idxtab[update] = idx
            solved[update] = True

        deadzone = np.where((Sigma > Sigdz) & (T < 800))[0]
        
        if logLog:
            axdz.loglog(circ.r, Sigma, color=_colors[i%7])
            axdz.axvline(circ.r[deadzone[0]], color=_colors[i%7])
            axdz.axvline(circ.r[deadzone[-1]], color=_colors[i%7])
            axT.loglog(circ.r, T, color=_colors[i%7])
            axT.axvline(circ.r[deadzone[0]], color=_colors[i%7])
            axT.axvline(circ.r[deadzone[-1]], color=_colors[i%7])
        
        else:
            axdz.semilogx(circ.r, Sigma, color=_colors[i%7])
            axdz.axvline(circ.r[deadzone[0]], color=_colors[i%7])
            axdz.axvline(circ.r[deadzone[-1]], color=_colors[i%7])
            axT.loglog(circ.r, T, color=_colors[i%7])
            axT.axvline(circ.r[deadzone[0]], color=_colors[i%7])
            axT.axvline(circ.r[deadzone[-1]], color=_colors[i%7])
    
    return fig


def ploticeline(circ, xlim=None, logLog=True, sigMin=0.0001):
    """
    Plot iceline
    """
    fig = plt.figure()
    
    times = np.logspace(np.log10(circ.times[0]), np.log10(circ.times[-1]), len(circ.times))

    if xlim==None:
        xlim=(circ.times[0], circ.times[-1])

    axice = plt.subplot(3, 1, 1)
    axTice = plt.subplot(3, 1, 2)
    axSigmaice = plt.subplot(3, 1, 3)
    
    axice.set_ylabel("Iceline (AU)")
    axice.set_xlabel("Time (MY)")
    axTice.set_ylabel("T (K)")
    axTice.set_xlabel("Time (MY)")
    axSigmaice.set_ylabel("Sigma (g/cm^2)")
    axSigmaice.set_xlabel("Time (MY)")

    axice.set_xlim(xlim)
    axTice.set_xlim(xlim)
    axSigmaice.set_xlim(xlim)
    
    iceline = np.zeros(times.shape)
    Tice = np.zeros(times.shape)
    Sigmaice = np.zeros(times.shape)
    
    for i, t in enumerate(times):
        circ.loadTime(t)
        Sigma = circ.dimensionalSigma()
        Sigma = np.maximum(sigMin, Sigma)
        r = circ.r*circ.gamma*a # Dimensional radius
        T = circ.T.value
        kappa = np.zeros(T.shape)
        idxtab = np.zeros(T.shape)
        solved = np.zeros(T.shape, dtype=bool)
        for idx in range(1, 13):
            Tmin, Tmax = thm.getBracket(r, Sigma, idx)
            good = np.logical_and(True, T > Tmin)
            good = np.logical_and(good, T < Tmax)
            update = np.logical_and(good, np.logical_not(solved))
            kappa[update] = thm.op(T[update], r[update], Sigma[update], idx)
            idxtab[update] = idx
            solved[update] = True
        iceline[i] = np.where((idxtab < 3) & (idxtab > 1))[0][-1]

    rad = np.zeros(times.shape)

    for i, ind in enumerate(iceline):
        rad[i] = circ.r[ind]
        Tice[i] = circ.T.value[ind]
        Sigmaice[i] = circ.dimensionalSigma()[ind]

    if logLog:
        axice.loglog(circ.dimensionalTime(circ.times)/1.0e6, rad*20)
        axTice.loglog(circ.dimensionalTime(circ.times)/1.0e6, Tice)
        axSigmaice.loglog(circ.dimensionalTime(circ.times)/1.0e6, Sigmaice)
    else:
        axice.semilogx(circ.dimensionalTime(circ.times)/1.0e6, rad*20)
        axTice.semilogx(circ.dimensionalTime(circ.times)/1.0e6, Tice)
        axSigmaice.semilogx(circ.dimensionalTime(circ.times)/1.0e6, Sigmaice)

    return fig

def plottrunc(circ, xlim=None, logLog=True, sigMin=0.0001):
    """
    Plot truncation radius
    """
    fig = plt.figure()
    
    times = np.logspace(np.log10(circ.times[0]), np.log10(circ.times[-1]), len(circ.times))
    
    if xlim==None:
        xlim=(circ.times[0], circ.times[-1])
    
    axtrunc= plt.subplot(2, 1, 1)
    axplat = plt.subplot(2, 1, 2)

    
    axtrunc.set_ylabel("Truncation Radius (AU)")
    axtrunc.set_xlabel("Time (MY)")
    axplat.set_ylabel("Plateau FJ (AU)")
    axplat.set_xlabel("Time (MY)")


    axtrunc.set_xlim(xlim)
    axplat.set_xlim(xlim)

    trunc = np.zeros(times.shape)
    plat = np.zeros(times.shape)
    
    for i, t in enumerate(times):
        circ.loadTime(t)
        FJ = circ.dimensionalFJ()
        r = circ.r*circ.gamma*a # Dimensional radius
        trunc[i] = np.where(FJ > 0.1*np.max(FJ))[0][-1]
        plat[i] = np.max(FJ)

    rad = np.zeros(times.shape)
    
    for i, ind in enumerate(trunc):
        rad[i] = circ.r[ind]
    
    if logLog:
        axtrunc.loglog(circ.dimensionalTime(circ.times)/1.0e6, rad*20)
        axplat.loglog(circ.dimensionalTime(circ.times)/1.0e6, plat)
    else:
        axtrunc.semilogx(circ.dimensionalTime(circ.times)/1.0e6, rad*20)
        axplat.semilogx(circ.dimensionalTime(circ.times)/1.0e6, plat)

    
    return fig

def geticeline(circ, sigMin=0.0001):
    """
    Returns two arrays, the first contains the times in years and the
    second the value of the radius at the iceline for that snapshot.
    """
    times = circ.dimensionalTime(circ.times)
    iceline = np.zeros(times.shape)
    for i, t in enumerate(circ.times):
        circ.loadTime(t)
        Sigma = circ.dimensionalSigma()
        Sigma = np.maximum(sigMin, Sigma)
        r = circ.r*circ.gamma*a # Dimensional radius
        T = circ.T.value
        kappa = np.zeros(T.shape)
        idxtab = np.zeros(T.shape)
        solved = np.zeros(T.shape, dtype=bool)
        for idx in range(1, 13):
            Tmin, Tmax = thm.getBracket(r, Sigma, idx)
            good = np.logical_and(True, T > Tmin)
            good = np.logical_and(good, T < Tmax)
            update = np.logical_and(good, np.logical_not(solved))
            kappa[update] = thm.op(T[update], r[update], Sigma[update], idx)
            idxtab[update] = idx
            solved[update] = True
        iceline[i] = circ.r[np.where((idxtab < 3) & (idxtab > 1))[0][-1]]*a*circ.gamma/AU
    return times, iceline

def getrinfl(circ):
    """
    Returns two arrays, the first contains the times in years and the
    second the value of rinfl for that snapshot.
    """
    times = circ.dimensionalTime(circ.times)
    rinfl = np.zeros(times.shape)
    for i, t in enumerate(circ.times):
        circ.loadTime(t)
        l = circ.Sigma.value*circ.r**2
        idx = l.argmax()
        rinfl[i] = circ.r[idx]*a*circ.gamma/AU
    return times, rinfl

def getFJt(circ):
    """
    Returns two arrays, the first contains the times in years and the
    second the maximum value of FJ for that snapshot.
    """
    times = circ.dimensionalTime(circ.times)
    FJ = np.zeros(times.shape)
    for i, t in enumerate(circ.times):
        circ.loadTime(t)
        FJ[i] = circ.dimensionalFJ().max()
    return times, FJ

def getKappa(circ):
    Sigma = circ.dimensionalSigma()
    T = circ.T.value
    r = circ.r*circ.gamma*a # Dimensional radius
    kappa = np.zeros(T.shape)
    solved = np.zeros(T.shape, dtype=bool)
    index = np.zeros(T.shape)
    for idx in range(1, 13):
        Tmin, Tmax = thm.getBracket(r, Sigma, idx)
        good = np.logical_and(True, T > Tmin)
        good = np.logical_and(good, T < Tmax)
        update = np.logical_and(good, np.logical_not(solved))
        kappa[update] = thm.op(T[update], r[update], Sigma[update], idx)
        index[update] = idx
        solved[update] = True
    return kappa

def getTeff(circ, tau=None, tauMin=0.0001):
    """
    Return an array with the effective temperature as defined
    in the paper
    """
    if tau is None:
        kappa = getKappa(circ)
        tau = np.maximum(tauMin, 0.5*circ.dimensionalSigma()*kappa)
    Sigma = circ.dimensionalSigma()
    T = circ.T.value
    r = circ.r*a*circ.gamma
    Fnu = thm.ftid(r, Sigma, circ.q, circ.fudge) + thm.fv(r, T, Sigma)
    Firr = sigma*thm.Tirr(r, circ.q)**4
    Teff = np.power(((1.0+1.0/tau)*Fnu + Firr)/sigma, 0.25)
    return Teff

def getBnu(nu, T):
    """
    Get the flux at frequency nu, for a blackbody at temperature T.
    T can be an array, in which case an array of fluxes is returned.
    """
    Bnu = np.zeros(T.shape)
    for i in range(len(T)):
        if T[i] > 0.0:
           if h*nu/k/T[i] < 1.0e2:
               Bnu[i] = 2*h*nu**3/c**2\
                        /(np.exp(h*nu/k/T[i]) - 1.0)
    return Bnu

def extrap(circ, nextrap=40, rmin = 6):
    Rs = 6.955e10
    r = np.exp(np.linspace(np.log(rmin*Rs/(a*circ.gamma)),np.log(circ.r[0]**2/circ.r[1]), nextrap))*a*circ.gamma
    r2 = np.append(r, circ.r[0]*a*circ.gamma)
    r2rev = r2[::-1]
    FJ = np.array([circ.dimensionalFJ()[0]])
    for i in range(len(r)):
        x = FJ[i]*(r2rev[i+1]/r2rev[i])**0.5
        FJ = np.append(FJ, x)
    Sigextrap = np.array([])

    table = thm.buildTempTable(r2rev, q=0)
    temp = np.power(10, table[2])
    Sig = np.power(10, table[1])

    Sig0 = np.array([circ.dimensionalSigma()[0]])
    for j in range(len(r)):
        idx = np.where(Sig < Sig0[-1])[-1][-1]
        def FJSig(Sig):
            return Sig - FJ[j]*mu/3/np.pi/alpha/k/r2rev[j]**2/temp[j,idx]
        Signew = broyden1(FJSig, Sig0[-1])
        Sig0 = np.append(Sig0,Signew)
        Sigextrap = np.append(Sigextrap, Signew)
    
    return Sigextrap[::-1]

def getSED(circ, extrap=False, CG = False, power=1.0/0.95, RStar = 1, MStar = 1, TStar = 5780, LStar = 1, \
            Rmin = 1, Rmax = 270, SH = False, Flared = False, nextrap = 40, Q = 1,\
            Teff=None, Tsh=None, Tirr=None, tau=None, nLambda=1000, tauMin=0.0001):
    """
    Returns four arrays:
    lamb: Wavelength in microns
    Q: the emission efficiency
    fnuD: The contribution of the disk to the SED 
    fnuSh: The contribution of superheated grains to the SED
    fnuS: The contribution of the binary/star to the SED
    fnuT: The total SED
    """
    Ts = np.array([TStar]) # Temperature of the star
    Rs = RStar * 6.955e10 #Radius of the star
    nu = np.linspace(10.0, 15.0, num = nLambda)
    nu = np.power(10.0, nu)
    fnuD = np.zeros(nu.shape)
    fnuS = np.zeros(nu.shape)
    fnuSh = np.zeros(nu.shape)
    fnuT = np.zeros(nu.shape)
    
    if CG:
        rout = np.where(circ.r*a*circ.gamma/AU < Rmax)[0][-1]
        r = np.append(np.exp(np.linspace(np.log(Rmin*Rs/(a*circ.gamma)),np.log(circ.r[0]**2/circ.r[1]),nextrap)),\
            circ.r[:-(circ.ncell - rout+1)])*a*circ.gamma     
        if Flared:
            alpha = 0.005*AU/r + 0.05*(r/AU)**(2.0/7)
        else:
            alpha = 0.005*AU/r 
        if Teff is None:
            Teff = (alpha/2)**0.25*(Rs/r)**0.5 * Ts
        if SH:
            Tsh = np.power(L*LStar/16/np.pi/sigma/Q/(r)**2, 0.25)
        else: 
            Tsh = 0 
        if Tirr is None:
            Tirr = (((eta/7.0)*0.5*L*LStar/(4*np.pi*sigma))**2* k/(G*MStar*M*mu))**(1.0/7.0)*r**(-3.0/7.0)
        Firr = sigma*Tirr**4
        for i in range(len(nu)):
            x = r
            y = getBnu(nu[i], Teff)
            z = alpha*getBnu(nu[i], Tsh)
            y *= 8*np.pi*np.pi*x
            z *= 8*np.pi*np.pi*x
            fnuD[i] = nu[i]*trapz(y, x)
            fnuSh[i] = nu[i]*trapz(z,x)
            fnuS[i] = nu[i]*np.pi*getBnu(nu[i], Ts)*4*np.pi*Rs**2
            
    if circ.q == 0:
        r = circ.r*a*circ.gamma
        kappa = getKappa(circ)
        if tau is None:
            tau = np.maximum(tauMin, 0.5*circ.dimensionalSigma()*kappa)
        if Teff is None:
            Teff = getTeff(circ, tau=tau)
        if Tsh is None:
            Tsh = np.power(L/16/np.pi/sigma/0.1/(r)**2, 0.25)
        Firr = sigma*thm.Tirr(r, circ.q)**4
        for i in range(len(nu)):
            x = r
            y = tau/(1.0 + tau)*getBnu(nu[i], Teff)
            z = (2.0+tau)/(1.0+tau)*Firr/sigma/np.maximum(1.0e1, Tsh)**4*getBnu(nu[i], Tsh)
            y *= 2*np.pi*np.pi*x
            z *= 2*np.pi*np.pi*x
            fnuD[i] = nu[i]*trapz(y, x)
            fnuSh[i] = nu[i]*trapz(z,x)
            fnuS[i] = nu[i]*np.pi*getBnu(nu[i], Ts)*np.pi*Rs**2
        
    
    elif circ.q == 1.0:
    # We don't include the gap for circumbinary disks
        r = circ.r*a*circ.gamma
        kappa = getKappa(circ)
        if tau is None:
            tau = np.maximum(tauMin, 0.5*circ.dimensionalSigma()*kappa)
        if Teff is None:
            Teff = getTeff(circ, tau=tau)
        if Tsh is None:
            Tsh = np.power(L/16/np.pi/sigma/0.1/(r)**2, 0.25)
        Firr = sigma*thm.Tirr(r, circ.q)**4
        for i in range(len(nu)):
            x = r
            y = tau/(1.0 + tau)*getBnu(nu[i], Teff)
            z = (2.0+tau)/(1.0+tau)*Firr/sigma/np.maximum(1.0e1, Tsh)**4*getBnu(nu[i], Tsh)
            y *= 2*np.pi*np.pi*x
            z *= 2*np.pi*np.pi*x
            fnuD[i] = nu[i]*trapz(y, x)
            fnuSh[i] = nu[i]*trapz(z,x)
            fnuS[i] = nu[i]*np.pi*getBnu(nu[i], Ts)*np.pi*Rs**2
        
        Teff[np.where(circ.r < circ.rF[0]*2)] = 0.0
        Tsh[np.where(circ.r < circ.rF[0]*2)] = 0.0
        Firr[np.where(circ.r < circ.rF[0]*2)] = 0.0
        
    elif circ.q != 0.0:
        raise ValueError("I only compute SEDs for q=1 and q=0, you specified q={0}".format(circ.q))
        
    # Integrate the set of blackbodies at each frequency using the trapezoidal rule
    
    fnuT = fnuD + fnuS + fnuSh
    lamb = c/nu*1.0e4 # In microns
    return lamb, fnuD, fnuSh, fnuS, fnuT

_cBinaries = ['/u/dvartany/circumaster/circumbinary/scripts/outputzz012',
              '/u/dvartany/circumaster/circumbinary/scripts/outputzz052',
              '/u/dvartany/circumaster/circumbinary/scripts/outputzz12']

_cStellars = ['/u/dvartany/circumaster/circumbinary/scripts/outputzzcs012',
              '/u/dvartany/circumaster/circumbinary/scripts/outputzzcs052',
              '/u/dvartany/circumaster/circumbinary/scripts/outputzzcs12']

_times = [5.0e3, 5.0e4, 5.0e5, 5.0e6]

def genSMInputs(cBinaries=None, cStellars=None, times=None, Sigmin=0.01, Tmin=1.0,
                tauMin=0.01, FJMin=1.0e35, nLambda=1000):
    """
    Generate the files that supermongo takes as inputs to generate the paper's plots
    Keywords:
    cBinaries: A list of the directories with the data for the circumbinary disks from
               lowest to highest mass. If None the above predefined lists are used.
    cStellars: A list of directories with the data  for the circumstellar disks from
               lowest to highest mass. If None the above predefined lists are used.
    times: A list of times to use for the snapshots. If None the above predefined
           lists are used.
    """

    if cBinaries is None:
        cBinaries = _cBinaries
    if cStellars is None:
        cStellars = _cStellars
    if times is None:
        times = _times

    # Generate the files to plot the Sigma, T, tau and FJ snapshots for
    # the circumbinary disks.
    for disk in cBinaries:
        circ = conv.loadResults(disk)
        for i, time in enumerate(times):
            outputArr = np.zeros((circ.ncell, 6))
            t = circ.dimensionlessTime(time)
            circ.loadTime(t)
            outputArr[:,0] = circ.r
            outputArr[:,1] = circ.r*a*circ.gamma/AU
            outputArr[:,2] = np.maximum(Sigmin, circ.dimensionalSigma())
            outputArr[:,3] = np.maximum(Tmin, circ.T.value)
            kappa = getKappa(circ)
            outputArr[:,4] = np.maximum(tauMin, circ.dimensionalSigma()*kappa)
            outputArr[:,5] = np.maximum(FJMin, circ.dimensionalFJ())
            np.savetxt('m{0}_{1}.dat'.format(circ.mDisk, i+1), outputArr)

    # Generate the files to plot the Sigma, T, tau and FJ snapshots for
    # the cicumstellar disks.
    for disk in cStellars:
        circ = conv.loadResults(disk)
        for i, time in enumerate(times):
            outputArr = np.zeros((circ.ncell, 6))
            t = circ.dimensionlessTime(time)
            circ.loadTime(t)
            outputArr[:,0] = circ.r
            outputArr[:,1] = circ.r*a*circ.gamma/AU
            outputArr[:,2] = np.maximum(Sigmin, circ.dimensionalSigma())
            outputArr[:,3] = np.maximum(Tmin, circ.T.value)
            kappa = getKappa(circ)
            outputArr[:,4] = np.maximum(tauMin, circ.dimensionalSigma()*kappa)
            outputArr[:,5] = np.maximum(FJMin, circ.dimensionalFJ())
            np.savetxt('m{0}_circumstellar_{1}.dat'.format(circ.mDisk, i+1), outputArr)

    # Generate the files to plot the value of FJ at the plateau as a function of
    # time. We only do this for the circumbinary disks.
    for disk in cBinaries:
        circ = conv.loadResults(disk)
        outputArr = np.zeros((len(circ.times), 3))
        Times, FJ = getFJt(circ)
        outputArr[:,0] = Times
        outputArr[:,1] = FJ
        # We also need to store the analytic fit
        outputArr[:,2] = 5.2e37*circ.mDisk/0.01*np.power(Times/3.0e6, -6.0/13)
        np.savetxt('m{0}_ftime.dat'.format(circ.mDisk), outputArr)

    # Generate the files to plot rinfl as a function of time, we also only do this
    # for circumbinary disks.
    for disk in cBinaries:
        circ = conv.loadResults(disk)
        outputArr = np.zeros((len(circ.times), 3))
        Times, rinfl = getrinfl(circ)
        outputArr[:,0] = Times
        outputArr[:,1] = rinfl
        # We also need to store the analytic fit
        outputArr[:,2] = 380.0*np.power(Times/3.0e6, 14.0/13)
        np.savetxt('m{0}_rinfl.dat'.format(circ.mDisk), outputArr)
       
    # Generate the files to plot the iceline as a function of time, we do this for both
    # circumbinary and circumstellar disks.
    for disk in cBinaries:
        circ = conv.loadResults(disk)
        outputArr = np.zeros((len(circ.times), 2))
        Times, iceline = geticeline(circ)
        outputArr[:,0] = Times
        outputArr[:,1] = iceline
        np.savetxt('m{0}_iceline.dat'.format(circ.mDisk), outputArr)
    for disk in cStellars:
        circ = conv.loadResults(disk)
        outputArr = np.zeros((len(circ.times), 2))
        Times, iceline = geticeline(circ)
        outputArr[:,0] = Times
        outputArr[:,1] = iceline
        np.savetxt('m{0}_iceline_circumstellar.dat'.format(circ.mDisk), outputArr)

    # Generate the files to plot the SEDs, we only do this for the disks with mass
    # 0.05 M_c
    for disk in [cBinaries[1], cStellars[1]]:
        circ = conv.loadResults(disk)
        for i, time in enumerate(times):
            outputArr = np.zeros((nLambda, 4))
            t = circ.dimensionlessTime(time)
            circ.loadTime(t)
            nu, fnuD, fnuS, fnuT = getSED(circ, nLambda=nLambda)
            outputArr[:,0] = nu
            outputArr[:,1] = fnuD
            outputArr[:,2] = fnuS
            outputArr[:,3] = fnuT
            if circ.q == 1.0:
                np.savetxt('m{0}_spectrum_{1}.dat'.format(circ.mDisk, i+1), outputArr)
            elif circ.q == 0.0:
                np.savetxt('m{0}_spectrum_circumstellar_{1}.dat'.format(circ.mDisk, i+1), outputArr)

if __name__ == '__main__':
    genSMInputs()
