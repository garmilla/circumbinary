#This code was copied from the astroML package
#https://github.com/astroML/astroML/blob/master/astroML/decorators.py
import pickle
import numpy as np
import matplotlib.pyplot as plt

from constants import *
import thermopy as thm

def pickle_results(filename=None, verbose=True):
    """Generator for decorator which allows pickling the results of a funcion
        
        Pickle is python's built-in object serialization.  This decorator, when
        used on a function, saves the results of the computation in the function
        to a pickle file.  If the function is called a second time with the
        same inputs, then the computation will not be repeated and the previous
        results will be used.
        
        This functionality is useful for computations which take a long time,
        but will need to be repeated (such as the first step of a data analysis).
        
        Parameters
        ----------
        filename : string (optional)
        pickle file to which results will be saved.
        If not specified, then the file is '<funcname>_output.pkl'
        where '<funcname>' is replaced by the name of the decorated function.
        verbose : boolean (optional)
        if True, then print a message to standard out specifying when the
        pickle file is written or read.
        
        Examples
        --------
        >>> @pickle_results('tmp.pkl', verbose=True)
        ... def f(x):
        ...     return x * x
        >>> f(4)
        @pickle_results: computing results and saving to 'tmp.pkl'
        16
        >>> f(4)
        @pickle_results: using precomputed results from 'tmp.pkl'
        16
        >>> f(6)
        @pickle_results: computing results and saving to 'tmp.pkl'
        36
        >>> import os; os.remove('tmp.pkl')
        """
    def pickle_func(f, filename=filename, verbose=verbose):
        if filename is None:
            filename = '%s_output.pkl' % f.__name__
        
        def new_f(*args, **kwargs):
            try:
                D = pickle.load(open(filename, 'rb'))
                cache_exists = True
            except:
                D = {}
                cache_exists = False
        
            # simple comparison doesn't work in the case of numpy arrays
            Dargs = D.get('args')
            Dkwargs = D.get('kwargs')
            
            try:
                args_match = (args == Dargs)
            except:
                args_match = np.all([np.all(a1 == a2)
                             for (a1, a2) in zip(Dargs, args)])
            
            try:
                kwargs_match = (kwargs == Dkwargs)
            except:
                kwargs_match = ((sorted(Dkwargs.keys())
                                 == sorted(kwargs.keys()))
                                and (np.all([np.all(Dkwargs[key]
                                                    == kwargs[key])
                                             for key in kwargs])))
                    
            if (type(D) == dict and D.get('funcname') == f.__name__
                    and args_match and kwargs_match):
                if verbose:
                    print("@pickle_results: using precomputed "
                            "results from '%s'" % filename)
                retval = D['retval']
                                                 
            else:
                if verbose:
                    print("@pickle_results: computing results "
                            "and saving to '%s'" % filename)
                    if cache_exists:
                        print("  warning: cache file '%s' exists" % filename)
                        print("    - args match:   %s" % args_match)
                        print("    - kwargs match: %s" % kwargs_match)
                retval = f(*args, **kwargs)
                                                                                 
                funcdict = dict(funcname=f.__name__, retval=retval,
                                args=args, kwargs=kwargs)
                with open(filename, 'wb') as outfile:
                    pickle.dump(funcdict, outfile)
                                                                                             
            return retval
        return new_f
    return pickle_func

_colors=['b', 'g', 'r', 'c', 'm', 'y', 'k']

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
            axheat.loglog(circ.r,sigma*thm.Tirr(r)**4,color='g')
        else:
            axheat.semilogx(circ.r,sigma*thm.Tirr(r)**4,color ='r')
            axheat.semilogx(circ.r,thm.fv(r,T,Sigma),color='b')
            axheat.semilogx(circ.r,sigma*thm.Tirr(r)**4,color='g')

    return fig

def plotdz(circ, xlim=None, times=None, nTimes=4, logLog=True, sigMin=0.0001):
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

        deadzone = np.where((Sigma > 20) & (T < 800))
        
        if logLog:
            axdz.loglog(circ.r, Sigma, color=_colors[i%7])
            axdz.axvline(circ.r[deadzone][0], circ.r[deadzone][-1])
            axT.loglog(circ.r, T, color=_colors[i%7])
            axT.axvline(circ.r[deadzone][0], circ.r[deadzone][-1])
        
        else:
            axdz.semilogx(circ.r, Sigma, color=_colors[i%7])
            axdz.axvline(circ.r[deadzone][0], circ.r[deadzone][-1])
            axT.semilogx(circ.r, T, color=_colors[i%7])
            axT.axvline(circ.r[deadzone][0], circ.r[deadzone][-1])
                
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
