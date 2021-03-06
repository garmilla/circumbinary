#This code was copied from the astroML package
#https://github.com/astroML/astroML/blob/master/astroML/decorators.py 
import pickle
import numpy as np
import matplotlib.pyplot as plt

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
