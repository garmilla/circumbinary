import os
import re
import argparse
import pickle
import numpy as np

from fipy import CylindricalGrid1D, CellVariable, FaceVariable, TransientTerm, ExponentialConvectionTerm

#from thermopy import buildTempTable
import thermopy
from constants import *
from utils import pickle_results

class Circumbinary(object):
    def __init__(self, rmax=1.0e2, ncell=200, nstep=100, dt=1.0e-6, delta=1.0e-100,
                 nsweep=10, titer=10, fudge=1.0e-2, q=1.0, gamma=100, mDisk=0.1, odir='output',
                 bellLin=True, **kargs):
        self.rmax = rmax
        self.ncell = ncell
        self.nstep = nstep
        self.dt = dt
        self.delta = delta
        self.nsweep = nsweep
        self.titer = titer
        self.mDisk = mDisk
        Omega0 = (G*M/(gamma*a)**3)**0.5
        nu0 = alpha*cs**2/Omega0
        self.chi = 2*fudge*q**2*np.sqrt(G*M)/nu0/a*(gamma*a)**1.5
        self.T0 = mu*Omega0/alpha/k*nu0
        self.gamma = gamma
        self.fudge = fudge
        self.q = q
        self.nu0 = nu0
        self.t = 0.0
        self.odir = odir
        self.bellLin = bellLin
        self._genGrid()
        self.r = self.mesh.cellCenters.value[0]
        self.rF = self.mesh.faceCenters.value[0]
        self.gap = np.where(self.rF < 1.7/gamma)
        self._genSigma()
        self._genTorque()
        if bellLin:
            # Pass the radial grid in phsyical units
            # Get back interpolator in logarithmic space
            log10Interp = thermopy.buildInterpolator(self.r, self.gamma, self.q, self.fudge, self.mDisk, **kargs)
            # Define callable function T(Sigma)
            rGrid = np.log10(self.r)
            r = self.r*a*self.gamma #In physical units (cgs)
            self.Omega = np.sqrt(G*M/r**3)
            Ti = np.power(np.square(eta/7*L/4/np.pi/sigma)*k/mu/G/M*r**(-3), 1.0/7)
            self._T = np.zeros(Ti.shape)
            self.T = CellVariable(name='Temperature',
                                 mesh=self.mesh, hasOld=True, value=Ti)
            # Define wrapper function that uses the interpolator and stores the results
            # in an array given as a second argument. It can handle zero or negative
            # Sigma values.
            def func(Sigma, T):
                good = np.where(Sigma > 0.0); bad = np.where(Sigma <= 0.0)
                T[good] = np.power(10.0, log10Interp.ev(rGrid[good], np.log10(Sigma[good])))
                T[bad] = Ti[bad]
            # Store interpolator as an instance method
            self._bellLinT = func
            self._bellLinUpdate()
        else:
            self._genT()
        self._genVr()
        self._buildEq()

    def _genGrid(self, inB=1.0):
        """Generate a logarithmically spaced grid"""
        logFaces = np.linspace(-np.log(self.gamma/inB), np.log(self.rmax), num=self.ncell+1)
        logFacesLeft = logFaces[:-1]
        logFacesRight = logFaces[1:]
        dr = tuple(np.exp(logFacesRight) - np.exp(logFacesLeft))
        self.mesh = CylindricalGrid1D(dr=dr, origin=(inB/self.gamma,))

    def _genSigma(self, width=0.1):
        """Create dependent variable Sigma"""
        # Gaussian initial condition
        value = self.mDisk*M/np.sqrt(2*np.pi)/(self.gamma*a*width)*\
                np.exp(-0.5*np.square(self.r-1.0)/width**2)/(2*np.pi*self.gamma*self.r*a)
        value /= M/(self.gamma*a)**2
        idxs = np.where(self.r < 0.1)
        value[idxs] = 0.0
        value = tuple(value)

        # Create the dependent variable and set the boundary conditions
        # to zero
        self.Sigma = CellVariable(name='Surface density',
                                 mesh=self.mesh, hasOld=True, value=value)
        self.Sigma.constrain(0, self.mesh.facesLeft)
        self.Sigma.constrain(0, self.mesh.facesRight)

    def _genTorque(self):
        """Generate Torque"""
        self.Lambda = FaceVariable(name='Torque at cell faces', mesh=self.mesh, rank=1)
        self.LambdaCell = CellVariable(name='Torque at cell centers', mesh=self.mesh)
        LambdaArr = np.zeros(self.rF.shape)
        LambdaArr[1:] = self.chi*np.power(1.0/(self.rF[1:]*self.gamma-1.0), 4)
        #LambdaArr[self.gap] = 0.0; LambdaArr[self.gap] = LambdaArr.max()
        self.Lambda.setValue(LambdaArr)
        self.LambdaCell.setValue(self.chi*np.power(1.0/(self.r*self.gamma-1.0), 4))
        self.LambdaCell[np.where(self.LambdaCell > LambdaArr.max())] = LambdaArr.max()

    def _interpT(self):
        """
        Get an initial guess for T using an interpolation of the solutions for T
        in the various thermodynamic limits.
        """
        Lambda = self.Lambda/self.chi*self.fudge*self.q**2*G*M/a
        LambdaCell = self.LambdaCell/self.chi*self.fudge*self.q**2*G*M/a
        Sigma = self.Sigma*M/(self.gamma*a)**2
        r = self.r*a*self.gamma #In physical units (cgs)
        self.Omega = np.sqrt(G*M/r**3)
        self.TvThin = np.power(9.0/4*alpha*k/sigma/mu/kappa0*self.Omega, 1.0/(3.0+beta))
        self.TtiThin = np.power(1/sigma/kappa0*(OmegaIn-self.Omega)*LambdaCell, 1.0/(4.0+beta))
        self.Ti = np.power(np.square(eta/7*L/4/np.pi/sigma)*k/mu/G/M*r**(-3), 1.0/7)
        self.TvThick = np.power(27.0/64*kappa0*alpha*k/sigma/mu*self.Omega*Sigma**2, 1.0/(3.0-beta))
        self.TtiThick = np.power(3*kappa0/16/sigma*Sigma**2*(OmegaIn-self.Omega)*LambdaCell, 1.0/(4.0-beta))
        #return np.power(self.TvThin**4 + self.TvThick**4 + self.TtiThin**4 + self.TtiThick**4 + self.Ti**4, 1.0/4)/self.T0
        return np.power(self.TvThin**4 + self.TvThick**4 + self.Ti**4, 1.0/4)/self.T0

    def _genT(self):
        """Create a cell variable for temperature"""
        # Initialize T with the interpolation of the various thermodynamic limits
        #r = self.r*a #In physical units (cgs)
        #self.Omega = np.sqrt(G*M/r**3)
        #self.TvThin = np.power(9.0/4*alpha*k/sigma/mu/kappa0*self.Omega, 1.0/(3.0+beta))
        #self.Ti = np.power(np.square(eta/7*L/4/np.pi/sigma)*k/mu/G/M*r**(-3), 1.0/7)
        self.T = self._interpT()

    def _genVr(self):
        """Generate the face variable that stores the velocity values"""
        r = self.r #In dimensionless units (cgs)
        # viscosity at cell centers in cgs
        self.nu = alpha*k*self.T/mu/self.Omega/self.nu0
        self.visc = r**0.5*self.nu*self.Sigma
        # I add the delta to avoid divisions by zero
        self.vrVisc = -3/self.rF**(0.5)/(self.Sigma.faceValue + self.delta)*self.visc.faceGrad
        self.vrTid = self.Lambda*np.sqrt(self.rF)

    def _buildEq(self):
        """
        Build the equation to solve, we can change this method to impelement other
        schemes, e.g. Crank-Nicholson.
        """
        # The current scheme is an implicit-upwind
        self.eq = TransientTerm(var=self.Sigma) == - ExponentialConvectionTerm(coeff=self.vrVisc + self.vrTid, var=self.Sigma)

    def _bellLinUpdate(self):
        """
        Update the temperature using the Bell & Lin opacities
        """
        self._bellLinT(self.Sigma.value, self._T)
        self.T.setValue(self._T)

    def singleTimestep(self, dt=None, update=True, emptyDt=False):
        """
        Evolve the system for a single timestep of size `dt`
        """
        if dt:
            self.dt = dt
        if emptyDt:
            vr = self.vrVisc.value[0] + self.vrTid.value[0]
            #vr[np.where(self.Sigma.value)] = self.delta
            self.flux = self.rF[1:]*vr[1:]-self.rF[:-1]*vr[:-1]
            self.flux = np.maximum(self.flux, self.delta)
            self.dts = self.mesh.cellVolumes/(self.flux)
            self.dts[np.where(self.Sigma.value == 0.0)] = np.inf
            self.dts[self.gap] = np.inf
            self.dt = 0.1*np.amin(self.dts)
        try:
            for i in range(self.nsweep):
                res = self.eq.sweep(dt=self.dt)
                if self.bellLin:
                    self._bellLinUpdate()
            if update:
                self.Sigma.updateOld()
            self.t += self.dt
        except FloatingPointError:
            import ipdb; ipdb.set_trace()

    def evolve(self, **kargs):
        """
        Evolve the system according to the values in its initialization
        self.dt, self.nstep, and self.nsweep
        """
        for i in range(self.nstep):
            self.singleTimestep(**kargs)

    def revert(self):
        """
        Revert evolve method if update=False was used, otherwise
        it has no effect.
        """
        self.Sigma.setValue(self.Sigma.old.value)

    def writeToFile(self):
        fName = self.odir + '/t{0}.pkl'.format(self.t)
        with open(fName, 'wb') as f:
            pickle.dump((self.t, self.Sigma.getValue()), f)

    def readFromFile(self, fName):
        with open(fName, 'rb') as f:
            t, Sigma = pickle.load(f)
        self.t = t
        self.Sigma.setValue(Sigma)
        if self.bellLin:
            self._bellLinUpdate()

    def loadTimesList(self):
        path = self.odir
        files = os.listdir(path)
        self.times = np.zeros((len(files)-1,))
        for i, f in enumerate(files[1:]):
            match = re.match(r"^t((\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?)\.pkl", f)
            if match == None:
                print "WARNING: File {0} has an unexepected name".format(f)
                continue
            self.times[i] = float(match.group(1))
        self.files = files[1:]

    def loadTime(self, t):
        """
        Load the file with the time closest to `t`
        """
        idx = (np.abs(self.times-t)).argmin()
        fName = self.odir + '/t'+str(self.times[idx]) + '.pkl'
        self.readFromFile(fName)

def loadRestults(path):
    fName = os.path.join(path, 'init.pkl')
    with open(fName, 'rb') as f:
        kargs = pickle.load(f)
    kargs['odir'] = path # In case the path has been changed
    circ = Circumbinary(**kargs)
    circ.loadTimesList()
    iMax = circ.times.argmax()
    fMax = os.path.join(path, circ.files[iMax])
    circ.readFromFile(fMax)
    return circ


def run(**kargs):
    tmax = kargs.get('tmax')
    kargs.pop('tmax')
    circ = Circumbinary(**kargs) 
    with open(circ.odir+'/init.pkl', 'wb') as f:
        pickle.dump(kargs, f)
    while circ.t < tmax:
        circ.evolve(emptyDt=True)
        circ.writeToFile()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
             description="Script that solves the convection problem in a cylindrical grid")
    parser.add_argument('--rmax', default=1.0e3, type=float,
                        help='The outer boundary of the grid in dimensionless units (r/rMin)')
    parser.add_argument('--ncell', default=200, type=int,
                        help='The number of cells to use in the grid')
    parser.add_argument('--nstep', default=100, type=int,
                        help='The number of time steps to do')
    parser.add_argument('--nsweep', default=10, type=int,
                        help='The number of sweeps to do')
    parser.add_argument('--titer', default=10, type=int,
                        help='The number of temprature iterations')
    parser.add_argument('--dt', default=1.0e-6, type=float,
                        help='The time step size (Constant for the moment)')
    parser.add_argument('--delta', default=1.0e-100, type=float,
                        help='Small number to add to avoid divisions by zero')
    parser.add_argument('--odir', default='output', type=str,
                        help='Directory where results are saved to')
    parser.add_argument('--tmax', default=3.0, type=float,
                        help='Maximum time to evolve the model to')
    kargs = vars(parser.parse_args())
    run(**kargs)
