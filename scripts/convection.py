import os
import re
import argparse
import pickle
import numpy as np
from scipy.interpolate import RectBivariateSpline

from fipy import CylindricalGrid1D, CellVariable, FaceVariable, TransientTerm, ExplicitUpwindConvectionTerm,\
                 ExponentialConvectionTerm, ImplicitSourceTerm, UpwindConvectionTerm
from fipy.steppers import sweepMonotonic
from fipy.boundaryConditions import FixedFlux

import thermopy
from constants import *
from utils import pickle_results

class Circumbinary(object):
    def __init__(self, rmax=1.0e4, ncell=300, dt=1.0e-6, delta=1.0e-100,
                 fudge=1.0e-3, q=1.0, gamma=100, mdisk=0.1, odir='output',
                 bellLin=True, emptydt=0.01, **kargs):
        self.rmax = rmax
        self.ncell = ncell
        self.dt = dt
        self.delta = delta
        self.mDisk = mdisk
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
        self.emptydt = emptydt
        self._genGrid()
        self.r = self.mesh.cellCenters.value[0]
        self.rF = self.mesh.faceCenters.value[0]
        if self.q > 0.0:
            self.gap = np.where(self.rF < 1.7/gamma)
        else:
            self.gap = np.where(self.rF < 1.0/gamma)
        self._genSigma()
        self._genTorque()
        self._genT(bellLin=self.bellLin, **kargs)
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
        # Make it dimensionless
        value /= self.mDisk*M/(self.gamma*a)**2
        idxs = np.where(self.r < 0.1)
        value[idxs] = 0.0
        value = tuple(value)

        # Create the dependent variable and set the boundary conditions
        # to zero
        self.Sigma = CellVariable(name='Surface density',
                                 mesh=self.mesh, hasOld=True, value=value)
        #self.Sigma.constrain(0, self.mesh.facesLeft)
        #self.Sigma.constrain(0, self.mesh.facesRight)

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

    def _genT(self, bellLin=True, **kargs):
        """Create a cell variable for temperature"""
        if bellLin:
            @pickle_results(os.path.join(self.odir, "interpolator.pkl"))
            def buildInterpolator(r, gamma, q, fudge, mDisk, **kargs):
                # Keep in mind that buildTemopTable() returns the log10's of the values
                rGrid, SigmaGrid, temp = thermopy.buildTempTable(r*a*gamma, q=q, f=fudge, **kargs)
                # Go back to dimensionless units
                rGrid -= np.log10(a*gamma)
                SigmaGrid -= np.log10(mDisk*M/gamma**2/a**2)
                # Get the range of values for Sigma in the table
                rangeSigma = (np.power(10.0, SigmaGrid.min()), np.power(10.0, SigmaGrid.max()))
                # Interpolate in the log of dimensionless units
                return rangeSigma, RectBivariateSpline(rGrid, SigmaGrid, temp)
            # Pass the radial grid in phsyical units
            # Get back interpolator in logarithmic space
            rangeSigma, log10Interp = buildInterpolator(self.r, self.gamma, self.q, self.fudge, self.mDisk, **kargs)
            rGrid = np.log10(self.r)
            SigmaMin = np.ones(rGrid.shape)*rangeSigma[0]
            SigmaMax = np.ones(rGrid.shape)*rangeSigma[1]
            r = self.r*a*self.gamma #In physical units (cgs)
            self.Omega = np.sqrt(G*M/r**3)
            Ti = np.power(np.square(eta/7*L/4/np.pi/sigma)*k/mu/G/M*r**(-3), 1.0/7)
            T = np.zeros(Ti.shape)
            # Define wrapper function that uses the interpolator and stores the results
            # in an array given as a second argument. It can handle zero or negative
            # Sigma values.
            def func(Sigma):
                good = np.logical_and(Sigma > rangeSigma[0], Sigma < rangeSigma[1])
                badMin = np.logical_and(True, Sigma < rangeSigma[0])
                badMax = np.logical_and(True, Sigma > rangeSigma[1])
                if np.sum(good) > 0:
                    T[good] = np.power(10.0, log10Interp.ev(rGrid[good], np.log10(Sigma[good])))
                if np.sum(badMin) > 0:
                    T[badMin] = np.power(10.0, log10Interp.ev(rGrid[badMin], np.log10(SigmaMin[badMin])))
                if np.sum(badMax) > 0:
                    raise ValueError("Extrapolation to large values of Sigma is not allowed, build a table with a larger Sigmax")
                    T[badMax] = np.power(10.0, log10Interp.ev(rGrid[badMax], np.log10(SigmaMax[badMax])))
                return T
            # Store interpolator as an instance method
            self._bellLinT = func
            # Save the temperature as an operator variable
            self.T = self.Sigma._UnaryOperatorVariable(lambda x: self._bellLinT(x))

        # Initialize T with the interpolation of the various thermodynamic limits
        else:
            self.T = self._interpT()

    def _genVr(self):
        """Generate the face variable that stores the velocity values"""
        r = self.r #In dimensionless units (cgs)
        # viscosity at cell centers in cgs
        self.nu = alpha*k/mu/self.Omega/self.nu0*self.T
        self.visc = r**0.5*self.nu*self.Sigma
        #self.visc.grad.constrain([self.visc/2/self.r[0]], self.mesh.facesLeft)
        #self.Sigma.constrain(self.visc.grad/self.nu*2*self.r**0.5, where=self.mesh.facesLeft)
        # I add the delta to avoid divisions by zero
        self.vrVisc = -3/self.rF**(0.5)/(self.Sigma.faceValue + self.delta)*self.visc.faceGrad
        if self.q > 0.0:
            self.vrTid = self.Lambda*np.sqrt(self.rF)

    def _buildEq(self):
        """
        Build the equation to solve, we can change this method to impelement other
        schemes, e.g. Crank-Nicholson.
        """
        # The current scheme is an implicit-upwind
        if self.q > 0.0:
            self.vr = self.vrVisc + self.vrTid
            self.eq = TransientTerm(var=self.Sigma) == - ExplicitUpwindConvectionTerm(coeff=self.vr, var=self.Sigma)
        else:
            self.vr = self.vrVisc
            mask_coeff = (self.mesh.facesLeft * self.mesh.faceNormals).getDivergence()
            self.eq = TransientTerm(var=self.Sigma) == - ExplicitUpwindConvectionTerm(coeff=self.vr, var=self.Sigma)\
                                                       - mask_coeff*3.0/2*self.nu/self.mesh.x*self.Sigma.old

    def dimensionalSigma(self):
        """
        Return Sigma in dimensional form (cgs)
        """
        return self.Sigma.value*self.mDisk*M/(self.gamma*a)**2

    def dimensionalFJ(self):
        """
        Return the viscous torque in dimensional units (cgs)
        """
        return 3*np.pi*self.nu.value*self.nu0*self.dimensionalSigma()*np.sqrt(G*M*self.r*a*self.gamma)

    def dimensionalTime(self, mode='yr'):
        """
        Return current time in dimensional units (years or seconds)
        """
        if mode == 'yr' or mode == 'years' or mode == 'year':
            return self.t*(a*self.gamma)**2/self.nu0/(365*24*60*60)
        else:
            return self.t*(a*self.gamma)**2/self.nu0

    def singleTimestep(self, dtMax=0.001, dt=None, update=True, emptyDt=False):
        """
        Evolve the system for a single timestep of size `dt`
        """
        if dt:
            self.dt = dt
        if emptyDt:
            vr = self.vr.value[0]
            if self.q == 0.0:
                vr[0] = -3.0/2*self.nu.faceValue.value[0]/self.rF[0]
            #vr[np.where(self.Sigma.value)] = self.delta
            self.flux = self.rF[1:]*vr[1:]-self.rF[:-1]*vr[:-1]
            self.flux = np.maximum(self.flux, self.delta)
            self.dts = self.mesh.cellVolumes/(self.flux)
            self.dts[np.where(self.Sigma.value == 0.0)] = np.inf
            self.dts[self.gap] = np.inf
            self.dt = self.emptydt*np.amin(self.dts)
        self.dt = min(dtMax, self.dt)
        try:
            if self.q > 0.0:
                self.eq.sweep(dt=self.dt)
            elif self.q == 0.0:
                self.eq.sweep(dt=self.dt)
            if np.any(self.Sigma.value < 0.0):
                self.singleTimestep(dt=self.dt/2)
            if update:
                self.Sigma.updateOld()
            self.t += self.dt
        except FloatingPointError:
            import ipdb; ipdb.set_trace()

    def evolve(self, deltaTime, **kargs):
        """
        Evolve the system using the singleTimestep method
        """
        tEv = self.t + deltaTime
        while self.t < tEv:
            dtMax = tEv - self.t
            self.singleTimestep(dtMax=dtMax, **kargs)

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

    def loadTimesList(self):
        path = self.odir
        files = os.listdir(path)
        if '.DS_Store' in files:
            files.remove('.DS_Store')
        if 'interpolator.pkl' in files:
            files.remove('interpolator.pkl')
        if 'init.pkl' in files:
            files.remove('init.pkl')
        self.times = np.zeros((len(files),))
        for i, f in enumerate(files):
            match = re.match(r"^t((\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?)\.pkl", f)
            if match == None:
                print "WARNING: File {0} has an unexepected name".format(f)
                files.remove(f)
                continue
            self.times[i] = float(match.group(1))
        self.times.sort()
        self.files = files

    def loadTime(self, t):
        """
        Load the file with the time closest to `t`
        """
        idx = (np.abs(self.times-t)).argmin()
        fName = self.odir + '/t'+str(self.times[idx]) + '.pkl'
        self.readFromFile(fName)

def loadResults(path):
    fName = os.path.join(path, 'init.pkl')
    with open(fName, 'rb') as f:
        kargs = pickle.load(f)
    kargs['odir'] = path # In case the path has been changed
    circ = Circumbinary(**kargs)
    circ.loadTimesList()
    try:
        iMax = circ.times.argmax()
        fMax = os.path.join(path, circ.files[iMax])
        circ.readFromFile(fMax)
    except ValueError:
        pass
    return circ

def run(**kargs):
    tmax = kargs.get('tmax')
    dstep = kargs.get('dstep')
    kargs.pop('tmax')
    kargs.pop('dstep')
    fName = os.path.join(kargs['odir'], 'init.pkl')
    if os.path.isfile(fName):
        with open(fName, 'rb') as f:
            kargsPrev = pickle.load(f)
        equal = True
        for k in kargs:
            if not kargs[k] == kargsPrev[k]:
                equal = False
        for k in kargsPrev:
            if not kargs[k] == kargsPrev[k]:
                equal = False
        if not equal:
            raise ValueError("The parameters you typed are different from the parameters in the previous run\n The parameters were {0}".format(kargsPrev))
        else:
           print "I found data in this folder, I'll resume from the last snapshot" 
           circ = loadResults(kargs['odir'])
    else:
        circ = Circumbinary(**kargs)
        with open(circ.odir+'/init.pkl', 'wb') as f:
            pickle.dump(kargs, f)
    tmax = tmax*365*24*60*60 # Convert time to seconds
    dstep = dstep*365*24*60*60 # Convert time to seconds
    tmax = tmax/((a*circ.gamma)**2/circ.nu0) # Go to dimensionless time
    dstep = dstep/((a*circ.gamma)**2/circ.nu0) # Go to dimensionless time
    #import ipdb; ipdb.set_trace()
    while circ.t < tmax:
        circ.evolve(dstep, emptyDt=True)
        circ.writeToFile()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
             description="Script that solves the convection problem in a cylindrical grid")
    parser.add_argument('--rmax', default=1.0e4, type=float,
                        help='The outer boundary of the grid in dimensionless units (r/rMin)')
    parser.add_argument('--ncell', default=300, type=int,
                        help='The number of cells to use in the grid')
    parser.add_argument('--dt', default=1.0e-6, type=float,
                        help='The time step size (Constant for the moment)')
    parser.add_argument('--q', default=1.0, type=float,
                        help='Binary mass ratio (smaller/bigger).')
    parser.add_argument('--fudge', default=0.002, type=float,
                        help='Fudge factor that the torque term is proportional to')
    parser.add_argument('--mdisk', default=0.1, type=float,
                        help='Total mass of the disk in units of central binary mass.')
    parser.add_argument('--emptydt', default=0.01, type=float,
                        help='Factor to use when using the emptyDt=True option')
    parser.add_argument('--delta', default=1.0e-100, type=float,
                        help='Small number to add to avoid divisions by zero')
    parser.add_argument('--odir', default='output', type=str,
                        help='Directory where results are saved to')
    parser.add_argument('--tmax', default=5.0e6, type=float,
                        help='Maximum time to evolve the model to in years')
    parser.add_argument('--dstep', default=5.0e3, type=float,
                        help='Intervals at which to save snapshots in years')
    kargs = vars(parser.parse_args())
    run(**kargs)
