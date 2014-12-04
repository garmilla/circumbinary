import argparse
import numpy as np
from scipy.optimize import root

from fipy import CylindricalGrid1D, CellVariable, FaceVariable, TransientTerm, UpwindConvectionTerm

from constants import *

# We can get the exact solution for T by finding the root of this function
# I doesn't work!
def fun(T, self, delta = 1.0e-20):
    r = self.mesh.cellCenters.value[0]*r0
    low = np.where(T <= 166.81)
    mid = np.where(np.logical_and(T > 166.81, T < 202.677))
    high = np.where(T >= 202.677)
    self.kappa[low] = 2.0e-4*np.square(T[low])
    self.kappa[mid] = 2.0e16*np.power(T[mid], -7)
    self.kappa[high] = 0.1*np.sqrt(T[high])
    tau = 0.5*self.kappa*self.Sigma
    nu = alpha*k/mu/self.Omega*T
    Fnu = 9.0/8*self.Omega**2*nu*self.Sigma
    self.hr.setValue(eta*np.power(k*T/G/M/mu, 1.5)*np.sqrt(r))
    Firr = 0.5*L/4/np.pi/r*np.maximum(self.hr.grad.value[0], 0.0)
    if np.count_nonzero(T) == 0:
        return fun(T+delta, self, delta=delta)
    else:
        return (3*tau/4 + 1.0/(tau+delta))*Fnu + Firr - sigma*T**4

class circumbinary(object):
    def __init__(self, rmax=1.0e3, ncell=100, nstep=100, dt=1.0e7, delta=1.0e-10,
                 nsweep=10, titer=10, fudge=1.0e-3, q=1.0):
        self.rmax = rmax
        self.ncell = ncell
        self.nstep = nstep
        self.dt = dt
        self.delta = delta
        self.nsweep = nsweep
        self.titer = titer
        self.fudge = fudge
        self.q = q
        self._genGrid()
        self.r = self.mesh.cellCenters.value[0]
        self._genSigma()
        self._genTorque()
        self._genT()
        self._genVr()
        self._buildEq()
        self.kappa = np.zeros((ncell,))

    def _genGrid(self):
        """Generate a logarithmically spaced grid"""
        logFaces = np.linspace(0, np.log(self.rmax), num=self.ncell+1)
        logFacesLeft = logFaces[:-1]
        logFacesRight = logFaces[1:]
        dr = tuple(np.exp(logFacesRight) - np.exp(logFacesLeft))
        self.mesh = CylindricalGrid1D(dr=dr, origin=(1,))

    def _genSigma(self):
        """Create dependent variable Sigma"""
        # Simplest possible initial condition,
        # zero everywhere except for a single cell at the middle.
        r = self.mesh.cellCenters.value[0]
        value = np.zeros(r.shape)
        value[len(value)/2] = 10000.0
        value = tuple(value)

        # Create the dependent variable and set the boundary conditions
        # to zero
        self.Sigma = CellVariable(name='Surface density',
                                 mesh=self.mesh, hasOld=True, value=value)
        self.Sigma.constrain(0, self.mesh.facesLeft)
        self.Sigma.constrain(0, self.mesh.facesRight)

    def _genTorque(self):
        """Generate Torque"""
        rF = self.mesh.faceCenters.value*r0 # radii at the cell faces
        self.Lambda = np.zeros(rF.shape)
        self.Lambda[0][1:] = self.fudge*self.q**2*M/r0*np.power(r0/(rF[0][1:]-r0), 4)

    def _genT(self):
        """Create a cell variable for temperature"""
        self.T = CellVariable(name='Temperature',
                                 mesh=self.mesh, hasOld=True)
        # Initialize T with the interpolation of the various thermodynamic limits
        r = self.r*r0 #In physical units (cgs)
        self.Omega = np.sqrt(G*M/r**3)
        self.TvThin = np.power(9.0/4*alpha*k/sigma/mu/kappa0*self.Omega, 1.0/(3.0+beta))
        self.Ti = np.power(np.square(eta/7*L/4/np.pi/sigma)*k/mu/G/M*r**(-3), 1.0/7)
        self.T.setValue(self._interpT())
        self.T.updateOld()
        # Create a binary operator variable for h/r
        self.hr = CellVariable(name='Ratio of thickness to radius, h/r',
                                 mesh=self.mesh, hasOld=True)
        self.hr.setValue(eta*np.power(k*self.T/G/M/mu, 1.5)*np.sqrt(r))

    def _genVr(self):
        """Generate the face variable that stores the velocity values"""
        r = self.r #In dimensionless units (cgs)
        # viscosity at cell centers in cgs
        nu = alpha*k*self.T/mu/self.Omega
        self.visc = CellVariable(name='Viscous Term', mesh=self.mesh)
        self.visc.setValue(r**0.5*nu*self.Sigma)
        # Create a face variable to store the values of the radial velocity
        self.vr = FaceVariable(name='Radial Velocity', mesh=self.mesh, rank=1)
        rF = self.mesh.faceCenters.value # radii at the cell faces
        # I add the delta to avoid divisions by zero
        self.vr.setValue(-3/r0**2/rF**(0.5)/(self.Sigma.faceValue + self.delta)*self.visc.faceGrad()
                         + 2/np.sqrt(r0)*self.Lambda*np.sqrt(rF/G/M))

    def _buildEq(self):
        """
        Build the equation to solve, we can change this method to impelement other
        schemes, e.g. Crank-Nicholson.
        """
        # The current scheme is an implicit-upwind
        self.eq = TransientTerm() == - UpwindConvectionTerm(coeff=self.vr)

    def _updateVr(self, mode='interp'):
        """
        Update the value of the velocity. This function assumes the value of
        Sigma has been updated by the solve or sweep method.

        Keywords
        mode: If equal to `interp`, use the simple interpolation scheme to update
              the temperature. If equal to `exact`, attempt an exact solutuion of
              the temperature equation using a nonlinear solver. If equal to
              `iterate`, use an iteration scheme on the temperature equation to
              approximate the solution.
        """
        r = self.r
        if mode == 'interp':
            self.T.setValue(self._interpT())
            self.T.updateOld()
        elif mode == 'exact':
            self._exactT()
        elif mode == 'iterate':
            self._iterT()
        #nu = alpha*k*self.T/mu/self.Omega
        nu = 6.0e14
        self.visc.setValue(r**0.5*nu*self.Sigma)
        rF = self.mesh.faceCenters.value # radii at the cell faces
        self.vr.setValue(-3/r0**2/rF**(0.5)/(self.Sigma.faceValue + self.delta)*self.visc.faceGrad())

    def _interpT(self):
        """
        Get an initial guess for T using an interpolation of the solutions for T
        in the various thermodynamic limits.
        """
        r = self.r*r0 #In physical units (cgs)
        self.Omega = np.sqrt(G*M/r**3)
        self.TvThin = np.power(9.0/4*alpha*k/sigma/mu/kappa0*self.Omega, 1.0/(3.0+beta))
        self.Ti = np.power(np.square(eta/7*L/4/np.pi/sigma)*k/mu/G/M*r**(-3), 1.0/7)
        TvThick = np.power(27.0/64*kappa0*alpha*k/sigma/mu*self.Omega*self.Sigma**2, 1.0/(3.0-beta))
        return np.power(self.TvThin**4 + TvThick**4 + self.Ti**4, 1.0/4)

    def _iterT(self, delta=1.0e-20):
        """
        Do an iteration on temperature
        """
        r = self.r*r0
        T = self.T.value
        low = np.where(T <= 166.81)
        mid = np.where(np.logical_and(T > 166.81, T < 202.677))
        high = np.where(T >= 202.677)
        self.kappa[low] = 2.0e-4*np.square(T[low])
        self.kappa[mid] = 2.0e16*np.power(T[mid], -7)
        self.kappa[high] = 0.1*np.sqrt(T[high])
        tau = 0.5*self.kappa*self.Sigma
        nu = alpha*k/mu/self.Omega*T
        Fnu = 9.0/8*self.Omega**2*nu*self.Sigma
        self.hr.setValue(eta*np.power(k*T/G/M/mu, 1.5)*np.sqrt(r))
        Firr = 0.5*L/4/np.pi/r*np.maximum(self.hr.grad.value[0], 0.0)
        self.T.setValue(np.power((3*tau/4 + 1.0/(tau+delta))*Fnu + Firr, 0.25)/sigma)

    def iterT(self):
        """
        Iterate over temperature to get a more reliable solution
        """
        for i in range(self.titer):
            self._iterT()

    def _exactT(self, init='old'):
        """
        Get the exact solution for T
        keywords:
        init: If equal to `old` then the old value of T is used as an initial guess, if
              equal to `inter` then self._interpT is used to get the initial guess.
        """
        # I have not been able to get this to work
        if init == 'old':
            T0 = self.T.old.value
        elif init == 'interp':
            T0 = self._interpT()

        result = root(fun, T0, args=(self,))

        if result.success:
            self.T.setValue(result.x)
            self.T.updateOld()
        else:
            print result.message
            print 'T was not updated'
            print result.x

    def singleTimestep(self, dt=None, update=True):
        """
        Evolve the system for a single timestep of size `dt`
        """
        if not dt:
            dt = self.dt
        try:
            for i in range(self.nsweep):
                self.eq.sweep(var=self.Sigma, dt=dt)
                self._updateVr()
            if update:
                self.Sigma.updateOld()
        except FloatingPointError:
            import ipdb; ipdb.set_trace()

    def evolve(self, dt=None, update=True):
        """
        Evolve the system according to the values in its initialization
        self.dt, self.nstep, and self.nsweep
        """
        for i in range(self.nstep):
            self.singleTimestep(dt=dt, update=update)

    def revert(self):
        """
        Revert evolve method if update=False was used, otherwise
        it has no effect.
        """
        self.Sigma.setValue(self.Sigma.old.value)

def run(**kargs):
    circ = circumbinary(**kargs) 
    circ.evolve()
    print circ.Sigma
    return circ.Sigma

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
             description="Script that solves the convection problem in a cylindrical grid")
    parser.add_argument('--rmax', default=1.0e3, type=float,
                        help='The outer boundary of the grid in dimensionless units (r/rMin)')
    parser.add_argument('--ncell', default=100, type=int,
                        help='The number of cells to use in the grid')
    parser.add_argument('--nstep', default=100, type=int,
                        help='The number of time steps to do')
    parser.add_argument('--nsweep', default=10, type=int,
                        help='The number of sweeps to do')
    parser.add_argument('--titer', default=10, type=int,
                        help='The number of temprature iterations')
    parser.add_argument('--dt', default=1.0e7, type=float,
                        help='The time step size (Constant for the moment)')
    parser.add_argument('--delta', default=1.0e-10, type=float,
                        help='Small number to add to avoid divisions by zero')
    kargs = vars(parser.parse_args())
    run(**kargs)
