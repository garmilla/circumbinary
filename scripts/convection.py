import argparse
import numpy as np

from fipy import CylindricalGrid1D, CellVariable, FaceVariable, TransientTerm, UpwindConvectionTerm

from constants import *

class circumbinary(object):
    def __init__(self, rmax=1.0e3, ncell=100, nstep=100, dt=1.0e7, delta=1.0e-10, nsweep=10):
        self.rmax = rmax
        self.ncell = ncell
        self.nstep = nstep
        self.dt = dt
        self.delta = delta
        self.nsweep = nsweep
        self._genGrid()
        self._genSigma()
        self._genVr()
        self._buildEq()

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

    def _genVr(self):
        """Generate the face variable that stores the velocity values"""
        r = self.mesh.cellCenters.value[0]
        r *= r0 #In physical units (cgs)
        self.Omega = np.sqrt(G*M/r**3)
        self.TvThin = np.power(9.0/4*alpha*k/sigma/mu/kappa0*self.Omega, 1.0/(3.0+beta))
        self.Ti = np.power(np.square(eta/7*L/4/np.pi/sigma)*k/mu/G/M*r**(-3), 1.0/7)
        TvThick = np.power(27.0/64*kappa0*alpha*k/sigma/mu*self.Omega*self.Sigma**2, 1.0/(3.0-beta))
        T = np.power(self.TvThin**4 + TvThick**4 + self.Ti**4, 1.0/4)
        r /= r0
        # viscosity at cell centers in cgs
        nu = alpha*k*T/mu/self.Omega
        self.visc = CellVariable(name='Viscous Term', mesh=self.mesh)
        self.visc.setValue(r**0.5*nu*self.Sigma)
        # Create a face variable to store the values of the radial velocity
        self.vr = FaceVariable(name='Radial Velocity', mesh=self.mesh, rank=1)
        rF = self.mesh.faceCenters.value # radii at the cell faces
        # I add the delta to avoid divisions by zero
        self.vr.setValue(-3/r0**2/rF**(0.5)/(self.Sigma.faceValue + self.delta)*self.visc.faceGrad())

    def _buildEq(self):
        """
        Build the equation to solve, we can change this method to impelement other
        schemes, e.g. Crank-Nicholson.
        """
        # The current scheme is an implicit-upwind
        self.eq = TransientTerm() == - UpwindConvectionTerm(coeff=self.vr)

    def _updateVr(self):
        """
        Update the value of the velocity. This function assumes the value of
        Sigma has been updated by the solve or sweep method.
        """
        r = self.mesh.cellCenters.value[0]
        TvThick = np.power(27.0/64*kappa0*alpha*k/sigma/mu*self.Omega*self.Sigma**2, 1.0/(3.0-beta))
        T = np.power(self.TvThin**4 + TvThick**4 + self.Ti**4, 1.0/4)
        nu = alpha*k*T/mu/self.Omega
        self.visc.setValue(r**0.5*nu*self.Sigma)
        rF = self.mesh.faceCenters.value # radii at the cell faces
        self.vr.setValue(-3/r0**2/rF**(0.5)/(self.Sigma.faceValue + self.delta)*self.visc.faceGrad())

    def singleTimestep(self, dt=None):
        """
        Evolve the system for a single timestep of size `dt`
        """
        if not dt:
            dt = self.dt
        try:
            for i in range(self.nsweep):
                self.eq.sweep(var=self.Sigma, dt=dt)
                self._updateVr()
            self.Sigma.updateOld()
        except FloatingPointError:
            import ipdb; ipdb.set_trace()

    def evolve(self):
        """
        Evolve the system according to the values in its inisialization
        self.dt, self.nstep, and self.nsweep
        """
        for i in range(self.nstep):
            self.singleTimestep()

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
    parser.add_argument('--dt', default=1.0e7, type=float,
                        help='The time step size (Constant for the moment)')
    parser.add_argument('--delta', default=1.0e-10, type=float,
                        help='Small number to add to avoid divisions by zero')
    kargs = vars(parser.parse_args())
    run(**kargs)
