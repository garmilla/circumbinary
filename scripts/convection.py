import argparse
import numpy as np

from fipy import CylindricalGrid1D, CellVariable, FaceVariable, TransientTerm, ExponentialConvectionTerm

from constants import *

def runScript(rmax = 1.0e3, ncell = 100):
    # Generate a logarithmically spaced grid
    logFaces = np.linspace(0, np.log(rmax), num=ncell+1)
    logFacesLeft = logFaces[:-1]
    logFacesRight = logFaces[1:]
    dr = tuple(np.exp(logFacesRight) - np.exp(logFacesLeft))
    mesh = CylindricalGrid1D(dr=dr, origin=(1,))

    # Simplest possible initial condition,
    # zero everywhere except for a single cell at the middle.
    r = mesh.cellCenters.value[0]
    value = np.zeros(r.shape)
    value[len(value)/2] = 10000.0
    value = tuple(value)

    # Create the dependent variable and set the boundary conditions
    # to zero
    Sigma = CellVariable(name='Dimensionless surface density',
                         mesh=mesh, value=value)
    Sigma.constrain(0, mesh.facesLeft)
    Sigma.constrain(0, mesh.facesRight)

    # Compute viscosity at cell centers in cgs
    r = r0*mesh.cellCenters.value[0] #In physical units (cgs)
    Omega = np.sqrt(G*M/r**3)
    TvThin = np.power(9.0/4*alpha*k/sigma/mu/kappa0*Omega, 1.0/(3.0+beta))
    TvThick = np.power(27.0/64*kappa0*alpha*k/sigma/mu*Omega*Sigma**2, 1.0/(3.0-beta))
    Ti = np.power(np.square(eta/7*L/4/np.pi/sigma)*k/mu/G/M*r**(-3), 1.0/7)
    T = np.power(TvThin**4 + TvThick**4 + Ti**4, 1.0/4)
    nu = alpha*k*T/mu/Omega
    r /= r0 # Get back the dimensionless values of the radius

    # Create a cell variable to store the values inside the derivative
    visc = CellVariable(name='Viscous Term', mesh=mesh)
    visc.setValue(r**0.5*nu*Sigma)

    # Create a face variable to store the values of the radial velocity
    vr = FaceVariable(name='Radial Velocity', mesh=mesh, rank=1)
    rF = mesh.faceCenters.value # radii at the cell faces
    # I add the 1.0e-10 to avoid divisions by zero
    vr.setValue(3/r0**2/rF**(0.5)/(Sigma.faceValue + 1.0e-10)*visc.faceGrad())

    # Build the equation
    eq = TransientTerm() ==ExponentialConvectionTerm(coeff=vr)

    # Solve the equation for a time step of dt=1.0e7
    eq.solve(var=Sigma, dt=1.0e7)
    print Sigma

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
             description="Script that solves the diffusion problem in a cylindrical grid")
    parser.add_argument('--rmax', default=1.0e3, type=float,
                        help='The outer boundary of the grid in dimensionless units (r/rMin)')
    parser.add_argument('--ncell', default=100, type=int,
                        help='The number of cells to use in the grid')
    kargs = vars(parser.parse_args())
    runScript(**kargs)
