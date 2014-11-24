import argparse
import numpy as np

from fipy import CylindricalGrid1D, CellVariable

def runScript(rmax = 1.0e3, ncell = 100):
    # Generate a logarithmically spaced grid
    etaFaces = np.linsapce(0, np.log(rmax), num=ncell+1)
    etaFacesLeft = etaFaces[:-1]
    etaFacesRight = etaFaces[1:]
    dr = tuple(np.exp(etaFacesRight) - np.exp(etaFacesLeft))
    mesh = CylindricalGrid1D(dr=dr, origin=(1,))

    # Simplest possible initial condition,
    # zero everywhere except for a single cell at the middle.
    cellCenters = mesh.cellCenters.value
    value = np.zeros(cellCenters.shape)
    value[len(value)/2] = 10.0
    value = tuple(value)

    sigma = CellVariable(name='Dimensionless surface density',
                         mesh=mesh,
                         value=value)
    # Build the equation
    eq = TransientTerm() == DiffusionTerm(coeff= (1 - phi[0]))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
             description="Script that solves the diffusion problem in a cylindrical grid")
    parser.add_argument('--rmax', default=1.0e3, type=float,
                        help='The outer boundary of the grid in dimensionless units (r/rMin)')
    parser.add_argument('--ncell', default=100, type=int,
                        help='The number of cells to use in the grid')
    kargs = vars(parser.parse_args())
    runScript(**kargs)
