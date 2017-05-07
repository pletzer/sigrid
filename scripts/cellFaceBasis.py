"""
Cell face basis functions
"""

import numpy
import argparse
from numpy import pi, sin, cos

class FaceBasis:

    def __init__(self, dims, xExpr, yExpr):
        
        # build the mesh
        dXi = 1.0 / float(dims[0] - 1)
        dEta = 1.0 / float(dims[1] - 1)
        xi1d = numpy.linspace(0., 1., dims[0])
        eta1d = numpy.linspace(0., 1., dims[1])
        xi = numpy.outer(xi1d, numpy.ones(dims[1], numpy.float64))
        eta = numpy.outer(numpy.ones(dims[0], numpy.float64), eta1d)

        # apply coordinate transformations
        x = eval(xExpr)
        y = eval(yExpr)

        # compute averages
        x1Avg = 0.5*(x[:-1, :] + x[1:, :]) # x average in the first coordinate
        x2Avg = 0.5*(x[:, :-1] + x[:, 1:]) # x average in the second coordinate
        y1Avg = 0.5*(y[:-1, :] + y[:1, :])
        y2Avg = 0.5*(y[:, :-1] + y[:, 1:])

        # partial derivatives
        dxdxi = (x2Avg[1:, :] - x2Avg[:-1, :]) / dXi
        dxdeta = (x1Avg[:, 1:] - x1Avg[:, :-1]) / dEta
        dydxi = (y2Avg[1:, :] - y2Avg[:-1, :]) / dXi
        dydeta = (y1Avg[:, 1:] - y1Avg[:, :-1]) / dEta

        # make xi and eta cell centred
        xi = 0.25*(xi[:-1, :-1] + xi[:-1, 1:] + xi[1:, 1:] + xi[1:, :-1])
        eta = 0.25*(eta[:-1, :-1] + eta[:-1, 1:] + eta[1:, 1:] + eta[1:, :-1])

        # Jacobian
        jac = dxdxi * dydeta - dydxi * dxdeta

        # basis functions
        self.bases = {(-1, 0): [(1 - eta)*dxdeta/jac, (1 - eta)*dydeta/jac], # ok
                      (-1, 1): [(0 + eta)*dxdeta/jac, (0 + eta)*dydeta/jac], # ok
                      (0, -1): [(1 - xi)*dxdxi/jac, (1 - xi)*dydxi/jac], # looks wrong
                      (1, -1): [(0 + xi)*dxdxi/jac, (0 + xi)*dydxi/jac], # looks wrong
                      }

        self.x = x
        self.y = y

    def saveVtk(self, filename):

        dims = self.x.shape

        f = open(filename, 'w')

        f.write("# vtk DataFile Version 2.0\n")
        f.write("bases\n")
        f.write("ASCII\n")
        f.write("DATASET STRUCTURED_GRID\n")
        f.write("DIMENSIONS 1 {} {}\n".format(dims[1], dims[0])) # inverse order
        numPoints = dims[0] * dims[1]
        f.write("POINTS {} float\n".format(numPoints))
        for i in range(dims[0]):
            for j in range(dims[1]):
                f.write("{} {} 0.0\n".format(self.x[i, j], self.y[i, j]))

        numCells = (dims[0] - 1) * (dims[1] - 1)

        f.write("CELL_DATA {}\n".format(numCells))

        # iterate over the basis functions
        for side, basis in self.bases.items():
            basisName = '_{}_{}'.format(side[0], side[1])
            f.write("SCALARS {} float 3\n".format(basisName))
            f.write("LOOKUP_TABLE default\n")
            vx, vy = basis[0], basis[1]
            for i in range(dims[0] - 1):
                for j in range(dims[1] - 1):
                    f.write("{} {} 0.0".format(vx[i, j], vy[i, j]))

        f.close()

###############################################################################

def main():
    
    parser = argparse.ArgumentParser(description='Generate VTK plots of the face basis functions in 2D.')
    parser.add_argument('--dims', dest='dims', default='11,11',
                    help='Number of nodes in xi and eta')
    parser.add_argument('--xExpr', dest='xExpr', default='xi*cos(eta) - 0.2*eta',
                    help='X coordinate transformation as a function of xi and eta')
    parser.add_argument('--yExpr', dest='yExpr', default='eta/2.0 + 0.4*sin(0.3*pi*xi*eta)',
                    help='Y coordinate transformation as a function of xi and eta')
    parser.add_argument('--filename', dest='filename', default='cellFaceBasis.vtk',
                    help='Filename for storing the cell centred vector field of each basis function')
    args = parser.parse_args()

    dims = eval(args.dims)
    fb = FaceBasis(dims, xExpr=args.xExpr, yExpr=args.yExpr)
    print('Now saving the basis functions in file {}'.format(args.filename))
    fb.saveVtk(args.filename)


if __name__ == '__main__': main()
