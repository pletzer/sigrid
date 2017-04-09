import numpy
from ctypes import c_int, c_double

def createRectilinearGrid(dims, xmin, xmax):
	x0 = numpy.linspace(xmin[0], xmax[0], dims[0])
	x1 = numpy.linspace(xmin[1], xmax[1], dims[1])
	xx0 = numpy.outer(x0, numpy.ones((dims[1],), numpy.float64))
	xx1 = numpy.outer(numpy.ones((dims[0],), numpy.float64), x1)
	return (xx0, xx1)

def saveStreamlinesVtk(filename, coords, psi):

    zeros = (0, 0)
    dims = psi.shape

    f = open(filename, 'w')

    f.write("# vtk DataFile Version 2.0\n")
    f.write("streamFunctionInterp\n")
    f.write("ASCII\n")
    f.write("DATASET STRUCTURED_GRID\n")
    f.write("DIMENSIONS 1 {} {}\n".format(dims[1], dims[0])) # inverse order
    numPoints = dims[0] * dims[1]
    f.write("POINTS {} float\n".format(numPoints))
    for i in range(dims[0]):
    	for j in range(dims[1]):
        	f.write("{} {} 0.0\n".format(coords[0][i, j], coords[1][i, j]))
 
    f.write("POINT_DATA {}\n".format(numPoints))
    f.write("SCALARS psi float\n")
    f.write("LOOKUP_TABLE default\n")
    for i in range(dims[0]):
    	for j in range(dims[1]):
        	f.write("{}\n".format(psi[i, j]))
        
    numCells = (dims[0] - 1) * (dims[1] - 1)

    f.write("CELL_DATA {}\n".format(numCells))
    f.write("SCALARS velocity float 3\n")
    f.write("LOOKUP_TABLE default\n")

    vx = -(psi[:, 1:] - psi[:, :-1])/(coords[1][:, 1:] - coords[1][:, :-1]) #- dpsi/dy
    vy = +(psi[1:, :] - psi[:-1, :])/(coords[0][1:, :] - coords[0][:-1, :]) #+ dpsi/dx

    # average the velocity field to cell centres
    vxAvg = 0.5*(vx[:-1, :] + vx[1:, :])
    vyAvg = 0.5*(vy[:, :-1] + vy[:, 1:])

    for i in range(dims[0] - 1):
    	for j in range(dims[1] - 1):
        	f.write("{} {} 0.0\n".format(vxAvg[i, j], vyAvg[i, j]))

    f.close();

