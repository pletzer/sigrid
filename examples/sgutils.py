import numpy
from ctypes import c_int, c_double

def createRectilinearGrid(dims, xmin, xmax):
	x0 = numpy.linspace(xmin[0], xmax[0], dims[0])
	x1 = numpy.linspace(xmin[1], xmax[1], dims[1])
	xx0 = numpy.outer(x0, numpy.ones((dims[1],), numpy.float64))
	xx1 = numpy.outer(numpy.ones((dims[0],), numpy.float64), x1)
	return (xx0, xx1)

def createPolarGrid(dims, radius):
    x0 = numpy.linspace(0.0, radius, dims[0])     # rho
    x1 = numpy.linspace(0.0, 2*numpy.pi, dims[1]) # theta
    rho = numpy.outer(x0, numpy.ones((dims[1],), numpy.float64))
    the = numpy.outer(numpy.ones((dims[0],), numpy.float64), x1)
    xx0 = rho*numpy.cos(the)
    xx1 = rho*numpy.sin(the)
    return (xx0, xx1)

def saveLineVtk(filename, coords, data):

    numPoints = len(coords[0])
    f = open(filename, 'w')

    f.write("# vtk DataFile Version 2.0\n")
    f.write("streamFunctionInterp\n")
    f.write("ASCII\n")
    f.write("DATASET POLYDATA\n")
    f.write("POINTS {} float\n".format(numPoints))
    for i in range(numPoints):
    	f.write("{} {} 0.0\n".format(coords[0][i], coords[1][i]))
    f.write("LINES 1 {}\n".format(numPoints + 1))
    f.write("{} ".format(numPoints))
    for i in range(numPoints):
    	f.write("{} ".format(i))
    f.write("\n")
    numCells = numPoints - 1
    f.write("CELL_DATA {}\n".format(numCells))
    f.write("SCALARS integrated_flux float\n")
    f.write("LOOKUP_TABLE default\n")
    for i in range(numCells):
        print i, data[i]
        f.write("{} ".format(data[i]))
    f.write("\n")

    f.close()

def saveFlowVtk(filename, coords, flux0, flux1):

    dims = coords[0].shape

    f = open(filename, 'w')

    f.write("# vtk DataFile Version 2.0\n")
    f.write("flow\n")
    f.write("ASCII\n")
    f.write("DATASET STRUCTURED_GRID\n")
    f.write("DIMENSIONS 1 {} {}\n".format(dims[1], dims[0])) # inverse order
    numPoints = dims[0] * dims[1]
    f.write("POINTS {} float\n".format(numPoints))
    for i in range(dims[0]):
        for j in range(dims[1]):
            f.write("{} {} 0.0\n".format(coords[0][i, j], coords[1][i, j]))
 
    numCells = (dims[0] - 1) * (dims[1] - 1)

    f.write("CELL_DATA {}\n".format(numCells))
    f.write("SCALARS velocity float 3\n")
    f.write("LOOKUP_TABLE default\n")

    # fluxes at cell centres
    # flux0 is the flux integrated along the 1st coordinate
    # flux1 is the flux integrated along the 2nd coordinate
    flux0Avg = 0.5*(flux0[:, :-1] + flux0[:, 1:])
    flux1Avg = 0.5*(flux1[:-1, :] + flux1[1:, :])

    xx, yy = coords[0], coords[1]

    # compute the mid face positions
    #
    # 0Lo is the (i,j), (i+1,j) face
    # 0Hi is the (i,j+1), (i+1,j+1) face
    # 1Lo is the (i,j), (i,j+1) face
    # 1Hi is the (i+1,j), (i+1,j+1) face
    xx0Lo = 0.5*(xx[:-1, :-1] + xx[1:, :-1])
    yy0Lo = 0.5*(yy[:-1, :-1] + yy[1:, :-1])

    xx0Hi = 0.5*(xx[:-1, 1:] + xx[1:, 1:])
    yy0Hi = 0.5*(yy[:-1, 1:] + yy[1:, 1:])

    xx1Lo = 0.5*(xx[:-1, :-1] + xx[:-1, 1:])
    yy1Lo = 0.5*(yy[:-1, :-1] + yy[:-1, 1:])

    xx1Hi = 0.5*(xx[1:, :-1] + xx[1:, 1:])
    yy1Hi = 0.5*(yy[1:, :-1] + yy[1:, 1:])

    # cell widths at the mid cell positions
    dx0 = xx1Hi - xx1Lo
    dy0 = yy1Hi - yy1Lo
    dx1 = xx0Hi - xx0Lo
    dy1 = yy0Hi - yy0Lo
    length0 = numpy.sqrt(dx0*dx0 + dy0*dy0)
    length1 = numpy.sqrt(dx1*dx1 + dy1*dy1)

    # flux densities
    ff0 = flux0Avg/length0
    ff1 = flux1Avg/length1

    # unit vectors perpendicular to each face
    u0x = + dy0/length0
    u0y = - dx0/length0
    u1x = + dy1/length1
    u1y = - dx1/length1

    # the velocity in the x, y directions
    vx = ff0*u0x + ff1*u1x
    vy = ff0*u0y + ff1*u1y

    for i in range(dims[0] - 1):
        for j in range(dims[1] - 1):
            f.write("{} {} 0.0\n".format(vx[i, j], vy[i, j]))

    f.close();


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

