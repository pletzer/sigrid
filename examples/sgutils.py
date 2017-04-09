import numpy

def createRectilinearGrid(dims, xmin, xmax):
	x0 = numpy.linspace(xmin[0], xmax[0], dims[0])
	x1 = numpy.linspace(xmin[1], xmax[1], dims[1])
	xx0 = numpy.outer(x0, numpy.ones((dims[1],), numpy.float64))
	xx1 = numpy.outer(numpy.ones((dims[1],), numpy.float64), x1)
	return (xx0, xx1)
