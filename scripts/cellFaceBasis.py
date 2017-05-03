"""
Cell face basis functions
"""

import numpy
import argparse
from numpy import pi, sin, cos

class FaceBasis(dims):

	def __init__(self, dims, xExpr, yExpr):
		
		# build the mesh
		xsi1d = numpy.linspace(0., 1., dims[0])
		est1d = numpy.linspace(0., 1., dims[1])
		self.xsi = numpy.outer(xsi1d, numpy.ones(dims[1], numpy.float64))
		self.eta = numpy.outer(numpy.ones(dims[0], numpy.float64), eta1d)

		# apply coordinate transformations
		self.x = eval(xExpr)
		self.y = eval(yExpr)

		# take finite differences



	def saveVtk(self, filename, side):
		pass
