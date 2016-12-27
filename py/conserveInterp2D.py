import ctypes
import numpy
import sigrid

NDIMS = 2

class ConserveInterp2D:

	def __init__(self,):
		"""
		Constructor
		"""

		self.handle = ctypes.c_void_p(0)
		sigrid.so.SgConserveInterp2D_new(ctypes.byref(self.handle))
		self.hasSrcCoords = False
		self.hasInterpWeights = False
		self.dstData = None

	def __del__(self):
		"""
		Destructor
		"""
		sigrid.so.SgConserveInterp2D_del(ctypes.byref(self.handle))

	def setDstGrid(self, xcoords, ycoords):
		"""
		Set destination grid
		@param xcoords array of point coordinates
		@param ycoords array of point coordinates
		"""
		n0, n1 = xcoords.shape
		dims = (ctypes.c_int * NDIMS)(n0, n1)
		coords = (ctypes.POINTER(ctypes.c_double) * NDIMS)(xcoords.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
		                                                   ycoords.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))

		sigrid.so.SgConserveInterp2D_setDstGrid(ctypes.byref(self.handle), dims, coords)
		self.dstData = numpy.zeros((n0, n1), numpy.float64)

	def setSrcGrid(self, periodicity, xcoords, ycoords):
		"""
		Set source grid
		@param periodicity array of True (periodic)/False (not periodic) values
		@param xcoords array of point coordinates
		@param ycoords array of point coordinates 
		"""
		n0, n1 = xcoords.shape
		dims = (ctypes.c_int * NDIMS)(n0, n1)
		p0, p1 = int(periodicity[0]), int(periodicity[1])
		periods = (ctypes.c_int * NDIMS)(p0, p1)
		coords = (ctypes.POINTER(ctypes.c_double) * NDIMS)(xcoords.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
		                                                   ycoords.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
		sigrid.so.SgConserveInterp2D_setSrcGrid(ctypes.byref(self.handle), dims, periods, coords)
		self.hasSrcCoords = True

	def apply(self, srcData):
		"""
		Apply the interpolation weights to the source grid dstData
		@param srcData cell centred data on the source grid
		@return interpolated data on the destination grid
		"""
		if not self.hasInterpWeights:
			raise RuntimeError('ERROR: must call "computeWeights" prior to invoking "apply"')
		self.dstData *= 0 # initialize to zero
		sigrid.so.SgConserveInterp2D_apply(ctypes.byref(self.handle),
			                               srcData.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
			                               self.dstData.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
		return self.dstData

	def computeWeights(self):
		"""
		Compute the interpolation weights
		"""
		if not self.hasSrcCoords:
			raise RuntimeError('ERROR: must call "setSrcGrid" prior to invoking "computeWeights"')
		if self.dstData is None:
			raise RuntimeError('ERROR: must call "setDstGrid" prior to invoking "computeWeights"')
		sigrid.so.SgConserveInterp2D_computeWeights(ctypes.byref(self.handle))
		self.hasInterpWeights = True



