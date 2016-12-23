import ctypes

class ConserveInterp2D:

	def __init__(self,):
		"""
		Constructor
		"""
		pass

	def setDstGrid(self, xcoords, ycoords):
		"""
		Set destination grid
		@param xcoords array of point coordinates
		@param ycoords array of point coordinates
		"""
		pass

	def setSrcGrid(self, periodicity, xcoords, ycoords):
		"""
		Set source grid
		@param periodicity array of True/False values, set to True if periodic
		@param xcoords array of point coordinates
		@param ycoords array of point coordinates 
		"""
		pass

	def apply(self, srcData):
		"""
		Apply the interpolation weights to the source grid dstData
		@param srcData cell centred data on the source grid
		@return interpolated data on the destination grid
		"""
		# TO DO
		return dstData

	def computeWeights(self):
		"""
		Compute the interpolation weights
		"""
		pass



