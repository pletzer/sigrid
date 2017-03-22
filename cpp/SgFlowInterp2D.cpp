/**
 * A class that computes the interpolation weights for a face centered field in 2D
 */
 
#include "SgFlowInterp2D.h"
#include <cmath>
#include <iostream>

extern "C"
int SgFlowInterp2D_new(SgFlowInterp2D_type** self) {
 	*self = new SgFlowInterp2D_type();
	return 0;
}
      
extern "C"                   
int SgFlowInterp2D_del(SgFlowInterp2D_type** self) {
 	delete *self;
 	return 0;
}

extern "C"
int SgFlowInterp2D_setDstGrid(SgFlowInterp2D_type** self, 
 	                              const int dims[],
 	                              const double** coords) {
	(*self)->setDstGrid(dims, coords);
	return 0;
}

extern "C"
int SgFlowInterp2D_setSrcGrid(SgFlowInterp2D_type** self, 
 	                              const int dims[], 
 	                              const double** coords) {
	(*self)->setSrcGrid(dims, coords);
	return 0;
}

extern "C"
int SgFlowInterp2D_computeWeights(SgFlowInterp2D_type** self) {
	(*self)->computeWeights();
	return 0;
}

extern "C"
int SgFlowInterp2D_apply(SgFlowInterp2D_type** self,
 	                     const double* srcData[], double dstData[]) {
	(*self)->apply(srcData, dstData);
	return 0;
}
