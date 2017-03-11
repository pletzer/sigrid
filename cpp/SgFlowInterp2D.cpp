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
 	                              const int periodicity[],
 	                              const double** coords) {
	(*self)->setSrcGrid(dims, periodicity, coords);
	return 0;
}

extern "C"
int SgFlowInterp2D_computeWeights(SgFlowInterp2D_type** self) {
	(*self)->computeWeights();
	return 0;
}

extern "C"
int SgFlowInterp2D_apply(SgFlowInterp2D_type** self,
 	                         const double srcData[], double dstData[]) {
	(*self)->apply(srcData, dstData);
	return 0;
}

extern "C"
int SgFlowInterp2D_reset(SgFlowInterp2D_type** self) {
	(*self)->reset();
	return 0;
}

extern "C"
int SgFlowInterp2D_next(SgFlowInterp2D_type** self) {
	return (*self)->next();
}

extern "C"
int SgFlowInterp2D_get(SgFlowInterp2D_type** self,
                           int* srcIndx, int* dstIndx, double* weight) {
	(*self)->get(srcIndx, dstIndx, weight);
	return 0;
}
