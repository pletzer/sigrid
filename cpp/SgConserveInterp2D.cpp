/**
 * A class that computes the interpolation weights for a cell centered field
 */
 
#include "SgConserveInterp2D.h"
#include <cmath>
#include <iostream>

extern "C"
int SgConserveInterp2D_new(SgConserveInterp2D_type** self) {
 	*self = new SgConserveInterp2D_type();

 	// reset
 	SgConserveInterp2D_reset(self);
 	return 0;
}
      
extern "C"                   
int SgConserveInterp2D_del(SgConserveInterp2D_type** self) {
 	delete *self;
 	return 0;
}

extern "C"
int SgConserveInterp2D_setDstGrid(SgConserveInterp2D_type** self, 
 	                              const int dims[], const double** coords) {
	(*self)->setDstGrid(dims, coords);
	return 0;
}

extern "C"
int SgConserveInterp2D_setSrcGrid(SgConserveInterp2D_type** self, 
 	                              const int dims[], const double** coords) {
	(*self)->setSrcGrid(dims, coords);
	return 0;
}

extern "C"
int SgConserveInterp2D_computeWeights(SgConserveInterp2D_type** self) {
	(*self)->computeWeights();
	return 0;
}

extern "C"
int SgConserveInterp2D_apply(SgConserveInterp2D_type** self,
 	                         const double srcData[], double dstData[]) {
	(*self)->apply(srcData, dstData);
	return 0;
}

extern "C"
int SgConserveInterp2D_reset(SgConserveInterp2D_type** self) {
	(*self)->weightIt = (*self)->weights.begin();
	(*self)->weightIt2 = (*self)->weightIt->second.begin();
	return 0;
}

extern "C"
int SgConserveInterp2D_next(SgConserveInterp2D_type** self) {

	if ((*self)->weightIt2 != (*self)->weightIt->second.end()) {
		// increment the src index counter
		(*self)->weightIt2++;
		return 0;
	}
	else if ((*self)->weightIt != (*self)->weights.end()) {
		// increment the dst index counter
		(*self)->weightIt++;
		(*self)->weightIt2 = (*self)->weightIt->second.begin();
		return 0;
	}

	// reached end (error)
	return 1;
}

extern "C"
int SgConserveInterp2D_get(SgConserveInterp2D_type** self,
                           int* srcIndx, int* dstIndx, double* weight) {
	*dstIndx = (*self)->weightIt->first;
	*srcIndx = (*self)->weightIt2->first;
	*weight = (*self)->weightIt2->second;
	return 0;
}
