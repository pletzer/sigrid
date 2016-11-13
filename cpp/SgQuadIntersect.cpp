/**
 * A class that finds all intersection points between two quads
 */
 
#include "SgQuadIntersect.h"
#include <iostream>

extern "C"
int SgQuadIntersect_new(SgQuadIntersect_type** self) {

 	*self = new SgQuadIntersect_type();

 	// tolerance for floating point comparisons
 	(*self)->tol = 1.e-12;

 	(*self)->slvr = NULL;
 	SgLinearSolve_new(&(*self)->slvr, 2, 2);

 	(*self)->quad1Coords = NULL;
 	(*self)->quad2Coords = NULL;

 	return 0;
}
      
extern "C"                   
int SgQuadIntersect_del(SgQuadIntersect_type** self) {

 	if ((*self)->slvr) SgLinearSolve_del(&(*self)->slvr);
 	delete *self;

 	return 0;
}

extern "C"
int SgQuadIntersect_setQuadPoints(SgQuadIntersect_type** self,
	                              const double* quad1Coords, const double* quad2Coords) {
 	(*self)->setQuadPoints(quad1Coords, quad2Coords);
 	return 0;
}


extern "C"
int SgQuadIntersect_getIntersectPoints(SgQuadIntersect_type** self,
 	                                   int* numPoints, double** points) {
	(*self)->collectIntersectPoints(numPoints, points);
 	return 0;
}
