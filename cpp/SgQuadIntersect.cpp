/**
 * A class that finds all intersection points between two quads
 */
 
#include "SgQuadIntersect.h"
#include <iostream>

extern "C"
int SgQuadIntersect_new(SgQuadIntersect_type** self) {

 	*self = new SgQuadIntersect_type();

 	return 0;
}
      
extern "C"                   
int SgQuadIntersect_del(SgQuadIntersect_type** self) {

 	if ((*self)->slvr) delete (*self)->slvr;
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
