/**
 * A class that finds all intersection points between a quad and a line
 */
 
#include "SgQuadLineIntersect.h"
#include <iostream>

extern "C"
int SgQuadLineIntersect_new(SgQuadLineIntersect_type** self) {

 	*self = new SgQuadLineIntersect_type();

 	return 0;
}
      
extern "C"                   
int SgQuadLineIntersect_del(SgQuadLineIntersect_type** self) {

 	if ((*self)->slvr) delete (*self)->slvr;
 	delete *self;

 	return 0;
}

extern "C"
int SgQuadLineIntersect_setQuadPoints(SgQuadLineIntersect_type** self,
	                                  const double* quadCoords) {
 	(*self)->setQuadPoints(quadCoords);
 	return 0;
}

extern "C"
int SgQuadLineIntersect_setLinePoints(SgQuadLineIntersect_type** self,
	                                  const double* lineCoords) {
 	(*self)->setLinePoints(lineCoords);
 	return 0;
}

extern "C"
int SgQuadLineIntersect_collectIntersectPoints(SgQuadLineIntersect_type** self,
 	                                           int* numPoints, double** points) {
	(*self)->collectIntersectPoints(numPoints, points);
 	return 0;
}
