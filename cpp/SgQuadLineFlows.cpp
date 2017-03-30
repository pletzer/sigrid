/**
 * A class that computes the projection of flux basis functions on a line in 2D
 */
 
#include "SgQuadLineFlows.h"
#include <iostream>

extern "C"
int SgQuadLineFlows_new(SgQuadLineFlows_type** self) {

 	*self = new SgQuadLineFlows_type();

 	return 0;
}
      
extern "C"                   
int SgQuadLineFlows_del(SgQuadLineFlows_type** self) {

 	delete *self;
 	
 	return 0;
}

extern "C"
int SgQuadLineFlows_setQuadPoints(SgQuadLineFlows_type** self,
	                                  const double* quadCoords) {
 	(*self)->setQuadPoints(quadCoords);
 	return 0;
}

extern "C"
int SgQuadLineFlows_setLinePoints(SgQuadLineFlows_type** self,
	                                  const double* lineCoords) {
 	(*self)->setLinePoints(lineCoords);
 	return 0;
}

extern "C"
int SgQuadLineFlows_computeProjections(SgQuadLineFlows_type** self) {
	(*self)->computeProjections();
 	return 0;
}
