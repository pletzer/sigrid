/**
 * A class that finds all intersection points between two quads
 */
 
#include "SgStructGrid2DIntersect.h"
#include <iostream>

extern "C"
int SgStructGrid2DIntersect_new(SgStructGrid2DIntersect_type** self) {

 	*self = new SgStructGrid2DIntersect_type();

 	return 0;
}
      
extern "C"                   
int SgStructGrid2DIntersect_del(SgStructGrid2DIntersect_type** self) {

 	if ((*self)->slvr) delete (*self)->slvr;
 	delete *self;

 	return 0;
}

extern "C"
int SgStructGrid2DIntersect_setDstGrid(SgStructGrid2DIntersect_type** self,
	                                   const int dims[],
                                       const double** coords) {
 	(*self)->setDstGrid(dims, coords);
 	return 0;
}


extern "C"
int SgStructGrid2DIntersect_getSrcGrid(SgStructGrid2DIntersect_type** self,
	                                   const int dims[], const int periodicity[],
                                       const double** coords) {
	(*self)->setSrcGrid(dims, periodicity, coords);
 	return 0;
}
