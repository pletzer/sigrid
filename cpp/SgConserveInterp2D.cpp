/**
 * A class that finds the location of a point in index space
 */
 
#include "SgConserveInterp2D.h"
#include <cmath>
#include <iostream>

extern "C"
int SgConserveInterp2D_new(SgConserveInterp2D_type** self) {
 	*self = new SgConserveInterp2D_type();
 	// more stuff

 	return 0;
}
      
extern "C"                   
int SgConserveInterp2D_del(SgConserveInterp2D_type** self) {

 	delete *self;

 	return 0;
}
