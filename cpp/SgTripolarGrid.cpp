/**
 * A class that constructs a tripolar grid
 */
 
#include "SgTripolarGrid.h"
#include <cmath>
#include <iostream>

extern "C"
int SgTripolarGrid_new(SgTripolarGrid_type** self,
                          const int dims[], int capLatIndex) {
    *self = new SgTripolarGrid_type(dims, capLatIndex);
    return 0;
}
      
extern "C"                   
int SgTripolarGrid_del(SgTripolarGrid_type** self) {
    delete *self;
    return 0;
}

extern "C"
int SgTripolarGrid_getGrid(SgTripolarGrid_type** self,
                           int dims[], double*** coords) {
    (*self)->getGrid(dims, coords);
    return 0;
}

