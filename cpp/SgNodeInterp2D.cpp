/**
 * A class that computes the interpolation weights for a cell centered field
 */
 
#include "SgNodeInterp2D.h"
#include <cmath>
#include <iostream>

extern "C"
int SgNodeInterp2D_new(SgNodeInterp2D_type** self, int nitermax, double tolpos) {
     *self = new SgNodeInterp2D_type(nitermax, tolpos);
    return 0;
}
      
extern "C"                   
int SgNodeInterp2D_del(SgNodeInterp2D_type** self) {
     delete *self;
     return 0;
}

extern "C"
int SgNodeInterp2D_setDstGrid(SgNodeInterp2D_type** self, 
                                   const int dims[],
                                   const double** coords) {
    (*self)->setDstGrid(dims, coords);
    return 0;
}

extern "C"
int SgNodeInterp2D_setSrcGrid(SgNodeInterp2D_type** self, 
                              const int dims[], 
                              const int periodicity[],
                              const double** coords) {
    (*self)->setSrcGrid(dims, periodicity, coords);
    return 0;
}

extern "C"
int SgNodeInterp2D_computeWeights(SgNodeInterp2D_type** self) {
    (*self)->computeWeights();
    return 0;
}

extern "C"
int SgNodeInterp2D_apply(SgNodeInterp2D_type** self,
                              const double srcData[], double dstData[]) {
    (*self)->apply(srcData, dstData);
    return 0;
}
