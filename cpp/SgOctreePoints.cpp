/**
 * A class that finds the location of a point in index space
 */
 
#include "SgOctreePoints.h"
#include <cmath>
#include <iostream>

extern "C"
int SgOctreePoints_new(SgOctreePoints_type** self,
                          int nitermax, double tolpos) {
    *self = new SgOctreePoints_type(nitermax, tolpos);
    return 0;
}
      
extern "C"                   
int SgOctreePoints_del(SgOctreePoints_type** self) {
    delete *self;
    return 0;
}

extern "C"
int SgOctreePoints_setGrid(SgOctreePoints_type** self, 
                              int ndims, const int dims[],
                              const int periodicity[], 
                              const double** coords) {
    (*self)->setGrid(ndims, dims, periodicity, coords);
    return 0;
}

extern "C"
int SgOctreePoints_getPosition(SgOctreePoints_type** self,
                                  double pos[]) {
    std::vector<double> position = (*self)->getPosition();
    for (size_t i = 0; i < position.size(); ++i) {
        pos[i] = position[i];
    }
    return 0;
}

extern "C"
int SgOctreePoints_getError(SgOctreePoints_type** self,
                               double* error) {
    *error = (*self)->getError();
    return 0;
}

extern "C"
int SgOctreePoints_reset(SgOctreePoints_type** self, 
                            const double dIndices[],
                            const double targetPoint[]) {
    (*self)->reset(dIndices, targetPoint);
    return 0;
}

extern "C"
int SgOctreePoints_next(SgOctreePoints_type** self) {
    return (*self)->next();
}

extern "C"
int SgOctreePoints_getIndices(SgOctreePoints_type** self,
                                 double dIndices[]) {
    std::vector<double> dInds = (*self)->getIndices();
    for (size_t i = 0; i < dInds.size(); ++i) {
        dIndices[i] = dInds[i];
    }
    return 0;
}

extern "C"
int SgOctreePoints_getErrorHistory(SgOctreePoints_type** self,
                                      int* niter, double** errorHistory) {
  (*self)->getErrorHistory(niter, errorHistory);
  return 0;
}
