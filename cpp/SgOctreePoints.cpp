/**
 * A class that finds the location of a point in index space
 */
 
#include "SgOctreePoints.h"
#include <cmath>
#include <iostream>

extern "C"
int SgOctreePoints_new(SgOctreePoints_type** self, 
                       int numLevels, int ndims, int npoints, const double* points) {
    std::vector<double> pts(points, points + npoints*ndims);
    *self = new SgOctreePoints_type((size_t) numLevels, (size_t) ndims, pts);
    return 0;
}
      
extern "C"                   
int SgOctreePoints_del(SgOctreePoints_type** self) {
    delete *self;
    return 0;
}
