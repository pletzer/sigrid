/**
 * A class that triangulates a set of points
 */
 
#include "SgTriangulate.h"
#include <iostream>

extern "C"
int SgTriangulate_new(SgTriangulate_type** self, int numPoints, const double points[]) {

    *self = new SgTriangulate_type(numPoints, points);

    return 0;
}
      
extern "C"                   
int SgTriangulate_del(SgTriangulate_type** self) {

     delete *self;

     return 0;
}


extern "C"
int SgTriangulate_getConvexHullArea(SgTriangulate_type** self,
                                    double* area) {
    *area = (*self)->getConvexHullArea();
    return 0;
}
