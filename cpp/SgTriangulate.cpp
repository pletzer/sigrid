/**
 * A class that triangulates a set of points
 */
 
#include "SgTriangulate.h"
#include <iostream>
#include <algorithm>

extern "C"
int SgTriangulate_new(SgTriangulate_type** self, int numPoints, const double** points) {

    *self = new SgTriangulate_type();

    (*self)->NDIMS = 2;

    // tolerance for floating point comparisons
    (*self)->eps = 1.e-12;

    (*self)->points.resize(2 * numPoints);

    SgSortByDistanceSquareFunctor sortFunc(numPoints, points);
 
    // set the indices before the sort
    std::vector<size_t> sortedInds(numPoints);
    for (int i = 0; i < numPoints; ++i) {
        sortedInds[i] = i;
    }
    // sort the point indices by increasing distance from the centre of gravity
    std::sort(sortedInds.begin(), sortedInds.end(), sortFunc);

    // set the points in increasing distance from the centre of gravity
    (*self)->points.resize(2 * numPoints);
    for (size_t i = 0; i < sortedInds.size(); ++i) {
        size_t indx = sortedInds[i];
        for (size_t j = 0; j < (*self)->NDIMS; ++j) {
            (*self)->points[2*i + j] = points[indx][j];
        }        
    }

    // remove degenerate points
    (*self)->removeDegenerateSegments();

    // create the first, non-degenerate triangle
    size_t lastIndexOfFirstTriangle = (*self)->makeFirstTriangle();

    // add the remaining points
    for (int ip = lastIndexOfFirstTriangle + 1; ip < numPoints; ++ip) {
        (*self)->addPoint(ip);
    }

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

    *area = 0;
    for (std::set< std::vector<size_t> >::const_iterator it = (*self)->triIndices.begin();
         it != (*self)->triIndices.end(); ++it) {
        *area += (*self)->getParallelipipedArea((*it)[0], (*it)[1], (*it)[2]);
    }
    *area *= 0.5;

     return 0;
}
