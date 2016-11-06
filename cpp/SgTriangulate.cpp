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
 
 	// store the points
 	(*self)->points.resize(2 * numPoints);

 	// set the indices before the sort
 	(*self)->sortedInds.resize(numPoints);
 	for (int i = 0; i < numPoints; ++i) {
 		(*self)->sortedInds[i] = i;
        (*self)->points[2*i + 0] = points[i][0];
        (*self)->points[2*i + 1] = points[i][1];
 	}
 	// sort the point indices by increasing distance from the centre of gravity
 	std::sort((*self)->sortedInds.begin(), (*self)->sortedInds.end(),
 		      sortFunc);

    if (numPoints < 3) return 0;

    // create the first triangle
    size_t i0 = (*self)->sortedInds[0];
    size_t i1 = (*self)->sortedInds[1];
    size_t i2 = (*self)->sortedInds[2];
    double area = (*self)->getParallelipipedArea(i0, i1, i2);
    if (std::abs(area) < (*self)->eps) {
        // degenerate triangle NEED TO FIX
        std::cerr << "*** degenerate first triangle\n";
    }
    else if (area < 0) {
        // change the ordering of the indices
        (*self)->makeCounterClockwise(&(*self)->sortedInds[0]);
    }
    i0 = (*self)->sortedInds[0];
    i1 = (*self)->sortedInds[1];
    i2 = (*self)->sortedInds[2];
    size_t edge[] = {i0, i1};
    (*self)->boundaryEdges.insert(std::vector<size_t>(edge, edge + 2));
    edge[0] = i1; edge[1] = i2;
    (*self)->boundaryEdges.insert(std::vector<size_t>(edge, edge + 2));
    edge[0] = i2; edge[1] = i0;

    // add remaining points
    for (int i = 3; i < numPoints; ++i) {
        (*self)->addPoint(i);
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
		*area += 0.5 * (*self)->getParallelipipedArea((*it)[0], (*it)[1], (*it)[2]);
	}

 	return 0;
}
