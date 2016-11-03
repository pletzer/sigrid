/**
 * A class that triangulates a set of points
 */
 
#include "SgTriangulate.h"
#include <iostream>
#include <algorithm>

extern "C"
int SgTriangulate_new(SgTriangulate_type** self, int numPoints, const double** points) {

 	*self = new SgTriangulate_type();

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
 	}
 	// sort the point indices by increasing distance from the centre of gravity
 	std::sort((*self)->sortedInds.begin(), (*self)->sortedInds.end(),
 		      sortFunc);


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
	size_t numTri = (*self)->triIndices.size() / 3;
	for (size_t i = 0; i < numTri; ++i) {
		size_t ia = (*self)->triIndices[3*i + 0];
		size_t ib = (*self)->triIndices[3*i + 1];
		size_t ic = (*self)->triIndices[3*i + 2];
		*area += 0.5 * (*self)->getParallelipipedArea(ia, ib, ic);
	}

 	return 0;
}
