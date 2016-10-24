/**
 * A class that finds the location of a point in index space
 */
 
#include "SgFindPointInCell.h"
#include <cmath>
#include <iostream>

extern "C"
int SgFindPointInCell_new(SgFindPointInCell_type** self,
                           int nitermax, double tolpos) {
 	*self = new SgFindPointInCell_type();
 	(*self)->tolpos = tolpos;
 	(*self)->nitermax = nitermax;
 	(*self)->slvr = NULL;
 	(*self)->iter = 0;

 	return 0;
}
      
extern "C"                   
int SgFindPointInCell_del(SgFindPointInCell_type** self) {

 	if ((*self)->slvr) SgLinearSolve_del(&(*self)->slvr);
 	delete *self;

 	return 0;
}

extern "C"
int SgFindPointInCell_setGrid(SgFindPointInCell_type** self, 
 	                           int ndims, const int dims[],
 	                           const double** coords) {

 	if (ndims <= 0) return 0;

 	(*self)->coords.resize(ndims);
 	(*self)->prodDims.resize(ndims);
 	(*self)->nodeProdDims.resize(ndims);
 	(*self)->dims.resize(ndims);
 	(*self)->jacMatrix.resize(ndims * ndims);
 	(*self)->rhs.resize(ndims);
 	(*self)->dIndices.resize(ndims);
 	(*self)->targetPoint.resize(ndims);

 	// must have at least one dimension
 	(*self)->prodDims[ndims - 1] = 1;
 	(*self)->nodeProdDims[ndims - 1] = 1;
 	for (int i = ndims - 2; i >= 0; --i) {
 		(*self)->prodDims[i] = (*self)->prodDims[i + 1] * dims[i + 1];
 		(*self)->nodeProdDims[i] = (*self)->nodeProdDims[i + 1] * 2;
 	}

 	// total number of nodes
 	int ntot = 1;
 	for (int i = 0; i < ndims; ++i) {
 		ntot *= dims[i];
 		(*self)->dims[i] = dims[i];
 	}

 	// set the coordinates
 	for (int i = 0; i < ndims; ++i) {
 		(*self)->coords[i].resize(ntot);
 		for (int j = 0; j < ntot; ++j) {
 			(*self)->coords[i][j] = coords[i][j];
 		}
 	}

 	// create a solver
 	SgLinearSolve_new(&(*self)->slvr, ndims, ndims);

 	return 0;
}

extern "C"
int SgFindPointInCell_getPosition(SgFindPointInCell_type** self,
 	                              double pos[]) {

	size_t ndims = (*self)->dims.size();
	for (size_t i = 0; i < ndims; ++i) {
		pos[i] = (*self)->interp((*self)->dIndices, (*self)->coords[i]);
	}

	return 0;
}


extern "C"
int SgFindPointInCell_reset(SgFindPointInCell_type** self, 
	                        const double dIndices[],
	                        const double targetPoint[]) {

	size_t ndims = (*self)->dims.size();
	for (size_t i = 0; i < ndims; ++i) {
		(*self)->dIndices[i] = dIndices[0];
		(*self)->targetPoint[i] = targetPoint[i];
	}
	(*self)->iter = 0;

	return 0;
}

extern "C"
int SgFindPointInCell_next(SgFindPointInCell_type** self) {

	(*self)->computeJacobianAndRHS();
	SgLinearSolve_setMatrix(&(*self)->slvr, &(*self)->jacMatrix[0]);
	SgLinearSolve_setRightHandSide(&(*self)->slvr, &(*self)->rhs[0]);
	SgLinearSolve_solve(&(*self)->slvr);

	double* sol;
	SgLinearSolve_getSolution(&(*self)->slvr, &sol);

	size_t ndims = (*self)->dims.size();
	// update the indices
	for (size_t i = 0; i < ndims; ++i) {
		(*self)->dIndices[i] += sol[i];
	}

	// check if the next iterator is still valid
	(*self)->iter++;
	std::vector<double> pos(ndims);
	SgFindPointInCell_getPosition(self, &pos[0]);

	// use Eulerian distance as error measure
	double posError = 0;
	for (size_t i = 0; i < ndims; ++i) {
		double dp = pos[i] - (*self)->targetPoint[i];
		posError += dp * dp;
	}
	posError = sqrt(posError);

	if ((*self)->iter >= (*self)->nitermax) {
		// reached max number of iterations
		return -1; 
	}
	if (posError < (*self)->tolpos) {
		// done!
		return 1;
	}

	// has not yet converged
	return 0;
}

extern "C"
int SgFindPointInCell_getIndices(SgFindPointInCell_type** self,
 	                             double dIndices[]) {
 	size_t ndims = (*self)->coords.size();
 	for (size_t i = 0; i < ndims; ++i) {
 		dIndices[i] = (*self)->dIndices[i];
 	}

 	return 0;
}
