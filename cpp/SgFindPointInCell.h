/**
 * A class that finds the location of a point in the index space of a structured grid
 */
 
#ifndef SG_FIND_POINT_IN_CELL_H
#define SG_FIND_POINT_IN_CELL_H
 
#include <vector>
#include <cstdio> // size_t
#include <cmath>
#include "SgLinearSolve.h"
 
struct SgFindPointInCell_type {

 	// linear solver
 	SgLinearSolve_type* slvr;

 	// the curvilinear coordinate points as ndims flat arrays
 	std::vector<std::vector<double> > coords;

 	// used to get the flat index
 	std::vector<size_t> prodDims;

 	// the number of grid points for each dimension
 	std::vector<size_t> dims;

 	// the target point
 	std::vector<double> targetPoint;

 	// max number of Newton iterations
 	int nitermax;

 	// tolerance in coordinate space
 	double tolpos;

 	// current iteration
 	int iter;

 	// current index location
 	std::vector<double> dIndices;

 	// current Jacobian matrix
 	std::vector<double> jacMatrix;

 	// target point minus current position (initially)
 	std::vector<double> rhs;

void getWeightsAndFlatIndices(const std::vector<double>& dInds,
 		                      std::vector<double>& weights, 
		                      std::vector<size_t>& flatInds) {

 	size_t ndims = this->dims.size();
 	size_t nnodes = weights.size();
 	for (size_t j = 0; j < nnodes; ++j) {
 		flatInds[j] = 0;
 		weights[j] = 1;
 		for (size_t i = 0; i < ndims; ++i) {
 			int loCornerIndx = (int) floor(dInds[i]);
 			int indx = loCornerIndx + (j / this->prodDims[i] % 2);
 			flatInds[j] += (size_t) this->prodDims[i] * indx;
 			weights[j] *= (double) indx - dInds[i];
 		}
 	}
}

double interp(const std::vector<double>& dInds,
 		      const std::vector<double>& nodalField) {

 	double res = 0;
 	size_t ndims = this->dims.size();
 	size_t nnodes = pow(2, ndims);
 	std::vector<double> weights(nnodes);
 	std::vector<size_t> flatInds(nnodes);
 	this->getWeightsAndFlatIndices(dInds, weights, flatInds);
 	for (size_t j = 0; j < nnodes; ++j) {
 		size_t bindx = flatInds[j];
 		res += weights[j] * nodalField[bindx];
 	}
 	return res;
}

void computeJacobianAndRHS() {

 	size_t ndims = this->dims.size();
 	std::vector<double> dInds(ndims);

 	// iterate over the physical space dimensions
 	for (size_t i = 0; i < ndims; ++i) {

 		this->rhs[i] = this->targetPoint[i] - this->interp(this->dIndices, this->coords[i]);

 		for (size_t j = 0; j < ndims; ++j) {
 			// mid cell
 			dInds[j] = (int) floor(this->dIndices[j]) + 0.5;
 		}

 		// iterate ove the index space dimensions
 		for (size_t j = 0; j < ndims; ++j) {

 			// high end of the cell
 			dInds[j] += 0.5;
 			double xHi = this->interp(dInds, this->coords[i]);

 			// low end of the cell
 			dInds[j] -= 1.0;
 			double xLo = this->interp(dInds, this->coords[i]);

 			// reset to mid cell
 			dInds[j] += 0.5;

			size_t k = i + ndims * j; // Fortran ordering

			// average difference of the i-th coordinate anlong the j-th topo direction
 			this->jacMatrix[k] = xHi - xLo; 			
 		}
 	}
}

};
 
#ifdef __cplusplus
extern "C" {
#endif

 int SgFindPointInCell_new(SgFindPointInCell_type** self,
                       int nitermax, double tolpos);
                       
 int SgFindPointInCell_del(SgFindPointInCell_type** self);

 int SgFindPointInCell_setGrid(SgFindPointInCell_type** self, 
 	                           int ndims, const int dims[], 
 	                           const double** coords);

 int SgFindPointInCell_getPosition(SgFindPointInCell_type** self, 
 	                               double pos[]);

 int SgFindPointInCell_reset(SgFindPointInCell_type** self, 
	                         const double dIndices[]);

 int SgFindPointInCell_next(SgFindPointInCell_type** self);

 int SgFindPointInCell_getIndices(SgFindPointInCell_type** self,
 	                              double dIndices[]);
 

#ifdef __cplusplus
}
#endif

#endif // SG_FIND_POINT_IN_CELL_H