/**
 * A class that finds the location of a point in index space
 */
 
 #include "SgFindPointInCell.h"
 #include <cmath>

 void SgGetWeightsAndFlatIndices(SgFindPointInCell_type** self,
		                         std::vector<double>& weights, 
		                         std::vector<size_t>& flatInds) {
 	size_t ndims = (*self)->dims.size();
 	size_t nnodes = weights.size();
 	size_t bif;
 	for (size_t j = 0; j < nnodes; ++j) {
 		for (size_t i = 0; i < ndims; ++i) {
 			int loCornerIndx = (int) floor((*self)->dIndices[i]);
 			size_t indx = loCornerIndx + (j / (*self)->prodDims[i] % 2);
 			bif += (*self)->prodDims[i] * indx;
 		}
 		flatInds[j] = bif;
 	}

 }

 void SgComputeJacobianAndRHS(SgFindPointInCell_type** self, 
 	                          std::vector<double>& jacMatrix,
 	                          std::vector<double>& rhs) {
 	
 }

 int SgFindPointInCell_new(SgFindPointInCell_type** self,
                           int nitermax, double tolpos) {
 	*self = new SgFindPointInCell_type();
 	(*self)->tolpos = tolpos;
 	(*self)->nitermax = nitermax;
 	(*self)->slvr = NULL;
 	return 0;
 }
                       
 int SgFindPointInCell_del(SgFindPointInCell_type** self) {
 	if ((*self)->slvr) SgLinearSolve_del(&(*self)->slvr);
 	delete *self;
 }

 int SgFindPointInCell_setGrid(SgFindPointInCell_type** self, 
 	                           int ndims, const int dims[],
 	                           const double** coords) {
 	(*self)->coords.resize(ndims);
 	(*self)->prodDims.resize(ndims);
 	(*self)->dims.resize(ndims);

 	// must have at least one dimension
 	(*self)->prodDims[ndims - 1] = 1;
 	for (int i = ndims - 2; i >= 0; --i) {
 		(*self)->prodDims[i] = (*self)->prodDims[i + 1] * 2;
 	}

 	int ntot = 1;
 	for (int i = 0; i < ndims; ++i) {
 		ntot *= dims[i];
 		(*self)->dims[i] = dims[i];
 	}

 	for (int i = 0; i < ndims; ++i) {
 		(*self)->coords[i].resize(ntot);
 		for (int j = 0; j < ntot; ++j) {
 			(*self)->coords[i][j] = coords[i][j];
 		}
 	}

 	SgLinearSolve_new(&(*self)->slvr, ndims, ndims);
 }

int SgFindPointInCell_getPosition(SgFindPointInCell_type** self,
 	                              double pos[]) {

	size_t ndims = (*self)->coords.size();
	size_t nnodes = pow(2, ndims);
	std::vector<double> weights(nnodes);
	std::vector<size_t> flatInds(nnodes);
	SgGetWeightsAndFlatIndices(self,
		                       weights, flatInds);
	for (size_t i = 0; i < ndims; ++i) {
		pos[i] = 0;
		for (size_t j = 0; j < nnodes; ++j) {
			pos[i] += weights[j] * (*self)->coords[i][flatInds[j]];
		}
	}
	return 0;
}

int SgFindPointInCell_next(SgFindPointInCell_type** self) {
	std::vector<double> jacMatrix((*self)->ndims * (*self)->ndims);
	std::vector<double> rhs((*self)->ndims);
	std::vector<double> sol((*self)->ndims);
	std::vector<double> pos((*self)->ndims);
	SgComputeJacobianAndRHS(self, jacMatrix, rhs);
	SgLinearSolve_setMatrix(&(*self)->slvr, &jacMatrix[0]);
	SgLinearSolve_setRightHandSide(&(*self)->slvr, &rhs[0]);
	SgLinearSolve_solve(&(*self)->slvr);
	SgLinearSolve_getSolution(&(*self)->slvr, &sol[0]);
	for (size_t i = 0; i < (*self)->ndims; ++i) {
		(*self)->dIndices[i] += sol[i];
	}
	(*self)->iter++;
	SgFindPointInCell_getPosition(self, &pos[0]);
	double posError = 0;
	for (size_t i = 0; i < (*self)->ndims; ++i) {
		double dp = pos[i] - (*self)->targetPoint[i];
		posError += dp * dp;
	}
	posError = sqrt(posError);
	if ((*self)->iter >= (*self)->nitermax || posError >= (*self)->tolpos) {
		return 1;
	}
	return 0;
}

int SgFindPointInCell_setIndices(SgFindPointInCell_type** self,
 	                             double dIndices[]) {
 	size_t ndims = (*self)->coords.size();
 	for (size_t i = 0; i < ndims; ++i) {
 		(*self)->dIndices[i] = dIndices[i];
 	}
 	return 0;
}

int SgFindPointInCell_getIndices(SgFindPointInCell_type** self,
 	                              double dIndices[]) {
 	size_t ndims = (*self)->coords.size();
 	for (size_t i = 0; i < ndims; ++i) {
 		dIndices[i] = (*self)->dIndices[i];
 	}
 	return 0;
}
