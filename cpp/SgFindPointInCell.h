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

    // used when iterating over cell nodes
    std::vector<size_t> nodeProdDims;

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
/**
 * Constructor
 * @param nitermax max number of Newton iterations
 * @param tolpos tolerance in physical space
 */
SgFindPointInCell_type(int nitermax, double tolpos) {
    this->tolpos = tolpos;
    this->nitermax = nitermax;
    this->slvr = NULL;
    this->iter = 0;
}

~SgFindPointInCell_type() {
    if (this->slvr) delete this->slvr;
}

/**
 * Get the interpolation weights and flat node indices
 * @param dInds position in index space
 * @param weights interpolation weights to be filled in
 * @param flatInds flat node indices to be filled in
 */
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
            int indx = loCornerIndx + (j / this->nodeProdDims[i] % 2);
            flatInds[j] += (size_t) this->prodDims[i] * indx;
            double dindx = (double) indx;
            double w = (dInds[i] >= dindx? dindx + 1 - dInds[i]: dInds[i] - dindx + 1);
            weights[j] *= w;
        }
    }
}

/**
 * Interpolate a nodal field
 * @param dInds position in index space
 * @param nodalField values of the nodal fields in the order returned by getWeightsAndFlatIndices
 */
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

/**
 * Compute the current Jacobian and right hand side vector using finite differences
 */
void computeJacobianAndRHS() {

    size_t ndims = this->dims.size();

    std::vector<double> dIndsHi(dIndices);
    std::vector<double> dIndsLo(dIndices);

    // iterate over the physical space dimensions
    size_t k = 0;
    for (size_t i = 0; i < ndims; ++i) {

        double pos = this->interp(this->dIndices, this->coords[i]);
        this->rhs[i] = this->targetPoint[i] - pos;

        // iterate over the index space dimensions
        for (size_t j = 0; j < ndims; ++j) {

            dIndsHi[j] += 0.5;
            dIndsHi[j] = (dIndsHi[j] < this->dims[j] - 1? dIndsHi[j]: this->dims[j] - 2);
            double xHi = this->interp(dIndsHi, this->coords[i]);

            dIndsLo[j] -= 0.5;
            dIndsLo[j] = (dIndsLo[j] >= 0? dIndsLo[j]: 0);
            double xLo = this->interp(dIndsLo, this->coords[i]);

            // average difference of the i-th coordinate anlong the j-th topo direction
            this->jacMatrix[k] = (xHi - xLo)/(dIndsHi[j] - dIndsLo[j]);

            // reset
            dIndsHi[j] = dIndices[j];
            dIndsLo[j] = dIndices[j];

            k++;            
        }
    }
}

/**
 * Set the grid 
 * @param ndims number of physical and topological dimensions
 * @param dims number of nodes along each dimensions
 * @param coords arrays of flat coordinates (component, coordinates)
 */
void setGrid(int ndims, const int dims[],
              const double** coords) {

    if (ndims <= 0) return;

    this->coords.resize(ndims);
    this->prodDims.resize(ndims);
    this->nodeProdDims.resize(ndims);
    this->dims.resize(ndims);
    this->jacMatrix.resize(ndims * ndims);
    this->rhs.resize(ndims);
    this->dIndices.resize(ndims);
    this->targetPoint.resize(ndims);

    // must have at least one dimension
    this->prodDims[ndims - 1] = 1;
    this->nodeProdDims[ndims - 1] = 1;
    for (int i = ndims - 2; i >= 0; --i) {
        this->prodDims[i] = this->prodDims[i + 1] * dims[i + 1];
        this->nodeProdDims[i] = this->nodeProdDims[i + 1] * 2;
    }

    // total number of nodes
    int ntot = 1;
    for (int i = 0; i < ndims; ++i) {
        ntot *= dims[i];
        this->dims[i] = dims[i];
    }

    // set the coordinates
    for (int i = 0; i < ndims; ++i) {
        this->coords[i].resize(ntot);
        for (int j = 0; j < ntot; ++j) {
            this->coords[i][j] = coords[i][j];
        }
    }

    // create a solver
    this->slvr = new SgLinearSolve_type(ndims, ndims);
}

/**
 * Get the current position
 * @return psoition
 */
std::vector<double> getPosition() {
    size_t ndims = this->dims.size();
    std::vector<double> pos(ndims);
    for (size_t i = 0; i < ndims; ++i) {
        pos[i] = this->interp(this->dIndices, this->coords[i]);
    }
    return pos;
}

/** 
 * Get the current error
 * @return error
 */
double getError() {
    double res = 0;
    std::vector<double> pos = this->getPosition();
    size_t ndims = this->dims.size();
    for (size_t i = 0; i < ndims; ++i) {
        double dp = pos[i] - this->targetPoint[i];
        res += dp * dp;
    }
    return sqrt(res);
}

void reset(const double dIndices[], const double targetPoint[]) {
    size_t ndims = this->dims.size();
    for (size_t i = 0; i < ndims; ++i) {
        this->dIndices[i] = dIndices[i];
        this->targetPoint[i] = targetPoint[i];
    }
    this->iter = 0;
}

/** 
 * Ensure that the indices fall inside the domain
 */
void truncateIndices() {
    for (size_t i = 0; i < this->dims.size(); ++i) {
        // make sure the indices are within the domain
        this->dIndices[i] = (this->dIndices[i] < 0? 
                             0: this->dIndices[i]);
        this->dIndices[i] = (this->dIndices[i] > this->dims[i] - 1? 
                             this->dims[i] - 1: this->dIndices[i]);
    }
}

/** 
 * Perform one Newton iteration
 * @return 0 if not yet reached target
 *         1 converged
 *         -1 hit max number of iterations
 */
int next() {

    this->computeJacobianAndRHS();
    this->slvr->setMatrix(&this->jacMatrix[0]);
    this->slvr->setRightHandSide(&this->rhs[0]);
    this->slvr->solve();

    double* sol;
    this->slvr->getSolution(&sol);

    size_t ndims = this->dims.size();

    // update the indices
    for (size_t i = 0; i < ndims; ++i) {
        this->dIndices[i] += sol[i];
    }

    // indices must be in valid range
    this->truncateIndices();

    // check if the next iterator is still valid
    this->iter++;

    // update the error
    double posError = this->getError();

    if (this->iter >= this->nitermax) {
        // reached max number of iterations
        return -1; 
    }
    if (posError < this->tolpos) {
        // done!
        return 1;
    }

    // has not yet converged
    return 0;
}

/**
 * Get the current index space position
 * @return index position
 */
std::vector<double> getIndices() const {
    return this->dIndices;
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

int SgFindPointInCell_getError(SgFindPointInCell_type** self, 
                               double* error);

int SgFindPointInCell_reset(SgFindPointInCell_type** self, 
                            const double dIndices[], 
                            const double targetPoint[]);

int SgFindPointInCell_next(SgFindPointInCell_type** self);

int SgFindPointInCell_getIndices(SgFindPointInCell_type** self,
                                 double dIndices[]);
 
#ifdef __cplusplus
}
#endif

#endif // SG_FIND_POINT_IN_CELL_H
