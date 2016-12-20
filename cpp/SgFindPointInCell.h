/**
 * A class that finds the location of a point in the index space of a structured grid
 */
 
#ifndef SG_FIND_POINT_IN_CELL_H
#define SG_FIND_POINT_IN_CELL_H
 
#include <vector>
#include <cstdio> // size_t
#include <cmath>
#include <limits>
#include "SgLinearSolve.h"
 
struct SgFindPointInCell_type {

    // linear solver
    SgLinearSolve_type* slvr;

    // the curvilinear coordinate points as ndims flat arrays
    std::vector<std::vector<double> > coords;

    // the min coordinate values of the src grid
    std::vector<double> xmins;

    // the max coordinate values of the src grid
    std::vector<double> xmaxs;

    // used to get the flat index
    std::vector<size_t> prodDims;

    // used when iterating over cell nodes
    std::vector<size_t> nodeProdDims;

    // periodicity length in index space, either 0 (no periodicity) or dims
    std::vector<double> indexPeriodicity;

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

    // history of the errors
    std::vector<double> errorHistory;

    // old iteration error
    double oldError;

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
    this->oldError = std::numeric_limits<double>::max();
    this->errorHistory.reserve(nitermax);
}

/**
 * Destructor
 */
~SgFindPointInCell_type() {
    if (this->slvr) delete this->slvr;
}

/**
 * Set the grid 
 * @param ndims number of physical and topological dimensions
 * @param dims number of nodes along each dimensions
 * @param periodicity flags indicating peridocity (1) or lack or it (0) along each dimension
 * @param coords arrays of flat coordinates (component, coordinates)
 */
void setGrid(int ndims, const int dims[], const int periodicity[],
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
    this->xmins.resize(ndims, std::numeric_limits<double>::max());
    this->xmaxs.resize(ndims, -std::numeric_limits<double>::max());
    // default is no periodicity
    this->indexPeriodicity.resize(ndims, 0.0);
    for (size_t i = 0; i < ndims; ++i) {
        if (periodicity[i] != 0) {
            this->indexPeriodicity[i] = (double) dims[i];
        }
    }

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
            this->xmins[i] = (this->xmins[i] < coords[i][j]? this->xmins[i]: coords[i][j]);
            this->xmaxs[i] = (this->xmaxs[i] > coords[i][j]? this->xmaxs[i]: coords[i][j]);
        }
    }

    // create a solver
    this->slvr = new SgLinearSolve_type(ndims, ndims);
}

/**
 * Get min domain corner 
 * @return min coordinates
 */
 const double* getSrcDomainLo() const {
    return &this->xmins.front();
 }

/**
 * Get max domain corner 
 * @return max coordinates
 */
 const double* getSrcDomainHi() const {
    return &this->xmaxs.front();
 }

/**
 * Get the current position
 * @return position
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
    this->errorHistory.resize(0);
}


/** 
 * Perform one Newton iteration
 * @return 0 if not yet reached target
 *         1 converged
 *         -1 hit max number of iterations
 *         -2 hit a fixed point?
 *         -3 outside of domain
 */
int next() {

    size_t ndims = this->dims.size();
    const double eps = std::numeric_limits<double>::epsilon();

    // check if inside domain
    bool inside = true;
    for (size_t i = 0; i < ndims; ++i) {
        inside &= this->targetPoint[i] >= this->xmins[i] - eps;
        inside &= this->targetPoint[i] <= this->xmaxs[i] + eps;
    }
    if (!inside) {
        return -3;
    }

    this->computeJacobianAndRHS();
    this->slvr->setMatrix(&this->jacMatrix[0]);
    this->slvr->setRightHandSide(&this->rhs[0]);
    this->slvr->solve();

    double* sol;
    this->slvr->getSolution(&sol);

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
    this->errorHistory.push_back(posError);
    if (posError > this->oldError + eps) {
        // we're not improving, fixed point?
        // stop if the error has not improved in the last 3 iterations
        size_t niter = this->errorHistory.size();
        //if (niter >= 3) {
        //    if (posError > this->errorHistory[niter - 3] + eps) {
        //        // not improving
        //        return -2;
        //    }
        //}
    }

    if (this->iter >= this->nitermax) {
        // reached max number of iterations
        return -1; 
    }
    if (posError < this->tolpos) {
        // done!
        return 1;
    }

    // has not yet converged
    this->oldError = posError;
    return 0;
}

/**
 * Get the current index space position
 * @return index position
 */
std::vector<double> getIndices() const {
    return this->dIndices;
}

/**
 * Get a pointer to the error history
 * @param niter number of iterations so far (output)
 * @param errorHistory pointer to the history of errors in physical space
 */
void getErrorHistory(int* niter, double* errorHistory[]) {
    *niter = (int) this->errorHistory.size();
    *errorHistory = &this->errorHistory.front();
}

/**
 * Get the flat source indices and interpolation weights 
 * @param dIndices index location in source grid
 * @param array of src flat index and weight pairs
 */
const std::vector< std::pair<size_t, double> > 
getSrcIndicesAndWeights(const std::vector<double>& dInds) const {

    size_t ndims = this->dims.size();
    size_t nnodes = pow(2, ndims);
    std::vector< std::pair<size_t, double> > res(nnodes);

    for (size_t j = 0; j < nnodes; ++j) {

        size_t flatIndx = 0;
        double wght = 1;

        // itereate over the dimensions
        for (size_t i = 0; i < ndims; ++i) {

            int loCornerIndx = (int) floor(dInds[i]);

            // need to a full cell to the right of loCornerIndx
            loCornerIndx = (loCornerIndx >= 0? loCornerIndx: 0);
            loCornerIndx = (loCornerIndx <= this->dims[i] - 2? loCornerIndx: this->dims[i] - 2);

            // index along this dimension
            int indx = loCornerIndx + (j / this->nodeProdDims[i] % 2);

            // update flat index
            flatIndx += (size_t) this->prodDims[i] * indx;

            double dindx = (double) indx;

            // weight depends on whether the node is on the low or high side
            wght *= (dInds[i] >= dindx? dindx + 1 - dInds[i]: dInds[i] - dindx + 1);
        }

        res[j] = std::pair<size_t, double>(flatIndx, wght);
    }

    return res;
}

private:

/**
 * Interpolate a nodal field
 * @param dInds position in index space
 * @param nodalField values of the nodal fields in the order returned by computeWeightsAndFlatIndices
 */
double interp(const std::vector<double>& dInds,
              const std::vector<double>& nodalField)  const {

    const std::vector< std::pair<size_t, double> > iw = this->getSrcIndicesAndWeights(dInds);

    double res = 0;
    for (std::vector< std::pair<size_t, double> >::const_iterator it = iw.begin();
                      it != iw.end(); ++it) {
        size_t srcNodeIndex = it->first;
        double wght = it->second;
        res += wght * nodalField[srcNodeIndex];
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
 * Ensure that the indices fall inside the index domain
 */
void truncateIndices() {

    // Float indices are truncated when they fall outside the topological domain 
    // (ie beyond the min/max dimensions). This can happen whenever one hits a 
    // pole in which case the Jacobian is singular and the index increment infinite.
    // Retreating by eps from the boundary of the index domain seems to help in this
    // case.
    const double eps = 1.234567890e-4;

    for (size_t i = 0; i < this->dims.size(); ++i) {
        // take into account periodicity 
        if (this->indexPeriodicity[i] != 0) {
            this->dIndices[i] = fmod(
                fmod(dIndices[i], this->indexPeriodicity[i]) + this->indexPeriodicity[i], 
                this->indexPeriodicity[i]);
        }
        // make sure the indices are within the domain
        this->dIndices[i] = (this->dIndices[i] < 0? 
                             eps: this->dIndices[i]);
        this->dIndices[i] = (this->dIndices[i] > this->dims[i] - 1? 
                             this->dims[i] - 1 - eps: this->dIndices[i]);
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
                              const int periodicity[], 
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
