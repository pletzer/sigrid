/**
 * A class that finds all intersection points between two grids
 */
 
#ifndef SG_STRUCT_GRID_2D_INTERSECT_H
#define SG_STRUCT_GRID_2D_INTERSECT_H
 
#include <vector>
#include <cstdio> // size_t
#include <cmath>
#include <limits>
#include <algorithm>
#include "SgLinearSolve.h"
#include "SgFindPointInCell.h"
#include "SgNdims.h"
 
struct SgStructGrid2DIntersect_type {

    // linear solver
    SgLinearSolve_type* slvr;

    // source grid
    int srcNodeDims[NDIMS_2D_TOPO];
    int periodicity[NDIMS_2D_TOPO];
    std::vector<double> srcGrdCoords; // flat array (node, components)

    // destination grid
    int dstNumPoints; 
    int dstNodeDims[];NDIMS_2D_TOPO
    int dstCellDims[NDIMS_2D_TOPO];
    int dstNodeDimProd[NDIMS_2D_TOPO];
    int dstCellDimProd[NDIMS_2D_TOPO];
    std::vector<double> dstGrdCoords; // flat array (node, components)

    // the dst grid nodes in src index space
    std::vector<double> dstNodesInSrcIndexSpace;
    std::vector<int> dstNodesValidMask;

    // maps dst edges to intersecting src grid line points
    std::map<int[2], std::vector<double> > dstEdges2SrcIndexPoints;

    double eps;

    // a tolerance to determine whether a two edges intersect
    // or a point is within a triangle
    double tol;

    // tolerance of point in index space
    double tolpos;

    // max number of Newton iterations
    int nitermax;

    /**
     * Constructor
     */
    SgStructGrid2DIntersect_type() {

        // tolerance to attribute a node to a given cell
        this->eps = 1.e-8;

        // tolerance for floating point comparisons
        this->tol = 1.e-12;

        // Newton solver
        this->tolpos = 1.e-6;
        this->nitermax = 1000;

        // we're in 2D
        this->slvr = new SgLinearSolve_type(2, 2);

        this->locator = new SgFindPointInCell_type(this->nitermax, this->tolpos);
    }

    /**
     * Destructor
     */
    ~SgStructGrid2DIntersect_type() {
        delete this->slvr;
        delete this->locator;
    }

    /**
     * Set the destination grid 
     * @param dims number of nodes in each direction
     * @param coords coordinates (component, node)
     */
    void setDstGrid(const int dims[], const double** coords) {

        this->dstNumPoints = 1;
        this->dstNumCells = 1;
        for (size_t j = 0; j < NDIMS_2D_TOPO; ++j) {
            this->dstNodeDims[j] = dims[j];
            this->dstCellDims[j] = dims[j] - 1;
            this->dstNumPoints *= dims[j];
            this->dstNumCells *= dims[j] - 1;
        }
        this->dstCellDimProd[NDIMS_2D_TOPO - 1] = 1;
        this->dstNodeDimProd[NDIMS_2D_TOPO - 1] = 1;
        for (int j = NDIMS_2D_TOPO - 2; j >= 0; --j) {
            // last index varies fastest
            this->dstCellDimProd[j] = this->dstCellDimProd[j + 1] * this->dstCellDims[j + 1];
            this->dstNodeDimProd[j] = this->dstNodeDimProd[j + 1] * this->dstNodeDims[j + 1];
        }

        this->dstGrdCoords.resize(NDIMS_2D_PHYS * this->dstNumPoints);
        for (size_t j = 0; j < NDIMS_2D_PHYS; ++j) {
            size_t k = 0;
            for (size_t i = 0; i < this->dstNumPoints; ++i) {
                this->dstGrdCoords[k*NDIMS_2D_PHYS + j] = coords[j][i];
                k++;
            }
        }
    }

    /**
     * Set the source grid 
     * @param dims number of nodes in each direction
     * @param coords coordinates (component, node)
     */
    void setSrcGrid(const int dims[], 
                    const int periodicity[], 
                    const double** coords) {

        int srcNumPoints = 1;
        for (size_t j = 0; j < NDIMS_2D_TOPO; ++j) {
            this->srcNodeDims[j] = dims[j];
            srcNumPoints *= dims[j];
            this->periodicity[j] = periodicity[j];
        }

        this->srcGrdCoords.resize(NDIMS_2D_PHYS * srcNumPoints);
        for (size_t j = 0; j < NDIMS_2D_PHYS; ++j) {
            size_t k = 0;
            for (size_t i = 0; i < srcNumPoints; ++i) {
                this->srcGrdCoords[k*NDIMS_2D_PHYS + j] = coords[j][i];
                k++;
            }
        }
    }

    void computeDstNodesInSrcIndexSpace() {

        this->dstNodesInSrcIndexSpace.reserve(NDIMS_2D_TOPO * this->dstNumPoints);

        // point in grid locator
        SgFindPointInCell_type locator(this->nitermax, this->tolpos);
        locator.setGrid(2, this->srcNodeDims, this->periodicity, this->srcGrdCoords);

        // initial guess, somewhere in the middle of the domain
        double dIndices[NDIMS_2D_TOPO];
        dIndices[0] = (this->srcNodeDims[0] - 1) / 2.11432432;
        dIndices[1] = (this->srcNodeDims[1] - 1) / 1.92568776;

        // iterate over the dst grid nodes
        for (int i = 0; i < this->dstNumPoints; ++i) { // consider using the snake iterator here!!!!
            const double* point = &this->dstGrdCoords[i * NDIMS_2D_PHYS];
            locator.reset(dIndices, point);
            int finished = 0;
            while(finished == 0) {
                finished = locator.next();
            }
            // 1 == success
            // -1 means hit max number of iterations (still ok)
            // -2 hit a fixed point? 
            // -3 outside the domain
            this->dstNodesValidMask[i] = 0;
            if (finished == 1 || finished == -1) {
                this->dstNodesValidMask[i] = 1;
                std::vector<double>  di = locator.getPosition();
                // store the result
                this->dstNodesInSrcIndexSpace.push_back(di[0]);
                this->dstNodesInSrcIndexSpace.push_back(di[1]);
                // use this as the next iteration guess
                dIndices[0] = di[0];
                dIndices[1] = di[1];
            }
        }
    }

    /**
     * 
     */
    void computeSrcGridLineWithDstCellEdgeIntersections() {

        // iterate over the horizontal src grid lines 
        for (int j = 0; j < this->srcNodeDims[0]) {
            
        }

        // iterate over the dst cell edges

        // find the start/end points of the dst edge in src grid index space

        // iterate over the src grid lines that intersect with this dst edge

        // add the points to the list

    }

    void findSrcPointsInsideDstCell() {

        // iterate over the dst grid cells

        // find the src grid nodes inside the dst cell

        // add the points to the list

    }

    void computeIntersectionPoints() {
        this->computeSrcGridLineWithDstCellEdgeIntersections();
        this->findSrcPointsInsideDstCell();
    }

};
 
#ifdef __cplusplus
extern "C" {
#endif

int SgStructGrid2DIntersect_new(SgStructGrid2DIntersect_type** self);
                       
int SgStructGrid2DIntersect_del(SgStructGrid2DIntersect_type** self);

int SgStructGrid2DIntersect_setSrcGrid(SgStructGrid2DIntersect_type** self,
                                       const int dims[], const int periodicity[], 
                                       const double** coords);

int SgStructGrid2DIntersect_setDstGrid(SgStructGrid2DIntersect_type** self,
                                       const int dims[],
                                       const double** coords);

#ifdef __cplusplus
}
#endif

#endif // SG_STRUCT_GRID_2D_INTERSECT_H
