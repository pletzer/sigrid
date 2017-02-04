/**
 * A class that finds all the structured grid cells that overlap with a polygon using simple box heuristics
 */
 
#ifndef SG_FIND_OVERLAPPING_CELLS_2D_H
#define SG_FIND_OVERLAPPING_CELLS_2D_H
 
#include "SgNdims.h"
#include "SgFindPointInCell.h"
#include "SgBoxIterator.h"
#include <vector>
#include <cstdio> // size_t

struct SgFindOverlappingCells2D_type {

    // low corner index set 
    std::vector<double> loIndxCorner;

    // high corner index set
    std::vector<double> hiIndxCorner;

    // the dimensions of the source grid
    std::vector<int> srcNodeDims;

    // to go from an index set to a flat index
    std::vector<int> srcCellProdDims;

    // collection of cell flat indices
    std::vector<int> srcCellFlatInds;

    // source grid coordinates (component, node)
    double** srcCoords;

    // Newton scheme to find the index set corresponding to a position 
    SgFindPointInCell_type* pointFinder;

    // polygon points as a flat array
    std::vector<double> polyPoints;

    // fuzzy region in index space so we're sure to include enough cells
    double eps;

    size_t numPolyPoints;

    /**
     * Constructor
     */
    SgFindOverlappingCells2D_type() {

        this->loIndxCorner.resize(NDIMS_2D_TOPO); // 2D
        this->hiIndxCorner.resize(NDIMS_2D_TOPO); // 2D
        const int nitermax = 100;
        const double tolpos = 1.e-6;
        this->pointFinder = new SgFindPointInCell_type(nitermax, tolpos);
        this->eps = 0.01;
        this->srcCellFlatInds.reserve(100); // rough estimate
    }

    /** 
     * Destructor
     */
    ~SgFindOverlappingCells2D_type() {
        delete this->pointFinder;
    }

    /**
     * Set the source grid 
     * @param dims number of nodes in each direction
     * @param periodicity array of 0 (not periodic) and 1 (periodic)
     * @param coords coordinates (component, node)
     */
    void setSrcGrid(const int dims[], const int periodicity[], const double** coords) {
        this->srcNodeDims.resize(NDIMS_2D_TOPO);
        this->srcCellProdDims.resize(NDIMS_2D_TOPO);
        for (size_t j = 0; j < NDIMS_2D_TOPO; ++j) {
            this->srcNodeDims[j] = dims[j];
        }
        this->srcCellProdDims[NDIMS_2D_TOPO - 1] = 1;
        for (int j = NDIMS_2D_TOPO - 2; j >= 0; --j) {
            // last index varies fastest
            this->srcCellProdDims[j] = this->srcCellProdDims[j + 1] * dims[j + 1];
          }
        this->pointFinder->setGrid(2, dims, periodicity, coords);
    }

    /** 
     * Set the polygon's nodes
     * @param numPoints number of points
     * @param coords point coordinates as a flat array (node, components)
     */
    void setPolygonPoints(int numPoints, const double coords[]) {
        this->numPolyPoints = numPoints;
        this->polyPoints.resize(this->numPolyPoints * NDIMS_2D_PHYS);
        for (size_t i = 0; i < this->numPolyPoints; ++i) {
            for (size_t j = 0; j < NDIMS_2D_PHYS; ++j) {
                this->polyPoints[NDIMS_2D_PHYS*i + j] = coords[NDIMS_2D_PHYS*i + j];
            }
        }
    }

    /**
     * Find the low/high corner indices of the polygon 
     */
    int findFloatIndexBoxOfPolygon() {

        // initial guess, somewhere in the middle
        double dIndices[] = {this->srcNodeDims[0]/2.123353, this->srcNodeDims[1]/1.9647354};

        // set the lo/hi corner values to the max/min values allowed in the domain
        for (size_t j = 0; j < NDIMS_2D_TOPO; ++j) {
            this->loIndxCorner[j] = this->srcNodeDims[j] - 1;
            this->hiIndxCorner[j] = 0;
        }

        // Newton iterator's error flag (0 = continue, 1 = finished, -1 = finished with an error)
        int ier;

        // total number of Newton search errors
        int totError = 0;

        // iterate over the polygon's points
        for (size_t i = 0; i < this->numPolyPoints; ++i) {

            bool iterFlag = true;

            // make this polygon's point th etarget point
            const double* targetPoint = &this->polyPoints[i*NDIMS_2D_PHYS];

            // start search
            this->pointFinder->reset(dIndices, targetPoint);
            while (iterFlag) {
                ier = this->pointFinder->next();
                if (ier != 0) {
                    // reached end of iterations
                    iterFlag = false;
                }
            }

            // error handling. In the case where the target point is outside the src grid, 
            // we expect pointFinder to return the index location just inside the src grid
            // corresponding to the physical position that is closest to the target point. 
            // I don't think we need to do anything special... TO CHECK!!!
            //if (ier < 0) {
            //    std::cerr << "ERROR: failed to find the index point for " << 
            //         this->polyPoints[i*NDIMS_2D_PHYS] << ", " << this->polyPoints[i*NDIMS_2D_PHYS + 1] << '\n';
            //    totError += 1;
            //}

            // update the low/high box corners
            std::vector<double> dInds = this->pointFinder->getIndices();
            for (size_t j = 0; j < this->loIndxCorner.size(); ++j) {
                this->loIndxCorner[j] = (dInds[j] < this->loIndxCorner[j]?
                                         dInds[j]: this->loIndxCorner[j]);
                this->hiIndxCorner[j] = (dInds[j] > this->hiIndxCorner[j]?
                                         dInds[j]: this->hiIndxCorner[j]);

                // make sure we don't go out of bounds
                this->loIndxCorner[j] = (this->loIndxCorner[j] >= 0?
                                         this->loIndxCorner[j]: 0);
                this->hiIndxCorner[j] = (this->hiIndxCorner[j] <= this->srcNodeDims[j] - 1?
                                         this->hiIndxCorner[j]: this->srcNodeDims[j] - 1);

            }
        }

        return totError;
    }

    /**
     * Find the source grid cells under the polygon
     */
    void findSrcCellIndices() {

        this->findFloatIndexBoxOfPolygon();

        // set the low/high of the index box
        int loInds[NDIMS_2D_TOPO];
        int hiInds[NDIMS_2D_TOPO];
        for (size_t j = 0; j < NDIMS_2D_TOPO; ++j) {
            // eps ensures that ceil and floor never give the same value
            loInds[j] = floor(this->loIndxCorner[j]);
            hiInds[j] = ceil(this->hiIndxCorner[j] + this->eps);
            // make sure the lo/hi index corners are in the domain
            loInds[j] = (loInds[j] >= 0? loInds[j]: 0);
            hiInds[j] = (hiInds[j] < this->srcNodeDims[j] - 1?
                         hiInds[j]: this->srcNodeDims[j] - 1);
        }

        // iterate over the cells inside the index box
        int inds[NDIMS_2D_TOPO];
        this->srcCellFlatInds.resize(0);
        SgBoxIterator_type boxIter(NDIMS_2D_TOPO, loInds, hiInds);
        int numCells = boxIter.getNumberOfElements();
        for (int i = 0; i < numCells; ++i) {
            boxIter.getElement(i, inds);
            // compute the flat index of the cell
            size_t cellIndx = 0;
            for (size_t j = 0; j < NDIMS_2D_TOPO; ++j) {
                cellIndx += this->srcCellProdDims[j] * inds[j];
            }
            // add the index to the list
            this->srcCellFlatInds.push_back(cellIndx);
        }
    }

    /**
     * Get the collection of source cells under this polygon
     * @return array of flat indices
     */
    std::vector<int> getSrcCellIndices() const {
        return this->srcCellFlatInds;
    }
};
 
#ifdef __cplusplus
extern "C" {
#endif

    int SgFindOverlappingCells2D_new(SgFindOverlappingCells2D_type** self);
                       
    int SgFindOverlappingCells2D_del(SgFindOverlappingCells2D_type** self);

    int SgFindOverlappingCells2D_setSrcGrid(SgFindOverlappingCells2D_type** self,
                                            const int dims[], const int periodicity[],
                                            const double** coords);

    int SgFindOverlappingCells2D_setPolygonPoints(SgFindOverlappingCells2D_type** self, 
                                                   int numPoints, const double coords[]);

    int SgFindOverlappingCells2D_findSrcCellIndices(SgFindOverlappingCells2D_type** self);

    int SgFindOverlappingCells2D_getNumberOfSrcCellIndices(SgFindOverlappingCells2D_type** self,
                                                           int* numSrcCells);

    int SgFindOverlappingCells2D_getSrcCellIndices(SgFindOverlappingCells2D_type** self,
                                                   int** srcCellInds);

#ifdef __cplusplus
}
#endif


#endif // SG_FIND_OVERLAPPING_CELLS_2D_H
