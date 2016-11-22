/**
 * A class that finds all the structured grid cells that overlap with a polygon
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

		this->loIndxCorner.resize(2); // 2D
		this->hiIndxCorner.resize(2); // 2D
		const int nitermax = 100;
		const double tolpos = 1.e-6;
		this->pointFinder = new SgFindPointInCell_type(nitermax, tolpos);
		this->eps = 0.01;
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
	 * @param coords coordinates (component, node)
	 */
	void setSrcGrid(const int dims[], const double** coords) {

		this->pointFinder->setGrid(2, dims, coords);
	}

	/** 
	 * Set the polygon's nodes
	 * @param numPoints number of points
	 * @param coords point coordinates as a flat array (node, components)
	 */
	void setPolygonPoints(int numPoints, const double coords[]) {
		this->numPolyPoints = numPoints;
		this->polyPoints.resize(this->numPolyPoints);
		for (size_t i = 0; i < this->numPolyPoints; ++i) {
			for (size_t j = 0; j < NDIMS_PHYS_2D; ++j) {
				this->polyPoints[NDIMS_PHYS_2D*i + j] = coords[NDIMS_PHYS_2D*i + j];
			}
		}
	}

	void findIndexBox() {
		// TO DO!!!
	}

	std::vector<size_t> getSrcCellIndices() {

		this->findIndexBox();

		// set the low/high of the index box
		int loInds[NDIMS_TOPO_2D];
		int hiInds[NDIMS_TOPO_2D];
		for (size_t j = 0; j < NDIMS_TOPO_2D; ++j) {
			loInds[j] = floor(this->loIndxCorner[j] - this->eps);
			hiInds[j] = floor(this->hiIndxCorner[j] + this->eps);
			// make sure the lo/hi index corners are in the domain
			loInds[j] = (loInds[j] >= 0? loInds[j]: 0);
			hiInds[j] = (hiInds[j] < this->srcCellDims[j]? hiInds[j]: this->srcCellDims[j] - 1);
		}

		// iterate over the cells inside the index box
		int inds[NDIMS_TOPO_2D];
		std::vector<size_t> cellFlatInds;
		cellFlatInds.capacity(100); // rough guess
		SgBoxIterator_type boxIter(loInds, hiInds);
		int numCells = boxIter.getNumberOfElements();
		for (int i = 0; i < numCells; ++i) {
			boxIter.getElement(i, inds);
			// compute the flat index of the cell
			size_t cellIndx = 0;
			for (size_t j = 0; j < NDIMS_TOPO_2D; ++j) {
				cellIndx += this->srcCellProdDims[j] * inds[j];
			}
			// add the index to the list
			cellFlatInds.push_back(cellIndx);
		}

		return cellFlatInds;
	}

};
 
#ifdef __cplusplus
extern "C" {
#endif

    int SgFindOverlappingCells2D_new(SgFindOverlappingCells2D_type** self);
                       
    int SgFindOverlappingCells2D_del(SgFindOverlappingCells2D_type** self);

    int SgFindOverlappingCells2D_setDstGrid(SgFindOverlappingCells2D_type** self, 
 	                                  const int dims[], const double** coords);

    int SgFindOverlappingCells2D_setSrcGrid(SgFindOverlappingCells2D_type** self, 
 	                                  const int dims[], const double** coords);

    int SgFindOverlappingCells2D_computeWeights(SgFindOverlappingCells2D_type** self);

    int SgFindOverlappingCells2D_apply(SgFindOverlappingCells2D_type** self,
 	                             const double srcData[], double dstData[]);

    int SgFindOverlappingCells2D_reset(SgFindOverlappingCells2D_type** self);

    int SgFindOverlappingCells2D_next(SgFindOverlappingCells2D_type** self);

    int SgFindOverlappingCells2D_get(SgFindOverlappingCells2D_type** self,
 	                           int* srcIndx, int* dstIndx, double* weight);

#ifdef __cplusplus
}
#endif


#endif // SG_CONSERVE_INTERP_2D_H