/**
 * A class that computes the interpolation weights for a cell centered field
 */
 
#ifndef SG_CONSERVE_INTERP_2D_H
#define SG_CONSERVE_INTERP_2D_H
 
#include "SgNdims.h"
#include "SgQuadIntersect.h"
#include "SgTriangulate.h"
#include <vector>
#include <map>
#include <cstdio> // size_t
#include <cmath>
 
struct SgConserveInterp2D_type {

	// the source grid
	int srcNodeDims[NDIMS_2D_TOPO];
	int srcCellDims[NDIMS_2D_TOPO];
	int srcNodeDimProd[NDIMS_2D_TOPO];
	int srcCellDimProd[NDIMS_2D_TOPO];
	std::vector<double> srcGrdCoords; // flat array (node, components)
	size_t srcNumPoints;
	size_t srcNumCells;

	// the destination grid 
	int dstNodeDims[NDIMS_2D_TOPO];
	int dstCellDims[NDIMS_2D_TOPO];
	int dstCellDimProd[NDIMS_2D_TOPO];
	int dstNodeDimProd[NDIMS_2D_TOPO];
	std::vector<double> dstGrdCoords; // flat array (node, components)
	size_t dstNumPoints;
	size_t dstNumCells;

	// the interpolation weights
	// dst flat index -> {src flat index: weight}
	std::map<size_t, std::vector< std::pair<size_t, double> > > weights;

	// iterators
    std::map<size_t, std::vector<std::pair<size_t, double> > >::const_iterator weightIt;
    std::vector<std::pair<size_t, double> >::const_iterator weightIt2;

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
	void setSrcGrid(const int dims[], const double** coords) {

		this->srcNumPoints = 1;
		this->srcNumCells = 1;
		for (size_t j = 0; j < NDIMS_2D_TOPO; ++j) {
			this->srcNodeDims[j] = dims[j];
			this->srcCellDims[j] = dims[j] - 1;
			this->srcNumPoints *= dims[j];
			this->srcNumCells *= dims[j] - 1;
		}
		this->srcCellDimProd[NDIMS_2D_TOPO - 1] = 1;
		this->srcNodeDimProd[NDIMS_2D_TOPO - 1] = 1;
		for (int j = NDIMS_2D_TOPO - 2; j >= 0; --j) {
    		// last index varies fastest
    		this->srcCellDimProd[j] = this->srcCellDimProd[j + 1] * this->srcCellDims[j + 1];
    		this->srcNodeDimProd[j] = this->srcNodeDimProd[j + 1] * this->srcNodeDims[j + 1];
  		}
 
  		this->srcGrdCoords.resize(NDIMS_2D_PHYS * this->srcNumPoints);
  		for (size_t j = 0; j < NDIMS_2D_PHYS; ++j) {
  			size_t k = 0;
  			for (size_t i = 0; i < this->srcNumPoints; ++i) {
  				this->srcGrdCoords[k*NDIMS_2D_PHYS + j] = coords[j][i];
  				k++;
  			}
  		}
	}

	/** 
	 * Apply the interpolation weights to the cell centered field
	 * @param srcData source cell data (input)
	 * @param dstData destination cell data (output)
	 */
	void apply(const double srcData[], double dstData[]) {

		for (std::map<size_t, std::vector<std::pair<size_t, double> > >::const_iterator 
			it = this->weights.begin();
			it != this->weights.end(); ++it) {
			size_t dstIndx = it->first;
			dstData[dstIndx] = 0;
			for (size_t i = 0; i < it->second.size(); ++i) {
				size_t srcIndx = it->second[i].first;
				double wght = it->second[i].second;
				dstData[dstIndx] += wght * srcData[srcIndx];
			}
		}
	}

	/** 
	 * Extract the destination cell coordinates from the grid
	 * @param indx nodal flat index
	 * @param offset displacement from the above node
	 * @param coords array of asize NDIMS_2D_PHYS to be filled in 
	 */
	void getDstQuadCoord(size_t indx, int offset[], double coords[]) const {
		for (size_t j = 0; j < NDIMS_2D_TOPO; ++j) {
			indx += this->dstNodeDimProd[j] * offset[j];
		}
		for (size_t j = 0; j < NDIMS_2D_PHYS; ++j) {
			coords[j] = this->dstGrdCoords[indx*NDIMS_2D_PHYS + j];
		}
	}

	/** 
	 * Extract the destination cell coordinates from the grid
	 * @param indx nodal flat index
	 * @param offset displacement from the above node
	 * @param coords array of asize NDIMS_2D_PHYS to be filled in 
	 */
	void getSrcQuadCoord(size_t indx, int offset[], double coords[]) const {
		for (size_t j = 0; j < NDIMS_2D_TOPO; ++j) {
			indx += this->srcNodeDimProd[j] * offset[j];
		}
		for (size_t j = 0; j < NDIMS_2D_PHYS; ++j) {
			coords[j] = this->srcGrdCoords[indx*NDIMS_2D_PHYS + j];
		}
	}

	/** 
	 * Compute the interpolation weights
	 */
	void computeWeights() {

		int numIntersectPoints;
		double* intersectPoints;
		SgQuadIntersect_type intersector;
		double dstQuadCoords[NDIMS_2D_PHYS*4]; // four nodes
		double srcQuadCoords[NDIMS_2D_PHYS*4]; // four nodes
		int offset[] = {0, 0};

		// iterate over the dst cells
		for (size_t dstIndx = 0; dstIndx < this->dstNumCells; ++dstIndx) {

			offset[0] = 0; offset[1] = 0;
			this->getDstQuadCoord(dstIndx, offset, &dstQuadCoords[0*NDIMS_2D_PHYS]);
			offset[0] = 1; offset[1] = 0;
			this->getDstQuadCoord(dstIndx, offset, &dstQuadCoords[1*NDIMS_2D_PHYS]);
			offset[0] = 1; offset[1] = 1;
			this->getDstQuadCoord(dstIndx, offset, &dstQuadCoords[2*NDIMS_2D_PHYS]);
			offset[0] = 0; offset[1] = 1;
			this->getDstQuadCoord(dstIndx, offset, &dstQuadCoords[3*NDIMS_2D_PHYS]);

			// compute the dst cell area
			SgTriangulate_type dstTriangulator(4, dstQuadCoords);
			double dstArea = dstTriangulator.getConvexHullArea();

			// iterate over the src cells
			std::vector< std::pair<size_t, double> > indWght;
			for (size_t srcIndx = 0; srcIndx < this->srcNumCells; ++srcIndx) {

				offset[0] = 0; offset[1] = 0;
				this->getSrcQuadCoord(srcIndx, offset, &srcQuadCoords[0*NDIMS_2D_PHYS]);
				offset[0] = 1; offset[1] = 0;
				this->getSrcQuadCoord(srcIndx, offset, &srcQuadCoords[1*NDIMS_2D_PHYS]);
				offset[0] = 1; offset[1] = 1;
				this->getSrcQuadCoord(srcIndx, offset, &srcQuadCoords[2*NDIMS_2D_PHYS]);
				offset[0] = 0; offset[1] = 1;
				this->getSrcQuadCoord(srcIndx, offset, &srcQuadCoords[3*NDIMS_2D_PHYS]);

				intersector.reset();
				intersector.setQuadPoints(dstQuadCoords, srcQuadCoords);
				intersector.collectIntersectPoints(&numIntersectPoints, &intersectPoints);
				if (numIntersectPoints >= 3) {
					// must be able to build at least one triangle
					SgTriangulate_type triangulator(numIntersectPoints, intersectPoints);
					double area = triangulator.getConvexHullArea();
					indWght.push_back(std::pair<size_t, double>(srcIndx, area/dstArea));
				}
			}
			if (indWght.size() > 0) {
				// add the overlap contributions
				std::pair<size_t, std::vector< std::pair<size_t, double> > > p(dstIndx, indWght);
				this->weights.insert(p);
			}
		}
	}
};
 
#ifdef __cplusplus
extern "C" {
#endif

    int SgConserveInterp2D_new(SgConserveInterp2D_type** self);
                       
    int SgConserveInterp2D_del(SgConserveInterp2D_type** self);

    int SgConserveInterp2D_setDstGrid(SgConserveInterp2D_type** self, 
 	                                  const int dims[], const double** coords);

    int SgConserveInterp2D_setSrcGrid(SgConserveInterp2D_type** self, 
 	                                  const int dims[], const double** coords);

    int SgConserveInterp2D_computeWeights(SgConserveInterp2D_type** self);

    int SgConserveInterp2D_apply(SgConserveInterp2D_type** self,
 	                             const double srcData[], double dstData[]);

    int SgConserveInterp2D_reset(SgConserveInterp2D_type** self);

    int SgConserveInterp2D_next(SgConserveInterp2D_type** self);

    int SgConserveInterp2D_get(SgConserveInterp2D_type** self,
 	                           int* srcIndx, int* dstIndx, double* weight);

#ifdef __cplusplus
}
#endif


#endif // SG_CONSERVE_INTERP_2D_H