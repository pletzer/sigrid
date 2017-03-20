/**
 * A class that computes the interpolation weights for a face centered field in 2D
 */
 
#ifndef SG_FLOW_INTERP_2D_H
#define SG_FLOW_INTERP_2D_H
 
#include "SgNdims.h"
#include "SgQuadLineIntersect.h"
#include "SgTriangulate.h"
#include <vector>
#include <algorithm>
#include <map>
#include <cstdio> // size_t
#include <cmath>
 
struct SgFlowInterp2D_type {

	// the source grid
	int srcNodeDims[NDIMS_2D_TOPO];
	int srcCellDims[NDIMS_2D_TOPO];
	int srcNodeDimProd[NDIMS_2D_TOPO];
	int srcCellDimProd[NDIMS_2D_TOPO];
	std::vector<double> srcGrdCoords; // flat array (node, components)
	size_t srcNumPoints;
	size_t srcNumCells;

	// the destination grid 
	int dstNodeDims[NDIMS_1D_TOPO];
	int dstCellDims[NDIMS_1D_TOPO];
	int dstCellDimProd[NDIMS_1D_TOPO];
	int dstNodeDimProd[NDIMS_1D_TOPO];
	std::vector<double> dstGrdCoords; // flat array (node, components)
	size_t dstNumPoints;
	size_t dstNumCells;

	// the flux integrals of x and y 2-forms
	// {dst cell index: {src node index: [wx, wy]}}
	std::map<size_t, std::vector< std::pair<size_t, std::vector<double> > > > weights;

    /** 
     * Constructor
     */
    SgFlowInterp2D_type() {
    }

    /**
     * Destructor
     */
    ~SgFlowInterp2D_type() {
    }

	/**
	 * Set the destination 1D grid 
	 * @param dims number of nodes in each in the single direction
	 * @param coords coordinates (component, node)
	 */
	void setDstGrid(const int dims[], const double** coords) {

		this->dstNumPoints = 1;
		this->dstNumCells = 1;
		for (size_t j = 0; j < NDIMS_1D_TOPO; ++j) {
			this->dstNodeDims[j] = dims[j];
			this->dstCellDims[j] = dims[j] - 1;
			this->dstNumPoints *= dims[j];
			this->dstNumCells *= dims[j] - 1;
		}
		this->dstCellDimProd[NDIMS_1D_TOPO - 1] = 1;
		this->dstNodeDimProd[NDIMS_1D_TOPO - 1] = 1;
		for (int j = NDIMS_1D_TOPO - 2; j >= 0; --j) {
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
		            const double** coords) {

		this->srcNumPoints = dims[0] * dims[1];
		this->srcNumCells = (dims[0] - 1) * (dims[1] - 1);

		this->srcNodeDims[0] = dims[0];
		this->srcNodeDims[1] = dims[1];
		this->srcCellDims[0] = dims[0] - 1;
		this->srcCellDims[1] = dims[1] - 1;

		this->srcCellDimProd[1] = 1;
		this->srcCellDimProd[0] = (dims[1] - 1);
		this->srcNodeDimProd[1] = 1;
		this->srcNodeDimProd[0] = dims[1];
 
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
	 * Apply the interpolation weights to the src edge field
     * @param index 0 or 1 (0 = x flux, 1 = y flux)
	 * @param srcData source edge data with number of nodes dimension (input)
	 * @param dstData destination cell data (output)
	 */
	void apply(size_t index, const double srcData[], double dstData[]) {

		for (std::map<size_t, std::vector<std::pair<size_t, std::vector<double> > > >::const_iterator 
			it = this->weights.begin();
			it != this->weights.end(); ++it) {
			size_t dstIndx = it->first;
			dstData[dstIndx] = 0;
			for (size_t i = 0; i < it->second.size(); ++i) {
				size_t srcIndx = it->second[i].first;
				double wght = it->second[i].second[index];
				dstData[dstIndx] += wght * srcData[srcIndx];
			}
		}
	}


	/** 
	 * Compute the interpolation weights
	 */
	void computeWeights() {

		double dstLineCoords[2*NDIMS_2D_PHYS];
		double srcQuadCoords[4*NDIMS_2D_PHYS];
		size_t srcNodeInds[4]; // 4 nodes
		int offset[2];

		SgQuadLineIntersect_type intersector;

		// iterate over the dst segments
		for (size_t dstIndx = 0; dstIndx < this->dstNumCells; ++dstIndx) {

			// get the start/end points
			offset[0] = 0;
			this->getDstLineCoord(dstIndx, offset, &dstLineCoords[0*NDIMS_2D_PHYS]);
			offset[0] = 1;
			this->getDstLineCoord(dstIndx, offset, &dstLineCoords[1*NDIMS_2D_PHYS]);

			intersector.setLinePoints(dstLineCoords);

			// iterate over the src quads
			for (size_t srcIndx = 0; srcIndx < this->srcNumCells; ++srcIndx) {

				// get the src quad's vertices
				offset[0] = 0; offset[1] = 0;
				this->getSrcQuadCoord(srcIndx, offset, &srcNodeInds[0], &srcQuadCoords[0*NDIMS_2D_PHYS]);
				offset[0] = 1; offset[1] = 0;
				this->getSrcQuadCoord(srcIndx, offset, &srcNodeInds[1], &srcQuadCoords[1*NDIMS_2D_PHYS]);
				offset[0] = 1; offset[1] = 1;
				this->getSrcQuadCoord(srcIndx, offset, &srcNodeInds[2], &srcQuadCoords[2*NDIMS_2D_PHYS]);
				offset[0] = 0; offset[1] = 1;
				this->getSrcQuadCoord(srcIndx, offset, &srcNodeInds[3], &srcQuadCoords[3*NDIMS_2D_PHYS]);

				intersector.setQuadPoints(srcQuadCoords);

				if (!intersector.checkIfOverlap()) {
					// no chance, skip
					continue;
				}

				double* points = NULL;
				int numPoints = 0;
				intersector.collectIntersectPoints(&numPoints, &points);

				if (numPoints == 2) {
					std::map<size_t, std::vector<std::pair<size_t, std::vector<double> > > >::iterator it;
					it = this->weights.find(dstIndx);
					if (it == this->weights.end()) {
						std::vector< std::pair<size_t, std::vector<double> > > v;
						std::pair<size_t, std::vector<std::pair<size_t, std::vector<double> > > > p(dstIndx, v);
						this->weights.insert(p);
					}
					it = this->weights.find(dstIndx);

					double* pA = &points[0*NDIMS_2D_PHYS];
					double* pB = &points[1*NDIMS_2D_PHYS];

					std::vector<double> w(2);
					w[0] = 0.5*(pB[0] - pA[0])*(pB[1] + pA[1]);
					w[1] = 0.5*(pB[0] + pA[0])*(pB[1] - pA[1]);

					std::pair<size_t, std::vector<double> > p(srcNodeInds[0], w);
					it->second.push_back(p);
				}
			}
		}

	}

private:

    /** 
	 * Extract the destination cell coordinates from the grid
	 * @param indx cell flat index
	 * @param offset displacement from the above node
	 * @param coords array of size NDIMS_2D_PHYS to be filled in 
	 */
	void getDstLineCoord(size_t indx, const int offset[], double coords[]) const {

		// compute the index set of the cell and add the offset
		size_t cellIndsOffset[NDIMS_1D_TOPO];
		for (size_t j = 0; j < NDIMS_1D_TOPO; ++j) {
    		cellIndsOffset[j] = indx / this->dstCellDimProd[j] % this->dstCellDims[j];
    		cellIndsOffset[j] += offset[j];
  		}

  		// compute the low-corner flat index of the node coorresponding to this cell
  		size_t nodeIndx = this->dstNodeDimProd[0] * cellIndsOffset[0];

		// fill in the node's coordinates
		for (size_t j = 0; j < NDIMS_1D_PHYS; ++j) {
			coords[j] = this->dstGrdCoords[nodeIndx * NDIMS_2D_PHYS + j];
		}
	}

	/** 
	 * Extract the source cell coordinates from the grid
	 * @param indx nodal flat index
	 * @param offset displacement from the above node
	 * @param nodeIndx grid node index (output)
	 * @param coords array of size NDIMS_2D_PHYS to be filled in 
	 */
	void getSrcQuadCoord(size_t indx, const int offset[], size_t* nodeIndx, double coords[]) const {

		// compute the index set of the cell and add the offset
		size_t cellIndsOffset[NDIMS_2D_TOPO];
		for (size_t j = 0; j < NDIMS_2D_TOPO; ++j) {
    		cellIndsOffset[j] = indx / this->srcCellDimProd[j] % this->srcCellDims[j];
    		cellIndsOffset[j] += offset[j];
  		}

  		// compute the low-corner flat index of the node coorresponding to this cell
  		*nodeIndx = 0;
		for (size_t j = 0; j < NDIMS_2D_TOPO; ++j) {
			*nodeIndx += this->srcNodeDimProd[j] * cellIndsOffset[j];
		}

		// fill in the node's coordinates
		for (size_t j = 0; j < NDIMS_2D_PHYS; ++j) {
			coords[j] = this->srcGrdCoords[*nodeIndx * NDIMS_2D_PHYS + j];
		}
	}

};
 
#ifdef __cplusplus
extern "C" {
#endif

    int SgFlowInterp2D_new(SgFlowInterp2D_type** self);
                       
    int SgFlowInterp2D_del(SgFlowInterp2D_type** self);

    int SgFlowInterp2D_setDstGrid(SgFlowInterp2D_type** self, 
 	                                  const int dims[], const double** coords);

    int SgFlowInterp2D_setSrcGrid(SgFlowInterp2D_type** self, 
 	                                  const int dims[], 
 	                                  const double** coords);

    int SgFlowInterp2D_computeWeights(SgFlowInterp2D_type** self);

    int SgFlowInterp2D_apply(SgFlowInterp2D_type** self, int index,
 	                         const double srcData[], double dstData[]);

#ifdef __cplusplus
}
#endif


#endif // SG_FLOW_INTERP_2D_H