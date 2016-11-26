/**
 * A class that computes the interpolation weights for a nodal field
 */
 
#ifndef SG_NODE_INTERP_2D_H
#define SG_NODE_INTERP_2D_H
 
#include "SgNdims.h"
#include "SgFindPointInCell.h"
#include <vector>
#include <map>
#include <cstdio> // size_t
#include <cmath>
#include <iostream>
 
struct SgNodeInterp2D_type {

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

	// dstNodeIndex -> [(srcNodeIndex, weight), ...]
	std::map<size_t, std::vector<std::pair<size_t, double> > > weights;

	SgFindPointInCell_type* pointFinder;

    /** 
     * Constructor
     */
    SgNodeInterp2D_type() {
    	const int nitermax = 100;
    	const double tolpos = 1.e-8;
    	this->pointFinder = new SgFindPointInCell_type(nitermax, tolpos);
    }

    /**
     * Destructor
     */
    ~SgNodeInterp2D_type() {
    	delete this->pointFinder;
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

		this->pointFinder->setGrid(NDIMS_2D_PHYS, dims, periodicity, coords);

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

		std::map<size_t, std::vector< std::pair<size_t, double> > >::const_iterator it;
		std::vector< std::pair<size_t, double> >::const_iterator it2;
		for (it = this->weights.begin(); it != this->weights.end(); ++it) {
			size_t dstNodeIndex = it->first;
			const std::vector< std::pair<size_t, double> >& iv = it->second;
			dstData[dstNodeIndex] = 0.0;
			for (it2 = iv.begin(); it2 != iv.end(); ++it2) {
				size_t srcNodeIndex = it2->first;
				double wght = it2->second;
				dstData[dstNodeIndex] += wght * srcData[srcNodeIndex];
			}
		}
	}

	/** 
	 * Compute the interpolation weights
	 */
	void computeWeights() {

		// search error code
		int ier;

		this->weights.clear();

		// start somewhere in the middle of the src grid
		std::vector<double> dIndices(NDIMS_2D_TOPO);
		dIndices[0] = this->srcCellDims[0]/2.14423;
		dIndices[1] = this->srcCellDims[1]/1.93265;

		for (size_t dstNodeIndex = 0; dstNodeIndex < this->dstNumPoints; ++dstNodeIndex) {

			// find the index position 
            bool iterFlag = true;
            const double* targetPoint = &this->dstGrdCoords[dstNodeIndex*NDIMS_2D_PHYS];

            // start search
            this->pointFinder->reset(&dIndices[0], targetPoint);
            while (iterFlag) {
                ier = this->pointFinder->next();
                if (ier != 0) {
                    // reached end of iterations
                    iterFlag = false;
                }
            }

            if (ier == 1) {
            	// normal termination, found the index location
            	const std::vector<double>& dInds = this->pointFinder->getIndices();
            	const std::vector< std::pair<size_t, double> > 
            	    iw = this->pointFinder->getSrcIndicesAndWeights(dInds);
            	std::pair<size_t, std::vector< std::pair<size_t, double> > > p(dstNodeIndex, iw);
            	this->weights.insert(p);

            }
            else {
            	// likely, dst point is outside the src domain
            	// (or search did not converge)
            	std::cerr << "ERROR: failed to find the index point for " << 
                          targetPoint[0] << ", " << targetPoint[1] << '\n';
                std::vector<double> pos = this->pointFinder->getPosition();
                std::cerr << "        best position so far is: " << pos[0] << ", " << pos[1] << '\n';
                std::cerr << "        error in phys space: " << this->pointFinder->getError() << '\n';
                std::cerr << "        error code: " << ier << '\n';
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
	void getDstQuadCoord(size_t indx, const int offset[], double coords[]) const {

		// compute the index set of the cell and add the offset
		size_t cellIndsOffset[NDIMS_2D_TOPO];
		for (size_t j = 0; j < NDIMS_2D_TOPO; ++j) {
    		cellIndsOffset[j] = indx / this->dstCellDimProd[j] % this->dstCellDims[j];
    		cellIndsOffset[j] += offset[j];
  		}

  		// compute the low-corner flat index of the node coorresponding to this cell
  		size_t nodeIndx = 0;
		for (size_t j = 0; j < NDIMS_2D_TOPO; ++j) {
			nodeIndx += this->dstNodeDimProd[j] * cellIndsOffset[j];
		}

		// fill in the node's coordinates
		for (size_t j = 0; j < NDIMS_2D_PHYS; ++j) {
			coords[j] = this->dstGrdCoords[nodeIndx * NDIMS_2D_PHYS + j];
		}
	}

	/** 
	 * Extract the source cell coordinates from the grid
	 * @param indx nodal flat index
	 * @param offset displacement from the above node
	 * @param coords array of size NDIMS_2D_PHYS to be filled in 
	 */
	void getSrcQuadCoord(size_t indx, const int offset[], double coords[]) const {

		// compute the index set of the cell and add the offset
		size_t cellIndsOffset[NDIMS_2D_TOPO];
		for (size_t j = 0; j < NDIMS_2D_TOPO; ++j) {
    		cellIndsOffset[j] = indx / this->srcCellDimProd[j] % this->srcCellDims[j];
    		cellIndsOffset[j] += offset[j];
  		}

  		// compute the low-corner flat index of the node coorresponding to this cell
  		size_t nodeIndx = 0;
		for (size_t j = 0; j < NDIMS_2D_TOPO; ++j) {
			nodeIndx += this->srcNodeDimProd[j] * cellIndsOffset[j];
		}

		// fill in the node's coordinates
		for (size_t j = 0; j < NDIMS_2D_PHYS; ++j) {
			coords[j] = this->srcGrdCoords[nodeIndx * NDIMS_2D_PHYS + j];
		}
	}

};
 
#ifdef __cplusplus
extern "C" {
#endif

    int SgNodeInterp2D_new(SgNodeInterp2D_type** self);
                       
    int SgNodeInterp2D_del(SgNodeInterp2D_type** self);

    int SgNodeInterp2D_setDstGrid(SgNodeInterp2D_type** self, 
 	                                  const int dims[], const double** coords);

    int SgNodeInterp2D_setSrcGrid(SgNodeInterp2D_type** self, 
 	                                  const int dims[], 
 	                                  const int periodicity[],
 	                                  const double** coords);

    int SgNodeInterp2D_computeWeights(SgNodeInterp2D_type** self);

    int SgNodeInterp2D_apply(SgNodeInterp2D_type** self,
 	                             const double srcData[], double dstData[]);

#ifdef __cplusplus
}
#endif


#endif // SG_NODE_INTERP_2D_H