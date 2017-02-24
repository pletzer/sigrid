/**
 * A class that computes the interpolation weights for a cell centered field
 */
 
#ifndef SG_CONSERVE_INTERP_2D_H
#define SG_CONSERVE_INTERP_2D_H
 
#include "SgNdims.h"
#include "SgQuadIntersect.h"
#include "SgTriangulate.h"
#include "SgFindOverlappingCells2D.h"
#include "SgOctreePoints.h"
#include <vector>
#include <algorithm>
#include <map>
#include <set>
#include <cstdio> // size_t
#include <cmath>

 
struct SgConserveInterp2D_type {

	// the source grid
	int srcNodeDims[NDIMS_2D_TOPO];
	int srcCellDims[NDIMS_2D_TOPO];
	int srcNodeDimProd[NDIMS_2D_TOPO];
	int srcCellDimProd[NDIMS_2D_TOPO];
	std::vector<double> srcGrdCoords; // flat array (node, components)
	std::vector<int> periodicity; // 0=not periodic, 1=periodic
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

	// the interpolation weights (ratios of overlap areas to dst cell areas)
	// dst flat index -> {src flat index: weight}
	std::map<size_t, std::vector< std::pair<size_t, double> > > weights;

	// finds all the src cell under a dst cell
	SgFindOverlappingCells2D_type* srcCover;

	// iterators
    std::map<size_t, std::vector<std::pair<size_t, double> > >::const_iterator weightIt;
    std::vector<std::pair<size_t, double> >::const_iterator weightIt2;

    SgOctreePoints_type* srcOctreePtr;
    std::map<std::vector<size_t>, std::vector<size_t> > partition2SrcCellInds;
    size_t numLevels;

    /** 
     * Constructor
     */
    SgConserveInterp2D_type() {
    	this->srcCover = new SgFindOverlappingCells2D_type;
    	this->reset();
    	this->srcOctreePtr = 0;
    	this->numLevels = 0;
    }

    /**
     * Destructor
     */
    ~SgConserveInterp2D_type() {
    	delete this->srcCover;
    	if (this->srcOctreePtr) delete this->srcOctreePtr;
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

		this->srcCover->setSrcGrid(dims, periodicity, coords);

		this->srcNumPoints = 1;
		this->srcNumCells = 1;
		this->periodicity.resize(NDIMS_2D_TOPO);
		for (size_t j = 0; j < NDIMS_2D_TOPO; ++j) {
			this->srcNodeDims[j] = dims[j];
			this->srcCellDims[j] = dims[j] - 1;
			this->srcNumPoints *= dims[j];
			this->srcNumCells *= dims[j] - 1;
			this->periodicity[j] = periodicity[j];
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

  		// octree classification
  		int numCellsPerPartition = 100; // ideal number of cells per partition
  		size_t minNumCells = std::min(this->srcNumCells, this->dstNumCells);
  		this->numLevels = std::max(1, (int) floor(log((double) minNumCells/(double) numCellsPerPartition)/log(2.)/2.));
  		this->srcOctreePtr = new SgOctreePoints_type(this->numLevels, 2, this->srcGrdCoords);

  		// iterate over the src grid cells and attach each cell to a partition. A cell belongs to an octree partition iff
  		// at least on node belongs to the partition
		double point[NDIMS_2D_PHYS];
		int offset[] = {0, 0};
		size_t srcNodeIndx;
  		for (size_t srcIndx = 0; srcIndx < this->srcNumCells; ++srcIndx) {

  			// quadtree partitions to which this src cell belongs to
  			std::set< std::vector<size_t> > partitions;
  			std::vector<size_t> part(this->numLevels);

  			for (size_t j = 0; j < 4; ++j) {
				offset[0] = j % 2;
				offset[1] = j / 2;
				this->getSrcQuadCoord(srcIndx, offset, &srcNodeIndx, point);
				//bool indomain = this->srcOctreePtr->getKey(point, this->numLevels, part);
				// should always belong to the domain
				partitions.insert(part);
			}

			// iterate over the partitions and build the partitions to src cell indices map
			for (std::set< std::vector<size_t> >::const_iterator it = partitions.begin();
				 it != partitions.end(); ++it) {

				// do we have an entry for this partition?
				std::map<std::vector<size_t>, std::vector<size_t> >::iterator it2 = this->partition2SrcCellInds.find(*it);
				if (it2 != this->partition2SrcCellInds.end()) {
					// yes, append the src cell index
					it2->second.push_back(srcIndx);
				}
				else {
					// new entry 
					std::vector<size_t> v;
					v.reserve(numCellsPerPartition);
					v.push_back(srcIndx);
					this->partition2SrcCellInds.insert(std::pair< std::vector<size_t>, std::vector<size_t> >(*it, v));
				}
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
	 * Compute the interpolation weights
	 */
	void computeWeights();

	/**
	 * reset the iterators
	 */
	void reset() {
		this->weightIt = this->weights.begin();
		this->weightIt2 = this->weightIt->second.begin();		
	}

	/** 
	 * Do one step
	 * @return 0 (= continue) or 1 (= stop)
	 */
	int next() {
		// increment source counter
		this->weightIt2++;
		if (this->weightIt2 == this->weightIt->second.end()) {
			// reached end of src counter, increment dst counter
			this->weightIt++;
			// reset src counter
			this->weightIt2 = this->weightIt->second.begin();
			if (this->weightIt == this->weights.end()) {
				// reached end of dst counter
				// reset dst and src counters
				this->reset();
				// end of iteration 
				return 1;
			}
		}
		return 0;
	}

	/**
	 * Get the current src grid index, dst grid index and weight
	 * @param srcIndx src grid flat index (output)
	 * @param dstIndx dst grid flat index (output)
	 * @param weight interpolation weight (output)
	 */
	void get(int* srcIndx, int* dstIndx, double* weight) {
		*dstIndx = this->weightIt->first;
		*srcIndx = this->weightIt2->first;
		*weight = this->weightIt2->second;
	}

private:

    /** 
	 * Extract the destination cell coordinates from the grid
	 * @param indx cell flat index
	 * @param offset displacement from the above node
	 * @param nodeIndx grid node index (output)
	 * @param coords array of size NDIMS_2D_PHYS to be filled in 
	 */
	void getDstQuadCoord(size_t indx, const int offset[], size_t* nodeIndx, double coords[]) const {

		// compute the index set of the cell and add the offset
		size_t cellIndsOffset[NDIMS_2D_TOPO];
		for (size_t j = 0; j < NDIMS_2D_TOPO; ++j) {
    		cellIndsOffset[j] = indx / this->dstCellDimProd[j] % this->dstCellDims[j];
    		cellIndsOffset[j] += offset[j];
  		}

  		// compute the low-corner flat index of the node coorresponding to this cell
  		*nodeIndx = 0;
		for (size_t j = 0; j < NDIMS_2D_TOPO; ++j) {
			*nodeIndx += this->dstNodeDimProd[j] * cellIndsOffset[j];
		}

		// fill in the node's coordinates
		for (size_t j = 0; j < NDIMS_2D_PHYS; ++j) {
			coords[j] = this->dstGrdCoords[*nodeIndx * NDIMS_2D_PHYS + j];
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

    int SgConserveInterp2D_new(SgConserveInterp2D_type** self);
                       
    int SgConserveInterp2D_del(SgConserveInterp2D_type** self);

    int SgConserveInterp2D_setDstGrid(SgConserveInterp2D_type** self, 
 	                                  const int dims[], const double** coords);

    int SgConserveInterp2D_setSrcGrid(SgConserveInterp2D_type** self, 
 	                                  const int dims[], 
 	                                  const int periodicity[],
 	                                  const double** coords);

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
