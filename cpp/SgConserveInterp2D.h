/**
 * A class that computes the interpolation weights for a cell centered field
 */
 
#ifndef SG_CONSERVE_INTERP_2D_H
#define SG_CONSERVE_INTERP_2D_H
 
#include "SgNdims.h"
#include "SgQuadIntersect.h"
#include "SgTriangulate.h"
#include "SgFindOverlappingCells2D.h"
#include <vector>
#include <algorithm>
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

    /** 
     * Constructor
     */
    SgConserveInterp2D_type() {
    	this->srcCover = new SgFindOverlappingCells2D_type;
    	this->reset();
    }

    /**
     * Destructor
     */
    ~SgConserveInterp2D_type() {
    	delete this->srcCover;
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
	void computeWeights() {

		// cache edges that are known not to intersect
		std::set< std::pair< std::vector<size_t>, std::vector<size_t> > > cacheEdgeNoX;

		int numIntersectPoints;
		std::vector<double> point(NDIMS_2D_PHYS); // edge to edge intersection point
		SgQuadIntersect_type intersector;
		double dstQuadCoords[NDIMS_2D_PHYS*4]; // four nodes
		double srcQuadCoords[NDIMS_2D_PHYS*4]; // four nodes
		int offset[] = {0, 0};
		size_t dstNodeInds[4];
		size_t srcNodeInds[4];

		// iterate over the dst cells
		for (size_t dstIndx = 0; dstIndx < this->dstNumCells; ++dstIndx) {

			std::vector<size_t> dstCellIndxSrcNodeIndx(2);
			dstCellIndxSrcNodeIndx[0] = dstIndx;

			offset[0] = 0; offset[1] = 0;
			this->getDstQuadCoord(dstIndx, offset, &dstNodeInds[0], &dstQuadCoords[0*NDIMS_2D_PHYS]);
			offset[0] = 1; offset[1] = 0;
			this->getDstQuadCoord(dstIndx, offset, &dstNodeInds[1], &dstQuadCoords[1*NDIMS_2D_PHYS]);
			offset[0] = 1; offset[1] = 1;
			this->getDstQuadCoord(dstIndx, offset, &dstNodeInds[2], &dstQuadCoords[2*NDIMS_2D_PHYS]);
			offset[0] = 0; offset[1] = 1;
			this->getDstQuadCoord(dstIndx, offset, &dstNodeInds[3], &dstQuadCoords[3*NDIMS_2D_PHYS]);

			// compute the dst cell area
			SgTriangulate_type dstTriangulator(4, dstQuadCoords);
			double dstArea = dstTriangulator.getConvexHullArea();

			// src index to fractional area
			std::vector< std::pair<size_t, double> > indWght;

			indWght.reserve(100);
			for (size_t srcIndx = 0; srcIndx < this->srcNumCells; ++srcIndx) {

				std::vector<size_t> dstNodeIndxSrcCellIndx(2);
				std::vector<size_t> dstCellIndxSrcNodeIndx(2);

				offset[0] = 0; offset[1] = 0;
				this->getSrcQuadCoord(srcIndx, offset, &srcNodeInds[0], &srcQuadCoords[0*NDIMS_2D_PHYS]);
				offset[0] = 1; offset[1] = 0;
				this->getSrcQuadCoord(srcIndx, offset, &srcNodeInds[1], &srcQuadCoords[1*NDIMS_2D_PHYS]);
				offset[0] = 1; offset[1] = 1;
				this->getSrcQuadCoord(srcIndx, offset, &srcNodeInds[2], &srcQuadCoords[2*NDIMS_2D_PHYS]);
				offset[0] = 0; offset[1] = 1;
				this->getSrcQuadCoord(srcIndx, offset, &srcNodeInds[3], &srcQuadCoords[3*NDIMS_2D_PHYS]);

				intersector.reset();
				intersector.setQuadPoints(dstQuadCoords, srcQuadCoords);

				if(!intersector.checkIfBoxesOverlap()) {
					// no chance
					continue;
				}
        
        		// iterate over dst cell edges and nodes
        		for (size_t i = 0; i < 4; ++i) {
            		// the edge of the first quad
            		size_t iA = i;
            		size_t iB = (i + 1) % 4;
            		size_t dstNodeA = dstNodeInds[iA];
            		size_t dstNodeB = dstNodeInds[iB];
            		std::vector<size_t> dstE(2); 
            		dstE[0] = dstNodeA; 
            		dstE[1] = dstNodeB;
            		std::sort(dstE.begin(), dstE.end());
            		double* dstCoordA = &dstQuadCoords[iA*NDIMS_2D_PHYS];
            		double* dstCoordB = &dstQuadCoords[iB*NDIMS_2D_PHYS];
            		std::vector<double> pA(dstCoordA, dstCoordA + NDIMS_2D_PHYS);

            		dstNodeIndxSrcCellIndx[0] = dstNodeA;
            		dstNodeIndxSrcCellIndx[1] = srcIndx;
            		bool ret = intersector.isPointInCell(dstCoordA, srcQuadCoords);
            		if (ret) {
            			intersector.addPoint(pA);
            		}            		

            		// iterate over the src cell edges and nodes
            		for (size_t j = 0; j < 4; ++j) {
                		// the edges of the second quad
                		size_t jA = j;
                		size_t jB = (j + 1) % 4;
                		size_t srcNodeA = srcNodeInds[jA];
                		size_t srcNodeB = srcNodeInds[jB];
                		std::vector<size_t> srcE(2); 
                		srcE[0] = srcNodeA; 
                		srcE[1] = srcNodeB;
                		std::sort(srcE.begin(), srcE.end());
            			double* srcCoordA = &srcQuadCoords[jA*NDIMS_2D_PHYS];
            			double* srcCoordB = &srcQuadCoords[jB*NDIMS_2D_PHYS];
                		std::vector<double> qA(srcCoordA, srcCoordA + NDIMS_2D_PHYS);

                		dstCellIndxSrcNodeIndx[0] = dstIndx;
                		dstCellIndxSrcNodeIndx[1] = srcNodeA;

            			bool ret = intersector.isPointInCell(srcCoordA, dstQuadCoords);
            			if (ret) {
            				intersector.addPoint(qA);
            			}

                		std::pair< std::vector<size_t>, std::vector<size_t> > dstSrcEdges(dstE, srcE);

                		{
                			// let's see if there is an intersection
                			int ret = intersector.collectEdgeToEdgeIntersectionPoints(dstCoordA, dstCoordB,
                                                                                      srcCoordA, srcCoordB,
                                                                                      &point[0]);

                			if (ret == 0) {
                				// no intersection
                				cacheEdgeNoX.insert(dstSrcEdges);
                			}
                		}

          			}
        		}

        		// ready to collect all the interesction points
				const std::vector<double>& pts = intersector.getIntersectionPoints();

				numIntersectPoints = pts.size() / NDIMS_2D_PHYS;
				if (numIntersectPoints >= 3) {
					// must be able to build at least one triangle
					SgTriangulate_type triangulator(numIntersectPoints, &pts[0]);
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