/**
 * A class that computes the interpolation weights for a cell centered field
 */
 
#include "SgConserveInterp2D.h"
#include <cmath>
#include <iostream>

void SgConserveInterp2D_type::computeWeights() {

	// cache the edge to edge intersections
	std::map< std::pair< std::vector<size_t>, std::vector<size_t> >, std::vector<double> > cacheEdgeX;

	// cache edges that are known not to intersect
	std::set< std::pair< std::vector<size_t>, std::vector<size_t> > > cacheEdgeNoX;

	// cache the dst node in src cell points
	std::map< std::vector<size_t>, std::vector<double> > cacheDstNodeInSrcCell;

	// cache the src node in dst cell points
	std::map< std::vector<size_t>, std::vector<double> > cacheSrcNodeInDstCell;

	int numIntersectPoints;
	std::vector<double> point(NDIMS_2D_PHYS); // edge to edge intersection point
	SgQuadIntersect_type intersector;
	double dstQuadCoords[NDIMS_2D_PHYS*4]; // four nodes
	double srcQuadCoords[NDIMS_2D_PHYS*4]; // four nodes
	int offsets[] = {0, 0, 
	                 1, 0,
	                 1, 1,
	                 0, 1};
	size_t dstNodeInds[4];
	size_t srcNodeInds[4];

	// iterate over the dst cells
	for (size_t dstIndx = 0; dstIndx < this->dstNumCells; ++dstIndx) {

		std::vector<size_t> dstCellIndxSrcNodeIndx(2);
		dstCellIndxSrcNodeIndx[0] = dstIndx;

		// fetch the dst cell quad coordinates, counterclockwise
		for (size_t k = 0; k < 4; ++k) {
			this->getDstQuadCoord(dstIndx, &offsets[k*2], &dstNodeInds[k], &dstQuadCoords[k*NDIMS_2D_PHYS]);
		}

		// compute the dst cell area
		SgTriangulate_type dstTriangulator(4, dstQuadCoords);
		double dstArea = dstTriangulator.getConvexHullArea();

		// src index to fractional area
		std::vector< std::pair<size_t, double> > indWght;

		indWght.reserve(100);

		// find the partitions of the dst cell with respect to the src grid
		std::set<std::vector<size_t> > partitions;
		std::vector<size_t> part(this->numLevels, 0);
		for (size_t i = 0; i < 4; ++i) {
			this->srcOctreePtr->getKey(&dstQuadCoords[i*NDIMS_2D_PHYS], this->numLevels, part);
			// always add since dst cell may contain src grid
			partitions.insert(part);
			// MIGHT NEED TO ADD MORE PARTITIONS WHEN THE DST CELL >> SRC CELL!!! (MIGHT NEED TO ADD 
			// ALL THE PARTITIONS INBETWEEN NODES)
		}

		// collect all the src cells in the dst cell partitions
		std::set<size_t> srcCells;
		for (std::set< std::vector<size_t> >::const_iterator it = partitions.begin();
			it != partitions.end(); ++it) {
			const std::vector<size_t>& part = *it;
			// all the src cells in this partition
			std::map< std::vector<size_t>, std::vector<size_t> >::const_iterator it2 = this->partition2SrcCellInds.find(part);
			if (it2 != this->partition2SrcCellInds.end()) {
				// collect the src cell indices
				for (size_t i = 0; i < it2->second.size(); ++i) {
					srcCells.insert(it2->second[i]);
				}
			}
		}

		// iterate over the src cells
		for (std::set<size_t>::const_iterator it = srcCells.begin(); it != srcCells.end(); ++it) {

			size_t srcIndx = *it;

			std::vector<size_t> dstNodeIndxSrcCellIndx(2);
			std::vector<size_t> dstCellIndxSrcNodeIndx(2);

			// fetch the src cell quad coordinates
			for (size_t k = 0; k < 4; ++k) {
				this->getSrcQuadCoord(srcIndx, &offsets[k*2], &srcNodeInds[k], &srcQuadCoords[k*NDIMS_2D_PHYS]);
			}

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
            	std::map< std::vector<size_t>, std::vector<double> >::const_iterator 
            		    it1 = cacheDstNodeInSrcCell.find(dstNodeIndxSrcCellIndx);
            	if (it1 != cacheDstNodeInSrcCell.end()) {
            		// we've already checked this
            		intersector.addPoint(pA);
            	}
            	else {
            		// first time
            		bool ret = intersector.isPointInCell(dstCoordA, srcQuadCoords);
            		if (ret) {
            			intersector.addPoint(pA);
            			std::pair< std::vector<size_t>, std::vector<double> > p(dstNodeIndxSrcCellIndx, pA);
            			cacheDstNodeInSrcCell.insert(p);
            		}
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
            		std::map< std::vector<size_t>, std::vector<double> >::const_iterator 
            		    	it2 = cacheSrcNodeInDstCell.find(dstCellIndxSrcNodeIndx);
            		if (it2 != cacheSrcNodeInDstCell.end()) {
            			// we've already checked this
            			intersector.addPoint(qA);
            		}
            		else {
            			// first time
            			bool ret = intersector.isPointInCell(srcCoordA, dstQuadCoords);
            			if (ret) {
            				intersector.addPoint(qA);
            				std::pair< std::vector<size_t>, std::vector<double> > p(dstCellIndxSrcNodeIndx, qA);
            				cacheSrcNodeInDstCell.insert(p);
            			}
            		}

                	std::pair< std::vector<size_t>, std::vector<size_t> > dstSrcEdges(dstE, srcE);

                	if (cacheEdgeNoX.find(dstSrcEdges) != cacheEdgeNoX.end()) {
                		// we already know there is no intersection, move on
                		continue;
                	}

                	std::map< std::pair< std::vector<size_t>, std::vector<size_t> >, std::vector<double> >::const_iterator 
                		    it = cacheEdgeX.find(dstSrcEdges);
                	if (cacheEdgeX.find(dstSrcEdges) != cacheEdgeX.end()) {
                		// we already know what the intersection point is
                		intersector.addPoint(it->second);
                	}
                	else {
                		// let's see if there is an intersection
                		int ret = intersector.collectEdgeToEdgeIntersectionPoints(dstCoordA, dstCoordB,
                                                                                  srcCoordA, srcCoordB,
                                                                                  &point[0]);

                		if (ret == 0) {
                			// no intersection
                			cacheEdgeNoX.insert(dstSrcEdges);
                		}
                		else if (ret == 1) {
                			// intersection, cache result for subsequent use
                			std::pair< std::pair<std::vector<size_t>, std::vector<size_t> >, std::vector<double> > p(dstSrcEdges, point);
                			cacheEdgeX.insert(p);
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


extern "C"
int SgConserveInterp2D_new(SgConserveInterp2D_type** self) {
 	*self = new SgConserveInterp2D_type();
	return 0;
}
      
extern "C"                   
int SgConserveInterp2D_del(SgConserveInterp2D_type** self) {
 	delete *self;
 	return 0;
}

extern "C"
int SgConserveInterp2D_setDstGrid(SgConserveInterp2D_type** self, 
 	                              const int dims[],
 	                              const double** coords) {
	(*self)->setDstGrid(dims, coords);
	return 0;
}

extern "C"
int SgConserveInterp2D_setSrcGrid(SgConserveInterp2D_type** self, 
 	                              const int dims[], 
 	                              const int periodicity[],
 	                              const double** coords) {
	(*self)->setSrcGrid(dims, periodicity, coords);
	return 0;
}

extern "C"
int SgConserveInterp2D_computeWeights(SgConserveInterp2D_type** self) {
	(*self)->computeWeights();
	return 0;
}

extern "C"
int SgConserveInterp2D_apply(SgConserveInterp2D_type** self,
 	                         const double srcData[], double dstData[]) {
	(*self)->apply(srcData, dstData);
	return 0;
}

extern "C"
int SgConserveInterp2D_reset(SgConserveInterp2D_type** self) {
	(*self)->reset();
	return 0;
}

extern "C"
int SgConserveInterp2D_next(SgConserveInterp2D_type** self) {
	return (*self)->next();
}

extern "C"
int SgConserveInterp2D_get(SgConserveInterp2D_type** self,
                           int* srcIndx, int* dstIndx, double* weight) {
	(*self)->get(srcIndx, dstIndx, weight);
	return 0;
}
