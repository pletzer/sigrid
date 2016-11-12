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

void setDstGrid(const int dims[], const double** coords) {

	this->dstNumPoints = 1;
	this->dstNumCells = 1;
	for (size_t j = 0; j < NDIMS_2D_TOPO; ++j) {
		this->dstNodeDims[j] = dims[j];
		this->dstCellDims[j] = dims[j] - 1;
		this->dstNumPoints *= dims[j];
		this->dstNumCells *= dims[j] - 1;
	}
	this->dstCellDimProd[0] = 1;
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

void setSrcGrid(const int dims[], const double** coords) {

	this->srcNumPoints = 1;
	this->srcNumCells = 1;
	for (size_t j = 0; j < NDIMS_2D_TOPO; ++j) {
		this->srcNodeDims[j] = dims[j];
		this->srcCellDims[j] = dims[j] - 1;
		this->srcNumPoints *= dims[j];
		this->srcNumCells *= dims[j] - 1;
	}
	this->srcCellDimProd[0] = 1;
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

void apply(double* dstData, const double* srcData) {

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

const double* getDstQuadCoord(size_t indx, int offset[]) const {
	for (size_t j = 0; j < NDIMS_2D_TOPO; ++j) {
		indx += this->dstNodeDimProd[j] * offset[j];
	}
	return &this->dstGrdCoords[indx];
}

const double* getSrcQuadCoord(size_t indx, int offset[]) const {
	for (size_t j = 0; j < NDIMS_2D_TOPO; ++j) {
		indx += this->srcNodeDimProd[j] * offset[j];
	}
	return &this->srcGrdCoords[indx];
}

void computeWeights() {
	int numIntersectPoints;
	double** intersectPoints;
	SgQuadIntersect_type intersector;
	const double* dstQuadCoords[] = {NULL, NULL, NULL, NULL};
	const double* srcQuadCoords[] = {NULL, NULL, NULL, NULL};
	int offset[] = {0, 0};
	// iterate over the dst cells
	for (size_t dstIndx = 0; dstIndx < this->dstNumCells; ++dstIndx) {

		offset[0] = 0; offset[1] = 0;
		dstQuadCoords[0] = this->getDstQuadCoord(dstIndx, offset);
		offset[0] = 1; offset[1] = 0;
		dstQuadCoords[1] = this->getDstQuadCoord(dstIndx, offset);
		offset[0] = 1; offset[1] = 1;
		dstQuadCoords[2] = this->getDstQuadCoord(dstIndx, offset);
		offset[0] = 0; offset[1] = 1;
		dstQuadCoords[3] = this->getDstQuadCoord(dstIndx, offset);

		// iterate over the src cells
		std::vector< std::pair<size_t, double> > indWght;
		for (size_t srcIndx = 0; srcIndx < this->srcNumCells; ++srcIndx) {

			offset[0] = 0; offset[1] = 0;
			srcQuadCoords[0] = this->getSrcQuadCoord(dstIndx, offset);
			offset[0] = 1; offset[1] = 0;
			srcQuadCoords[1] = this->getSrcQuadCoord(dstIndx, offset);
			offset[0] = 1; offset[1] = 1;
			srcQuadCoords[2] = this->getSrcQuadCoord(dstIndx, offset);
			offset[0] = 0; offset[1] = 1;
			srcQuadCoords[3] = this->getSrcQuadCoord(dstIndx, offset);

			intersector.setQuadPoints(dstQuadCoords, srcQuadCoords);
			intersector.collectIntersectPoints(&numIntersectPoints, intersectPoints);
			if (numIntersectPoints >= 3) {
				// must be able to build at least one triangle
				SgTriangulate_type triangulator(numIntersectPoints, (const double**) intersectPoints);
				double area = triangulator.getConvexHullArea();
				indWght.push_back(std::pair<size_t, double>(srcIndx, area));
			}
		}
		if (indWght.size() > 0) {
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

#ifdef __cplusplus
}
#endif


#endif // SG_CONSERVE_INTERP_2D_H