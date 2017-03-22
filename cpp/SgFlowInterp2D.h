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
    int srcEdgeXDims[NDIMS_2D_TOPO];
    int srcEdgeYDims[NDIMS_2D_TOPO];
    int srcNodeDimProd[NDIMS_2D_TOPO];
    int srcCellDimProd[NDIMS_2D_TOPO];
    int srcEdgeXDimProd[NDIMS_2D_TOPO];
    int srcEdgeYDimProd[NDIMS_2D_TOPO];
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
    // {dst cell index: {src cell index: [wXLo, wXHi, wYLo, wYHi]}}
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
        for (size_t k = 0; k < NDIMS_2D_PHYS; ++k) {
            for (size_t i = 0; i < this->dstNumPoints; ++i) {
                this->dstGrdCoords[i*NDIMS_2D_PHYS + k] = coords[k][i];
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

        this->srcCellDims[0] = dims[0] - 1;
        this->srcCellDims[1] = dims[1] - 1;

        this->srcEdgeXDims[0] = dims[0];     // y direction
        this->srcEdgeXDims[1] = dims[1] - 1; // x dimension

        this->srcEdgeYDims[0] = dims[0] - 1; // y direction
        this->srcEdgeYDims[1] = dims[1];     // x dimension

        this->srcNodeDims[0] = dims[0];
        this->srcNodeDims[1] = dims[1];

        this->srcCellDimProd[1] = 1;
        this->srcCellDimProd[0] = this->srcCellDims[1];
 
        this->srcEdgeXDimProd[1] = 1;
        this->srcEdgeXDimProd[0] = this->srcEdgeXDims[1];

        this->srcEdgeYDimProd[1] = 1;
        this->srcEdgeYDimProd[0] = this->srcEdgeYDims[1];

        this->srcNodeDimProd[1] = 1;
        this->srcNodeDimProd[0] = this->srcNodeDims[1];

        this->srcGrdCoords.resize(NDIMS_2D_PHYS * this->srcNumPoints);
        // iterate over components
        for (size_t k = 0; k < NDIMS_2D_PHYS; ++k) {
            // iterate over nodes
            for (size_t i = 0; i < this->srcNumPoints; ++i) {
                this->srcGrdCoords[i*NDIMS_2D_PHYS + k] = coords[k][i];
            }
        }
    }

    /** 
     * Apply the interpolation weights to the src edge field
     * @param srcData source edge data with dimensions numCells x numNodes (x) or numNodes x numCells (y)
     * @param dstData destination cell data (output)
     */
    void apply(const double* srcData[], double dstData[]) {

        size_t inds[2];
        int* edgeDimProd;
        int* edgeDimProds[] = {this->srcEdgeXDimProd, this->srcEdgeYDimProd};

        // iterate over the dst segments
        for (std::map<size_t, std::vector<std::pair<size_t, std::vector<double> > > >::const_iterator 
            it = this->weights.begin(); it != this->weights.end(); ++it) {

            // dst cell (flat) index 
            size_t dstIndx = it->first;

            // initialize
            dstData[dstIndx] = 0;

            // iterate over the src cells intersected by this segment
            for (size_t i = 0; i < it->second.size(); ++i) {

                size_t srcIndx = it->second[i].first;

                // iterate over the two directions/components
                for (size_t k = 0; k < 2; ++k) {

                    // select the right dimensions (number of edges along x and y are different)
                    edgeDimProd = edgeDimProds[k];

                    // compute the src cell index set
                    inds[0] = srcIndx / edgeDimProd[0] % edgeDimProd[0];
                    inds[1] = srcIndx / edgeDimProd[1] % edgeDimProd[1];
                
                    // iterate over the lo/hi edges
                    for (size_t loHi = 0; loHi < 2; ++loHi) {
                    
                        // aply offset (or not)
                        inds[k] += loHi;

                        // get the interpolation weight
                        double weight = it->second[i].second[2*k + loHi];

                        // compute the flat index for the src edge
                        size_t srcEdgeIndx = edgeDimProd[0]*inds[0] + edgeDimProd[1]*inds[1];

                        // update the flux
                        dstData[dstIndx] += weight * srcData[k][srcEdgeIndx];
                    }
                }
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
        const int offset1D[] = {0, 1};
        const int offset2D[] = {0, 0,
                                1, 0,
                                1, 1,
                                0, 1};

        SgQuadLineIntersect_type intersector;

        // iterate over the dst segments
        for (size_t dstIndx = 0; dstIndx < this->dstNumCells; ++dstIndx) {

            // get the start/end points
            this->getDstLineCoord(dstIndx, offset1D, dstLineCoords);

            intersector.setLinePoints(dstLineCoords);

            // iterate over the src quads
            for (size_t srcIndx = 0; srcIndx < this->srcNumCells; ++srcIndx) {

                // get the src quad's vertices
                this->getSrcQuadCoord(srcIndx, offset2D, srcNodeInds, srcQuadCoords);

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
                        it = this->weights.find(dstIndx);
                    }

                    double* pA = &points[0*NDIMS_2D_PHYS];
                    double* pB = &points[1*NDIMS_2D_PHYS];

                    std::vector<double> w(4);
                    w[0] = (pB[0] - pA[0])*(1.0 - 0.5*(pB[1] + pA[1])); // x low side
                    w[1] = (pB[0] - pA[0])*(0.0 + 0.5*(pB[1] + pA[1])); // x high side
                    w[2] = (1.0 - 0.5*(pB[0] + pA[0]))*(pB[1] - pA[1]); // y low side
                    w[3] = (0.0 + 0.5*(pB[0] + pA[0]))*(pB[1] - pA[1]); // y high side

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
     * @param offset displacement from the above node, either 0 or 1
     * @param coords array of size NDIMS_2D_PHYS to be filled in 
     */
    void getDstLineCoord(size_t indx, const int offset[], double coords[]) const {

        // iterate over the 2 nodes
        for (size_t i = 0; i < 2; ++i) {
            size_t nodeIndx = indx + offset[i];
            // fill in the node's coordinates
            for (size_t j = 0; j < NDIMS_2D_PHYS; ++j) {
                coords[i*NDIMS_2D_PHYS + j] = this->dstGrdCoords[nodeIndx*NDIMS_2D_PHYS + j];
            }
        }
    }

    /** 
     * Extract the source cell coordinates from the grid
     * @param indx cell flat index
     * @param offset list of displacements from the above node, one per node, going counterclockwise
     * @param nodeIndx array of size 4 for the grid node indices (output)
     * @param coords array of size 4*NDIMS_2D_PHYS to be filled in 
     */
    void getSrcQuadCoord(size_t indx, const int offset[], size_t nodeIndx[], double coords[]) const {

        size_t cellInds[NDIMS_2D_TOPO];
        for (size_t j = 0; j < NDIMS_2D_TOPO; ++j) {
            cellInds[j] = indx / this->srcCellDimProd[j] % this->srcCellDims[j];
        }

        // iterate over the quad's nodes
        for (size_t i = 0; i < 4; ++i) {

            // compute the low-corner flat index of the node coorresponding to this cell
            nodeIndx[i] = 0;
            for (size_t j = 0; j < NDIMS_2D_TOPO; ++j) {
                nodeIndx[i] += this->srcNodeDimProd[j]*(cellInds[j] + offset[i*NDIMS_2D_TOPO + j]);
            }

            // fill in the node's coordinates
            for (size_t j = 0; j < NDIMS_2D_PHYS; ++j) {
                coords[i*NDIMS_2D_PHYS + j] = this->srcGrdCoords[nodeIndx[i]*NDIMS_2D_PHYS + j];
            }
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

    int SgFlowInterp2D_apply(SgFlowInterp2D_type** self,
                              const double* srcData[], double dstData[]);

#ifdef __cplusplus
}
#endif


#endif // SG_FLOW_INTERP_2D_H
