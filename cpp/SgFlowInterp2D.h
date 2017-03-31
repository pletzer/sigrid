/**
 * A class that computes the interpolation weights for a face centered field in 2D
 */
 
#ifndef SG_FLOW_INTERP_2D_H
#define SG_FLOW_INTERP_2D_H
 
#include "SgNdims.h"
#include "SgQuadLineIntersect.h"
#include "SgQuadLineFlows.h"
#include "SgBoxIterator.h"
#include <vector>
#include <algorithm>
#include <map>
#include <cstdio> // size_t
#include <cmath>
#include <iostream>
 
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
    int dstCellDims[NDIMS_1D_TOPO];
    std::vector<double> dstGrdCoords; // flat array (node, components)
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

        size_t dstNumPoints = 1;
        this->dstNumCells = 1;
        for (size_t j = 0; j < NDIMS_1D_TOPO; ++j) {
            this->dstCellDims[j] = dims[j] - 1;
            dstNumPoints *= dims[j];
            this->dstNumCells *= dims[j] - 1;
        }

        this->dstGrdCoords.resize(NDIMS_2D_PHYS * dstNumPoints);
        for (size_t k = 0; k < NDIMS_2D_PHYS; ++k) {
            for (size_t i = 0; i < dstNumPoints; ++i) {
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

        int inds[2];
        const int zeros[] = {0, 0};

        SgBoxIterator_type cellIt(2, zeros, this->srcCellDims);
        SgBoxIterator_type edgeXIt(2, zeros, this->srcEdgeXDims);
        SgBoxIterator_type edgeYIt(2, zeros, this->srcEdgeYDims);
        SgBoxIterator_type edgeIts[] = {edgeXIt, edgeYIt};

        // iterate over components
        for (size_t k = 0; k < 2; ++k) {

            // iterate over the dst segments
            for (std::map<size_t, std::vector<std::pair<size_t, std::vector<double> > > >::const_iterator 
                 it = this->weights.begin(); it != this->weights.end(); ++it) {

                // dst cell (flat) index 
                size_t dstIndx = it->first;

                // initialize
                dstData[dstIndx] = 0.0;

                // iterate over the src cells intersected by this segment
                for (size_t i = 0; i < it->second.size(); ++i) {

                    // weights: xLo, xHi, yLo, yHi
                    const std::vector<double>& weightsXLoXHiYLoYHi = it->second[i].second;

                    // src cell index
                    size_t srcIndx = it->second[i].first;

                    // cell indices
                    cellIt.getElement(srcIndx, inds);

                    // iterate over the low/high sides
                    for (size_t loHi = 0; loHi < 2; ++loHi) {

                        // note that we're using y, z order for indexing
                        // so xHi is in the y direction (index 0) and 
                        // yHi is in the x direction (index 1)
                        inds[k] += loHi;

                        // convert the index set to a flat edge iterator index
                        int index = edgeIts[k].getIndex(inds);

                        // get the weight 
                        double weight = weightsXLoXHiYLoYHi[2*k + loHi];

                        // update
                        dstData[dstIndx] += weight*srcData[k][index];
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
        const int offset1D[] = {0, 1};
        const int offset2D[] = {0, 0,
                                1, 0,
                                1, 1,
                                0, 1};

        SgQuadLineIntersect_type intersector;
        SgQuadLineFlows_type weightCalc;

        // iterate over the dst segments
        for (size_t dstIndx = 0; dstIndx < this->dstNumCells; ++dstIndx) {

            // get the start/end points
            this->getDstLineCoord(dstIndx, offset1D, dstLineCoords);

            intersector.setLinePoints(dstLineCoords);

            // iterate over the src quads
            for (size_t srcIndx = 0; srcIndx < this->srcNumCells; ++srcIndx) {

                // get the src quad's vertices
                this->getSrcQuadCoord(srcIndx, offset2D, srcQuadCoords);

                // reset the number of intersection points to zero
                intersector.reset();

                // set the quad's vertices in counterclock ordering
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

                    weightCalc.setQuadPoints(srcQuadCoords);
                    weightCalc.setLinePoints(points);
                    weightCalc.computeProjections();

                    std::vector<double> w(4);
                    w[0] = weightCalc.getProjection(EDGE_LO_0);
                    w[1] = weightCalc.getProjection(EDGE_HI_0);
                    w[2] = weightCalc.getProjection(EDGE_LO_1);
                    w[3] = weightCalc.getProjection(EDGE_HI_1);

                    std::pair<size_t, std::vector<double> > p(srcIndx, w);
                    it->second.push_back(p);
                }
            }
        }

    }

    /**
     * Write debug information
     */
    void debug() const {
        std::cerr << "SgFlowInterp2D:\n";
        for (std::map<size_t, std::vector<std::pair<size_t, std::vector<double> > > >::const_iterator 
            it = this->weights.begin(); it != this->weights.end(); ++it) {

            size_t dstIndx = it->first;
            std::cerr << "dst segment index: " << dstIndx << ":\n";
            const std::vector< std::pair<size_t, std::vector<double> > >& srcIndx2Weights = it->second;

            for (size_t i = 0; i < srcIndx2Weights.size(); ++i) {
                std::cerr << "\tsrc cell index: " << srcIndx2Weights[i].first << " -> weights (xLo, xHi, yLo, yHi): ";
                for (size_t j = 0; j < srcIndx2Weights[i].second.size(); ++j) {
                    std::cerr << srcIndx2Weights[i].second[j] << ", ";
                }
                std::cerr << '\n';
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
    void getSrcQuadCoord(size_t indx, const int offset[], double coords[]) const {

        size_t cellInds[NDIMS_2D_TOPO];
        for (size_t j = 0; j < NDIMS_2D_TOPO; ++j) {
            cellInds[j] = indx / this->srcCellDimProd[j] % this->srcCellDims[j];
        }

        // iterate over the quad's nodes
        for (size_t i = 0; i < 4; ++i) {

            // compute the low-corner flat index of the node coorresponding to this cell
            size_t nodeIndx = 0;
            for (size_t j = 0; j < NDIMS_2D_TOPO; ++j) {
                nodeIndx += this->srcNodeDimProd[j]*(cellInds[j] + offset[i*NDIMS_2D_TOPO + j]);
            }

            // fill in the node's coordinates
            for (size_t j = 0; j < NDIMS_2D_PHYS; ++j) {
                coords[i*NDIMS_2D_PHYS + j] = this->srcGrdCoords[nodeIndx*NDIMS_2D_PHYS + j];
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

    int SgFlowInterp2D_debug(SgFlowInterp2D_type** self);

#ifdef __cplusplus
}
#endif


#endif // SG_FLOW_INTERP_2D_H
