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
#include <limits>
 
struct SgFlowInterp2D_type {

    // the source grid
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

    // iterators
    SgBoxIterator_type* srcNodeIt;
    SgBoxIterator_type* srcCellIt;
    std::vector<SgBoxIterator_type*> srcEdgeIts;

    /** 
     * Constructor
     */
    SgFlowInterp2D_type() {
        this->srcNodeIt = NULL;
        this->srcCellIt = NULL;
        this->srcEdgeIts.resize(2);
        for (size_t i = 0; i < this->srcEdgeIts.size(); ++i) {
            this->srcEdgeIts[i] = NULL;
        }
    }

    /**
     * Destructor
     */
    ~SgFlowInterp2D_type() {
        if (this->srcNodeIt) {
            delete this->srcNodeIt;
        }
        for (size_t i = 0; i < this->srcEdgeIts.size(); ++i) {
            delete this->srcEdgeIts[i];
        }
    }

    /**
     * Set the destination 1D grid 
     * @param dims number of nodes in each in the single direction
     * @param coords coordinates (component, node)
     */
    void setDstGrid(const int dims[], const double** coords) {

        const double eps = 100 * std::numeric_limits<double>::epsilon();

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
                // apply small perturbation
                //double pertAngle = 2.0 * M_PI * double(i)/double(dstNumPoints - 1);
                //double pertAmpl = eps * (k + 1);
                //this->dstGrdCoords[i*NDIMS_2D_PHYS + k] = coords[k][i] + pertAmpl*cos(pertAngle);
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

        int srcCellDims[] = {dims[0] - 1, dims[1] - 1};

        int srcEdge0Dims[] = {dims[0] - 1, dims[1]};
        int srcEdge1Dims[] = {dims[0], dims[1] - 1};

        const int zeros[] = {0, 0};
        this->srcNodeIt = new SgBoxIterator_type(2, zeros, dims);
        this->srcCellIt = new SgBoxIterator_type(2, zeros, srcCellDims);
        this->srcEdgeIts[0] = new SgBoxIterator_type(2, zeros, srcEdge0Dims);
        this->srcEdgeIts[1] = new SgBoxIterator_type(2, zeros, srcEdge1Dims);

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
        int indsOffset[2];
        const int offset[] = {
            0, 0,
            0, 1,
            0, 0,
            1, 0
        };

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
                this->srcCellIt->getElement(srcIndx, inds);

                // iterate over the sides of the cell
                for (size_t side = 0; side < 4; ++side) {

                    indsOffset[0] = inds[0] + offset[side*2 + 0];
                    indsOffset[1] = inds[1] + offset[side*2 + 1];

                    // component; 0 = x edge, 1 = y edge
                    size_t k = side / 2;

                    // flat index of edge
                    size_t index = this->srcEdgeIts[k]->getIndex(indsOffset);

                    // update flux
                    dstData[dstIndx] += weightsXLoXHiYLoYHi[side]*srcData[k][index];
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

                // reset the number of intersection points to zero
                intersector.reset();

                // get the src quad's vertices
                this->getSrcQuadCoord(srcIndx, offset2D, srcQuadCoords);

                // set the quad's vertices in counterclock ordering
                intersector.setQuadPoints(srcQuadCoords);

                if (!intersector.checkIfOverlap()) {
                    // no chance, skip
                    continue;
                }

                double* points = NULL;
                int numPoints = 0;
                intersector.collectIntersectPoints(&numPoints, &points);

                if (numPoints >= 2) {
                    // must have at least 2 points
                    
                    // find the weight entries for this dstIndx
                    std::map<size_t, std::vector<std::pair<size_t, std::vector<double> > > >::iterator it;
                    it = this->weights.find(dstIndx);
                    
                    // add an entry if not present
                    if (it == this->weights.end()) {
                        std::vector< std::pair<size_t, std::vector<double> > > v;
                        std::pair<size_t, std::vector<std::pair<size_t, std::vector<double> > > > p(dstIndx, v);

                        this->weights.insert(p);
                        it = this->weights.find(dstIndx);
                    }

                    // initialize the weights
                    std::vector<double> w(4, 0.);
                    
                    weightCalc.setQuadPoints(srcQuadCoords);

                    // iterate over the segments
                    for (int i = 0; i < numPoints - 1; ++i) {
                        weightCalc.setLinePoints(&points[i*NDIMS_2D_PHYS]);
                        weightCalc.computeProjections();
                        
                        w[0] += weightCalc.getProjection(EDGE_LO_0);
                        w[1] += weightCalc.getProjection(EDGE_HI_0);
                        w[2] += weightCalc.getProjection(EDGE_LO_1);
                        w[3] += weightCalc.getProjection(EDGE_HI_1);
                    }
                    // done

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

        int cellInds[NDIMS_2D_TOPO];
        int inds[NDIMS_2D_TOPO];
        this->srcCellIt->getElement(indx, cellInds);

        // iterate over the quad's nodes
        for (size_t i = 0; i < 4; ++i) {

            // apply offset
            for (size_t j = 0; j < NDIMS_2D_TOPO; ++j) {
                inds[j] = cellInds[j] + offset[i*NDIMS_2D_TOPO + j];
            }

            // get the nodal flat index
            int nodeIndx = this->srcNodeIt->getIndex(inds);

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
