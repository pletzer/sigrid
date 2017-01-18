/**
 * A class that finds all intersection points between two grids
 */
 
#ifndef SG_STRUCT_GRID_2D_INTERSECT_H
#define SG_STRUCT_GRID_2D_INTERSECT_H
 
#include <vector>
#include <cstdio> // size_t
#include <cmath>
#include <limits>
#include <algorithm>
#include "SgLinearSolve.h"
#include "SgNdims.h"
 
struct SgStructGrid2DIntersect_type {

    // linear solver
    SgLinearSolve_type* slvr;

    // grid 1
    int nodeDims1[NDIMS_2D_TOPO];
    int nodeDimsProd1[NDIMS_2D_TOPO]; // to go from node flat index to index set
    int cellDims1[NDIMS_2D_TOPO];
    int cellDimsProd1[NDIMS_2D_TOPO]; // to go from cell flat index to index set
    int numNodes1; 
    int numCells1;
    std::vector<double> coords1; // flat array (node, components)

    // grid 2
    int nodeDims2[NDIMS_2D_TOPO];
    int nodeDimsProd2[NDIMS_2D_TOPO]; // to go from node flat index to index set
    int cellDims2[NDIMS_2D_TOPO];
    int cellDimsProd2[NDIMS_2D_TOPO]; // to go from cell flat index to index set
    int numNodes2; 
    int numCells2;
    std::vector<double> coords2; // flat array (node, components)

    // (cellIndx1, cellIndx2) -> points
    std::map<int[2], std::vector<double> >  cellIndx1Indx2_points;

    // a tolerance to determine whether a two edges intersect
    // or a point is within a triangle
    double tol;

    /**
     * Constructor
     */
    SgStructGrid2DIntersect_type() {

        // tolerance for floating point comparisons
        this->tol = 1.e-12;

        // we're in 2D
        this->slvr = new SgLinearSolve_type(2, 2);
    }

    /**
     * Destructor
     */
    ~SgStructGrid2DIntersect_type() {
        delete this->slvr;
    }

    /**
     * Set the grids
     * @param nodeDims1 number of nodes in each direction (grid 1)
     * @param coords1 coordinates (grid 1)
     * @param nodeDims2 number of nodes in each direction (grid 2)
     * @param coords2 coordinates (grid 2)
     */
    void setGrids(int nodeDims1, const double** coords1,
                  int nodeDims2, const double** coords2) {

        // grid 1
        this->numNodes1 = 1;
        this->numCells1 = 1;
        this->numNodes2 = 1;
        this->numCells2 = 1;
        this->nodeDimsProd[NDIMS_2D_TOPO - 1] = 1;
        this->cellDimsProd[NDIMS_2D_TOPO - 1] = 1;
        for (size_t j = 0; j < NDIMS_2D_TOPO; ++j) {
            this->numNodes1 *= nodeDims1[j];
            this->nodeDims1[j] = this->nodeDims1[j];
            this->numCells1 *= nodeDims1[j] - 1;
            this->cellDims1[j] = this->nodeDims1[j] - 1;

            this->numNodes2 *= nodeDims2[j];
            this->nodeDims2[j] = this->nodeDims2[j];
            this->numCells2 *= nodeDims2[j] - 1;
            this->cellDims2[j] = this->nodeDims1[j] - 1;
        }

        this->nodeDimsProd1[NDIMS_2D_TOPO - 1] = 1;
        this->cellDimsProd1[NDIMS_2D_TOPO - 1] = 1;  
        this->nodeDimsProd2[NDIMS_2D_TOPO - 1] = 1;
        this->cellDimsProd2[NDIMS_2D_TOPO - 1] = 1;  
        for (int j = NDIMS_2D_TOPO - 2; j >= 0; --j) {
            // last index varies fastest
            this->nodeDimsProd1[j] = this->nodeDimsProd1[j + 1] * this->nodeDims1[j + 1];
            this->cellDimsProd1[j] = this->cellDimsProd1[j + 1] * this->cellDims1[j + 1];
            this->nodeDimsProd2[j] = this->nodeDimsProd2[j + 1] * this->nodeDims2[j + 1];
            this->cellDimsProd2[j] = this->cellDimsProd2[j + 1] * this->cellDims2[j + 1];
        }
    }

    void collectIntersectionPoints() {

        // edge 1, egde 2 -> points
        //std::map<int[4], std::vector<double> > cacheEdgeIntersectPoints;

        int numIntersectionPoints;
        std::vector<int> cellIndx12(2);

        double quadCoords1[NDIMS_2D_PHYS*4]; // four nodes
        double quadCoords2[NDIMS_2D_PHYS*4]; // four nodes
        int offset[] = {0, 0};

        for (int cellIndx1 = 0; cellIndx1 < this->numCells1; ++cellIndx1) {

            std::vector<int> nodeInds1(4);
            // build quad 1
            offset[0] = 0; offset[1] = 0;
            this->getQuadCoords1(cellIndx1, offset, &nodeInds1[0], &quadCoords1[0*NDIMS_2D_PHYS]);
            offset[0] = 1; offset[1] = 0;
            this->getQuadCoords1(cellIndx1, offset, &nodeInds1[1], &quadCoords1[1*NDIMS_2D_PHYS]);
            offset[0] = 1; offset[1] = 1;
            this->getQuadCoords1(cellIndx1, offset, &nodeInds1[2], &quadCoords1[2*NDIMS_2D_PHYS]);
            offset[0] = 0; offset[1] = 1;
            this->getQuadCoords1(cellIndx1, offset, &nodeInds1[3], &quadCoords1[3*NDIMS_2D_PHYS]);

            for (int cellIndx2 = 0; cellIndx2 < this->numCells2; ++cellIndx2) {

                std::vector<int> nodeInds2(4);
                // build quad 2
                offset[0] = 0; offset[1] = 0;
                this->getQuadCoords2(cellIndx2, offset, &nodeInds2[0], &quadCoords2[0*NDIMS_2D_PHYS]);
                offset[0] = 1; offset[1] = 0;
                this->getQuadCoords2(cellIndx2, offset, &nodeInds2[1], &quadCoords2[1*NDIMS_2D_PHYS]);
                offset[0] = 1; offset[1] = 1;
                this->getQuadCoords2(cellIndx2, offset, &nodeInds2[2], &quadCoords2[2*NDIMS_2D_PHYS]);
                offset[0] = 0; offset[1] = 1;
                this->getQuadCoords2(cellIndx2, offset, &nodeInds2[3], &quadCoords2[3*NDIMS_2D_PHYS]);

                std::vector<double> intersectPoints;
                intersectPoints.reserve(100);; // rough guess

                //
                // find the quad 1 nodes inside quad2
                //

                //
                // find the quad2 nodes inside quad 1
                //

                //
                // find the intersection points between the quad 1 and quad 2 edges
                //

                // iterate over the quad 1 edges
                for (size_t j = 0; j < 4; ++j) {
                    size_t jA = j;
                    size_t jB = (j + 1) % 4;
                    double* coords1A = &quadCoords1[jA*NDIMS_2D_PHYS];
                    double* coords1B = &quadCoords1[jB*NDIMS_2D_PHYS];
                    // iterate over the quad 2 edges
                    for (size_t i = 0; i < 4; ++i) {
                        size_t iA = i;
                        size_t iB = (i + 1) % 4;
                        double* coords2A = &quadCoords2[iA*NDIMS_2D_PHYS];
                        double* coords2B = &quadCoords2[iB*NDIMS_2D_PHYS];

                        // find intersections
                        this->collectEdgeToEdgeIntersectionPoints(coords1A, coords1B,
                                                                  coords2A, coords2B,
                                                                  intersectionPoints);
                    }
                }

                if (intersectPoints.size() > 0) {
                    cellIndx12[0] = cellIndx1;
                    cellIndx12[1] = cellindx2;
                    std::pair< std::vector<int>, std::vector<double> > p(cellIndx12, intersectPoints);
                    // add the intersection points
                    this->cellIndx1Indx2_points.insert(p);
                }

            }
        }
    }

private:

    void collectEdgeToEdgeIntersectionPoints(const double edge1Point0[],
                                             const double edge1Point1[],
                                             const double edge2Point0[],
                                             const double edge2Point1[]
                                             std::vector<double>& intersectionPoints) {
        std::vector<double> mat(2 * 2);
        std::vector<double> rhs(2);
        // edge1Point0 + xi*(edge1Point1 - edge1Point0) = edge2Point0 + eta*(edge2Point1 - edge2Point0)
        for (size_t i = 0; i < 2; ++i) {
            rhs[i] = edge2Point0[i] - edge1Point0[i];
            mat[2*i + 0] = edge1Point1[i] - edge1Point0[i];
            mat[2*i + 1] = edge2Point0[i] - edge2Point1[i];
        }
        this->slvr->setMatrix(&mat[0]);
        this->slvr->setRightHandSide(&rhs[0]);
        int ier = this->slvr->solve();
        if (ier > 0) {
            // singular system, likely because the two edges are parallel 
            // not adding any point, even if the edges are degenerate
            return;
        }
        double* xis;
        this->slvr->getSolution(&xis);
        // make sure the parametric coordinates are within the (0+, 1-) range
        // no need to include the end points since they are already taken into 
        // account when looking for nodes inside cell
        if (xis[0] > this->tol && xis[0] <= 1 - this->tol && 
            xis[1] > this->tol && xis[1] <= 1 - this->tol) {
            // the two edges intersect
            for (size_t j = 0; j < NDIMS_2D_PHYS; ++j) {
                double p = edge1Point0[j] + xis[0]*(edge1Point1[j] - edge1Point0[j]);
                intersectionPoints.push_back(p);
            }
        }
    }

    /** 
     * Extract the coordinates from grid 1
     * @param indx cell flat index
     * @param offset displacement from the above node
     * @param nodeIndx flat node index (output)
     * @param coords array of size NDIMS_2D_PHYS to be filled in 
     */
    void getQuadCoords1(size_t indx, const int offset[], int* nodeIndx, double coords[]) const {

        // compute the index set of the cell and add the offset
        size_t cellIndsOffset[NDIMS_2D_TOPO];
        for (size_t j = 0; j < NDIMS_2D_TOPO; ++j) {
            cellIndsOffset[j] = indx / this->cellDimsProd1[j] % this->cellDims1[j];
            cellIndsOffset[j] += offset[j];
        }

        // compute the low-corner flat index of the node coorresponding to this cell
        nodeIndx = 0;
        for (size_t j = 0; j < NDIMS_2D_TOPO; ++j) {
            nodeIndx += this->nodeDimsProd1[j] * cellIndsOffset[j];
        }

        // fill in the node's coordinates
        for (size_t j = 0; j < NDIMS_2D_PHYS; ++j) {
            coords[j] = this->coords1[nodeIndx * NDIMS_2D_PHYS + j];
        }
    }

    /** 
     * Extract the coordinates from grid 2
     * @param indx cell flat index
     * @param offset displacement from the above node
     * @param nodeIndx flat node index (output)
     * @param coords array of size NDIMS_2D_PHYS to be filled in 
     */
    void getQuadCoords2(size_t indx, const int offset[], int* nodeIndx, double coords[]) const {

        // compute the index set of the cell and add the offset
        size_t cellIndsOffset[NDIMS_2D_TOPO];
        for (size_t j = 0; j < NDIMS_2D_TOPO; ++j) {
            cellIndsOffset[j] = indx / this->cellDimsProd2[j] % this->cellDims2[j];
            cellIndsOffset[j] += offset[j];
        }

        // compute the low-corner flat index of the node coorresponding to this cell
        nodeIndx = 0;
        for (size_t j = 0; j < NDIMS_2D_TOPO; ++j) {
            nodeIndx += this->nodeDimsProd2[j] * cellIndsOffset[j];
        }

        // fill in the node's coordinates
        for (size_t j = 0; j < NDIMS_2D_PHYS; ++j) {
            coords[j] = this->coords2[nodeIndx * NDIMS_2D_PHYS + j];
        }
    }

};
 
#ifdef __cplusplus
extern "C" {
#endif

int SgStructGrid2DIntersect_new(SgStructGrid2DIntersect_type** self);
                       
int SgStructGrid2DIntersect_del(SgStructGrid2DIntersect_type** self);


#ifdef __cplusplus
}
#endif

#endif // SG_STRUCT_GRID_2D_INTERSECT_H
