/**
 * A class that finds all intersection points between two quads
 */
 
#ifndef SG_QUAD_INTERSECT_H
#define SG_QUAD_INTERSECT_H
 
#include <vector>
#include <cstdio> // size_t
#include <cmath>
#include <limits>
#include <algorithm>
#include "SgLinearSolve.h"
#include "SgNdims.h"
 
struct SgQuadIntersect_type {

    // linear solver
    SgLinearSolve_type* slvr;

    // flat array (node, component) of the nodal coordinates
    double* quad1Coords;
    double* quad2Coords;

    // the collection of intersection points as a flat array
    std::vector<double> intersectionPoints;

    // a tolerance to determine whether a two edges intersect
    // or a point is within a triangle
    double tol;

    /**
     * Check if the boxes containg the two quads overlap
     * @return true if they overlap, false otherwise
     */
    bool checkIfBoxesOverlap() {

        // initialize the box min/max coordinates
        std::vector<double> quad1Min(2, std::numeric_limits<double>::infinity());
        std::vector<double> quad1Max(2, -std::numeric_limits<double>::infinity());
        std::vector<double> quad2Min(2, std::numeric_limits<double>::infinity());
        std::vector<double> quad2Max(2, -std::numeric_limits<double>::infinity());

         // find each quad's box corner points
        for (size_t i = 0; i < 4; ++i) { // 4 nodes
            for (size_t j = 0; j < NDIMS_2D_PHYS; ++j) {  // 2 dims

                size_t k = i*NDIMS_2D_PHYS + j;

                // quad 1
                quad1Min[j] = std::min(this->quad1Coords[k], quad1Min[j]);
                quad1Max[j] = std::max(this->quad1Coords[k], quad1Max[j]);

                // quad 2
                quad2Min[j] = std::min(this->quad2Coords[k], quad2Min[j]);
                quad2Max[j] = std::max(this->quad2Coords[k], quad2Max[j]);
            }
        }

        // overlap is when the min of the high corner points is higher than the max
        // of the low corner points
        bool overlap = true;
        for (size_t j = 0; j < NDIMS_2D_PHYS; ++j) {
            overlap &= std::min(quad1Max[j], quad2Max[j]) >= std::max(quad1Min[j], quad2Min[j]);
        }
        return overlap;
    }

    /**
     * Check if a point is inside a triangle
     * @param pt the target point
     * @param node0 1st point of the triangle
     * @param node1 2nd point of the triangle
     * @param node2 3rd point of the triangle
     * @return true if the point in the triangle or on its boundary
     */
    bool checkIfPointIsInsideTriangle(const double* pt, 
                                       const double* node0,
                                       const double* node1,
                                       const double* node2) {

        std::vector<double> mat(2 * 2);
        std::vector<double> rhs(2);
        // x0 + xi*(x1 - x0) + eta*(x2 - x0) = pt
        for (size_t i = 0; i < 2; ++i) {
            rhs[i] = pt[i] - node0[i];
            mat[2*i + 0] = node1[i] - node0[i];
            mat[2*i + 1] = node2[i] - node0[i];
        }
        SgLinearSolve_setMatrix(&this->slvr, &mat[0]);
        SgLinearSolve_setRightHandSide(&this->slvr, &rhs[0]);
        int ier = SgLinearSolve_solve(&this->slvr);
        if (ier > 0) {
            // singular system, likely degenerate triangle
            return false;
        }
        double* xis;
        SgLinearSolve_getSolution(&this->slvr, &xis);
        // make sure the parametric coordinates are within triangle
        // (0 <= xi < 1 and 0 <= xi + eta < 1). Allow for a small tolerance
        // in case the point is right on the triangle's edges
        return xis[0] > -this->tol && xis[0] <= 1 + this->tol && 
                xis[1] > -this->tol && xis[0] + xis[1] <= 1 + this->tol;

    }

    /**
     * Collect all the quad's nodes that are inside the other quad
     * @param nodes flat array of the 1st quad's coordinates
     * @param quad flat array of the 2nd quad's cooridnates
     */
    void collectNodesInsideQuad(const double* nodes, const double* quad) {

        for (size_t i = 0; i < 4; ++i) {

            if (this->checkIfPointIsInsideTriangle(&nodes[i*NDIMS_2D_PHYS], 
                                                   &quad[0*NDIMS_2D_PHYS],
                                                   &quad[1*NDIMS_2D_PHYS],
                                                   &quad[3*NDIMS_2D_PHYS])) {
                // pt is in inside triangle 0, 1, 3 of quad
                for (size_t j = 0; j < NDIMS_2D_PHYS; ++j) {
                    this->intersectionPoints.push_back(nodes[i*NDIMS_2D_PHYS + j]);
                }
            }
            else if (this->checkIfPointIsInsideTriangle(&nodes[i*NDIMS_2D_PHYS], 
                                                        &quad[2*NDIMS_2D_PHYS],
                                                        &quad[3*NDIMS_2D_PHYS],
                                                        &quad[1*NDIMS_2D_PHYS])) {
                // pt is inside triangle 2, 3, 1 of quad
                for (size_t j = 0; j < NDIMS_2D_PHYS; ++j) {
                    this->intersectionPoints.push_back(nodes[i*NDIMS_2D_PHYS + j]);
                }
            }
        }
    }

    void collectEdgeToEdgeIntersectionPoints(const double edge1Point0[],
                                             const double edge1Point1[],
                                             const double edge2Point0[],
                                             const double edge2Point1[]) {
        std::vector<double> mat(2 * 2);
        std::vector<double> rhs(2);
        // edge1Point0 + xi*(edge1Point1 - edge1Point0) = edge2Point0 + eta*(edge2Point1 - edge2Point0)
        for (size_t i = 0; i < 2; ++i) {
            rhs[i] = edge2Point0[i] - edge1Point0[i];
            mat[2*i + 0] = edge1Point1[i] - edge1Point0[i];
            mat[2*i + 1] = edge2Point0[i] - edge2Point1[i];
        }
        SgLinearSolve_setMatrix(&this->slvr, &mat[0]);
        SgLinearSolve_setRightHandSide(&this->slvr, &rhs[0]);
        int ier = SgLinearSolve_solve(&this->slvr);
        if (ier > 0) {
            // singular system, likely because the two edges are parallel 
            // not adding any point, even if the edges are degenerate
            return;
        }
        double* xis;
        SgLinearSolve_getSolution(&this->slvr, &xis);
        // make sure the parametric coordinates are within the (0, 1) range
        if (xis[0] > -this->tol && xis[0] <= 1 + this->tol && 
            xis[1] > -this->tol && xis[1] <= 1 + this->tol) {
            // the two edges intersect
            for (size_t j = 0; j < NDIMS_2D_PHYS; ++j) {
                double p = edge1Point0[j] + xis[0]*(edge1Point1[j] - edge1Point0[j]);
                this->intersectionPoints.push_back(p);
            }
        }
    }

    void setQuadPoints(const double* quad1Coords, const double* quad2Coords) {
        this->quad1Coords = (double*) quad1Coords;
        this->quad2Coords = (double*) quad2Coords;
    }


    void collectIntersectPoints(int* numPoints, double** points) {

        *numPoints = 0;

        // quickly check if there is any chance of overlap
        if (!this->checkIfBoxesOverlap()) {
            // no chance to have an overlap
            return;
        }

        // seems like the quads are at least partially overlapping 
        this->collectNodesInsideQuad((const double*)this->quad1Coords, 
                                     (const double*)this->quad2Coords);
        this->collectNodesInsideQuad((const double*)this->quad2Coords,
                                     (const double*)this->quad1Coords);
        // iterate over edges
        for (size_t i = 0; i < 4; ++i) {
            // the edge of one of the first quad
            size_t quad1IndxA = i;
            size_t quad1IndxB = (i + 1) % 4;
            double* quad1CoordA = &this->quad1Coords[quad1IndxA*NDIMS_2D_PHYS];
            double* quad1CoordB = &this->quad1Coords[quad1IndxB*NDIMS_2D_PHYS];
            // iterate over the other quad's edges
            for (size_t j = 0; j < 4; ++j) {
                // the edges of the second quad
                size_t quad2IndxA = j;
                size_t quad2IndxB = (j + 1) % 4;
                double* quad2CoordA = &this->quad2Coords[quad2IndxA*NDIMS_2D_PHYS];
                double* quad2CoordB = &this->quad2Coords[quad2IndxB*NDIMS_2D_PHYS];
                this->collectEdgeToEdgeIntersectionPoints(quad1CoordA, quad1CoordB,
                                                          quad2CoordA, quad2CoordB);
            }
        }

        // set the return values
        *numPoints = this->intersectionPoints.size() / NDIMS_2D_PHYS;
        *points = NULL;
        if (*numPoints > 0) {
            *points = &this->intersectionPoints[0];
        }
    }

};
 
#ifdef __cplusplus
extern "C" {
#endif

int SgQuadIntersect_new(SgQuadIntersect_type** self);
                       
int SgQuadIntersect_del(SgQuadIntersect_type** self);

int SgQuadIntersect_setQuadPoints(SgQuadIntersect_type** self, 
                                  const double* quad1Points, const double* quad2Points);

int SgQuadIntersect_getIntersectPoints(SgQuadIntersect_type** self, 
                                       int *numPoints, double** points);

#ifdef __cplusplus
}
#endif

#endif // SG_QUAD_INTERSECT_H
