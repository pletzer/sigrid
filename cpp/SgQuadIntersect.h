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
 
struct SgQuadIntersect_type {

 	// linear solver
 	SgLinearSolve_type* slvr;

 	// the nodes of the two quads
 	double** quad1Coords;
 	double** quad2Coords;

 	// the collection of intersection points as a flat array
 	// 6 is the max number of intersection points
 	double intersectionPoints[6*2];
 	int numIntersectionPoints;

 	// a tolerance to determine whether a two edges intersect
 	// or a point is within a triangle
 	double tol;

 	bool checkIfBoxesOverlap() {

 		// initialize the box min/max coordinates
 		std::vector<double> quad1Min(2, std::numeric_limits<double>::infinity());
 		std::vector<double> quad1Max(2, -std::numeric_limits<double>::infinity());
 		std::vector<double> quad2Min(2, std::numeric_limits<double>::infinity());
 		std::vector<double> quad2Max(2, -std::numeric_limits<double>::infinity());

 		// find each quad's box corner points
 		for (size_t j = 0; j < 4; ++j) { // 4 nodes
 			for (size_t i = 0; i < 2; ++i) {  // 2 dims

 				// quad 1
 				quad1Min[i] = std::min(this->quad1Coords[j][i], quad1Min[i]);
 				quad1Max[i] = std::max(this->quad1Coords[j][i], quad1Max[i]);

 				// quad 2
 				quad2Min[i] = std::min(this->quad2Coords[j][i], quad2Min[i]);
 				quad2Max[i] = std::max(this->quad2Coords[j][i], quad2Max[i]);
 			}
		}

		// overlap is when the min of the high corner points is higher than the max
		// of the low corner points
 		bool overlap = true;
 		for (size_t i = 0; i < 2; ++i) {
 			overlap &= std::min(quad1Max[i], quad2Max[i]) >= std::max(quad1Min[i], quad2Min[i]);
 		}
 		return overlap;
 	}

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

    void collectNodesInsideQuad(const double** nodes, const double** quad) {

    	for (size_t i = 0; i < 4; ++i) {
    		if (checkIfPointIsInsideTriangle(nodes[i], quad[0], quad[1], quad[3])) {
    			// pt is in inside triangle 0, 1, 3 of quad
    			this->intersectionPoints[this->numIntersectionPoints*2 + 0] = nodes[i][0];
    			this->intersectionPoints[this->numIntersectionPoints*2 + 1] = nodes[i][1];
    			this->numIntersectionPoints++;
    		}
    		else if (checkIfPointIsInsideTriangle(nodes[i], quad[2], quad[3], quad[1])) {
    			// pt is inside triangle 2, 3, 1 of quad
    			this->intersectionPoints[this->numIntersectionPoints*2 + 0] = nodes[i][0];
    			this->intersectionPoints[this->numIntersectionPoints*2 + 1] = nodes[i][1];
    			this->numIntersectionPoints++;
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
 			double intersectPt[] = {edge1Point0[0] + xis[0]*(edge1Point1[0] - edge1Point0[0]),
 			                        edge1Point0[1] + xis[0]*(edge1Point1[1] - edge1Point0[1])};
 			this->intersectionPoints[this->numIntersectionPoints*2 + 0] = intersectPt[0];
 			this->intersectionPoints[this->numIntersectionPoints*2 + 1] = intersectPt[1];
 			this->numIntersectionPoints++;
 		}
	
    }

    void setQuadPoints(const double** quad1Coords, const double** quad2Coords) {
        this->quad1Coords = (double**) quad1Coords;
        this->quad2Coords = (double**) quad2Coords;
    }


    void collectIntersectPoints(int* numPoints, double** points) {

        // quickly check if there is any chance of overlap
        if (!this->checkIfBoxesOverlap()) {
            // no chance to have an overlap
            return;
        }

        // seems like the quads are at least partially overlapping 
        this->collectNodesInsideQuad((const double**)this->quad1Coords, 
                                     (const double**)this->quad2Coords);
        this->collectNodesInsideQuad((const double**)this->quad2Coords,
                                     (const double**)this->quad1Coords);
        // iterate over edges
        for (size_t i = 0; i < 4; ++i) {
            // the edge of one of the first quad
            size_t quad1Indx0 = i;
            size_t quad1Indx1 = (i + 1) % 4;
            const double* quad1Coord0 = this->quad1Coords[quad1Indx0];
            const double* quad1Coord1 = this->quad1Coords[quad1Indx1];
            // iterate over the other quad's edges
            for (size_t j = 0; j < 4; ++j) {
                // the edges of the second quad
                size_t quad2Indx0 = j;
                size_t quad2Indx1 = (j + 1) % 4;
                const double* quad2Coord0 = this->quad2Coords[quad2Indx0];
                const double* quad2Coord1 = this->quad2Coords[quad2Indx1];
                this->collectEdgeToEdgeIntersectionPoints(quad1Coord0, quad1Coord1,
                                                          quad2Coord0, quad2Coord1);
            }
        }

        *numPoints = this->numIntersectionPoints;
        *points = this->intersectionPoints;
    }

};
 
#ifdef __cplusplus
extern "C" {
#endif

int SgQuadIntersect_new(SgQuadIntersect_type** self);
                       
int SgQuadIntersect_del(SgQuadIntersect_type** self);

int SgQuadIntersect_setQuadPoints(SgQuadIntersect_type** self, 
                                  const double** quad1Points, const double** quad2Points);

int SgQuadIntersect_getIntersectPoints(SgQuadIntersect_type** self, 
 	                                    int *numPoints, double** points);

#ifdef __cplusplus
}
#endif

#endif // SG_QUAD_INTERSECT_H