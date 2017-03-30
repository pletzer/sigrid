/**
 * A class that finds all intersection points between a QuadLine and a line
 */
 
#ifndef SG_QUAD_LINE_INTERSECT_H
#define SG_QUAD_LINE_INTERSECT_H
 
#include <vector>
#include <cstdio> // size_t
#include <cmath>
#include <limits>
#include <algorithm>
#include "SgLinearSolve.h"
#include "SgTriangulate.h"
#include "SgNdims.h"
 
struct SgQuadLineIntersect_type {

    // linear solver
    SgLinearSolve_type* slvr;

    // flat array (node, component) of the nodal coordinates
    double* quadCoords;
    double* lineCoords;

    // intersection points as a flat array (array of zero, one or two elements)
    std::vector<double> intersectionPoints;

    // a tolerance to determine whether a two edges intersect
    // or a point is within a triangle
    double tol;

    // for finding intersections of edges
    double mat[2 * 2];
    double rhs[2];

    // the quad corner points
    double quadMin[NDIMS_2D_PHYS];
    double quadMax[NDIMS_2D_PHYS];

    // the straight line corner points
    double lineMin[NDIMS_2D_PHYS];
    double lineMax[NDIMS_2D_PHYS];

    /**
     * Constructor
     */
     SgQuadLineIntersect_type() {

        // tolerance for floating point comparisons
        this->tol = 1.e-12;

        this->slvr = new SgLinearSolve_type(2, 2);

        this->quadCoords = NULL;
        this->lineCoords = NULL;
    }

    void reset() {
        this->intersectionPoints.resize(0);
    }

    /**
     * Check if the boxes containg the quad and line overlap
     * @return true if they overlap, false otherwise
     */
    bool checkIfOverlap() {

      // overlap is when the min of the high corner points is higher than the max
      // of the low corner points
      bool overlap = true;
      for (size_t j = 0; j < NDIMS_2D_PHYS; ++j) {
        overlap &= std::min(this->quadMax[j], this->lineMax[j]) >= std::max(this->quadMin[j], this->lineMin[j]);
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

        // x0 + xi*(x1 - x0) + eta*(x2 - x0) = pt
        for (size_t i = 0; i < 2; ++i) {
            this->rhs[i] = pt[i] - node0[i];
            this->mat[2*i + 0] = node1[i] - node0[i];
            this->mat[2*i + 1] = node2[i] - node0[i];
        }
        this->slvr->setMatrix(this->mat);
        this->slvr->setRightHandSide(this->rhs);
        int ier = this->slvr->solve();
        if (ier > 0) {
            // singular system, likely degenerate triangle
            return false;
        }
        double* xis;
        this->slvr->getSolution(&xis);
        // make sure the parametric coordinates are within triangle
        // (0 <= xi < 1 and 0 <= xi + eta < 1). Allow for a small tolerance
        // in case the point is right on the triangle's edges
        return xis[0] > -this->tol && xis[0] <= 1 + this->tol && 
                xis[1] > -this->tol && xis[0] + xis[1] <= 1 + this->tol;
    }

    /**
     * Add a point to the list of intersections
     * @param point point
     */
    void addPoint(const std::vector<double>& point) {
      for (size_t i = 0; i < NDIMS_2D_PHYS; ++i) {
        this->intersectionPoints.push_back(point[i]);
      }
    }

    /**
     * Get the intersection points gathered so far
     * @return array
     */
    const std::vector<double>& getIntersectionPoints() const {
        return this->intersectionPoints;
    }

    /**
     * Is point inside quad?
     * @param point target coordinates
     * @param quad 4 vertex points as flat array
     */
    bool isPointInQuad(const double* point, const double* quad) {
      if (this->checkIfPointIsInsideTriangle(point, &quad[0*NDIMS_2D_PHYS],
                                                    &quad[1*NDIMS_2D_PHYS],
                                                    &quad[3*NDIMS_2D_PHYS])) {
        return true;
      }
      if (this->checkIfPointIsInsideTriangle(point, &quad[2*NDIMS_2D_PHYS],
                                                    &quad[3*NDIMS_2D_PHYS],
                                                    &quad[1*NDIMS_2D_PHYS])) {
        return true;
      }
      return false;
    }

    /**
     * Collect edge to line intersection points
     * @param edgePoint1 start point of edge 1
     * @param edgePoint1 end point of edge 1
     * @param edgePoint2 start point of edge 2
     * @param edgePoint2 end point of edge 2
     * @param pt intersection point (if present)
     * @return 1 if intersection, 0 otherwise
     */
    int collectEdgeToLineIntersectionPoints(const double edgePoint0[],
                                            const double edgePoint1[],
                                            const double linePoint0[],
                                            const double linePoint1[],
                                            double pt[]) {
        // edgePoint0 + xi*(edgePoint1 - edgePoint0) = linePoint0 + eta*(linePoint1 - linePoint0)
        for (size_t i = 0; i < 2; ++i) {
            this->rhs[i] = linePoint0[i] - edgePoint0[i];
            this->mat[2*i + 0] = edgePoint1[i] - edgePoint0[i];
            this->mat[2*i + 1] = linePoint0[i] - linePoint1[i];
        }
        this->slvr->setMatrix(this->mat);
        this->slvr->setRightHandSide(this->rhs);
        int ier = this->slvr->solve();
        if (ier > 0) {
            // singular system, likely because the two edges are parallel 
            // not adding any point, even if the edges are degenerate
            return 0;
        }
        double* xis;
        this->slvr->getSolution(&xis);
        // make sure the parametric coordinates are within the (0+, 1-) range
        // no need to include the end points since they are already taken into 
        // account when looking for nodes inside cell
        if (xis[0] > this->tol && xis[0] < 1.0 - this->tol && 
            xis[1] > this->tol && xis[1] < 1.0 - this->tol) {
            // the two edges intersect
            for (size_t j = 0; j < NDIMS_2D_PHYS; ++j) {
                pt[j] = edgePoint0[j] + xis[0]*(edgePoint1[j] - edgePoint0[j]);
                this->intersectionPoints.push_back(pt[j]);
            }
            return 1;
        }
        return 0;
    }

    /**
     * Set the quad's coordinates
     * @param quadCoords quad coordinates
     */
    void setQuadPoints(const double* quadCoords) {

      this->quadCoords = (double*)quadCoords;

      // compute the corner points
      // initialize the box min/max coordinates
      for (size_t j = 0; j < NDIMS_2D_PHYS; ++j) {
        this->quadMin[j] = std::numeric_limits<double>::infinity();
        this->quadMax[j] = -std::numeric_limits<double>::infinity();
        for (size_t k = 0; k < 4; ++k) { // 4 nodes
          size_t i = NDIMS_2D_PHYS*k + j;
          this->quadMin[j] = std::min(this->quadCoords[i], this->quadMin[j]);
          this->quadMax[j] = std::max(this->quadCoords[i], this->quadMax[j]);
        }
      }
    }

    /**
     * Set the line's coordinates
     * @param lineCoords line coordinates
     */
    void setLinePoints(const double* lineCoords) {

      this->lineCoords = (double*) lineCoords;

      // compute the corner points
      // initialize the box min/max coordinates
      for (size_t j = 0; j < NDIMS_2D_PHYS; ++j) {
        this->lineMin[j] = std::numeric_limits<double>::infinity();
        this->lineMax[j] = -std::numeric_limits<double>::infinity();
        for (size_t k = 0; k < 2; ++k) { // 2 nodes
          size_t i = NDIMS_2D_PHYS*k + j;
          // perturb lineCoords to avoid hitting quad nodes
          //this->lineCoords[i] += (double(i) - 1.5)*12.8645237120137646*this->tol;
          this->lineMin[j] = std::min(this->lineCoords[i], this->lineMin[j]);
          this->lineMax[j] = std::max(this->lineCoords[i], this->lineMax[j]);
        }
      }
    }

    /**
     * Collect all the intersection points
     * @param numPoints number of points (output)
     * @param points flat array of points (output)
     */
    void collectIntersectPoints(int* numPoints, double** points) {

        *numPoints = 0;
        double pt[NDIMS_2D_PHYS];

        // add the starting line point if inside the quad
        if (this->isPointInQuad(&this->lineCoords[0*NDIMS_2D_PHYS], this->quadCoords)) {
          std::vector<double> pt(&this->lineCoords[0*NDIMS_2D_PHYS], 
                                 &this->lineCoords[0*NDIMS_2D_PHYS] + 2);
          this->addPoint(pt);
        }

        // iterate over edges
        for (size_t i = 0; i < 4; ++i) {
            // the edge of one of the first QuadLine
            size_t quadIndxA = i;
            size_t quadIndxB = (i + 1) % 4;
            double* quadCoordA = &this->quadCoords[quadIndxA*NDIMS_2D_PHYS];
            double* quadCoordB = &this->quadCoords[quadIndxB*NDIMS_2D_PHYS];
            double* lineCoordA = &this->lineCoords[0*NDIMS_2D_PHYS];
            double* lineCoordB = &this->lineCoords[1*NDIMS_2D_PHYS];
            this->collectEdgeToLineIntersectionPoints(quadCoordA, quadCoordB,
                                                      lineCoordA, lineCoordB, pt);
        }

        // add the ending line point if inside the quad
        size_t npts = this->intersectionPoints.size() / NDIMS_2D_PHYS;
        if (npts < 2 && 
            this->isPointInQuad(&this->lineCoords[1*NDIMS_2D_PHYS], this->quadCoords)) {
          std::vector<double> pt(&this->lineCoords[1*NDIMS_2D_PHYS], &this->lineCoords[1*NDIMS_2D_PHYS] + 2);
          this->addPoint(pt);
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

int SgQuadLineIntersect_new(SgQuadLineIntersect_type** self);
                       
int SgQuadLineIntersect_del(SgQuadLineIntersect_type** self);

int SgQuadLineIntersect_setQuadPoints(SgQuadLineIntersect_type** self, 
                                      const double* quadPoints);

int SgQuadLineIntersect_setLinePoints(SgQuadLineIntersect_type** self, 
                                      const double* linePoints);

int SgQuadLineIntersect_collectIntersectPoints(SgQuadLineIntersect_type** self, 
                                               int *numPoints, double** points);

#ifdef __cplusplus
}
#endif

#endif // SG_QUAD_LINE_INTERSECT_H
