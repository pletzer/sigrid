/**
 * A class that computes the projection of flux basis functions on a line in 2D
 */
 
#ifndef SG_QUAD_LINE_FLOWS_H
#define SG_QUAD_LINE_FLOWS_H
 
#include <vector>
#include <cstdio> // size_t
#include <cmath>
#include <limits>
#include <algorithm>
#include "SgNdims.h"

const int EDGE_LO_X = 0;
const int EDGE_HI_X = 1;
const int EDGE_LO_Y = 2;
const int EDGE_HI_Y = 3;

struct SgQuadLineFlows_type {

    // flat array (node, component) of the nodal coordinates
    double* quadCoords;
    double* lineCoords;

    double projections[4]; // 4 edges/basis functions

    /**
     * Constructor
     */
     SgQuadLineFlows_type() {

        this->quadCoords = NULL;
        this->lineCoords = NULL;
    }

    /**
     * Set the quad's coordinates
     * @param quadCoords quad coordinates
     */
    void setQuadPoints(const double* quadCoords) {

      this->quadCoords = (double*)quadCoords;
    }

    /**
     * Set the line's coordinates
     * @param lineCoords line coordinates
     */
    void setLinePoints(const double* lineCoords) {

      this->lineCoords = (double*) lineCoords;
    }

    /**
     * Get the flux projection onto an edge
     * @param edgeIndex (0 = low x, 1 = high x, 2 = low y, 3 = high y)
     * @return flux 
     */
    double getProjection(int edgeIndex) {
      return this->projections[edgeIndex];
    }

    /**
     * Compute projections
     */
    void computeProjections() {

      double xa = this->lineCoords[0];
      double ya = this->lineCoords[1];
      double xb = this->lineCoords[2];
      double yb = this->lineCoords[3];
      double p0x = this->quadCoords[0];
      double p0y = this->quadCoords[1];
      double p1x = this->quadCoords[2];
      double p1y = this->quadCoords[3];
      double p2x = this->quadCoords[4];
      double p2y = this->quadCoords[5];
      double p3x = this->quadCoords[6];
      double p3y = this->quadCoords[7];

      double integrals[4];

      integrals[0] = (-((xa - xb)*(ya - yb)*pow(-(p0x*p2y*xa) + p1x*p2y*xa + p0x*p3y*xa - p1x*p3y*xa + 
          p0x*p2y*xb - p1x*p2y*xb - p0x*p3y*xb + p1x*p3y*xb + p0x*p2y*ya - 
          p2y*p3x*ya - p1x*p3y*ya + p2x*p3y*ya - p0x*p2y*yb + p2y*p3x*yb + 
          p1x*p3y*yb - p2x*p3y*yb + 
          p1y*(p2x*(-xa + xb) + p3x*(xa - xb + ya - yb) + p0x*(-ya + yb)) + 
          p0y*(p3x*(-xa + xb) + p1x*(ya - yb) + p2x*(xa - xb - ya + yb)),2)) + 
     2*(-(p0x*p2y*xa) + p1x*p2y*xa + p0x*p3y*xa - p1x*p3y*xa + p0x*p2y*xb - 
        p1x*p2y*xb - p0x*p3y*xb + p1x*p3y*xb + p0x*p2y*ya - p2y*p3x*ya - 
        p1x*p3y*ya + p2x*p3y*ya - p0x*p2y*yb + p2y*p3x*yb + p1x*p3y*yb - 
        p2x*p3y*yb + p1y*(p2x*(-xa + xb) + p3x*(xa - xb + ya - yb) + 
           p0x*(-ya + yb)) + p0y*(p3x*(-xa + xb) + p1x*(ya - yb) + 
           p2x*(xa - xb - ya + yb)))*
      (p0x*p2y*pow(xa,2) - p1x*p2y*pow(xa,2) - p0x*p3y*pow(xa,2) + 
        p1x*p3y*pow(xa,2) - 2*p0x*p2y*xa*xb + 2*p1x*p2y*xa*xb + 2*p0x*p3y*xa*xb - 
        2*p1x*p3y*xa*xb + p0x*p2y*pow(xb,2) - p1x*p2y*pow(xb,2) - 
        p0x*p3y*pow(xb,2) + p1x*p3y*pow(xb,2) - p1x*p2y*xa*ya + p2y*p3x*xa*ya + 
        p1x*p3y*xa*ya - p2x*p3y*xa*ya - p0x*p2y*pow(xa,2)*ya + 
        p1x*p2y*pow(xa,2)*ya + p0x*p3y*pow(xa,2)*ya - p1x*p3y*pow(xa,2)*ya + 
        p1x*p2y*xb*ya - p2y*p3x*xb*ya - p1x*p3y*xb*ya + p2x*p3y*xb*ya + 
        2*p0x*p2y*xa*xb*ya - 2*p1x*p2y*xa*xb*ya - 2*p0x*p3y*xa*xb*ya + 
        2*p1x*p3y*xa*xb*ya - p0x*p2y*pow(xb,2)*ya + p1x*p2y*pow(xb,2)*ya + 
        p0x*p3y*pow(xb,2)*ya - p1x*p3y*pow(xb,2)*ya - p0x*p2y*pow(ya,2) + 
        p2y*p3x*pow(ya,2) + p1x*p3y*pow(ya,2) - p2x*p3y*pow(ya,2) + 
        p0x*p2y*xa*pow(ya,2) - p2y*p3x*xa*pow(ya,2) - p1x*p3y*xa*pow(ya,2) + 
        p2x*p3y*xa*pow(ya,2) + p0y*
         (-(p3x*pow(xa - xb,2)*(-1 + ya)) + 
           p2x*(pow(xa,2)*(-1 + ya) + pow(xb,2)*(-1 + ya) - 
              xa*(2*xb*(-1 + ya) + pow(ya - yb,2)) + pow(ya - yb,2)) + 
           p1x*(-1 + xa)*pow(ya - yb,2)) + p1x*p2y*xa*yb - p2y*p3x*xa*yb - 
        p1x*p3y*xa*yb + p2x*p3y*xa*yb - p1x*p2y*xb*yb + p2y*p3x*xb*yb + 
        p1x*p3y*xb*yb - p2x*p3y*xb*yb + 2*p0x*p2y*ya*yb - 2*p2y*p3x*ya*yb - 
        2*p1x*p3y*ya*yb + 2*p2x*p3y*ya*yb - 2*p0x*p2y*xa*ya*yb + 
        2*p2y*p3x*xa*ya*yb + 2*p1x*p3y*xa*ya*yb - 2*p2x*p3y*xa*ya*yb - 
        p0x*p2y*pow(yb,2) + p2y*p3x*pow(yb,2) + p1x*p3y*pow(yb,2) - 
        p2x*p3y*pow(yb,2) + p0x*p2y*xa*pow(yb,2) - p2y*p3x*xa*pow(yb,2) - 
        p1x*p3y*xa*pow(yb,2) + p2x*p3y*xa*pow(yb,2) + 
        p1y*(-(p0x*(-1 + xa)*pow(ya - yb,2)) - 
           p2x*(xa - xb)*(xb + xa*(-1 + ya) - ya - xb*ya + yb) + 
           p3x*(pow(xa,2)*(-1 + ya) + pow(xb,2)*(-1 + ya) + xb*(ya - yb) - 
              pow(ya - yb,2) + xa*
               (-2*xb*(-1 + ya) - ya + pow(ya,2) + yb - 2*ya*yb + pow(yb,2)))))\
      - 2*(-(pow(p1y,2)*pow(p2x,2)*pow(xa,2)) + 
        2*p1x*p1y*p2x*p2y*pow(xa,2) - pow(p1x,2)*pow(p2y,2)*pow(xa,2) + 
        pow(p1y,2)*p2x*p3x*pow(xa,2) - p1x*p1y*p2y*p3x*pow(xa,2) - 
        p1x*p1y*p2x*p3y*pow(xa,2) + pow(p1x,2)*p2y*p3y*pow(xa,2) + 
        2*pow(p1y,2)*pow(p2x,2)*xa*xb - 4*p1x*p1y*p2x*p2y*xa*xb + 
        2*pow(p1x,2)*pow(p2y,2)*xa*xb - 2*pow(p1y,2)*p2x*p3x*xa*xb + 
        2*p1x*p1y*p2y*p3x*xa*xb + 2*p1x*p1y*p2x*p3y*xa*xb - 
        2*pow(p1x,2)*p2y*p3y*xa*xb - pow(p1y,2)*pow(p2x,2)*pow(xb,2) + 
        2*p1x*p1y*p2x*p2y*pow(xb,2) - pow(p1x,2)*pow(p2y,2)*pow(xb,2) + 
        pow(p1y,2)*p2x*p3x*pow(xb,2) - p1x*p1y*p2y*p3x*pow(xb,2) - 
        p1x*p1y*p2x*p3y*pow(xb,2) + pow(p1x,2)*p2y*p3y*pow(xb,2) + 
        pow(p1y,2)*p2x*p3x*xa*ya - p1x*p1y*p2y*p3x*xa*ya - 
        2*p1y*p2x*p2y*p3x*xa*ya + 2*p1x*pow(p2y,2)*p3x*xa*ya - 
        pow(p1y,2)*pow(p3x,2)*xa*ya + p1y*p2y*pow(p3x,2)*xa*ya - 
        p1x*p1y*p2x*p3y*xa*ya + 2*p1y*pow(p2x,2)*p3y*xa*ya + 
        pow(p1x,2)*p2y*p3y*xa*ya - 2*p1x*p2x*p2y*p3y*xa*ya + 
        2*p1x*p1y*p3x*p3y*xa*ya - p1y*p2x*p3x*p3y*xa*ya - p1x*p2y*p3x*p3y*xa*ya - 
        pow(p1x,2)*pow(p3y,2)*xa*ya + p1x*p2x*pow(p3y,2)*xa*ya - 
        pow(p1y,2)*p2x*p3x*xb*ya + p1x*p1y*p2y*p3x*xb*ya + 
        2*p1y*p2x*p2y*p3x*xb*ya - 2*p1x*pow(p2y,2)*p3x*xb*ya + 
        pow(p1y,2)*pow(p3x,2)*xb*ya - p1y*p2y*pow(p3x,2)*xb*ya + 
        p1x*p1y*p2x*p3y*xb*ya - 2*p1y*pow(p2x,2)*p3y*xb*ya - 
        pow(p1x,2)*p2y*p3y*xb*ya + 2*p1x*p2x*p2y*p3y*xb*ya - 
        2*p1x*p1y*p3x*p3y*xb*ya + p1y*p2x*p3x*p3y*xb*ya + p1x*p2y*p3x*p3y*xb*ya + 
        pow(p1x,2)*pow(p3y,2)*xb*ya - p1x*p2x*pow(p3y,2)*xb*ya - 
        pow(p1y,2)*pow(p2x,2)*xa*xb*ya + 2*p1x*p1y*p2x*p2y*xa*xb*ya - 
        pow(p1x,2)*pow(p2y,2)*xa*xb*ya + p1y*p2x*p2y*p3x*xa*xb*ya - 
        p1x*pow(p2y,2)*p3x*xa*xb*ya + pow(p1y,2)*pow(p3x,2)*xa*xb*ya - 
        p1y*p2y*pow(p3x,2)*xa*xb*ya - p1y*pow(p2x,2)*p3y*xa*xb*ya + 
        p1x*p2x*p2y*p3y*xa*xb*ya - 2*p1x*p1y*p3x*p3y*xa*xb*ya + 
        p1y*p2x*p3x*p3y*xa*xb*ya + p1x*p2y*p3x*p3y*xa*xb*ya + 
        pow(p1x,2)*pow(p3y,2)*xa*xb*ya - p1x*p2x*pow(p3y,2)*xa*xb*ya + 
        pow(p1y,2)*pow(p2x,2)*pow(xb,2)*ya - 
        2*p1x*p1y*p2x*p2y*pow(xb,2)*ya + 
        pow(p1x,2)*pow(p2y,2)*pow(xb,2)*ya - p1y*p2x*p2y*p3x*pow(xb,2)*ya + 
        p1x*pow(p2y,2)*p3x*pow(xb,2)*ya - 
        pow(p1y,2)*pow(p3x,2)*pow(xb,2)*ya + 
        p1y*p2y*pow(p3x,2)*pow(xb,2)*ya + p1y*pow(p2x,2)*p3y*pow(xb,2)*ya - 
        p1x*p2x*p2y*p3y*pow(xb,2)*ya + 2*p1x*p1y*p3x*p3y*pow(xb,2)*ya - 
        p1y*p2x*p3x*p3y*pow(xb,2)*ya - p1x*p2y*p3x*p3y*pow(xb,2)*ya - 
        pow(p1x,2)*pow(p3y,2)*pow(xb,2)*ya + 
        p1x*p2x*pow(p3y,2)*pow(xb,2)*ya + p1y*p2y*pow(p3x,2)*pow(ya,2) - 
        pow(p2y,2)*pow(p3x,2)*pow(ya,2) - p1y*p2x*p3x*p3y*pow(ya,2) - 
        p1x*p2y*p3x*p3y*pow(ya,2) + 2*p2x*p2y*p3x*p3y*pow(ya,2) + 
        p1x*p2x*pow(p3y,2)*pow(ya,2) - pow(p2x,2)*pow(p3y,2)*pow(ya,2) + 
        pow(p1y,2)*p2x*p3x*xb*pow(ya,2) - p1x*p1y*p2y*p3x*xb*pow(ya,2) - 
        p1y*p2x*p2y*p3x*xb*pow(ya,2) + p1x*pow(p2y,2)*p3x*xb*pow(ya,2) - 
        pow(p1y,2)*pow(p3x,2)*xb*pow(ya,2) + 
        pow(p2y,2)*pow(p3x,2)*xb*pow(ya,2) - p1x*p1y*p2x*p3y*xb*pow(ya,2) + 
        p1y*pow(p2x,2)*p3y*xb*pow(ya,2) + pow(p1x,2)*p2y*p3y*xb*pow(ya,2) - 
        p1x*p2x*p2y*p3y*xb*pow(ya,2) + 2*p1x*p1y*p3x*p3y*xb*pow(ya,2) - 
        2*p2x*p2y*p3x*p3y*xb*pow(ya,2) - 
        pow(p1x,2)*pow(p3y,2)*xb*pow(ya,2) + 
        pow(p2x,2)*pow(p3y,2)*xb*pow(ya,2) - 
        pow(p1y,2)*p2x*p3x*pow(xb,2)*pow(ya,2) + 
        p1x*p1y*p2y*p3x*pow(xb,2)*pow(ya,2) + 
        p1y*p2x*p2y*p3x*pow(xb,2)*pow(ya,2) - 
        p1x*pow(p2y,2)*p3x*pow(xb,2)*pow(ya,2) + 
        pow(p1y,2)*pow(p3x,2)*pow(xb,2)*pow(ya,2) - 
        p1y*p2y*pow(p3x,2)*pow(xb,2)*pow(ya,2) + 
        p1x*p1y*p2x*p3y*pow(xb,2)*pow(ya,2) - 
        p1y*pow(p2x,2)*p3y*pow(xb,2)*pow(ya,2) - 
        pow(p1x,2)*p2y*p3y*pow(xb,2)*pow(ya,2) + 
        p1x*p2x*p2y*p3y*pow(xb,2)*pow(ya,2) - 
        2*p1x*p1y*p3x*p3y*pow(xb,2)*pow(ya,2) + 
        p1y*p2x*p3x*p3y*pow(xb,2)*pow(ya,2) + 
        p1x*p2y*p3x*p3y*pow(xb,2)*pow(ya,2) + 
        pow(p1x,2)*pow(p3y,2)*pow(xb,2)*pow(ya,2) - 
        p1x*p2x*pow(p3y,2)*pow(xb,2)*pow(ya,2) + 
        pow(p0y,2)*(p1x - p2x)*(p2x - p3x)*
         pow(xb + ya - xb*ya + xa*(-1 + yb) - yb,2) + 
        pow(p0x,2)*(p1y - p2y)*(p2y - p3y)*
         pow(xb + ya - xb*ya + xa*(-1 + yb) - yb,2) - pow(p1y,2)*p2x*p3x*xa*yb + 
        p1x*p1y*p2y*p3x*xa*yb + 2*p1y*p2x*p2y*p3x*xa*yb - 
        2*p1x*pow(p2y,2)*p3x*xa*yb + pow(p1y,2)*pow(p3x,2)*xa*yb - 
        p1y*p2y*pow(p3x,2)*xa*yb + p1x*p1y*p2x*p3y*xa*yb - 
        2*p1y*pow(p2x,2)*p3y*xa*yb - pow(p1x,2)*p2y*p3y*xa*yb + 
        2*p1x*p2x*p2y*p3y*xa*yb - 2*p1x*p1y*p3x*p3y*xa*yb + p1y*p2x*p3x*p3y*xa*yb + 
        p1x*p2y*p3x*p3y*xa*yb + pow(p1x,2)*pow(p3y,2)*xa*yb - 
        p1x*p2x*pow(p3y,2)*xa*yb + pow(p1y,2)*pow(p2x,2)*pow(xa,2)*yb - 
        2*p1x*p1y*p2x*p2y*pow(xa,2)*yb + 
        pow(p1x,2)*pow(p2y,2)*pow(xa,2)*yb - p1y*p2x*p2y*p3x*pow(xa,2)*yb + 
        p1x*pow(p2y,2)*p3x*pow(xa,2)*yb - 
        pow(p1y,2)*pow(p3x,2)*pow(xa,2)*yb + 
        p1y*p2y*pow(p3x,2)*pow(xa,2)*yb + p1y*pow(p2x,2)*p3y*pow(xa,2)*yb - 
        p1x*p2x*p2y*p3y*pow(xa,2)*yb + 2*p1x*p1y*p3x*p3y*pow(xa,2)*yb - 
        p1y*p2x*p3x*p3y*pow(xa,2)*yb - p1x*p2y*p3x*p3y*pow(xa,2)*yb - 
        pow(p1x,2)*pow(p3y,2)*pow(xa,2)*yb + 
        p1x*p2x*pow(p3y,2)*pow(xa,2)*yb + pow(p1y,2)*p2x*p3x*xb*yb - 
        p1x*p1y*p2y*p3x*xb*yb - 2*p1y*p2x*p2y*p3x*xb*yb + 
        2*p1x*pow(p2y,2)*p3x*xb*yb - pow(p1y,2)*pow(p3x,2)*xb*yb + 
        p1y*p2y*pow(p3x,2)*xb*yb - p1x*p1y*p2x*p3y*xb*yb + 
        2*p1y*pow(p2x,2)*p3y*xb*yb + pow(p1x,2)*p2y*p3y*xb*yb - 
        2*p1x*p2x*p2y*p3y*xb*yb + 2*p1x*p1y*p3x*p3y*xb*yb - p1y*p2x*p3x*p3y*xb*yb - 
        p1x*p2y*p3x*p3y*xb*yb - pow(p1x,2)*pow(p3y,2)*xb*yb + 
        p1x*p2x*pow(p3y,2)*xb*yb - pow(p1y,2)*pow(p2x,2)*xa*xb*yb + 
        2*p1x*p1y*p2x*p2y*xa*xb*yb - pow(p1x,2)*pow(p2y,2)*xa*xb*yb + 
        p1y*p2x*p2y*p3x*xa*xb*yb - p1x*pow(p2y,2)*p3x*xa*xb*yb + 
        pow(p1y,2)*pow(p3x,2)*xa*xb*yb - p1y*p2y*pow(p3x,2)*xa*xb*yb - 
        p1y*pow(p2x,2)*p3y*xa*xb*yb + p1x*p2x*p2y*p3y*xa*xb*yb - 
        2*p1x*p1y*p3x*p3y*xa*xb*yb + p1y*p2x*p3x*p3y*xa*xb*yb + 
        p1x*p2y*p3x*p3y*xa*xb*yb + pow(p1x,2)*pow(p3y,2)*xa*xb*yb - 
        p1x*p2x*pow(p3y,2)*xa*xb*yb - 2*p1y*p2y*pow(p3x,2)*ya*yb + 
        2*pow(p2y,2)*pow(p3x,2)*ya*yb + 2*p1y*p2x*p3x*p3y*ya*yb + 
        2*p1x*p2y*p3x*p3y*ya*yb - 4*p2x*p2y*p3x*p3y*ya*yb - 
        2*p1x*p2x*pow(p3y,2)*ya*yb + 2*pow(p2x,2)*pow(p3y,2)*ya*yb - 
        pow(p1y,2)*p2x*p3x*xa*ya*yb + p1x*p1y*p2y*p3x*xa*ya*yb + 
        p1y*p2x*p2y*p3x*xa*ya*yb - p1x*pow(p2y,2)*p3x*xa*ya*yb + 
        pow(p1y,2)*pow(p3x,2)*xa*ya*yb - pow(p2y,2)*pow(p3x,2)*xa*ya*yb + 
        p1x*p1y*p2x*p3y*xa*ya*yb - p1y*pow(p2x,2)*p3y*xa*ya*yb - 
        pow(p1x,2)*p2y*p3y*xa*ya*yb + p1x*p2x*p2y*p3y*xa*ya*yb - 
        2*p1x*p1y*p3x*p3y*xa*ya*yb + 2*p2x*p2y*p3x*p3y*xa*ya*yb + 
        pow(p1x,2)*pow(p3y,2)*xa*ya*yb - pow(p2x,2)*pow(p3y,2)*xa*ya*yb - 
        pow(p1y,2)*p2x*p3x*xb*ya*yb + p1x*p1y*p2y*p3x*xb*ya*yb + 
        p1y*p2x*p2y*p3x*xb*ya*yb - p1x*pow(p2y,2)*p3x*xb*ya*yb + 
        pow(p1y,2)*pow(p3x,2)*xb*ya*yb - pow(p2y,2)*pow(p3x,2)*xb*ya*yb + 
        p1x*p1y*p2x*p3y*xb*ya*yb - p1y*pow(p2x,2)*p3y*xb*ya*yb - 
        pow(p1x,2)*p2y*p3y*xb*ya*yb + p1x*p2x*p2y*p3y*xb*ya*yb - 
        2*p1x*p1y*p3x*p3y*xb*ya*yb + 2*p2x*p2y*p3x*p3y*xb*ya*yb + 
        pow(p1x,2)*pow(p3y,2)*xb*ya*yb - pow(p2x,2)*pow(p3y,2)*xb*ya*yb + 
        2*pow(p1y,2)*p2x*p3x*xa*xb*ya*yb - 2*p1x*p1y*p2y*p3x*xa*xb*ya*yb - 
        2*p1y*p2x*p2y*p3x*xa*xb*ya*yb + 2*p1x*pow(p2y,2)*p3x*xa*xb*ya*yb - 
        2*pow(p1y,2)*pow(p3x,2)*xa*xb*ya*yb + 
        2*p1y*p2y*pow(p3x,2)*xa*xb*ya*yb - 2*p1x*p1y*p2x*p3y*xa*xb*ya*yb + 
        2*p1y*pow(p2x,2)*p3y*xa*xb*ya*yb + 2*pow(p1x,2)*p2y*p3y*xa*xb*ya*yb - 
        2*p1x*p2x*p2y*p3y*xa*xb*ya*yb + 4*p1x*p1y*p3x*p3y*xa*xb*ya*yb - 
        2*p1y*p2x*p3x*p3y*xa*xb*ya*yb - 2*p1x*p2y*p3x*p3y*xa*xb*ya*yb - 
        2*pow(p1x,2)*pow(p3y,2)*xa*xb*ya*yb + 
        2*p1x*p2x*pow(p3y,2)*xa*xb*ya*yb + p1y*p2y*pow(p3x,2)*pow(yb,2) - 
        pow(p2y,2)*pow(p3x,2)*pow(yb,2) - p1y*p2x*p3x*p3y*pow(yb,2) - 
        p1x*p2y*p3x*p3y*pow(yb,2) + 2*p2x*p2y*p3x*p3y*pow(yb,2) + 
        p1x*p2x*pow(p3y,2)*pow(yb,2) - pow(p2x,2)*pow(p3y,2)*pow(yb,2) + 
        pow(p1y,2)*p2x*p3x*xa*pow(yb,2) - p1x*p1y*p2y*p3x*xa*pow(yb,2) - 
        p1y*p2x*p2y*p3x*xa*pow(yb,2) + p1x*pow(p2y,2)*p3x*xa*pow(yb,2) - 
        pow(p1y,2)*pow(p3x,2)*xa*pow(yb,2) + 
        pow(p2y,2)*pow(p3x,2)*xa*pow(yb,2) - p1x*p1y*p2x*p3y*xa*pow(yb,2) + 
        p1y*pow(p2x,2)*p3y*xa*pow(yb,2) + pow(p1x,2)*p2y*p3y*xa*pow(yb,2) - 
        p1x*p2x*p2y*p3y*xa*pow(yb,2) + 2*p1x*p1y*p3x*p3y*xa*pow(yb,2) - 
        2*p2x*p2y*p3x*p3y*xa*pow(yb,2) - 
        pow(p1x,2)*pow(p3y,2)*xa*pow(yb,2) + 
        pow(p2x,2)*pow(p3y,2)*xa*pow(yb,2) - 
        pow(p1y,2)*p2x*p3x*pow(xa,2)*pow(yb,2) + 
        p1x*p1y*p2y*p3x*pow(xa,2)*pow(yb,2) + 
        p1y*p2x*p2y*p3x*pow(xa,2)*pow(yb,2) - 
        p1x*pow(p2y,2)*p3x*pow(xa,2)*pow(yb,2) + 
        pow(p1y,2)*pow(p3x,2)*pow(xa,2)*pow(yb,2) - 
        p1y*p2y*pow(p3x,2)*pow(xa,2)*pow(yb,2) + 
        p1x*p1y*p2x*p3y*pow(xa,2)*pow(yb,2) - 
        p1y*pow(p2x,2)*p3y*pow(xa,2)*pow(yb,2) - 
        pow(p1x,2)*p2y*p3y*pow(xa,2)*pow(yb,2) + 
        p1x*p2x*p2y*p3y*pow(xa,2)*pow(yb,2) - 
        2*p1x*p1y*p3x*p3y*pow(xa,2)*pow(yb,2) + 
        p1y*p2x*p3x*p3y*pow(xa,2)*pow(yb,2) + 
        p1x*p2y*p3x*p3y*pow(xa,2)*pow(yb,2) + 
        pow(p1x,2)*pow(p3y,2)*pow(xa,2)*pow(yb,2) - 
        p1x*p2x*pow(p3y,2)*pow(xa,2)*pow(yb,2) + 
        p0x*(xb + ya - xb*ya + xa*(-1 + yb) - yb)*
         (pow(p1y,2)*(p2x - p3x)*(xb - xb*ya + xa*(-1 + yb)) + 
           p1x*(p2y - p3y)*(p2y*(-(xb*(-2 + ya)) + xa*(-2 + yb)) + 
              p3y*(ya - xb*ya + (-1 + xa)*yb)) - 
           (p2y*p3x - p2x*p3y)*(p2y*((-2 + xb)*ya - (-2 + xa)*yb) + 
              p3y*(ya - xb*ya + (-1 + xa)*yb)) + 
           p1y*(p2x*p2y*xb*(-2 + ya) + p2x*p3y*(-xa + xb + ya - yb) - 
              p2x*p2y*xa*(-2 + yb) + p2y*p3x*(-xa + xb - 2*ya + 2*yb) + 
              p3x*p3y*(ya - xb*ya + (-1 + xa)*yb) + 
              p1x*(p2y - p3y)*(xa + xb*(-1 + ya) - xa*yb))) - 
        p0y*(xb + ya - xb*ya + xa*(-1 + yb) - yb)*
         (2*p1y*pow(p2x,2)*xa - 2*p1y*p2x*p3x*xa - 2*p1y*pow(p2x,2)*xb + 
           2*p1y*p2x*p3x*xb - p1y*p2x*p3x*ya + 2*p2x*p2y*p3x*ya + 
           p1y*pow(p3x,2)*ya - p2y*pow(p3x,2)*ya - 2*pow(p2x,2)*p3y*ya + 
           p2x*p3x*p3y*ya + p1y*pow(p2x,2)*xb*ya - p2x*p2y*p3x*xb*ya - 
           p1y*pow(p3x,2)*xb*ya + p2y*pow(p3x,2)*xb*ya + 
           pow(p2x,2)*p3y*xb*ya - p2x*p3x*p3y*xb*ya + 
           p0x*(-2*p2x*p2y + p1y*(p2x - p3x) + p2y*p3x + p1x*(p2y - p3y) + p2x*p3y)*
            (xb + ya - xb*ya + xa*(-1 + yb) - yb) + p1y*p2x*p3x*yb - 
           2*p2x*p2y*p3x*yb - p1y*pow(p3x,2)*yb + p2y*pow(p3x,2)*yb + 
           2*pow(p2x,2)*p3y*yb - p2x*p3x*p3y*yb - p1y*pow(p2x,2)*xa*yb + 
           p2x*p2y*p3x*xa*yb + p1y*pow(p3x,2)*xa*yb - p2y*pow(p3x,2)*xa*yb - 
           pow(p2x,2)*p3y*xa*yb + p2x*p3x*p3y*xa*yb + 
           pow(p1x,2)*(p2y - p3y)*(xa + xb*(-1 + ya) - xa*yb) + 
           p1x*(-(p2x*p2y*xb*(-2 + ya)) + 
              p1y*(p2x - p3x)*(xb - xb*ya + xa*(-1 + yb)) + 
              p2x*p3y*(xa - xb + 2*ya - 2*yb) + p2x*p2y*xa*(-2 + yb) + 
              p2y*p3x*(xa - xb - ya + yb) + p3x*p3y*((-1 + xb)*ya + yb - xa*yb))))*
      log(p1y*p3x - p1x*p3y + p1y*p2x*xa - p1x*p2y*xa - p1y*p3x*xa + p1x*p3y*xa + 
        p0x*(p3y - p3y*xa + p2y*(xa - ya) + p1y*(-1 + ya)) - p1y*p3x*ya + 
        p2y*p3x*ya + p1x*p3y*ya - p2x*p3y*ya + 
        p0y*(p1x + p3x*(-1 + xa) - p1x*ya + p2x*(-xa + ya))) + 
     2*(-(pow(p1y,2)*pow(p2x,2)*pow(xa,2)) + 2*p1x*p1y*p2x*p2y*pow(xa,2) - 
        pow(p1x,2)*pow(p2y,2)*pow(xa,2) + pow(p1y,2)*p2x*p3x*pow(xa,2) - 
        p1x*p1y*p2y*p3x*pow(xa,2) - p1x*p1y*p2x*p3y*pow(xa,2) + 
        pow(p1x,2)*p2y*p3y*pow(xa,2) + 2*pow(p1y,2)*pow(p2x,2)*xa*xb - 
        4*p1x*p1y*p2x*p2y*xa*xb + 2*pow(p1x,2)*pow(p2y,2)*xa*xb - 
        2*pow(p1y,2)*p2x*p3x*xa*xb + 2*p1x*p1y*p2y*p3x*xa*xb + 
        2*p1x*p1y*p2x*p3y*xa*xb - 2*pow(p1x,2)*p2y*p3y*xa*xb - 
        pow(p1y,2)*pow(p2x,2)*pow(xb,2) + 2*p1x*p1y*p2x*p2y*pow(xb,2) - 
        pow(p1x,2)*pow(p2y,2)*pow(xb,2) + pow(p1y,2)*p2x*p3x*pow(xb,2) - 
        p1x*p1y*p2y*p3x*pow(xb,2) - p1x*p1y*p2x*p3y*pow(xb,2) + 
        pow(p1x,2)*p2y*p3y*pow(xb,2) + pow(p1y,2)*p2x*p3x*xa*ya - 
        p1x*p1y*p2y*p3x*xa*ya - 2*p1y*p2x*p2y*p3x*xa*ya + 
        2*p1x*pow(p2y,2)*p3x*xa*ya - pow(p1y,2)*pow(p3x,2)*xa*ya + 
        p1y*p2y*pow(p3x,2)*xa*ya - p1x*p1y*p2x*p3y*xa*ya + 
        2*p1y*pow(p2x,2)*p3y*xa*ya + pow(p1x,2)*p2y*p3y*xa*ya - 
        2*p1x*p2x*p2y*p3y*xa*ya + 2*p1x*p1y*p3x*p3y*xa*ya - p1y*p2x*p3x*p3y*xa*ya - 
        p1x*p2y*p3x*p3y*xa*ya - pow(p1x,2)*pow(p3y,2)*xa*ya + 
        p1x*p2x*pow(p3y,2)*xa*ya - pow(p1y,2)*p2x*p3x*xb*ya + 
        p1x*p1y*p2y*p3x*xb*ya + 2*p1y*p2x*p2y*p3x*xb*ya - 
        2*p1x*pow(p2y,2)*p3x*xb*ya + pow(p1y,2)*pow(p3x,2)*xb*ya - 
        p1y*p2y*pow(p3x,2)*xb*ya + p1x*p1y*p2x*p3y*xb*ya - 
        2*p1y*pow(p2x,2)*p3y*xb*ya - pow(p1x,2)*p2y*p3y*xb*ya + 
        2*p1x*p2x*p2y*p3y*xb*ya - 2*p1x*p1y*p3x*p3y*xb*ya + p1y*p2x*p3x*p3y*xb*ya + 
        p1x*p2y*p3x*p3y*xb*ya + pow(p1x,2)*pow(p3y,2)*xb*ya - 
        p1x*p2x*pow(p3y,2)*xb*ya - pow(p1y,2)*pow(p2x,2)*xa*xb*ya + 
        2*p1x*p1y*p2x*p2y*xa*xb*ya - pow(p1x,2)*pow(p2y,2)*xa*xb*ya + 
        p1y*p2x*p2y*p3x*xa*xb*ya - p1x*pow(p2y,2)*p3x*xa*xb*ya + 
        pow(p1y,2)*pow(p3x,2)*xa*xb*ya - p1y*p2y*pow(p3x,2)*xa*xb*ya - 
        p1y*pow(p2x,2)*p3y*xa*xb*ya + p1x*p2x*p2y*p3y*xa*xb*ya - 
        2*p1x*p1y*p3x*p3y*xa*xb*ya + p1y*p2x*p3x*p3y*xa*xb*ya + 
        p1x*p2y*p3x*p3y*xa*xb*ya + pow(p1x,2)*pow(p3y,2)*xa*xb*ya - 
        p1x*p2x*pow(p3y,2)*xa*xb*ya + pow(p1y,2)*pow(p2x,2)*pow(xb,2)*ya - 
        2*p1x*p1y*p2x*p2y*pow(xb,2)*ya + 
        pow(p1x,2)*pow(p2y,2)*pow(xb,2)*ya - p1y*p2x*p2y*p3x*pow(xb,2)*ya + 
        p1x*pow(p2y,2)*p3x*pow(xb,2)*ya - 
        pow(p1y,2)*pow(p3x,2)*pow(xb,2)*ya + 
        p1y*p2y*pow(p3x,2)*pow(xb,2)*ya + p1y*pow(p2x,2)*p3y*pow(xb,2)*ya - 
        p1x*p2x*p2y*p3y*pow(xb,2)*ya + 2*p1x*p1y*p3x*p3y*pow(xb,2)*ya - 
        p1y*p2x*p3x*p3y*pow(xb,2)*ya - p1x*p2y*p3x*p3y*pow(xb,2)*ya - 
        pow(p1x,2)*pow(p3y,2)*pow(xb,2)*ya + 
        p1x*p2x*pow(p3y,2)*pow(xb,2)*ya + p1y*p2y*pow(p3x,2)*pow(ya,2) - 
        pow(p2y,2)*pow(p3x,2)*pow(ya,2) - p1y*p2x*p3x*p3y*pow(ya,2) - 
        p1x*p2y*p3x*p3y*pow(ya,2) + 2*p2x*p2y*p3x*p3y*pow(ya,2) + 
        p1x*p2x*pow(p3y,2)*pow(ya,2) - pow(p2x,2)*pow(p3y,2)*pow(ya,2) + 
        pow(p1y,2)*p2x*p3x*xb*pow(ya,2) - p1x*p1y*p2y*p3x*xb*pow(ya,2) - 
        p1y*p2x*p2y*p3x*xb*pow(ya,2) + p1x*pow(p2y,2)*p3x*xb*pow(ya,2) - 
        pow(p1y,2)*pow(p3x,2)*xb*pow(ya,2) + 
        pow(p2y,2)*pow(p3x,2)*xb*pow(ya,2) - p1x*p1y*p2x*p3y*xb*pow(ya,2) + 
        p1y*pow(p2x,2)*p3y*xb*pow(ya,2) + pow(p1x,2)*p2y*p3y*xb*pow(ya,2) - 
        p1x*p2x*p2y*p3y*xb*pow(ya,2) + 2*p1x*p1y*p3x*p3y*xb*pow(ya,2) - 
        2*p2x*p2y*p3x*p3y*xb*pow(ya,2) - 
        pow(p1x,2)*pow(p3y,2)*xb*pow(ya,2) + 
        pow(p2x,2)*pow(p3y,2)*xb*pow(ya,2) - 
        pow(p1y,2)*p2x*p3x*pow(xb,2)*pow(ya,2) + 
        p1x*p1y*p2y*p3x*pow(xb,2)*pow(ya,2) + 
        p1y*p2x*p2y*p3x*pow(xb,2)*pow(ya,2) - 
        p1x*pow(p2y,2)*p3x*pow(xb,2)*pow(ya,2) + 
        pow(p1y,2)*pow(p3x,2)*pow(xb,2)*pow(ya,2) - 
        p1y*p2y*pow(p3x,2)*pow(xb,2)*pow(ya,2) + 
        p1x*p1y*p2x*p3y*pow(xb,2)*pow(ya,2) - 
        p1y*pow(p2x,2)*p3y*pow(xb,2)*pow(ya,2) - 
        pow(p1x,2)*p2y*p3y*pow(xb,2)*pow(ya,2) + 
        p1x*p2x*p2y*p3y*pow(xb,2)*pow(ya,2) - 
        2*p1x*p1y*p3x*p3y*pow(xb,2)*pow(ya,2) + 
        p1y*p2x*p3x*p3y*pow(xb,2)*pow(ya,2) + 
        p1x*p2y*p3x*p3y*pow(xb,2)*pow(ya,2) + 
        pow(p1x,2)*pow(p3y,2)*pow(xb,2)*pow(ya,2) - 
        p1x*p2x*pow(p3y,2)*pow(xb,2)*pow(ya,2) + 
        pow(p0y,2)*(p1x - p2x)*(p2x - p3x)*
         pow(xb + ya - xb*ya + xa*(-1 + yb) - yb,2) + 
        pow(p0x,2)*(p1y - p2y)*(p2y - p3y)*
         pow(xb + ya - xb*ya + xa*(-1 + yb) - yb,2) - pow(p1y,2)*p2x*p3x*xa*yb + 
        p1x*p1y*p2y*p3x*xa*yb + 2*p1y*p2x*p2y*p3x*xa*yb - 
        2*p1x*pow(p2y,2)*p3x*xa*yb + pow(p1y,2)*pow(p3x,2)*xa*yb - 
        p1y*p2y*pow(p3x,2)*xa*yb + p1x*p1y*p2x*p3y*xa*yb - 
        2*p1y*pow(p2x,2)*p3y*xa*yb - pow(p1x,2)*p2y*p3y*xa*yb + 
        2*p1x*p2x*p2y*p3y*xa*yb - 2*p1x*p1y*p3x*p3y*xa*yb + p1y*p2x*p3x*p3y*xa*yb + 
        p1x*p2y*p3x*p3y*xa*yb + pow(p1x,2)*pow(p3y,2)*xa*yb - 
        p1x*p2x*pow(p3y,2)*xa*yb + pow(p1y,2)*pow(p2x,2)*pow(xa,2)*yb - 
        2*p1x*p1y*p2x*p2y*pow(xa,2)*yb + 
        pow(p1x,2)*pow(p2y,2)*pow(xa,2)*yb - p1y*p2x*p2y*p3x*pow(xa,2)*yb + 
        p1x*pow(p2y,2)*p3x*pow(xa,2)*yb - 
        pow(p1y,2)*pow(p3x,2)*pow(xa,2)*yb + 
        p1y*p2y*pow(p3x,2)*pow(xa,2)*yb + p1y*pow(p2x,2)*p3y*pow(xa,2)*yb - 
        p1x*p2x*p2y*p3y*pow(xa,2)*yb + 2*p1x*p1y*p3x*p3y*pow(xa,2)*yb - 
        p1y*p2x*p3x*p3y*pow(xa,2)*yb - p1x*p2y*p3x*p3y*pow(xa,2)*yb - 
        pow(p1x,2)*pow(p3y,2)*pow(xa,2)*yb + 
        p1x*p2x*pow(p3y,2)*pow(xa,2)*yb + pow(p1y,2)*p2x*p3x*xb*yb - 
        p1x*p1y*p2y*p3x*xb*yb - 2*p1y*p2x*p2y*p3x*xb*yb + 
        2*p1x*pow(p2y,2)*p3x*xb*yb - pow(p1y,2)*pow(p3x,2)*xb*yb + 
        p1y*p2y*pow(p3x,2)*xb*yb - p1x*p1y*p2x*p3y*xb*yb + 
        2*p1y*pow(p2x,2)*p3y*xb*yb + pow(p1x,2)*p2y*p3y*xb*yb - 
        2*p1x*p2x*p2y*p3y*xb*yb + 2*p1x*p1y*p3x*p3y*xb*yb - p1y*p2x*p3x*p3y*xb*yb - 
        p1x*p2y*p3x*p3y*xb*yb - pow(p1x,2)*pow(p3y,2)*xb*yb + 
        p1x*p2x*pow(p3y,2)*xb*yb - pow(p1y,2)*pow(p2x,2)*xa*xb*yb + 
        2*p1x*p1y*p2x*p2y*xa*xb*yb - pow(p1x,2)*pow(p2y,2)*xa*xb*yb + 
        p1y*p2x*p2y*p3x*xa*xb*yb - p1x*pow(p2y,2)*p3x*xa*xb*yb + 
        pow(p1y,2)*pow(p3x,2)*xa*xb*yb - p1y*p2y*pow(p3x,2)*xa*xb*yb - 
        p1y*pow(p2x,2)*p3y*xa*xb*yb + p1x*p2x*p2y*p3y*xa*xb*yb - 
        2*p1x*p1y*p3x*p3y*xa*xb*yb + p1y*p2x*p3x*p3y*xa*xb*yb + 
        p1x*p2y*p3x*p3y*xa*xb*yb + pow(p1x,2)*pow(p3y,2)*xa*xb*yb - 
        p1x*p2x*pow(p3y,2)*xa*xb*yb - 2*p1y*p2y*pow(p3x,2)*ya*yb + 
        2*pow(p2y,2)*pow(p3x,2)*ya*yb + 2*p1y*p2x*p3x*p3y*ya*yb + 
        2*p1x*p2y*p3x*p3y*ya*yb - 4*p2x*p2y*p3x*p3y*ya*yb - 
        2*p1x*p2x*pow(p3y,2)*ya*yb + 2*pow(p2x,2)*pow(p3y,2)*ya*yb - 
        pow(p1y,2)*p2x*p3x*xa*ya*yb + p1x*p1y*p2y*p3x*xa*ya*yb + 
        p1y*p2x*p2y*p3x*xa*ya*yb - p1x*pow(p2y,2)*p3x*xa*ya*yb + 
        pow(p1y,2)*pow(p3x,2)*xa*ya*yb - pow(p2y,2)*pow(p3x,2)*xa*ya*yb + 
        p1x*p1y*p2x*p3y*xa*ya*yb - p1y*pow(p2x,2)*p3y*xa*ya*yb - 
        pow(p1x,2)*p2y*p3y*xa*ya*yb + p1x*p2x*p2y*p3y*xa*ya*yb - 
        2*p1x*p1y*p3x*p3y*xa*ya*yb + 2*p2x*p2y*p3x*p3y*xa*ya*yb + 
        pow(p1x,2)*pow(p3y,2)*xa*ya*yb - pow(p2x,2)*pow(p3y,2)*xa*ya*yb - 
        pow(p1y,2)*p2x*p3x*xb*ya*yb + p1x*p1y*p2y*p3x*xb*ya*yb + 
        p1y*p2x*p2y*p3x*xb*ya*yb - p1x*pow(p2y,2)*p3x*xb*ya*yb + 
        pow(p1y,2)*pow(p3x,2)*xb*ya*yb - pow(p2y,2)*pow(p3x,2)*xb*ya*yb + 
        p1x*p1y*p2x*p3y*xb*ya*yb - p1y*pow(p2x,2)*p3y*xb*ya*yb - 
        pow(p1x,2)*p2y*p3y*xb*ya*yb + p1x*p2x*p2y*p3y*xb*ya*yb - 
        2*p1x*p1y*p3x*p3y*xb*ya*yb + 2*p2x*p2y*p3x*p3y*xb*ya*yb + 
        pow(p1x,2)*pow(p3y,2)*xb*ya*yb - pow(p2x,2)*pow(p3y,2)*xb*ya*yb + 
        2*pow(p1y,2)*p2x*p3x*xa*xb*ya*yb - 2*p1x*p1y*p2y*p3x*xa*xb*ya*yb - 
        2*p1y*p2x*p2y*p3x*xa*xb*ya*yb + 2*p1x*pow(p2y,2)*p3x*xa*xb*ya*yb - 
        2*pow(p1y,2)*pow(p3x,2)*xa*xb*ya*yb + 
        2*p1y*p2y*pow(p3x,2)*xa*xb*ya*yb - 2*p1x*p1y*p2x*p3y*xa*xb*ya*yb + 
        2*p1y*pow(p2x,2)*p3y*xa*xb*ya*yb + 2*pow(p1x,2)*p2y*p3y*xa*xb*ya*yb - 
        2*p1x*p2x*p2y*p3y*xa*xb*ya*yb + 4*p1x*p1y*p3x*p3y*xa*xb*ya*yb - 
        2*p1y*p2x*p3x*p3y*xa*xb*ya*yb - 2*p1x*p2y*p3x*p3y*xa*xb*ya*yb - 
        2*pow(p1x,2)*pow(p3y,2)*xa*xb*ya*yb + 
        2*p1x*p2x*pow(p3y,2)*xa*xb*ya*yb + p1y*p2y*pow(p3x,2)*pow(yb,2) - 
        pow(p2y,2)*pow(p3x,2)*pow(yb,2) - p1y*p2x*p3x*p3y*pow(yb,2) - 
        p1x*p2y*p3x*p3y*pow(yb,2) + 2*p2x*p2y*p3x*p3y*pow(yb,2) + 
        p1x*p2x*pow(p3y,2)*pow(yb,2) - pow(p2x,2)*pow(p3y,2)*pow(yb,2) + 
        pow(p1y,2)*p2x*p3x*xa*pow(yb,2) - p1x*p1y*p2y*p3x*xa*pow(yb,2) - 
        p1y*p2x*p2y*p3x*xa*pow(yb,2) + p1x*pow(p2y,2)*p3x*xa*pow(yb,2) - 
        pow(p1y,2)*pow(p3x,2)*xa*pow(yb,2) + 
        pow(p2y,2)*pow(p3x,2)*xa*pow(yb,2) - p1x*p1y*p2x*p3y*xa*pow(yb,2) + 
        p1y*pow(p2x,2)*p3y*xa*pow(yb,2) + pow(p1x,2)*p2y*p3y*xa*pow(yb,2) - 
        p1x*p2x*p2y*p3y*xa*pow(yb,2) + 2*p1x*p1y*p3x*p3y*xa*pow(yb,2) - 
        2*p2x*p2y*p3x*p3y*xa*pow(yb,2) - 
        pow(p1x,2)*pow(p3y,2)*xa*pow(yb,2) + 
        pow(p2x,2)*pow(p3y,2)*xa*pow(yb,2) - 
        pow(p1y,2)*p2x*p3x*pow(xa,2)*pow(yb,2) + 
        p1x*p1y*p2y*p3x*pow(xa,2)*pow(yb,2) + 
        p1y*p2x*p2y*p3x*pow(xa,2)*pow(yb,2) - 
        p1x*pow(p2y,2)*p3x*pow(xa,2)*pow(yb,2) + 
        pow(p1y,2)*pow(p3x,2)*pow(xa,2)*pow(yb,2) - 
        p1y*p2y*pow(p3x,2)*pow(xa,2)*pow(yb,2) + 
        p1x*p1y*p2x*p3y*pow(xa,2)*pow(yb,2) - 
        p1y*pow(p2x,2)*p3y*pow(xa,2)*pow(yb,2) - 
        pow(p1x,2)*p2y*p3y*pow(xa,2)*pow(yb,2) + 
        p1x*p2x*p2y*p3y*pow(xa,2)*pow(yb,2) - 
        2*p1x*p1y*p3x*p3y*pow(xa,2)*pow(yb,2) + 
        p1y*p2x*p3x*p3y*pow(xa,2)*pow(yb,2) + 
        p1x*p2y*p3x*p3y*pow(xa,2)*pow(yb,2) + 
        pow(p1x,2)*pow(p3y,2)*pow(xa,2)*pow(yb,2) - 
        p1x*p2x*pow(p3y,2)*pow(xa,2)*pow(yb,2) + 
        p0x*(xb + ya - xb*ya + xa*(-1 + yb) - yb)*
         (pow(p1y,2)*(p2x - p3x)*(xb - xb*ya + xa*(-1 + yb)) + 
           p1x*(p2y - p3y)*(p2y*(-(xb*(-2 + ya)) + xa*(-2 + yb)) + 
              p3y*(ya - xb*ya + (-1 + xa)*yb)) - 
           (p2y*p3x - p2x*p3y)*(p2y*((-2 + xb)*ya - (-2 + xa)*yb) + 
              p3y*(ya - xb*ya + (-1 + xa)*yb)) + 
           p1y*(p2x*p2y*xb*(-2 + ya) + p2x*p3y*(-xa + xb + ya - yb) - 
              p2x*p2y*xa*(-2 + yb) + p2y*p3x*(-xa + xb - 2*ya + 2*yb) + 
              p3x*p3y*(ya - xb*ya + (-1 + xa)*yb) + 
              p1x*(p2y - p3y)*(xa + xb*(-1 + ya) - xa*yb))) - 
        p0y*(xb + ya - xb*ya + xa*(-1 + yb) - yb)*
         (2*p1y*pow(p2x,2)*xa - 2*p1y*p2x*p3x*xa - 2*p1y*pow(p2x,2)*xb + 
           2*p1y*p2x*p3x*xb - p1y*p2x*p3x*ya + 2*p2x*p2y*p3x*ya + 
           p1y*pow(p3x,2)*ya - p2y*pow(p3x,2)*ya - 2*pow(p2x,2)*p3y*ya + 
           p2x*p3x*p3y*ya + p1y*pow(p2x,2)*xb*ya - p2x*p2y*p3x*xb*ya - 
           p1y*pow(p3x,2)*xb*ya + p2y*pow(p3x,2)*xb*ya + 
           pow(p2x,2)*p3y*xb*ya - p2x*p3x*p3y*xb*ya + 
           p0x*(-2*p2x*p2y + p1y*(p2x - p3x) + p2y*p3x + p1x*(p2y - p3y) + p2x*p3y)*
            (xb + ya - xb*ya + xa*(-1 + yb) - yb) + p1y*p2x*p3x*yb - 
           2*p2x*p2y*p3x*yb - p1y*pow(p3x,2)*yb + p2y*pow(p3x,2)*yb + 
           2*pow(p2x,2)*p3y*yb - p2x*p3x*p3y*yb - p1y*pow(p2x,2)*xa*yb + 
           p2x*p2y*p3x*xa*yb + p1y*pow(p3x,2)*xa*yb - p2y*pow(p3x,2)*xa*yb - 
           pow(p2x,2)*p3y*xa*yb + p2x*p3x*p3y*xa*yb + 
           pow(p1x,2)*(p2y - p3y)*(xa + xb*(-1 + ya) - xa*yb) + 
           p1x*(-(p2x*p2y*xb*(-2 + ya)) + 
              p1y*(p2x - p3x)*(xb - xb*ya + xa*(-1 + yb)) + 
              p2x*p3y*(xa - xb + 2*ya - 2*yb) + p2x*p2y*xa*(-2 + yb) + 
              p2y*p3x*(xa - xb - ya + yb) + p3x*p3y*((-1 + xb)*ya + yb - xa*yb))))*
      log(p1y*p3x - p1x*p3y + p1y*p2x*xb - p1x*p2y*xb - p1y*p3x*xb + p1x*p3y*xb + 
        p0x*(p3y - p3y*xb + p2y*(xb - yb) + p1y*(-1 + yb)) - p1y*p3x*yb + 
        p2y*p3x*yb + p1x*p3y*yb - p2x*p3y*yb + 
        p0y*(p1x + p3x*(-1 + xb) - p1x*yb + p2x*(-xb + yb))))/
   (2.*pow(-(p0x*p2y*xa) + p1x*p2y*xa + p0x*p3y*xa - p1x*p3y*xa + p0x*p2y*xb - 
       p1x*p2y*xb - p0x*p3y*xb + p1x*p3y*xb + p0x*p2y*ya - p2y*p3x*ya - p1x*p3y*ya + 
       p2x*p3y*ya - p0x*p2y*yb + p2y*p3x*yb + p1x*p3y*yb - p2x*p3y*yb + 
       p1y*(p2x*(-xa + xb) + p3x*(xa - xb + ya - yb) + p0x*(-ya + yb)) + 
       p0y*(p3x*(-xa + xb) + p1x*(ya - yb) + p2x*(xa - xb - ya + yb)),3));

    integrals[1] = ((xa - xb)*(ya - yb)*pow(-(p0x*p2y*xa) + p1x*p2y*xa + p0x*p3y*xa - p1x*p3y*xa + 
        p0x*p2y*xb - p1x*p2y*xb - p0x*p3y*xb + p1x*p3y*xb + p0x*p2y*ya - 
        p2y*p3x*ya - p1x*p3y*ya + p2x*p3y*ya - p0x*p2y*yb + p2y*p3x*yb + 
        p1x*p3y*yb - p2x*p3y*yb + p1y*
         (p2x*(-xa + xb) + p3x*(xa - xb + ya - yb) + p0x*(-ya + yb)) + 
        p0y*(p3x*(-xa + xb) + p1x*(ya - yb) + p2x*(xa - xb - ya + yb)),2) + 
     2*(-(p0x*p2y*xa) + p1x*p2y*xa + p0x*p3y*xa - p1x*p3y*xa + p0x*p2y*xb - 
        p1x*p2y*xb - p0x*p3y*xb + p1x*p3y*xb + p0x*p2y*ya - p2y*p3x*ya - 
        p1x*p3y*ya + p2x*p3y*ya - p0x*p2y*yb + p2y*p3x*yb + p1x*p3y*yb - 
        p2x*p3y*yb + p1y*(p2x*(-xa + xb) + p3x*(xa - xb + ya - yb) + 
           p0x*(-ya + yb)) + p0y*(p3x*(-xa + xb) + p1x*(ya - yb) + 
           p2x*(xa - xb - ya + yb)))*
      (-(p0x*p2y*pow(xa,2)) + p1x*p2y*pow(xa,2) + p0x*p3y*pow(xa,2) - 
        p1x*p3y*pow(xa,2) + 2*p0x*p2y*xa*xb - 2*p1x*p2y*xa*xb - 2*p0x*p3y*xa*xb + 
        2*p1x*p3y*xa*xb - p0x*p2y*pow(xb,2) + p1x*p2y*pow(xb,2) + 
        p0x*p3y*pow(xb,2) - p1x*p3y*pow(xb,2) + p0x*p2y*xa*ya - p2y*p3x*xa*ya - 
        p0x*p3y*xa*ya + p2x*p3y*xa*ya + p0x*p2y*pow(xa,2)*ya - 
        p1x*p2y*pow(xa,2)*ya - p0x*p3y*pow(xa,2)*ya + p1x*p3y*pow(xa,2)*ya - 
        p0x*p2y*xb*ya + p2y*p3x*xb*ya + p0x*p3y*xb*ya - p2x*p3y*xb*ya - 
        2*p0x*p2y*xa*xb*ya + 2*p1x*p2y*xa*xb*ya + 2*p0x*p3y*xa*xb*ya - 
        2*p1x*p3y*xa*xb*ya + p0x*p2y*pow(xb,2)*ya - p1x*p2y*pow(xb,2)*ya - 
        p0x*p3y*pow(xb,2)*ya + p1x*p3y*pow(xb,2)*ya - p0x*p2y*xa*pow(ya,2) + 
        p2y*p3x*xa*pow(ya,2) + p1x*p3y*xa*pow(ya,2) - p2x*p3y*xa*pow(ya,2) + 
        p1y*(p2x*pow(xa - xb,2)*(-1 + ya) + 
           p3x*(-(pow(xa,2)*(-1 + ya)) + 2*xa*xb*(-1 + ya) - 
              pow(xb,2)*(-1 + ya) - xa*pow(ya - yb,2)) + p0x*xa*pow(ya - yb,2))
          - p0x*p2y*xa*yb + p2y*p3x*xa*yb + p0x*p3y*xa*yb - p2x*p3y*xa*yb + 
        p0x*p2y*xb*yb - p2y*p3x*xb*yb - p0x*p3y*xb*yb + p2x*p3y*xb*yb + 
        2*p0x*p2y*xa*ya*yb - 2*p2y*p3x*xa*ya*yb - 2*p1x*p3y*xa*ya*yb + 
        2*p2x*p3y*xa*ya*yb - p0x*p2y*xa*pow(yb,2) + p2y*p3x*xa*pow(yb,2) + 
        p1x*p3y*xa*pow(yb,2) - p2x*p3y*xa*pow(yb,2) + 
        p0y*(-(p1x*xa*pow(ya - yb,2)) + 
           p3x*(xa - xb)*(xb + xa*(-1 + ya) + ya - xb*ya - yb) + 
           p2x*(-(pow(xa,2)*(-1 + ya)) + xb*(xb + ya - xb*ya - yb) + 
              xa*(2*xb*(-1 + ya) - ya + pow(ya,2) + yb - 2*ya*yb + pow(yb,2)))))\
      + 2*(-(pow(p1y,2)*p2x*p3x*pow(xa,2)) + p1x*p1y*p2y*p3x*pow(xa,2) + 
        pow(p1y,2)*pow(p3x,2)*pow(xa,2) + p1x*p1y*p2x*p3y*pow(xa,2) - 
        pow(p1x,2)*p2y*p3y*pow(xa,2) - 2*p1x*p1y*p3x*p3y*pow(xa,2) + 
        pow(p1x,2)*pow(p3y,2)*pow(xa,2) + 2*pow(p1y,2)*p2x*p3x*xa*xb - 
        2*p1x*p1y*p2y*p3x*xa*xb - 2*pow(p1y,2)*pow(p3x,2)*xa*xb - 
        2*p1x*p1y*p2x*p3y*xa*xb + 2*pow(p1x,2)*p2y*p3y*xa*xb + 
        4*p1x*p1y*p3x*p3y*xa*xb - 2*pow(p1x,2)*pow(p3y,2)*xa*xb - 
        pow(p1y,2)*p2x*p3x*pow(xb,2) + p1x*p1y*p2y*p3x*pow(xb,2) + 
        pow(p1y,2)*pow(p3x,2)*pow(xb,2) + p1x*p1y*p2x*p3y*pow(xb,2) - 
        pow(p1x,2)*p2y*p3y*pow(xb,2) - 2*p1x*p1y*p3x*p3y*pow(xb,2) + 
        pow(p1x,2)*pow(p3y,2)*pow(xb,2) - p1y*p2y*pow(p3x,2)*xa*ya + 
        p1y*p2x*p3x*p3y*xa*ya + p1x*p2y*p3x*p3y*xa*ya - p1x*p2x*pow(p3y,2)*xa*ya + 
        p1y*p2y*pow(p3x,2)*xb*ya - p1y*p2x*p3x*p3y*xb*ya - p1x*p2y*p3x*p3y*xb*ya + 
        p1x*p2x*pow(p3y,2)*xb*ya - 2*pow(p1y,2)*p2x*p3x*xa*xb*ya + 
        2*p1x*p1y*p2y*p3x*xa*xb*ya + p1y*p2x*p2y*p3x*xa*xb*ya - 
        p1x*pow(p2y,2)*p3x*xa*xb*ya + 2*pow(p1y,2)*pow(p3x,2)*xa*xb*ya - 
        p1y*p2y*pow(p3x,2)*xa*xb*ya + 2*p1x*p1y*p2x*p3y*xa*xb*ya - 
        p1y*pow(p2x,2)*p3y*xa*xb*ya - 2*pow(p1x,2)*p2y*p3y*xa*xb*ya + 
        p1x*p2x*p2y*p3y*xa*xb*ya - 4*p1x*p1y*p3x*p3y*xa*xb*ya + 
        p1y*p2x*p3x*p3y*xa*xb*ya + p1x*p2y*p3x*p3y*xa*xb*ya + 
        2*pow(p1x,2)*pow(p3y,2)*xa*xb*ya - p1x*p2x*pow(p3y,2)*xa*xb*ya + 
        2*pow(p1y,2)*p2x*p3x*pow(xb,2)*ya - 2*p1x*p1y*p2y*p3x*pow(xb,2)*ya - 
        p1y*p2x*p2y*p3x*pow(xb,2)*ya + p1x*pow(p2y,2)*p3x*pow(xb,2)*ya - 
        2*pow(p1y,2)*pow(p3x,2)*pow(xb,2)*ya + 
        p1y*p2y*pow(p3x,2)*pow(xb,2)*ya - 2*p1x*p1y*p2x*p3y*pow(xb,2)*ya + 
        p1y*pow(p2x,2)*p3y*pow(xb,2)*ya + 
        2*pow(p1x,2)*p2y*p3y*pow(xb,2)*ya - p1x*p2x*p2y*p3y*pow(xb,2)*ya + 
        4*p1x*p1y*p3x*p3y*pow(xb,2)*ya - p1y*p2x*p3x*p3y*pow(xb,2)*ya - 
        p1x*p2y*p3x*p3y*pow(xb,2)*ya - 
        2*pow(p1x,2)*pow(p3y,2)*pow(xb,2)*ya + 
        p1x*p2x*pow(p3y,2)*pow(xb,2)*ya - p1y*p2y*pow(p3x,2)*xb*pow(ya,2) + 
        pow(p2y,2)*pow(p3x,2)*xb*pow(ya,2) + p1y*p2x*p3x*p3y*xb*pow(ya,2) + 
        p1x*p2y*p3x*p3y*xb*pow(ya,2) - 2*p2x*p2y*p3x*p3y*xb*pow(ya,2) - 
        p1x*p2x*pow(p3y,2)*xb*pow(ya,2) + 
        pow(p2x,2)*pow(p3y,2)*xb*pow(ya,2) - 
        pow(p1y,2)*p2x*p3x*pow(xb,2)*pow(ya,2) + 
        p1x*p1y*p2y*p3x*pow(xb,2)*pow(ya,2) + 
        p1y*p2x*p2y*p3x*pow(xb,2)*pow(ya,2) - 
        p1x*pow(p2y,2)*p3x*pow(xb,2)*pow(ya,2) + 
        pow(p1y,2)*pow(p3x,2)*pow(xb,2)*pow(ya,2) - 
        p1y*p2y*pow(p3x,2)*pow(xb,2)*pow(ya,2) + 
        p1x*p1y*p2x*p3y*pow(xb,2)*pow(ya,2) - 
        p1y*pow(p2x,2)*p3y*pow(xb,2)*pow(ya,2) - 
        pow(p1x,2)*p2y*p3y*pow(xb,2)*pow(ya,2) + 
        p1x*p2x*p2y*p3y*pow(xb,2)*pow(ya,2) - 
        2*p1x*p1y*p3x*p3y*pow(xb,2)*pow(ya,2) + 
        p1y*p2x*p3x*p3y*pow(xb,2)*pow(ya,2) + 
        p1x*p2y*p3x*p3y*pow(xb,2)*pow(ya,2) + 
        pow(p1x,2)*pow(p3y,2)*pow(xb,2)*pow(ya,2) - 
        p1x*p2x*pow(p3y,2)*pow(xb,2)*pow(ya,2) + p1y*p2y*pow(p3x,2)*xa*yb - 
        p1y*p2x*p3x*p3y*xa*yb - p1x*p2y*p3x*p3y*xa*yb + p1x*p2x*pow(p3y,2)*xa*yb + 
        2*pow(p1y,2)*p2x*p3x*pow(xa,2)*yb - 2*p1x*p1y*p2y*p3x*pow(xa,2)*yb - 
        p1y*p2x*p2y*p3x*pow(xa,2)*yb + p1x*pow(p2y,2)*p3x*pow(xa,2)*yb - 
        2*pow(p1y,2)*pow(p3x,2)*pow(xa,2)*yb + 
        p1y*p2y*pow(p3x,2)*pow(xa,2)*yb - 2*p1x*p1y*p2x*p3y*pow(xa,2)*yb + 
        p1y*pow(p2x,2)*p3y*pow(xa,2)*yb + 
        2*pow(p1x,2)*p2y*p3y*pow(xa,2)*yb - p1x*p2x*p2y*p3y*pow(xa,2)*yb + 
        4*p1x*p1y*p3x*p3y*pow(xa,2)*yb - p1y*p2x*p3x*p3y*pow(xa,2)*yb - 
        p1x*p2y*p3x*p3y*pow(xa,2)*yb - 
        2*pow(p1x,2)*pow(p3y,2)*pow(xa,2)*yb + 
        p1x*p2x*pow(p3y,2)*pow(xa,2)*yb - p1y*p2y*pow(p3x,2)*xb*yb + 
        p1y*p2x*p3x*p3y*xb*yb + p1x*p2y*p3x*p3y*xb*yb - p1x*p2x*pow(p3y,2)*xb*yb - 
        2*pow(p1y,2)*p2x*p3x*xa*xb*yb + 2*p1x*p1y*p2y*p3x*xa*xb*yb + 
        p1y*p2x*p2y*p3x*xa*xb*yb - p1x*pow(p2y,2)*p3x*xa*xb*yb + 
        2*pow(p1y,2)*pow(p3x,2)*xa*xb*yb - p1y*p2y*pow(p3x,2)*xa*xb*yb + 
        2*p1x*p1y*p2x*p3y*xa*xb*yb - p1y*pow(p2x,2)*p3y*xa*xb*yb - 
        2*pow(p1x,2)*p2y*p3y*xa*xb*yb + p1x*p2x*p2y*p3y*xa*xb*yb - 
        4*p1x*p1y*p3x*p3y*xa*xb*yb + p1y*p2x*p3x*p3y*xa*xb*yb + 
        p1x*p2y*p3x*p3y*xa*xb*yb + 2*pow(p1x,2)*pow(p3y,2)*xa*xb*yb - 
        p1x*p2x*pow(p3y,2)*xa*xb*yb + p1y*p2y*pow(p3x,2)*xa*ya*yb - 
        pow(p2y,2)*pow(p3x,2)*xa*ya*yb - p1y*p2x*p3x*p3y*xa*ya*yb - 
        p1x*p2y*p3x*p3y*xa*ya*yb + 2*p2x*p2y*p3x*p3y*xa*ya*yb + 
        p1x*p2x*pow(p3y,2)*xa*ya*yb - pow(p2x,2)*pow(p3y,2)*xa*ya*yb + 
        p1y*p2y*pow(p3x,2)*xb*ya*yb - pow(p2y,2)*pow(p3x,2)*xb*ya*yb - 
        p1y*p2x*p3x*p3y*xb*ya*yb - p1x*p2y*p3x*p3y*xb*ya*yb + 
        2*p2x*p2y*p3x*p3y*xb*ya*yb + p1x*p2x*pow(p3y,2)*xb*ya*yb - 
        pow(p2x,2)*pow(p3y,2)*xb*ya*yb + 2*pow(p1y,2)*p2x*p3x*xa*xb*ya*yb - 
        2*p1x*p1y*p2y*p3x*xa*xb*ya*yb - 2*p1y*p2x*p2y*p3x*xa*xb*ya*yb + 
        2*p1x*pow(p2y,2)*p3x*xa*xb*ya*yb - 
        2*pow(p1y,2)*pow(p3x,2)*xa*xb*ya*yb + 
        2*p1y*p2y*pow(p3x,2)*xa*xb*ya*yb - 2*p1x*p1y*p2x*p3y*xa*xb*ya*yb + 
        2*p1y*pow(p2x,2)*p3y*xa*xb*ya*yb + 2*pow(p1x,2)*p2y*p3y*xa*xb*ya*yb - 
        2*p1x*p2x*p2y*p3y*xa*xb*ya*yb + 4*p1x*p1y*p3x*p3y*xa*xb*ya*yb - 
        2*p1y*p2x*p3x*p3y*xa*xb*ya*yb - 2*p1x*p2y*p3x*p3y*xa*xb*ya*yb - 
        2*pow(p1x,2)*pow(p3y,2)*xa*xb*ya*yb + 
        2*p1x*p2x*pow(p3y,2)*xa*xb*ya*yb - p1y*p2y*pow(p3x,2)*xa*pow(yb,2) + 
        pow(p2y,2)*pow(p3x,2)*xa*pow(yb,2) + p1y*p2x*p3x*p3y*xa*pow(yb,2) + 
        p1x*p2y*p3x*p3y*xa*pow(yb,2) - 2*p2x*p2y*p3x*p3y*xa*pow(yb,2) - 
        p1x*p2x*pow(p3y,2)*xa*pow(yb,2) + 
        pow(p2x,2)*pow(p3y,2)*xa*pow(yb,2) - 
        pow(p1y,2)*p2x*p3x*pow(xa,2)*pow(yb,2) + 
        p1x*p1y*p2y*p3x*pow(xa,2)*pow(yb,2) + 
        p1y*p2x*p2y*p3x*pow(xa,2)*pow(yb,2) - 
        p1x*pow(p2y,2)*p3x*pow(xa,2)*pow(yb,2) + 
        pow(p1y,2)*pow(p3x,2)*pow(xa,2)*pow(yb,2) - 
        p1y*p2y*pow(p3x,2)*pow(xa,2)*pow(yb,2) + 
        p1x*p1y*p2x*p3y*pow(xa,2)*pow(yb,2) - 
        p1y*pow(p2x,2)*p3y*pow(xa,2)*pow(yb,2) - 
        pow(p1x,2)*p2y*p3y*pow(xa,2)*pow(yb,2) + 
        p1x*p2x*p2y*p3y*pow(xa,2)*pow(yb,2) - 
        2*p1x*p1y*p3x*p3y*pow(xa,2)*pow(yb,2) + 
        p1y*p2x*p3x*p3y*pow(xa,2)*pow(yb,2) + 
        p1x*p2y*p3x*p3y*pow(xa,2)*pow(yb,2) + 
        pow(p1x,2)*pow(p3y,2)*pow(xa,2)*pow(yb,2) - 
        p1x*p2x*pow(p3y,2)*pow(xa,2)*pow(yb,2) + 
        pow(p0y,2)*(p2x - p3x)*(xb + ya - xb*ya + xa*(-1 + yb) - yb)*
         (p3x*xa - p3x*xb + p2x*xb*ya + p1x*(xb - xb*ya + xa*(-1 + yb)) - p2x*xa*yb)\
         + pow(p0x,2)*(p2y - p3y)*(xb + ya - xb*ya + xa*(-1 + yb) - yb)*
         (p3y*xa - p3y*xb + p2y*xb*ya + p1y*(xb - xb*ya + xa*(-1 + yb)) - p2y*xa*yb)\
         - p0y*(-2*p1y*p2x*p3x*pow(xa,2) + 2*p1y*pow(p3x,2)*pow(xa,2) + 
           4*p1y*p2x*p3x*xa*xb - 4*p1y*pow(p3x,2)*xa*xb - 
           2*p1y*p2x*p3x*pow(xb,2) + 2*p1y*pow(p3x,2)*pow(xb,2) + 
           p1y*p2x*p3x*xa*ya - p1y*pow(p3x,2)*xa*ya - p2y*pow(p3x,2)*xa*ya + 
           p2x*p3x*p3y*xa*ya - p1y*p2x*p3x*xb*ya + p1y*pow(p3x,2)*xb*ya + 
           p2y*pow(p3x,2)*xb*ya - p2x*p3x*p3y*xb*ya - p1y*pow(p2x,2)*xa*xb*ya - 
           2*p1y*p2x*p3x*xa*xb*ya + p2x*p2y*p3x*xa*xb*ya + 
           3*p1y*pow(p3x,2)*xa*xb*ya - p2y*pow(p3x,2)*xa*xb*ya - 
           pow(p2x,2)*p3y*xa*xb*ya + p2x*p3x*p3y*xa*xb*ya + 
           p1y*pow(p2x,2)*pow(xb,2)*ya + 2*p1y*p2x*p3x*pow(xb,2)*ya - 
           p2x*p2y*p3x*pow(xb,2)*ya - 3*p1y*pow(p3x,2)*pow(xb,2)*ya + 
           p2y*pow(p3x,2)*pow(xb,2)*ya + pow(p2x,2)*p3y*pow(xb,2)*ya - 
           p2x*p3x*p3y*pow(xb,2)*ya + p1y*p2x*p3x*xb*pow(ya,2) - 
           2*p2x*p2y*p3x*xb*pow(ya,2) - p1y*pow(p3x,2)*xb*pow(ya,2) + 
           p2y*pow(p3x,2)*xb*pow(ya,2) + 2*pow(p2x,2)*p3y*xb*pow(ya,2) - 
           p2x*p3x*p3y*xb*pow(ya,2) - p1y*pow(p2x,2)*pow(xb,2)*pow(ya,2) + 
           p2x*p2y*p3x*pow(xb,2)*pow(ya,2) + 
           p1y*pow(p3x,2)*pow(xb,2)*pow(ya,2) - 
           p2y*pow(p3x,2)*pow(xb,2)*pow(ya,2) - 
           pow(p2x,2)*p3y*pow(xb,2)*pow(ya,2) + 
           p2x*p3x*p3y*pow(xb,2)*pow(ya,2) - 
           pow(p1x,2)*(p2y - p3y)*pow(xb - xb*ya + xa*(-1 + yb),2) - 
           p1y*p2x*p3x*xa*yb + p1y*pow(p3x,2)*xa*yb + p2y*pow(p3x,2)*xa*yb - 
           p2x*p3x*p3y*xa*yb + p1y*pow(p2x,2)*pow(xa,2)*yb + 
           2*p1y*p2x*p3x*pow(xa,2)*yb - p2x*p2y*p3x*pow(xa,2)*yb - 
           3*p1y*pow(p3x,2)*pow(xa,2)*yb + p2y*pow(p3x,2)*pow(xa,2)*yb + 
           pow(p2x,2)*p3y*pow(xa,2)*yb - p2x*p3x*p3y*pow(xa,2)*yb + 
           p1y*p2x*p3x*xb*yb - p1y*pow(p3x,2)*xb*yb - p2y*pow(p3x,2)*xb*yb + 
           p2x*p3x*p3y*xb*yb - p1y*pow(p2x,2)*xa*xb*yb - 2*p1y*p2x*p3x*xa*xb*yb + 
           p2x*p2y*p3x*xa*xb*yb + 3*p1y*pow(p3x,2)*xa*xb*yb - 
           p2y*pow(p3x,2)*xa*xb*yb - pow(p2x,2)*p3y*xa*xb*yb + 
           p2x*p3x*p3y*xa*xb*yb - p1y*p2x*p3x*xa*ya*yb + 2*p2x*p2y*p3x*xa*ya*yb + 
           p1y*pow(p3x,2)*xa*ya*yb - p2y*pow(p3x,2)*xa*ya*yb - 
           2*pow(p2x,2)*p3y*xa*ya*yb + p2x*p3x*p3y*xa*ya*yb - 
           p1y*p2x*p3x*xb*ya*yb + 2*p2x*p2y*p3x*xb*ya*yb + 
           p1y*pow(p3x,2)*xb*ya*yb - p2y*pow(p3x,2)*xb*ya*yb - 
           2*pow(p2x,2)*p3y*xb*ya*yb + p2x*p3x*p3y*xb*ya*yb + 
           2*p1y*pow(p2x,2)*xa*xb*ya*yb - 2*p2x*p2y*p3x*xa*xb*ya*yb - 
           2*p1y*pow(p3x,2)*xa*xb*ya*yb + 2*p2y*pow(p3x,2)*xa*xb*ya*yb + 
           2*pow(p2x,2)*p3y*xa*xb*ya*yb - 2*p2x*p3x*p3y*xa*xb*ya*yb + 
           p1y*p2x*p3x*xa*pow(yb,2) - 2*p2x*p2y*p3x*xa*pow(yb,2) - 
           p1y*pow(p3x,2)*xa*pow(yb,2) + p2y*pow(p3x,2)*xa*pow(yb,2) + 
           2*pow(p2x,2)*p3y*xa*pow(yb,2) - p2x*p3x*p3y*xa*pow(yb,2) - 
           p1y*pow(p2x,2)*pow(xa,2)*pow(yb,2) + 
           p2x*p2y*p3x*pow(xa,2)*pow(yb,2) + 
           p1y*pow(p3x,2)*pow(xa,2)*pow(yb,2) - 
           p2y*pow(p3x,2)*pow(xa,2)*pow(yb,2) - 
           pow(p2x,2)*p3y*pow(xa,2)*pow(yb,2) + 
           p2x*p3x*p3y*pow(xa,2)*pow(yb,2) + 
           p0x*(xb + ya - xb*ya + xa*(-1 + yb) - yb)*
            (p2y*p3x*xa + p2x*p3y*xa - 2*p3x*p3y*xa - p2y*p3x*xb - p2x*p3y*xb + 
              2*p3x*p3y*xb + 2*p2x*p2y*xb*ya - p2y*p3x*xb*ya - p2x*p3y*xb*ya + 
              p1y*(p2x - p3x)*(xb - xb*ya + xa*(-1 + yb)) + 
              p1x*(p2y - p3y)*(xb - xb*ya + xa*(-1 + yb)) - 2*p2x*p2y*xa*yb + 
              p2y*p3x*xa*yb + p2x*p3y*xa*yb) + 
           p1x*(xb - xb*ya + xa*(-1 + yb))*
            (-(p2x*p2y*xb*ya) + p1y*(p2x - p3x)*(xb - xb*ya + xa*(-1 + yb)) + 
              p2x*p3y*(-xa + xb + 2*ya - 2*yb) + p2x*p2y*xa*yb + 
              p2y*p3x*(-xa + xb - ya + yb) + 
              p3x*p3y*(xb*(-2 + ya) - ya - xa*(-2 + yb) + yb))) + 
        p0x*(pow(p1y,2)*(p2x - p3x)*pow(xb - xb*ya + xa*(-1 + yb),2) + 
           p1x*(p2y - p3y)*(xb - xb*ya + xa*(-1 + yb))*
            (-(p2y*xb*ya) + p3y*(-(xb*(-2 + ya)) + ya + xa*(-2 + yb) - yb) + 
              p2y*xa*yb) - p1y*(xb - xb*ya + xa*(-1 + yb))*
            (-(p2x*p2y*xb*ya) + p1x*(p2y - p3y)*(xb - xb*ya + xa*(-1 + yb)) + 
              p2y*p3x*(-xa + xb + 2*ya - 2*yb) + p2x*p2y*xa*yb + 
              p2x*p3y*(-xa + xb - ya + yb) + 
              p3x*p3y*(xb*(-2 + ya) - ya - xa*(-2 + yb) + yb)) + 
           (p2y*p3x - p2x*p3y)*(-(p3y*
                 (pow(xa,2)*(-1 + yb)*yb + 
                   xb*((-1 + xb)*pow(ya,2) + yb + ya*(-1 - xb + yb)) + 
                   xa*((-1 + xb - yb)*yb + ya*(1 + xb + yb - 2*xb*yb)))) + 
              p2y*(pow(xa,2)*(-1 + yb)*yb + xb*ya*(xb*(-1 + ya) - 2*ya + 2*yb) + 
                 xa*(2*(ya - yb)*yb + xb*(ya + yb - 2*ya*yb))))))*
      log(p1y*p3x - p1x*p3y + p1y*p2x*xa - p1x*p2y*xa - p1y*p3x*xa + p1x*p3y*xa + 
        p0x*(p3y - p3y*xa + p2y*(xa - ya) + p1y*(-1 + ya)) - p1y*p3x*ya + 
        p2y*p3x*ya + p1x*p3y*ya - p2x*p3y*ya + 
        p0y*(p1x + p3x*(-1 + xa) - p1x*ya + p2x*(-xa + ya))) - 
     2*(-(pow(p1y,2)*p2x*p3x*pow(xa,2)) + p1x*p1y*p2y*p3x*pow(xa,2) + 
        pow(p1y,2)*pow(p3x,2)*pow(xa,2) + p1x*p1y*p2x*p3y*pow(xa,2) - 
        pow(p1x,2)*p2y*p3y*pow(xa,2) - 2*p1x*p1y*p3x*p3y*pow(xa,2) + 
        pow(p1x,2)*pow(p3y,2)*pow(xa,2) + 2*pow(p1y,2)*p2x*p3x*xa*xb - 
        2*p1x*p1y*p2y*p3x*xa*xb - 2*pow(p1y,2)*pow(p3x,2)*xa*xb - 
        2*p1x*p1y*p2x*p3y*xa*xb + 2*pow(p1x,2)*p2y*p3y*xa*xb + 
        4*p1x*p1y*p3x*p3y*xa*xb - 2*pow(p1x,2)*pow(p3y,2)*xa*xb - 
        pow(p1y,2)*p2x*p3x*pow(xb,2) + p1x*p1y*p2y*p3x*pow(xb,2) + 
        pow(p1y,2)*pow(p3x,2)*pow(xb,2) + p1x*p1y*p2x*p3y*pow(xb,2) - 
        pow(p1x,2)*p2y*p3y*pow(xb,2) - 2*p1x*p1y*p3x*p3y*pow(xb,2) + 
        pow(p1x,2)*pow(p3y,2)*pow(xb,2) - p1y*p2y*pow(p3x,2)*xa*ya + 
        p1y*p2x*p3x*p3y*xa*ya + p1x*p2y*p3x*p3y*xa*ya - p1x*p2x*pow(p3y,2)*xa*ya + 
        p1y*p2y*pow(p3x,2)*xb*ya - p1y*p2x*p3x*p3y*xb*ya - p1x*p2y*p3x*p3y*xb*ya + 
        p1x*p2x*pow(p3y,2)*xb*ya - 2*pow(p1y,2)*p2x*p3x*xa*xb*ya + 
        2*p1x*p1y*p2y*p3x*xa*xb*ya + p1y*p2x*p2y*p3x*xa*xb*ya - 
        p1x*pow(p2y,2)*p3x*xa*xb*ya + 2*pow(p1y,2)*pow(p3x,2)*xa*xb*ya - 
        p1y*p2y*pow(p3x,2)*xa*xb*ya + 2*p1x*p1y*p2x*p3y*xa*xb*ya - 
        p1y*pow(p2x,2)*p3y*xa*xb*ya - 2*pow(p1x,2)*p2y*p3y*xa*xb*ya + 
        p1x*p2x*p2y*p3y*xa*xb*ya - 4*p1x*p1y*p3x*p3y*xa*xb*ya + 
        p1y*p2x*p3x*p3y*xa*xb*ya + p1x*p2y*p3x*p3y*xa*xb*ya + 
        2*pow(p1x,2)*pow(p3y,2)*xa*xb*ya - p1x*p2x*pow(p3y,2)*xa*xb*ya + 
        2*pow(p1y,2)*p2x*p3x*pow(xb,2)*ya - 2*p1x*p1y*p2y*p3x*pow(xb,2)*ya - 
        p1y*p2x*p2y*p3x*pow(xb,2)*ya + p1x*pow(p2y,2)*p3x*pow(xb,2)*ya - 
        2*pow(p1y,2)*pow(p3x,2)*pow(xb,2)*ya + 
        p1y*p2y*pow(p3x,2)*pow(xb,2)*ya - 2*p1x*p1y*p2x*p3y*pow(xb,2)*ya + 
        p1y*pow(p2x,2)*p3y*pow(xb,2)*ya + 
        2*pow(p1x,2)*p2y*p3y*pow(xb,2)*ya - p1x*p2x*p2y*p3y*pow(xb,2)*ya + 
        4*p1x*p1y*p3x*p3y*pow(xb,2)*ya - p1y*p2x*p3x*p3y*pow(xb,2)*ya - 
        p1x*p2y*p3x*p3y*pow(xb,2)*ya - 
        2*pow(p1x,2)*pow(p3y,2)*pow(xb,2)*ya + 
        p1x*p2x*pow(p3y,2)*pow(xb,2)*ya - p1y*p2y*pow(p3x,2)*xb*pow(ya,2) + 
        pow(p2y,2)*pow(p3x,2)*xb*pow(ya,2) + p1y*p2x*p3x*p3y*xb*pow(ya,2) + 
        p1x*p2y*p3x*p3y*xb*pow(ya,2) - 2*p2x*p2y*p3x*p3y*xb*pow(ya,2) - 
        p1x*p2x*pow(p3y,2)*xb*pow(ya,2) + 
        pow(p2x,2)*pow(p3y,2)*xb*pow(ya,2) - 
        pow(p1y,2)*p2x*p3x*pow(xb,2)*pow(ya,2) + 
        p1x*p1y*p2y*p3x*pow(xb,2)*pow(ya,2) + 
        p1y*p2x*p2y*p3x*pow(xb,2)*pow(ya,2) - 
        p1x*pow(p2y,2)*p3x*pow(xb,2)*pow(ya,2) + 
        pow(p1y,2)*pow(p3x,2)*pow(xb,2)*pow(ya,2) - 
        p1y*p2y*pow(p3x,2)*pow(xb,2)*pow(ya,2) + 
        p1x*p1y*p2x*p3y*pow(xb,2)*pow(ya,2) - 
        p1y*pow(p2x,2)*p3y*pow(xb,2)*pow(ya,2) - 
        pow(p1x,2)*p2y*p3y*pow(xb,2)*pow(ya,2) + 
        p1x*p2x*p2y*p3y*pow(xb,2)*pow(ya,2) - 
        2*p1x*p1y*p3x*p3y*pow(xb,2)*pow(ya,2) + 
        p1y*p2x*p3x*p3y*pow(xb,2)*pow(ya,2) + 
        p1x*p2y*p3x*p3y*pow(xb,2)*pow(ya,2) + 
        pow(p1x,2)*pow(p3y,2)*pow(xb,2)*pow(ya,2) - 
        p1x*p2x*pow(p3y,2)*pow(xb,2)*pow(ya,2) + p1y*p2y*pow(p3x,2)*xa*yb - 
        p1y*p2x*p3x*p3y*xa*yb - p1x*p2y*p3x*p3y*xa*yb + p1x*p2x*pow(p3y,2)*xa*yb + 
        2*pow(p1y,2)*p2x*p3x*pow(xa,2)*yb - 2*p1x*p1y*p2y*p3x*pow(xa,2)*yb - 
        p1y*p2x*p2y*p3x*pow(xa,2)*yb + p1x*pow(p2y,2)*p3x*pow(xa,2)*yb - 
        2*pow(p1y,2)*pow(p3x,2)*pow(xa,2)*yb + 
        p1y*p2y*pow(p3x,2)*pow(xa,2)*yb - 2*p1x*p1y*p2x*p3y*pow(xa,2)*yb + 
        p1y*pow(p2x,2)*p3y*pow(xa,2)*yb + 
        2*pow(p1x,2)*p2y*p3y*pow(xa,2)*yb - p1x*p2x*p2y*p3y*pow(xa,2)*yb + 
        4*p1x*p1y*p3x*p3y*pow(xa,2)*yb - p1y*p2x*p3x*p3y*pow(xa,2)*yb - 
        p1x*p2y*p3x*p3y*pow(xa,2)*yb - 
        2*pow(p1x,2)*pow(p3y,2)*pow(xa,2)*yb + 
        p1x*p2x*pow(p3y,2)*pow(xa,2)*yb - p1y*p2y*pow(p3x,2)*xb*yb + 
        p1y*p2x*p3x*p3y*xb*yb + p1x*p2y*p3x*p3y*xb*yb - p1x*p2x*pow(p3y,2)*xb*yb - 
        2*pow(p1y,2)*p2x*p3x*xa*xb*yb + 2*p1x*p1y*p2y*p3x*xa*xb*yb + 
        p1y*p2x*p2y*p3x*xa*xb*yb - p1x*pow(p2y,2)*p3x*xa*xb*yb + 
        2*pow(p1y,2)*pow(p3x,2)*xa*xb*yb - p1y*p2y*pow(p3x,2)*xa*xb*yb + 
        2*p1x*p1y*p2x*p3y*xa*xb*yb - p1y*pow(p2x,2)*p3y*xa*xb*yb - 
        2*pow(p1x,2)*p2y*p3y*xa*xb*yb + p1x*p2x*p2y*p3y*xa*xb*yb - 
        4*p1x*p1y*p3x*p3y*xa*xb*yb + p1y*p2x*p3x*p3y*xa*xb*yb + 
        p1x*p2y*p3x*p3y*xa*xb*yb + 2*pow(p1x,2)*pow(p3y,2)*xa*xb*yb - 
        p1x*p2x*pow(p3y,2)*xa*xb*yb + p1y*p2y*pow(p3x,2)*xa*ya*yb - 
        pow(p2y,2)*pow(p3x,2)*xa*ya*yb - p1y*p2x*p3x*p3y*xa*ya*yb - 
        p1x*p2y*p3x*p3y*xa*ya*yb + 2*p2x*p2y*p3x*p3y*xa*ya*yb + 
        p1x*p2x*pow(p3y,2)*xa*ya*yb - pow(p2x,2)*pow(p3y,2)*xa*ya*yb + 
        p1y*p2y*pow(p3x,2)*xb*ya*yb - pow(p2y,2)*pow(p3x,2)*xb*ya*yb - 
        p1y*p2x*p3x*p3y*xb*ya*yb - p1x*p2y*p3x*p3y*xb*ya*yb + 
        2*p2x*p2y*p3x*p3y*xb*ya*yb + p1x*p2x*pow(p3y,2)*xb*ya*yb - 
        pow(p2x,2)*pow(p3y,2)*xb*ya*yb + 2*pow(p1y,2)*p2x*p3x*xa*xb*ya*yb - 
        2*p1x*p1y*p2y*p3x*xa*xb*ya*yb - 2*p1y*p2x*p2y*p3x*xa*xb*ya*yb + 
        2*p1x*pow(p2y,2)*p3x*xa*xb*ya*yb - 
        2*pow(p1y,2)*pow(p3x,2)*xa*xb*ya*yb + 
        2*p1y*p2y*pow(p3x,2)*xa*xb*ya*yb - 2*p1x*p1y*p2x*p3y*xa*xb*ya*yb + 
        2*p1y*pow(p2x,2)*p3y*xa*xb*ya*yb + 2*pow(p1x,2)*p2y*p3y*xa*xb*ya*yb - 
        2*p1x*p2x*p2y*p3y*xa*xb*ya*yb + 4*p1x*p1y*p3x*p3y*xa*xb*ya*yb - 
        2*p1y*p2x*p3x*p3y*xa*xb*ya*yb - 2*p1x*p2y*p3x*p3y*xa*xb*ya*yb - 
        2*pow(p1x,2)*pow(p3y,2)*xa*xb*ya*yb + 
        2*p1x*p2x*pow(p3y,2)*xa*xb*ya*yb - p1y*p2y*pow(p3x,2)*xa*pow(yb,2) + 
        pow(p2y,2)*pow(p3x,2)*xa*pow(yb,2) + p1y*p2x*p3x*p3y*xa*pow(yb,2) + 
        p1x*p2y*p3x*p3y*xa*pow(yb,2) - 2*p2x*p2y*p3x*p3y*xa*pow(yb,2) - 
        p1x*p2x*pow(p3y,2)*xa*pow(yb,2) + 
        pow(p2x,2)*pow(p3y,2)*xa*pow(yb,2) - 
        pow(p1y,2)*p2x*p3x*pow(xa,2)*pow(yb,2) + 
        p1x*p1y*p2y*p3x*pow(xa,2)*pow(yb,2) + 
        p1y*p2x*p2y*p3x*pow(xa,2)*pow(yb,2) - 
        p1x*pow(p2y,2)*p3x*pow(xa,2)*pow(yb,2) + 
        pow(p1y,2)*pow(p3x,2)*pow(xa,2)*pow(yb,2) - 
        p1y*p2y*pow(p3x,2)*pow(xa,2)*pow(yb,2) + 
        p1x*p1y*p2x*p3y*pow(xa,2)*pow(yb,2) - 
        p1y*pow(p2x,2)*p3y*pow(xa,2)*pow(yb,2) - 
        pow(p1x,2)*p2y*p3y*pow(xa,2)*pow(yb,2) + 
        p1x*p2x*p2y*p3y*pow(xa,2)*pow(yb,2) - 
        2*p1x*p1y*p3x*p3y*pow(xa,2)*pow(yb,2) + 
        p1y*p2x*p3x*p3y*pow(xa,2)*pow(yb,2) + 
        p1x*p2y*p3x*p3y*pow(xa,2)*pow(yb,2) + 
        pow(p1x,2)*pow(p3y,2)*pow(xa,2)*pow(yb,2) - 
        p1x*p2x*pow(p3y,2)*pow(xa,2)*pow(yb,2) + 
        pow(p0y,2)*(p2x - p3x)*(xb + ya - xb*ya + xa*(-1 + yb) - yb)*
         (p3x*xa - p3x*xb + p2x*xb*ya + p1x*(xb - xb*ya + xa*(-1 + yb)) - p2x*xa*yb)\
         + pow(p0x,2)*(p2y - p3y)*(xb + ya - xb*ya + xa*(-1 + yb) - yb)*
         (p3y*xa - p3y*xb + p2y*xb*ya + p1y*(xb - xb*ya + xa*(-1 + yb)) - p2y*xa*yb)\
         - p0y*(-2*p1y*p2x*p3x*pow(xa,2) + 2*p1y*pow(p3x,2)*pow(xa,2) + 
           4*p1y*p2x*p3x*xa*xb - 4*p1y*pow(p3x,2)*xa*xb - 
           2*p1y*p2x*p3x*pow(xb,2) + 2*p1y*pow(p3x,2)*pow(xb,2) + 
           p1y*p2x*p3x*xa*ya - p1y*pow(p3x,2)*xa*ya - p2y*pow(p3x,2)*xa*ya + 
           p2x*p3x*p3y*xa*ya - p1y*p2x*p3x*xb*ya + p1y*pow(p3x,2)*xb*ya + 
           p2y*pow(p3x,2)*xb*ya - p2x*p3x*p3y*xb*ya - p1y*pow(p2x,2)*xa*xb*ya - 
           2*p1y*p2x*p3x*xa*xb*ya + p2x*p2y*p3x*xa*xb*ya + 
           3*p1y*pow(p3x,2)*xa*xb*ya - p2y*pow(p3x,2)*xa*xb*ya - 
           pow(p2x,2)*p3y*xa*xb*ya + p2x*p3x*p3y*xa*xb*ya + 
           p1y*pow(p2x,2)*pow(xb,2)*ya + 2*p1y*p2x*p3x*pow(xb,2)*ya - 
           p2x*p2y*p3x*pow(xb,2)*ya - 3*p1y*pow(p3x,2)*pow(xb,2)*ya + 
           p2y*pow(p3x,2)*pow(xb,2)*ya + pow(p2x,2)*p3y*pow(xb,2)*ya - 
           p2x*p3x*p3y*pow(xb,2)*ya + p1y*p2x*p3x*xb*pow(ya,2) - 
           2*p2x*p2y*p3x*xb*pow(ya,2) - p1y*pow(p3x,2)*xb*pow(ya,2) + 
           p2y*pow(p3x,2)*xb*pow(ya,2) + 2*pow(p2x,2)*p3y*xb*pow(ya,2) - 
           p2x*p3x*p3y*xb*pow(ya,2) - p1y*pow(p2x,2)*pow(xb,2)*pow(ya,2) + 
           p2x*p2y*p3x*pow(xb,2)*pow(ya,2) + 
           p1y*pow(p3x,2)*pow(xb,2)*pow(ya,2) - 
           p2y*pow(p3x,2)*pow(xb,2)*pow(ya,2) - 
           pow(p2x,2)*p3y*pow(xb,2)*pow(ya,2) + 
           p2x*p3x*p3y*pow(xb,2)*pow(ya,2) - 
           pow(p1x,2)*(p2y - p3y)*pow(xb - xb*ya + xa*(-1 + yb),2) - 
           p1y*p2x*p3x*xa*yb + p1y*pow(p3x,2)*xa*yb + p2y*pow(p3x,2)*xa*yb - 
           p2x*p3x*p3y*xa*yb + p1y*pow(p2x,2)*pow(xa,2)*yb + 
           2*p1y*p2x*p3x*pow(xa,2)*yb - p2x*p2y*p3x*pow(xa,2)*yb - 
           3*p1y*pow(p3x,2)*pow(xa,2)*yb + p2y*pow(p3x,2)*pow(xa,2)*yb + 
           pow(p2x,2)*p3y*pow(xa,2)*yb - p2x*p3x*p3y*pow(xa,2)*yb + 
           p1y*p2x*p3x*xb*yb - p1y*pow(p3x,2)*xb*yb - p2y*pow(p3x,2)*xb*yb + 
           p2x*p3x*p3y*xb*yb - p1y*pow(p2x,2)*xa*xb*yb - 2*p1y*p2x*p3x*xa*xb*yb + 
           p2x*p2y*p3x*xa*xb*yb + 3*p1y*pow(p3x,2)*xa*xb*yb - 
           p2y*pow(p3x,2)*xa*xb*yb - pow(p2x,2)*p3y*xa*xb*yb + 
           p2x*p3x*p3y*xa*xb*yb - p1y*p2x*p3x*xa*ya*yb + 2*p2x*p2y*p3x*xa*ya*yb + 
           p1y*pow(p3x,2)*xa*ya*yb - p2y*pow(p3x,2)*xa*ya*yb - 
           2*pow(p2x,2)*p3y*xa*ya*yb + p2x*p3x*p3y*xa*ya*yb - 
           p1y*p2x*p3x*xb*ya*yb + 2*p2x*p2y*p3x*xb*ya*yb + 
           p1y*pow(p3x,2)*xb*ya*yb - p2y*pow(p3x,2)*xb*ya*yb - 
           2*pow(p2x,2)*p3y*xb*ya*yb + p2x*p3x*p3y*xb*ya*yb + 
           2*p1y*pow(p2x,2)*xa*xb*ya*yb - 2*p2x*p2y*p3x*xa*xb*ya*yb - 
           2*p1y*pow(p3x,2)*xa*xb*ya*yb + 2*p2y*pow(p3x,2)*xa*xb*ya*yb + 
           2*pow(p2x,2)*p3y*xa*xb*ya*yb - 2*p2x*p3x*p3y*xa*xb*ya*yb + 
           p1y*p2x*p3x*xa*pow(yb,2) - 2*p2x*p2y*p3x*xa*pow(yb,2) - 
           p1y*pow(p3x,2)*xa*pow(yb,2) + p2y*pow(p3x,2)*xa*pow(yb,2) + 
           2*pow(p2x,2)*p3y*xa*pow(yb,2) - p2x*p3x*p3y*xa*pow(yb,2) - 
           p1y*pow(p2x,2)*pow(xa,2)*pow(yb,2) + 
           p2x*p2y*p3x*pow(xa,2)*pow(yb,2) + 
           p1y*pow(p3x,2)*pow(xa,2)*pow(yb,2) - 
           p2y*pow(p3x,2)*pow(xa,2)*pow(yb,2) - 
           pow(p2x,2)*p3y*pow(xa,2)*pow(yb,2) + 
           p2x*p3x*p3y*pow(xa,2)*pow(yb,2) + 
           p0x*(xb + ya - xb*ya + xa*(-1 + yb) - yb)*
            (p2y*p3x*xa + p2x*p3y*xa - 2*p3x*p3y*xa - p2y*p3x*xb - p2x*p3y*xb + 
              2*p3x*p3y*xb + 2*p2x*p2y*xb*ya - p2y*p3x*xb*ya - p2x*p3y*xb*ya + 
              p1y*(p2x - p3x)*(xb - xb*ya + xa*(-1 + yb)) + 
              p1x*(p2y - p3y)*(xb - xb*ya + xa*(-1 + yb)) - 2*p2x*p2y*xa*yb + 
              p2y*p3x*xa*yb + p2x*p3y*xa*yb) + 
           p1x*(xb - xb*ya + xa*(-1 + yb))*
            (-(p2x*p2y*xb*ya) + p1y*(p2x - p3x)*(xb - xb*ya + xa*(-1 + yb)) + 
              p2x*p3y*(-xa + xb + 2*ya - 2*yb) + p2x*p2y*xa*yb + 
              p2y*p3x*(-xa + xb - ya + yb) + 
              p3x*p3y*(xb*(-2 + ya) - ya - xa*(-2 + yb) + yb))) + 
        p0x*(pow(p1y,2)*(p2x - p3x)*pow(xb - xb*ya + xa*(-1 + yb),2) + 
           p1x*(p2y - p3y)*(xb - xb*ya + xa*(-1 + yb))*
            (-(p2y*xb*ya) + p3y*(-(xb*(-2 + ya)) + ya + xa*(-2 + yb) - yb) + 
              p2y*xa*yb) - p1y*(xb - xb*ya + xa*(-1 + yb))*
            (-(p2x*p2y*xb*ya) + p1x*(p2y - p3y)*(xb - xb*ya + xa*(-1 + yb)) + 
              p2y*p3x*(-xa + xb + 2*ya - 2*yb) + p2x*p2y*xa*yb + 
              p2x*p3y*(-xa + xb - ya + yb) + 
              p3x*p3y*(xb*(-2 + ya) - ya - xa*(-2 + yb) + yb)) + 
           (p2y*p3x - p2x*p3y)*(-(p3y*
                 (pow(xa,2)*(-1 + yb)*yb + 
                   xb*((-1 + xb)*pow(ya,2) + yb + ya*(-1 - xb + yb)) + 
                   xa*((-1 + xb - yb)*yb + ya*(1 + xb + yb - 2*xb*yb)))) + 
              p2y*(pow(xa,2)*(-1 + yb)*yb + xb*ya*(xb*(-1 + ya) - 2*ya + 2*yb) + 
                 xa*(2*(ya - yb)*yb + xb*(ya + yb - 2*ya*yb))))))*
      log(p1y*p3x - p1x*p3y + p1y*p2x*xb - p1x*p2y*xb - p1y*p3x*xb + p1x*p3y*xb + 
        p0x*(p3y - p3y*xb + p2y*(xb - yb) + p1y*(-1 + yb)) - p1y*p3x*yb + 
        p2y*p3x*yb + p1x*p3y*yb - p2x*p3y*yb + 
        p0y*(p1x + p3x*(-1 + xb) - p1x*yb + p2x*(-xb + yb))))/
   (2.*pow(-(p0x*p2y*xa) + p1x*p2y*xa + p0x*p3y*xa - p1x*p3y*xa + p0x*p2y*xb - 
       p1x*p2y*xb - p0x*p3y*xb + p1x*p3y*xb + p0x*p2y*ya - p2y*p3x*ya - p1x*p3y*ya + 
       p2x*p3y*ya - p0x*p2y*yb + p2y*p3x*yb + p1x*p3y*yb - p2x*p3y*yb + 
       p1y*(p2x*(-xa + xb) + p3x*(xa - xb + ya - yb) + p0x*(-ya + yb)) + 
       p0y*(p3x*(-xa + xb) + p1x*(ya - yb) + p2x*(xa - xb - ya + yb)),3));

    integrals[2] = (-((xa - xb)*(ya - yb)*pow(-(p0x*p2y*xa) + p1x*p2y*xa + p0x*p3y*xa - p1x*p3y*xa + 
          p0x*p2y*xb - p1x*p2y*xb - p0x*p3y*xb + p1x*p3y*xb + p0x*p2y*ya - 
          p2y*p3x*ya - p1x*p3y*ya + p2x*p3y*ya - p0x*p2y*yb + p2y*p3x*yb + 
          p1x*p3y*yb - p2x*p3y*yb + 
          p1y*(p2x*(-xa + xb) + p3x*(xa - xb + ya - yb) + p0x*(-ya + yb)) + 
          p0y*(p3x*(-xa + xb) + p1x*(ya - yb) + p2x*(xa - xb - ya + yb)),2)) + 
     2*(-(p0x*p2y*xa) + p1x*p2y*xa + p0x*p3y*xa - p1x*p3y*xa + p0x*p2y*xb - 
        p1x*p2y*xb - p0x*p3y*xb + p1x*p3y*xb + p0x*p2y*ya - p2y*p3x*ya - 
        p1x*p3y*ya + p2x*p3y*ya - p0x*p2y*yb + p2y*p3x*yb + p1x*p3y*yb - 
        p2x*p3y*yb + p1y*(p2x*(-xa + xb) + p3x*(xa - xb + ya - yb) + 
           p0x*(-ya + yb)) + p0y*(p3x*(-xa + xb) + p1x*(ya - yb) + 
           p2x*(xa - xb - ya + yb)))*
      (p1y*p3x*xa*ya - p1x*p3y*xa*ya - p1y*p2x*pow(xa,2)*ya + 
        p1x*p2y*pow(xa,2)*ya + p1y*p3x*pow(xa,2)*ya - p1x*p3y*pow(xa,2)*ya - 
        p1y*p3x*xb*ya + p1x*p3y*xb*ya + 2*p1y*p2x*xa*xb*ya - 2*p1x*p2y*xa*xb*ya - 
        2*p1y*p3x*xa*xb*ya + 2*p1x*p3y*xa*xb*ya - p1y*p2x*pow(xb,2)*ya + 
        p1x*p2y*pow(xb,2)*ya + p1y*p3x*pow(xb,2)*ya - p1x*p3y*pow(xb,2)*ya + 
        p1y*p3x*xa*pow(ya,2) - p2y*p3x*xa*pow(ya,2) - p1x*p3y*xa*pow(ya,2) + 
        p2x*p3y*xa*pow(ya,2) + p0y*
         (p2x*(pow(xa,2)*ya + pow(xb,2)*ya - xa*(2*xb*ya + pow(ya - yb,2))) + 
           p1x*(-xb + xa*(1 + ya - yb))*(ya - yb) - 
           p3x*(xa - xb)*((1 + xa - xb)*ya - yb)) - p1y*p3x*xa*yb + p1x*p3y*xa*yb + 
        p1y*p3x*xb*yb - p1x*p3y*xb*yb - 2*p1y*p3x*xa*ya*yb + 2*p2y*p3x*xa*ya*yb + 
        2*p1x*p3y*xa*ya*yb - 2*p2x*p3y*xa*ya*yb + p1y*p3x*xa*pow(yb,2) - 
        p2y*p3x*xa*pow(yb,2) - p1x*p3y*xa*pow(yb,2) + p2x*p3y*xa*pow(yb,2) + 
        p0x*(p2y*(-(pow(xa,2)*ya) - pow(xb,2)*ya + 
              xa*(2*xb*ya + pow(ya - yb,2))) + 
           p3y*(xa - xb)*((1 + xa - xb)*ya - yb) + 
           p1y*(ya - yb)*(xb + xa*(-1 - ya + yb)))) + 
     2*(pow(p1y,2)*pow(p3x,2)*xa*ya - 2*p1x*p1y*p3x*p3y*xa*ya + 
        pow(p1x,2)*pow(p3y,2)*xa*ya - pow(p1y,2)*pow(p3x,2)*xb*ya + 
        2*p1x*p1y*p3x*p3y*xb*ya - pow(p1x,2)*pow(p3y,2)*xb*ya + 
        pow(p1y,2)*p2x*p3x*xa*xb*ya - p1x*p1y*p2y*p3x*xa*xb*ya - 
        pow(p1y,2)*pow(p3x,2)*xa*xb*ya - p1x*p1y*p2x*p3y*xa*xb*ya + 
        pow(p1x,2)*p2y*p3y*xa*xb*ya + 2*p1x*p1y*p3x*p3y*xa*xb*ya - 
        pow(p1x,2)*pow(p3y,2)*xa*xb*ya - pow(p1y,2)*p2x*p3x*pow(xb,2)*ya + 
        p1x*p1y*p2y*p3x*pow(xb,2)*ya + pow(p1y,2)*pow(p3x,2)*pow(xb,2)*ya + 
        p1x*p1y*p2x*p3y*pow(xb,2)*ya - pow(p1x,2)*p2y*p3y*pow(xb,2)*ya - 
        2*p1x*p1y*p3x*p3y*pow(xb,2)*ya + 
        pow(p1x,2)*pow(p3y,2)*pow(xb,2)*ya + 
        pow(p1y,2)*pow(p3x,2)*xb*pow(ya,2) - 
        p1y*p2y*pow(p3x,2)*xb*pow(ya,2) - 2*p1x*p1y*p3x*p3y*xb*pow(ya,2) + 
        p1y*p2x*p3x*p3y*xb*pow(ya,2) + p1x*p2y*p3x*p3y*xb*pow(ya,2) + 
        pow(p1x,2)*pow(p3y,2)*xb*pow(ya,2) - 
        p1x*p2x*pow(p3y,2)*xb*pow(ya,2) + 
        pow(p1y,2)*p2x*p3x*pow(xb,2)*pow(ya,2) - 
        p1x*p1y*p2y*p3x*pow(xb,2)*pow(ya,2) - 
        p1y*p2x*p2y*p3x*pow(xb,2)*pow(ya,2) + 
        p1x*pow(p2y,2)*p3x*pow(xb,2)*pow(ya,2) - 
        pow(p1y,2)*pow(p3x,2)*pow(xb,2)*pow(ya,2) + 
        p1y*p2y*pow(p3x,2)*pow(xb,2)*pow(ya,2) - 
        p1x*p1y*p2x*p3y*pow(xb,2)*pow(ya,2) + 
        p1y*pow(p2x,2)*p3y*pow(xb,2)*pow(ya,2) + 
        pow(p1x,2)*p2y*p3y*pow(xb,2)*pow(ya,2) - 
        p1x*p2x*p2y*p3y*pow(xb,2)*pow(ya,2) + 
        2*p1x*p1y*p3x*p3y*pow(xb,2)*pow(ya,2) - 
        p1y*p2x*p3x*p3y*pow(xb,2)*pow(ya,2) - 
        p1x*p2y*p3x*p3y*pow(xb,2)*pow(ya,2) - 
        pow(p1x,2)*pow(p3y,2)*pow(xb,2)*pow(ya,2) + 
        p1x*p2x*pow(p3y,2)*pow(xb,2)*pow(ya,2) - 
        pow(p1y,2)*pow(p3x,2)*xa*yb + 2*p1x*p1y*p3x*p3y*xa*yb - 
        pow(p1x,2)*pow(p3y,2)*xa*yb - pow(p1y,2)*p2x*p3x*pow(xa,2)*yb + 
        p1x*p1y*p2y*p3x*pow(xa,2)*yb + pow(p1y,2)*pow(p3x,2)*pow(xa,2)*yb + 
        p1x*p1y*p2x*p3y*pow(xa,2)*yb - pow(p1x,2)*p2y*p3y*pow(xa,2)*yb - 
        2*p1x*p1y*p3x*p3y*pow(xa,2)*yb + 
        pow(p1x,2)*pow(p3y,2)*pow(xa,2)*yb + pow(p1y,2)*pow(p3x,2)*xb*yb - 
        2*p1x*p1y*p3x*p3y*xb*yb + pow(p1x,2)*pow(p3y,2)*xb*yb + 
        pow(p1y,2)*p2x*p3x*xa*xb*yb - p1x*p1y*p2y*p3x*xa*xb*yb - 
        pow(p1y,2)*pow(p3x,2)*xa*xb*yb - p1x*p1y*p2x*p3y*xa*xb*yb + 
        pow(p1x,2)*p2y*p3y*xa*xb*yb + 2*p1x*p1y*p3x*p3y*xa*xb*yb - 
        pow(p1x,2)*pow(p3y,2)*xa*xb*yb - pow(p1y,2)*pow(p3x,2)*xa*ya*yb + 
        p1y*p2y*pow(p3x,2)*xa*ya*yb + 2*p1x*p1y*p3x*p3y*xa*ya*yb - 
        p1y*p2x*p3x*p3y*xa*ya*yb - p1x*p2y*p3x*p3y*xa*ya*yb - 
        pow(p1x,2)*pow(p3y,2)*xa*ya*yb + p1x*p2x*pow(p3y,2)*xa*ya*yb - 
        pow(p1y,2)*pow(p3x,2)*xb*ya*yb + p1y*p2y*pow(p3x,2)*xb*ya*yb + 
        2*p1x*p1y*p3x*p3y*xb*ya*yb - p1y*p2x*p3x*p3y*xb*ya*yb - 
        p1x*p2y*p3x*p3y*xb*ya*yb - pow(p1x,2)*pow(p3y,2)*xb*ya*yb + 
        p1x*p2x*pow(p3y,2)*xb*ya*yb - 2*pow(p1y,2)*p2x*p3x*xa*xb*ya*yb + 
        2*p1x*p1y*p2y*p3x*xa*xb*ya*yb + 2*p1y*p2x*p2y*p3x*xa*xb*ya*yb - 
        2*p1x*pow(p2y,2)*p3x*xa*xb*ya*yb + 
        2*pow(p1y,2)*pow(p3x,2)*xa*xb*ya*yb - 
        2*p1y*p2y*pow(p3x,2)*xa*xb*ya*yb + 2*p1x*p1y*p2x*p3y*xa*xb*ya*yb - 
        2*p1y*pow(p2x,2)*p3y*xa*xb*ya*yb - 2*pow(p1x,2)*p2y*p3y*xa*xb*ya*yb + 
        2*p1x*p2x*p2y*p3y*xa*xb*ya*yb - 4*p1x*p1y*p3x*p3y*xa*xb*ya*yb + 
        2*p1y*p2x*p3x*p3y*xa*xb*ya*yb + 2*p1x*p2y*p3x*p3y*xa*xb*ya*yb + 
        2*pow(p1x,2)*pow(p3y,2)*xa*xb*ya*yb - 
        2*p1x*p2x*pow(p3y,2)*xa*xb*ya*yb + 
        pow(p1y,2)*pow(p3x,2)*xa*pow(yb,2) - 
        p1y*p2y*pow(p3x,2)*xa*pow(yb,2) - 2*p1x*p1y*p3x*p3y*xa*pow(yb,2) + 
        p1y*p2x*p3x*p3y*xa*pow(yb,2) + p1x*p2y*p3x*p3y*xa*pow(yb,2) + 
        pow(p1x,2)*pow(p3y,2)*xa*pow(yb,2) - 
        p1x*p2x*pow(p3y,2)*xa*pow(yb,2) + 
        pow(p1y,2)*p2x*p3x*pow(xa,2)*pow(yb,2) - 
        p1x*p1y*p2y*p3x*pow(xa,2)*pow(yb,2) - 
        p1y*p2x*p2y*p3x*pow(xa,2)*pow(yb,2) + 
        p1x*pow(p2y,2)*p3x*pow(xa,2)*pow(yb,2) - 
        pow(p1y,2)*pow(p3x,2)*pow(xa,2)*pow(yb,2) + 
        p1y*p2y*pow(p3x,2)*pow(xa,2)*pow(yb,2) - 
        p1x*p1y*p2x*p3y*pow(xa,2)*pow(yb,2) + 
        p1y*pow(p2x,2)*p3y*pow(xa,2)*pow(yb,2) + 
        pow(p1x,2)*p2y*p3y*pow(xa,2)*pow(yb,2) - 
        p1x*p2x*p2y*p3y*pow(xa,2)*pow(yb,2) + 
        2*p1x*p1y*p3x*p3y*pow(xa,2)*pow(yb,2) - 
        p1y*p2x*p3x*p3y*pow(xa,2)*pow(yb,2) - 
        p1x*p2y*p3x*p3y*pow(xa,2)*pow(yb,2) - 
        pow(p1x,2)*pow(p3y,2)*pow(xa,2)*pow(yb,2) + 
        p1x*p2x*pow(p3y,2)*pow(xa,2)*pow(yb,2) + 
        pow(p0y,2)*(pow(p3x,2)*(xa - xb)*(ya - xb*ya + (-1 + xa)*yb) + 
           pow(p1x,2)*(ya - yb)*(xa + xb*(-1 + ya) - xa*yb) + 
           pow(p2x,2)*pow(xb*ya - xa*yb,2) - 
           p2x*p3x*(pow(xa,2)*yb*(1 + yb) + xb*ya*(xb - ya + xb*ya + yb) - 
              xa*(yb*(-ya + yb) + xb*(ya + yb + 2*ya*yb))) + 
           p1x*(p3x*((-1 + xb)*xb*pow(ya,2) - 2*xb*yb + xa*(2 + xb - yb)*yb + 
                 pow(xa,2)*(-1 + yb)*yb + xb*ya*(2 - xb + yb) + 
                 xa*ya*(-2 + xb + yb - 2*xb*yb)) - 
              p2x*(xb*ya*(xb*(-1 + ya) + ya - yb) + pow(xa,2)*(-1 + yb)*yb + 
                 xa*(yb*(-ya + yb) + xb*(ya + yb - 2*ya*yb))))) + 
        pow(p0x,2)*(pow(p3y,2)*(xa - xb)*(ya - xb*ya + (-1 + xa)*yb) + 
           pow(p1y,2)*(ya - yb)*(xa + xb*(-1 + ya) - xa*yb) + 
           pow(p2y,2)*pow(xb*ya - xa*yb,2) - 
           p2y*p3y*(pow(xa,2)*yb*(1 + yb) + xb*ya*(xb - ya + xb*ya + yb) - 
              xa*(yb*(-ya + yb) + xb*(ya + yb + 2*ya*yb))) + 
           p1y*(p3y*((-1 + xb)*xb*pow(ya,2) - 2*xb*yb + xa*(2 + xb - yb)*yb + 
                 pow(xa,2)*(-1 + yb)*yb + xb*ya*(2 - xb + yb) + 
                 xa*ya*(-2 + xb + yb - 2*xb*yb)) - 
              p2y*(xb*ya*(xb*(-1 + ya) + ya - yb) + pow(xa,2)*(-1 + yb)*yb + 
                 xa*(yb*(-ya + yb) + xb*(ya + yb - 2*ya*yb))))) + 
        p0x*(pow(p1y,2)*(xb - xb*ya + xa*(-1 + yb))*
            (-(p3x*(-2 + xb)*ya) + p3x*(-2 + xa)*yb + p2x*(xb*ya - xa*yb)) + 
           (-(p2y*xb*ya) + p2y*xa*yb + p3y*((-1 + xb)*ya + yb - xa*yb))*
            ((p2y*p3x - p2x*p3y)*(xb*ya - xa*yb) + 
              p1x*(p3y*xb*(-2 + ya) - p3y*xa*(-2 + yb) + p2y*(xb*ya - xa*yb))) + 
           p1y*(p1x*(xa + xb*(-1 + ya) - xa*yb)*
               (-(p3y*(-2 + xb)*ya) + p3y*(-2 + xa)*yb + p2y*(xb*ya - xa*yb)) + 
              p2x*(xb*ya - xa*yb)*(p3y*(xa - xb - ya + yb) + p2y*(xb*ya - xa*yb)) - 
              p3x*(p2y*(pow(xa,2)*yb + xb*ya*(xb - 2*ya + 2*yb) - 
                    xa*(2*yb*(-ya + yb) + xb*(ya + yb))) + 
                 p3y*(pow(xa,2)*(-2 + yb)*yb + xb*(-2 + ya)*((-1 + xb)*ya + yb) + 
                    xa*((2 + 2*xb - yb)*yb + ya*(-2 + 2*xb + yb - 2*xb*yb)))))) + 
        p0y*(pow(p1x,2)*(xb - xb*ya + xa*(-1 + yb))*
            (-(p3y*(-2 + xb)*ya) + p3y*(-2 + xa)*yb + p2y*(xb*ya - xa*yb)) + 
           (-(p2x*xb*ya) + p2x*xa*yb + p3x*((-1 + xb)*ya + yb - xa*yb))*
            (-((p2y*p3x - p2x*p3y)*(xb*ya - xa*yb)) + 
              p1y*(p3x*xb*(-2 + ya) - p3x*xa*(-2 + yb) + p2x*(xb*ya - xa*yb))) + 
           p1x*(p1y*(xa + xb*(-1 + ya) - xa*yb)*
               (-(p3x*(-2 + xb)*ya) + p3x*(-2 + xa)*yb + p2x*(xb*ya - xa*yb)) + 
              p2x*(xb*ya - xa*yb)*
               (p3y*(xa - xb + 2*ya - 2*yb) + p2y*(xb*ya - xa*yb)) - 
              p3x*(p2y*(xb*ya*(xb + ya - yb) + pow(xa,2)*yb - 
                    xa*((ya - yb)*yb + xb*(ya + yb))) + 
                 p3y*(pow(xa,2)*(-2 + yb)*yb + xb*(-2 + ya)*((-1 + xb)*ya + yb) + 
                    xa*((2 + 2*xb - yb)*yb + ya*(-2 + 2*xb + yb - 2*xb*yb))))) + 
           p0x*(-2*p3x*p3y*xa*ya + 2*p3x*p3y*xb*ya - p2y*p3x*xa*xb*ya - 
              p2x*p3y*xa*xb*ya + 2*p3x*p3y*xa*xb*ya + p2y*p3x*pow(xb,2)*ya + 
              p2x*p3y*pow(xb,2)*ya - 2*p3x*p3y*pow(xb,2)*ya - 
              p2y*p3x*xb*pow(ya,2) - p2x*p3y*xb*pow(ya,2) - 
              2*p2x*p2y*pow(xb,2)*pow(ya,2) + p2y*p3x*pow(xb,2)*pow(ya,2) + 
              p2x*p3y*pow(xb,2)*pow(ya,2) + 2*p3x*p3y*xa*yb + 
              p2y*p3x*pow(xa,2)*yb + p2x*p3y*pow(xa,2)*yb - 
              2*p3x*p3y*pow(xa,2)*yb - 2*p3x*p3y*xb*yb - p2y*p3x*xa*xb*yb - 
              p2x*p3y*xa*xb*yb + 2*p3x*p3y*xa*xb*yb + p2y*p3x*xa*ya*yb + 
              p2x*p3y*xa*ya*yb + p2y*p3x*xb*ya*yb + p2x*p3y*xb*ya*yb + 
              4*p2x*p2y*xa*xb*ya*yb - 2*p2y*p3x*xa*xb*ya*yb - 
              2*p2x*p3y*xa*xb*ya*yb - p2y*p3x*xa*pow(yb,2) - 
              p2x*p3y*xa*pow(yb,2) - 2*p2x*p2y*pow(xa,2)*pow(yb,2) + 
              p2y*p3x*pow(xa,2)*pow(yb,2) + p2x*p3y*pow(xa,2)*pow(yb,2) - 
              p1y*p3x*((-1 + xb)*xb*pow(ya,2) - 2*xb*yb + xa*(2 + xb - yb)*yb + 
                 pow(xa,2)*(-1 + yb)*yb + xb*ya*(2 - xb + yb) + 
                 xa*ya*(-2 + xb + yb - 2*xb*yb)) + 
              p1y*p2x*(xb*ya*(xb*(-1 + ya) + ya - yb) + pow(xa,2)*(-1 + yb)*yb + 
                 xa*(yb*(-ya + yb) + xb*(ya + yb - 2*ya*yb))) + 
              p1x*(2*p1y*(xb - xb*ya + xa*(-1 + yb))*(ya - yb) - 
                 p3y*((-1 + xb)*xb*pow(ya,2) - 2*xb*yb + xa*(2 + xb - yb)*yb + 
                    pow(xa,2)*(-1 + yb)*yb + xb*ya*(2 - xb + yb) + 
                    xa*ya*(-2 + xb + yb - 2*xb*yb)) + 
                 p2y*(xb*ya*(xb*(-1 + ya) + ya - yb) + pow(xa,2)*(-1 + yb)*yb + 
                    xa*(yb*(-ya + yb) + xb*(ya + yb - 2*ya*yb)))))))*
      log(p1y*p3x - p1x*p3y + p1y*p2x*xa - p1x*p2y*xa - p1y*p3x*xa + p1x*p3y*xa + 
        p0x*(p3y - p3y*xa + p2y*(xa - ya) + p1y*(-1 + ya)) - p1y*p3x*ya + 
        p2y*p3x*ya + p1x*p3y*ya - p2x*p3y*ya + 
        p0y*(p1x + p3x*(-1 + xa) - p1x*ya + p2x*(-xa + ya))) + 
     2*(-(pow(p1y,2)*pow(p3x,2)*xa*ya) + 2*p1x*p1y*p3x*p3y*xa*ya - 
        pow(p1x,2)*pow(p3y,2)*xa*ya + pow(p1y,2)*pow(p3x,2)*xb*ya - 
        2*p1x*p1y*p3x*p3y*xb*ya + pow(p1x,2)*pow(p3y,2)*xb*ya - 
        pow(p1y,2)*p2x*p3x*xa*xb*ya + p1x*p1y*p2y*p3x*xa*xb*ya + 
        pow(p1y,2)*pow(p3x,2)*xa*xb*ya + p1x*p1y*p2x*p3y*xa*xb*ya - 
        pow(p1x,2)*p2y*p3y*xa*xb*ya - 2*p1x*p1y*p3x*p3y*xa*xb*ya + 
        pow(p1x,2)*pow(p3y,2)*xa*xb*ya + pow(p1y,2)*p2x*p3x*pow(xb,2)*ya - 
        p1x*p1y*p2y*p3x*pow(xb,2)*ya - pow(p1y,2)*pow(p3x,2)*pow(xb,2)*ya - 
        p1x*p1y*p2x*p3y*pow(xb,2)*ya + pow(p1x,2)*p2y*p3y*pow(xb,2)*ya + 
        2*p1x*p1y*p3x*p3y*pow(xb,2)*ya - 
        pow(p1x,2)*pow(p3y,2)*pow(xb,2)*ya - 
        pow(p1y,2)*pow(p3x,2)*xb*pow(ya,2) + 
        p1y*p2y*pow(p3x,2)*xb*pow(ya,2) + 2*p1x*p1y*p3x*p3y*xb*pow(ya,2) - 
        p1y*p2x*p3x*p3y*xb*pow(ya,2) - p1x*p2y*p3x*p3y*xb*pow(ya,2) - 
        pow(p1x,2)*pow(p3y,2)*xb*pow(ya,2) + 
        p1x*p2x*pow(p3y,2)*xb*pow(ya,2) - 
        pow(p1y,2)*p2x*p3x*pow(xb,2)*pow(ya,2) + 
        p1x*p1y*p2y*p3x*pow(xb,2)*pow(ya,2) + 
        p1y*p2x*p2y*p3x*pow(xb,2)*pow(ya,2) - 
        p1x*pow(p2y,2)*p3x*pow(xb,2)*pow(ya,2) + 
        pow(p1y,2)*pow(p3x,2)*pow(xb,2)*pow(ya,2) - 
        p1y*p2y*pow(p3x,2)*pow(xb,2)*pow(ya,2) + 
        p1x*p1y*p2x*p3y*pow(xb,2)*pow(ya,2) - 
        p1y*pow(p2x,2)*p3y*pow(xb,2)*pow(ya,2) - 
        pow(p1x,2)*p2y*p3y*pow(xb,2)*pow(ya,2) + 
        p1x*p2x*p2y*p3y*pow(xb,2)*pow(ya,2) - 
        2*p1x*p1y*p3x*p3y*pow(xb,2)*pow(ya,2) + 
        p1y*p2x*p3x*p3y*pow(xb,2)*pow(ya,2) + 
        p1x*p2y*p3x*p3y*pow(xb,2)*pow(ya,2) + 
        pow(p1x,2)*pow(p3y,2)*pow(xb,2)*pow(ya,2) - 
        p1x*p2x*pow(p3y,2)*pow(xb,2)*pow(ya,2) + 
        pow(p1y,2)*pow(p3x,2)*xa*yb - 2*p1x*p1y*p3x*p3y*xa*yb + 
        pow(p1x,2)*pow(p3y,2)*xa*yb + pow(p1y,2)*p2x*p3x*pow(xa,2)*yb - 
        p1x*p1y*p2y*p3x*pow(xa,2)*yb - pow(p1y,2)*pow(p3x,2)*pow(xa,2)*yb - 
        p1x*p1y*p2x*p3y*pow(xa,2)*yb + pow(p1x,2)*p2y*p3y*pow(xa,2)*yb + 
        2*p1x*p1y*p3x*p3y*pow(xa,2)*yb - 
        pow(p1x,2)*pow(p3y,2)*pow(xa,2)*yb - pow(p1y,2)*pow(p3x,2)*xb*yb + 
        2*p1x*p1y*p3x*p3y*xb*yb - pow(p1x,2)*pow(p3y,2)*xb*yb - 
        pow(p1y,2)*p2x*p3x*xa*xb*yb + p1x*p1y*p2y*p3x*xa*xb*yb + 
        pow(p1y,2)*pow(p3x,2)*xa*xb*yb + p1x*p1y*p2x*p3y*xa*xb*yb - 
        pow(p1x,2)*p2y*p3y*xa*xb*yb - 2*p1x*p1y*p3x*p3y*xa*xb*yb + 
        pow(p1x,2)*pow(p3y,2)*xa*xb*yb + pow(p1y,2)*pow(p3x,2)*xa*ya*yb - 
        p1y*p2y*pow(p3x,2)*xa*ya*yb - 2*p1x*p1y*p3x*p3y*xa*ya*yb + 
        p1y*p2x*p3x*p3y*xa*ya*yb + p1x*p2y*p3x*p3y*xa*ya*yb + 
        pow(p1x,2)*pow(p3y,2)*xa*ya*yb - p1x*p2x*pow(p3y,2)*xa*ya*yb + 
        pow(p1y,2)*pow(p3x,2)*xb*ya*yb - p1y*p2y*pow(p3x,2)*xb*ya*yb - 
        2*p1x*p1y*p3x*p3y*xb*ya*yb + p1y*p2x*p3x*p3y*xb*ya*yb + 
        p1x*p2y*p3x*p3y*xb*ya*yb + pow(p1x,2)*pow(p3y,2)*xb*ya*yb - 
        p1x*p2x*pow(p3y,2)*xb*ya*yb + 2*pow(p1y,2)*p2x*p3x*xa*xb*ya*yb - 
        2*p1x*p1y*p2y*p3x*xa*xb*ya*yb - 2*p1y*p2x*p2y*p3x*xa*xb*ya*yb + 
        2*p1x*pow(p2y,2)*p3x*xa*xb*ya*yb - 
        2*pow(p1y,2)*pow(p3x,2)*xa*xb*ya*yb + 
        2*p1y*p2y*pow(p3x,2)*xa*xb*ya*yb - 2*p1x*p1y*p2x*p3y*xa*xb*ya*yb + 
        2*p1y*pow(p2x,2)*p3y*xa*xb*ya*yb + 2*pow(p1x,2)*p2y*p3y*xa*xb*ya*yb - 
        2*p1x*p2x*p2y*p3y*xa*xb*ya*yb + 4*p1x*p1y*p3x*p3y*xa*xb*ya*yb - 
        2*p1y*p2x*p3x*p3y*xa*xb*ya*yb - 2*p1x*p2y*p3x*p3y*xa*xb*ya*yb - 
        2*pow(p1x,2)*pow(p3y,2)*xa*xb*ya*yb + 
        2*p1x*p2x*pow(p3y,2)*xa*xb*ya*yb - 
        pow(p1y,2)*pow(p3x,2)*xa*pow(yb,2) + 
        p1y*p2y*pow(p3x,2)*xa*pow(yb,2) + 2*p1x*p1y*p3x*p3y*xa*pow(yb,2) - 
        p1y*p2x*p3x*p3y*xa*pow(yb,2) - p1x*p2y*p3x*p3y*xa*pow(yb,2) - 
        pow(p1x,2)*pow(p3y,2)*xa*pow(yb,2) + 
        p1x*p2x*pow(p3y,2)*xa*pow(yb,2) - 
        pow(p1y,2)*p2x*p3x*pow(xa,2)*pow(yb,2) + 
        p1x*p1y*p2y*p3x*pow(xa,2)*pow(yb,2) + 
        p1y*p2x*p2y*p3x*pow(xa,2)*pow(yb,2) - 
        p1x*pow(p2y,2)*p3x*pow(xa,2)*pow(yb,2) + 
        pow(p1y,2)*pow(p3x,2)*pow(xa,2)*pow(yb,2) - 
        p1y*p2y*pow(p3x,2)*pow(xa,2)*pow(yb,2) + 
        p1x*p1y*p2x*p3y*pow(xa,2)*pow(yb,2) - 
        p1y*pow(p2x,2)*p3y*pow(xa,2)*pow(yb,2) - 
        pow(p1x,2)*p2y*p3y*pow(xa,2)*pow(yb,2) + 
        p1x*p2x*p2y*p3y*pow(xa,2)*pow(yb,2) - 
        2*p1x*p1y*p3x*p3y*pow(xa,2)*pow(yb,2) + 
        p1y*p2x*p3x*p3y*pow(xa,2)*pow(yb,2) + 
        p1x*p2y*p3x*p3y*pow(xa,2)*pow(yb,2) + 
        pow(p1x,2)*pow(p3y,2)*pow(xa,2)*pow(yb,2) - 
        p1x*p2x*pow(p3y,2)*pow(xa,2)*pow(yb,2) + 
        pow(p0y,2)*(pow(p1x,2)*(xb - xb*ya + xa*(-1 + yb))*(ya - yb) - 
           pow(p3x,2)*(xa - xb)*(ya - xb*ya + (-1 + xa)*yb) - 
           pow(p2x,2)*pow(xb*ya - xa*yb,2) + 
           p2x*p3x*(pow(xa,2)*yb*(1 + yb) + xb*ya*(xb - ya + xb*ya + yb) - 
              xa*(yb*(-ya + yb) + xb*(ya + yb + 2*ya*yb))) + 
           p1x*(-(p3x*((-1 + xb)*xb*pow(ya,2) - 2*xb*yb + xa*(2 + xb - yb)*yb + 
                   pow(xa,2)*(-1 + yb)*yb + xb*ya*(2 - xb + yb) + 
                   xa*ya*(-2 + xb + yb - 2*xb*yb))) + 
              p2x*(xb*ya*(xb*(-1 + ya) + ya - yb) + pow(xa,2)*(-1 + yb)*yb + 
                 xa*(yb*(-ya + yb) + xb*(ya + yb - 2*ya*yb))))) + 
        pow(p0x,2)*(pow(p1y,2)*(xb - xb*ya + xa*(-1 + yb))*(ya - yb) - 
           pow(p3y,2)*(xa - xb)*(ya - xb*ya + (-1 + xa)*yb) - 
           pow(p2y,2)*pow(xb*ya - xa*yb,2) + 
           p2y*p3y*(pow(xa,2)*yb*(1 + yb) + xb*ya*(xb - ya + xb*ya + yb) - 
              xa*(yb*(-ya + yb) + xb*(ya + yb + 2*ya*yb))) + 
           p1y*(-(p3y*((-1 + xb)*xb*pow(ya,2) - 2*xb*yb + xa*(2 + xb - yb)*yb + 
                   pow(xa,2)*(-1 + yb)*yb + xb*ya*(2 - xb + yb) + 
                   xa*ya*(-2 + xb + yb - 2*xb*yb))) + 
              p2y*(xb*ya*(xb*(-1 + ya) + ya - yb) + pow(xa,2)*(-1 + yb)*yb + 
                 xa*(yb*(-ya + yb) + xb*(ya + yb - 2*ya*yb))))) + 
        p0x*(pow(p1y,2)*(xa + xb*(-1 + ya) - xa*yb)*
            (-(p3x*(-2 + xb)*ya) + p3x*(-2 + xa)*yb + p2x*(xb*ya - xa*yb)) - 
           (-(p2y*xb*ya) + p2y*xa*yb + p3y*((-1 + xb)*ya + yb - xa*yb))*
            ((p2y*p3x - p2x*p3y)*(xb*ya - xa*yb) + 
              p1x*(p3y*xb*(-2 + ya) - p3y*xa*(-2 + yb) + p2y*(xb*ya - xa*yb))) + 
           p1y*(p1x*(xb - xb*ya + xa*(-1 + yb))*
               (-(p3y*(-2 + xb)*ya) + p3y*(-2 + xa)*yb + p2y*(xb*ya - xa*yb)) + 
              p2x*(-(xb*ya) + xa*yb)*
               (p3y*(xa - xb - ya + yb) + p2y*(xb*ya - xa*yb)) + 
              p3x*(p2y*(pow(xa,2)*yb + xb*ya*(xb - 2*ya + 2*yb) - 
                    xa*(2*yb*(-ya + yb) + xb*(ya + yb))) + 
                 p3y*(pow(xa,2)*(-2 + yb)*yb + xb*(-2 + ya)*((-1 + xb)*ya + yb) + 
                    xa*((2 + 2*xb - yb)*yb + ya*(-2 + 2*xb + yb - 2*xb*yb)))))) - 
        p0y*(pow(p1x,2)*(xb - xb*ya + xa*(-1 + yb))*
            (-(p3y*(-2 + xb)*ya) + p3y*(-2 + xa)*yb + p2y*(xb*ya - xa*yb)) + 
           (-(p2x*xb*ya) + p2x*xa*yb + p3x*((-1 + xb)*ya + yb - xa*yb))*
            (-((p2y*p3x - p2x*p3y)*(xb*ya - xa*yb)) + 
              p1y*(p3x*xb*(-2 + ya) - p3x*xa*(-2 + yb) + p2x*(xb*ya - xa*yb))) + 
           p1x*(p1y*(xa + xb*(-1 + ya) - xa*yb)*
               (-(p3x*(-2 + xb)*ya) + p3x*(-2 + xa)*yb + p2x*(xb*ya - xa*yb)) + 
              p2x*(xb*ya - xa*yb)*
               (p3y*(xa - xb + 2*ya - 2*yb) + p2y*(xb*ya - xa*yb)) - 
              p3x*(p2y*(xb*ya*(xb + ya - yb) + pow(xa,2)*yb - 
                    xa*((ya - yb)*yb + xb*(ya + yb))) + 
                 p3y*(pow(xa,2)*(-2 + yb)*yb + xb*(-2 + ya)*((-1 + xb)*ya + yb) + 
                    xa*((2 + 2*xb - yb)*yb + ya*(-2 + 2*xb + yb - 2*xb*yb))))) + 
           p0x*(-2*p3x*p3y*xa*ya + 2*p3x*p3y*xb*ya - p2y*p3x*xa*xb*ya - 
              p2x*p3y*xa*xb*ya + 2*p3x*p3y*xa*xb*ya + p2y*p3x*pow(xb,2)*ya + 
              p2x*p3y*pow(xb,2)*ya - 2*p3x*p3y*pow(xb,2)*ya - 
              p2y*p3x*xb*pow(ya,2) - p2x*p3y*xb*pow(ya,2) - 
              2*p2x*p2y*pow(xb,2)*pow(ya,2) + p2y*p3x*pow(xb,2)*pow(ya,2) + 
              p2x*p3y*pow(xb,2)*pow(ya,2) + 2*p3x*p3y*xa*yb + 
              p2y*p3x*pow(xa,2)*yb + p2x*p3y*pow(xa,2)*yb - 
              2*p3x*p3y*pow(xa,2)*yb - 2*p3x*p3y*xb*yb - p2y*p3x*xa*xb*yb - 
              p2x*p3y*xa*xb*yb + 2*p3x*p3y*xa*xb*yb + p2y*p3x*xa*ya*yb + 
              p2x*p3y*xa*ya*yb + p2y*p3x*xb*ya*yb + p2x*p3y*xb*ya*yb + 
              4*p2x*p2y*xa*xb*ya*yb - 2*p2y*p3x*xa*xb*ya*yb - 
              2*p2x*p3y*xa*xb*ya*yb - p2y*p3x*xa*pow(yb,2) - 
              p2x*p3y*xa*pow(yb,2) - 2*p2x*p2y*pow(xa,2)*pow(yb,2) + 
              p2y*p3x*pow(xa,2)*pow(yb,2) + p2x*p3y*pow(xa,2)*pow(yb,2) - 
              p1y*p3x*((-1 + xb)*xb*pow(ya,2) - 2*xb*yb + xa*(2 + xb - yb)*yb + 
                 pow(xa,2)*(-1 + yb)*yb + xb*ya*(2 - xb + yb) + 
                 xa*ya*(-2 + xb + yb - 2*xb*yb)) + 
              p1y*p2x*(xb*ya*(xb*(-1 + ya) + ya - yb) + pow(xa,2)*(-1 + yb)*yb + 
                 xa*(yb*(-ya + yb) + xb*(ya + yb - 2*ya*yb))) + 
              p1x*(2*p1y*(xb - xb*ya + xa*(-1 + yb))*(ya - yb) - 
                 p3y*((-1 + xb)*xb*pow(ya,2) - 2*xb*yb + xa*(2 + xb - yb)*yb + 
                    pow(xa,2)*(-1 + yb)*yb + xb*ya*(2 - xb + yb) + 
                    xa*ya*(-2 + xb + yb - 2*xb*yb)) + 
                 p2y*(xb*ya*(xb*(-1 + ya) + ya - yb) + pow(xa,2)*(-1 + yb)*yb + 
                    xa*(yb*(-ya + yb) + xb*(ya + yb - 2*ya*yb)))))))*
      log(p1y*p3x - p1x*p3y + p1y*p2x*xb - p1x*p2y*xb - p1y*p3x*xb + p1x*p3y*xb + 
        p0x*(p3y - p3y*xb + p2y*(xb - yb) + p1y*(-1 + yb)) - p1y*p3x*yb + 
        p2y*p3x*yb + p1x*p3y*yb - p2x*p3y*yb + 
        p0y*(p1x + p3x*(-1 + xb) - p1x*yb + p2x*(-xb + yb))))/
   (2.*pow(-(p0x*p2y*xa) + p1x*p2y*xa + p0x*p3y*xa - p1x*p3y*xa + p0x*p2y*xb - 
       p1x*p2y*xb - p0x*p3y*xb + p1x*p3y*xb + p0x*p2y*ya - p2y*p3x*ya - p1x*p3y*ya + 
       p2x*p3y*ya - p0x*p2y*yb + p2y*p3x*yb + p1x*p3y*yb - p2x*p3y*yb + 
       p1y*(p2x*(-xa + xb) + p3x*(xa - xb + ya - yb) + p0x*(-ya + yb)) + 
       p0y*(p3x*(-xa + xb) + p1x*(ya - yb) + p2x*(xa - xb - ya + yb)),3));

    integrals[3] = ((xa - xb)*(ya - yb)*pow(-(p0x*p2y*xa) + p1x*p2y*xa + p0x*p3y*xa - p1x*p3y*xa + 
        p0x*p2y*xb - p1x*p2y*xb - p0x*p3y*xb + p1x*p3y*xb + p0x*p2y*ya - 
        p2y*p3x*ya - p1x*p3y*ya + p2x*p3y*ya - p0x*p2y*yb + p2y*p3x*yb + 
        p1x*p3y*yb - p2x*p3y*yb + p1y*
         (p2x*(-xa + xb) + p3x*(xa - xb + ya - yb) + p0x*(-ya + yb)) + 
        p0y*(p3x*(-xa + xb) + p1x*(ya - yb) + p2x*(xa - xb - ya + yb)),2) - 
     2*(-(p0x*p2y*xa) + p1x*p2y*xa + p0x*p3y*xa - p1x*p3y*xa + p0x*p2y*xb - 
        p1x*p2y*xb - p0x*p3y*xb + p1x*p3y*xb + p0x*p2y*ya - p2y*p3x*ya - 
        p1x*p3y*ya + p2x*p3y*ya - p0x*p2y*yb + p2y*p3x*yb + p1x*p3y*yb - 
        p2x*p3y*yb + p1y*(p2x*(-xa + xb) + p3x*(xa - xb + ya - yb) + 
           p0x*(-ya + yb)) + p0y*(p3x*(-xa + xb) + p1x*(ya - yb) + 
           p2x*(xa - xb - ya + yb)))*
      (p1y*p2x*xa*ya - p1x*p2y*xa*ya - p1y*p2x*pow(xa,2)*ya + 
        p1x*p2y*pow(xa,2)*ya + p1y*p3x*pow(xa,2)*ya - p1x*p3y*pow(xa,2)*ya - 
        p1y*p2x*xb*ya + p1x*p2y*xb*ya + 2*p1y*p2x*xa*xb*ya - 2*p1x*p2y*xa*xb*ya - 
        2*p1y*p3x*xa*xb*ya + 2*p1x*p3y*xa*xb*ya - p1y*p2x*pow(xb,2)*ya + 
        p1x*p2y*pow(xb,2)*ya + p1y*p3x*pow(xb,2)*ya - p1x*p3y*pow(xb,2)*ya - 
        p1y*p3x*pow(ya,2) + p2y*p3x*pow(ya,2) + p1x*p3y*pow(ya,2) - 
        p2x*p3y*pow(ya,2) + p1y*p3x*xa*pow(ya,2) - p2y*p3x*xa*pow(ya,2) - 
        p1x*p3y*xa*pow(ya,2) + p2x*p3y*xa*pow(ya,2) - p1y*p2x*xa*yb + 
        p1x*p2y*xa*yb + p1y*p2x*xb*yb - p1x*p2y*xb*yb + 2*p1y*p3x*ya*yb - 
        2*p2y*p3x*ya*yb - 2*p1x*p3y*ya*yb + 2*p2x*p3y*ya*yb - 2*p1y*p3x*xa*ya*yb + 
        2*p2y*p3x*xa*ya*yb + 2*p1x*p3y*xa*ya*yb - 2*p2x*p3y*xa*ya*yb - 
        p1y*p3x*pow(yb,2) + p2y*p3x*pow(yb,2) + p1x*p3y*pow(yb,2) - 
        p2x*p3y*pow(yb,2) + p1y*p3x*xa*pow(yb,2) - p2y*p3x*xa*pow(yb,2) - 
        p1x*p3y*xa*pow(yb,2) + p2x*p3y*xa*pow(yb,2) + 
        p0y*(-(p3x*pow(xa - xb,2)*ya) + 
           p1x*(ya - yb)*(-xb - ya + xa*(1 + ya - yb) + yb) + 
           p2x*(pow(xa,2)*ya + pow(xb,2)*ya + xb*(ya - yb) + pow(ya - yb,2) - 
              xa*(pow(ya,2) + ya*(1 + 2*xb - 2*yb) + (-1 + yb)*yb))) + 
        p0x*(p3y*pow(xa - xb,2)*ya - 
           p1y*(ya - yb)*(-xb - ya + xa*(1 + ya - yb) + yb) - 
           p2y*(pow(xa,2)*ya + pow(xb,2)*ya + xb*(ya - yb) + pow(ya - yb,2) - 
              xa*(pow(ya,2) + ya*(1 + 2*xb - 2*yb) + (-1 + yb)*yb)))) + 
     2*(-(pow(p1y,2)*p2x*p3x*xa*ya) + p1x*p1y*p2y*p3x*xa*ya + 
        p1x*p1y*p2x*p3y*xa*ya - pow(p1x,2)*p2y*p3y*xa*ya + 
        pow(p1y,2)*p2x*p3x*xb*ya - p1x*p1y*p2y*p3x*xb*ya - p1x*p1y*p2x*p3y*xb*ya + 
        pow(p1x,2)*p2y*p3y*xb*ya - pow(p1y,2)*pow(p2x,2)*xa*xb*ya + 
        2*p1x*p1y*p2x*p2y*xa*xb*ya - pow(p1x,2)*pow(p2y,2)*xa*xb*ya + 
        pow(p1y,2)*p2x*p3x*xa*xb*ya - p1x*p1y*p2y*p3x*xa*xb*ya - 
        p1x*p1y*p2x*p3y*xa*xb*ya + pow(p1x,2)*p2y*p3y*xa*xb*ya + 
        pow(p1y,2)*pow(p2x,2)*pow(xb,2)*ya - 
        2*p1x*p1y*p2x*p2y*pow(xb,2)*ya + 
        pow(p1x,2)*pow(p2y,2)*pow(xb,2)*ya - 
        pow(p1y,2)*p2x*p3x*pow(xb,2)*ya + p1x*p1y*p2y*p3x*pow(xb,2)*ya + 
        p1x*p1y*p2x*p3y*pow(xb,2)*ya - pow(p1x,2)*p2y*p3y*pow(xb,2)*ya + 
        pow(p1y,2)*pow(p3x,2)*pow(ya,2) - p1y*p2y*pow(p3x,2)*pow(ya,2) - 
        2*p1x*p1y*p3x*p3y*pow(ya,2) + p1y*p2x*p3x*p3y*pow(ya,2) + 
        p1x*p2y*p3x*p3y*pow(ya,2) + pow(p1x,2)*pow(p3y,2)*pow(ya,2) - 
        p1x*p2x*pow(p3y,2)*pow(ya,2) + pow(p1y,2)*p2x*p3x*xb*pow(ya,2) - 
        p1x*p1y*p2y*p3x*xb*pow(ya,2) - p1y*p2x*p2y*p3x*xb*pow(ya,2) + 
        p1x*pow(p2y,2)*p3x*xb*pow(ya,2) - 
        2*pow(p1y,2)*pow(p3x,2)*xb*pow(ya,2) + 
        2*p1y*p2y*pow(p3x,2)*xb*pow(ya,2) - p1x*p1y*p2x*p3y*xb*pow(ya,2) + 
        p1y*pow(p2x,2)*p3y*xb*pow(ya,2) + pow(p1x,2)*p2y*p3y*xb*pow(ya,2) - 
        p1x*p2x*p2y*p3y*xb*pow(ya,2) + 4*p1x*p1y*p3x*p3y*xb*pow(ya,2) - 
        2*p1y*p2x*p3x*p3y*xb*pow(ya,2) - 2*p1x*p2y*p3x*p3y*xb*pow(ya,2) - 
        2*pow(p1x,2)*pow(p3y,2)*xb*pow(ya,2) + 
        2*p1x*p2x*pow(p3y,2)*xb*pow(ya,2) - 
        pow(p1y,2)*p2x*p3x*pow(xb,2)*pow(ya,2) + 
        p1x*p1y*p2y*p3x*pow(xb,2)*pow(ya,2) + 
        p1y*p2x*p2y*p3x*pow(xb,2)*pow(ya,2) - 
        p1x*pow(p2y,2)*p3x*pow(xb,2)*pow(ya,2) + 
        pow(p1y,2)*pow(p3x,2)*pow(xb,2)*pow(ya,2) - 
        p1y*p2y*pow(p3x,2)*pow(xb,2)*pow(ya,2) + 
        p1x*p1y*p2x*p3y*pow(xb,2)*pow(ya,2) - 
        p1y*pow(p2x,2)*p3y*pow(xb,2)*pow(ya,2) - 
        pow(p1x,2)*p2y*p3y*pow(xb,2)*pow(ya,2) + 
        p1x*p2x*p2y*p3y*pow(xb,2)*pow(ya,2) - 
        2*p1x*p1y*p3x*p3y*pow(xb,2)*pow(ya,2) + 
        p1y*p2x*p3x*p3y*pow(xb,2)*pow(ya,2) + 
        p1x*p2y*p3x*p3y*pow(xb,2)*pow(ya,2) + 
        pow(p1x,2)*pow(p3y,2)*pow(xb,2)*pow(ya,2) - 
        p1x*p2x*pow(p3y,2)*pow(xb,2)*pow(ya,2) + pow(p1y,2)*p2x*p3x*xa*yb - 
        p1x*p1y*p2y*p3x*xa*yb - p1x*p1y*p2x*p3y*xa*yb + pow(p1x,2)*p2y*p3y*xa*yb + 
        pow(p1y,2)*pow(p2x,2)*pow(xa,2)*yb - 
        2*p1x*p1y*p2x*p2y*pow(xa,2)*yb + 
        pow(p1x,2)*pow(p2y,2)*pow(xa,2)*yb - 
        pow(p1y,2)*p2x*p3x*pow(xa,2)*yb + p1x*p1y*p2y*p3x*pow(xa,2)*yb + 
        p1x*p1y*p2x*p3y*pow(xa,2)*yb - pow(p1x,2)*p2y*p3y*pow(xa,2)*yb - 
        pow(p1y,2)*p2x*p3x*xb*yb + p1x*p1y*p2y*p3x*xb*yb + p1x*p1y*p2x*p3y*xb*yb - 
        pow(p1x,2)*p2y*p3y*xb*yb - pow(p1y,2)*pow(p2x,2)*xa*xb*yb + 
        2*p1x*p1y*p2x*p2y*xa*xb*yb - pow(p1x,2)*pow(p2y,2)*xa*xb*yb + 
        pow(p1y,2)*p2x*p3x*xa*xb*yb - p1x*p1y*p2y*p3x*xa*xb*yb - 
        p1x*p1y*p2x*p3y*xa*xb*yb + pow(p1x,2)*p2y*p3y*xa*xb*yb - 
        2*pow(p1y,2)*pow(p3x,2)*ya*yb + 2*p1y*p2y*pow(p3x,2)*ya*yb + 
        4*p1x*p1y*p3x*p3y*ya*yb - 2*p1y*p2x*p3x*p3y*ya*yb - 
        2*p1x*p2y*p3x*p3y*ya*yb - 2*pow(p1x,2)*pow(p3y,2)*ya*yb + 
        2*p1x*p2x*pow(p3y,2)*ya*yb - pow(p1y,2)*p2x*p3x*xa*ya*yb + 
        p1x*p1y*p2y*p3x*xa*ya*yb + p1y*p2x*p2y*p3x*xa*ya*yb - 
        p1x*pow(p2y,2)*p3x*xa*ya*yb + 2*pow(p1y,2)*pow(p3x,2)*xa*ya*yb - 
        2*p1y*p2y*pow(p3x,2)*xa*ya*yb + p1x*p1y*p2x*p3y*xa*ya*yb - 
        p1y*pow(p2x,2)*p3y*xa*ya*yb - pow(p1x,2)*p2y*p3y*xa*ya*yb + 
        p1x*p2x*p2y*p3y*xa*ya*yb - 4*p1x*p1y*p3x*p3y*xa*ya*yb + 
        2*p1y*p2x*p3x*p3y*xa*ya*yb + 2*p1x*p2y*p3x*p3y*xa*ya*yb + 
        2*pow(p1x,2)*pow(p3y,2)*xa*ya*yb - 2*p1x*p2x*pow(p3y,2)*xa*ya*yb - 
        pow(p1y,2)*p2x*p3x*xb*ya*yb + p1x*p1y*p2y*p3x*xb*ya*yb + 
        p1y*p2x*p2y*p3x*xb*ya*yb - p1x*pow(p2y,2)*p3x*xb*ya*yb + 
        2*pow(p1y,2)*pow(p3x,2)*xb*ya*yb - 2*p1y*p2y*pow(p3x,2)*xb*ya*yb + 
        p1x*p1y*p2x*p3y*xb*ya*yb - p1y*pow(p2x,2)*p3y*xb*ya*yb - 
        pow(p1x,2)*p2y*p3y*xb*ya*yb + p1x*p2x*p2y*p3y*xb*ya*yb - 
        4*p1x*p1y*p3x*p3y*xb*ya*yb + 2*p1y*p2x*p3x*p3y*xb*ya*yb + 
        2*p1x*p2y*p3x*p3y*xb*ya*yb + 2*pow(p1x,2)*pow(p3y,2)*xb*ya*yb - 
        2*p1x*p2x*pow(p3y,2)*xb*ya*yb + 2*pow(p1y,2)*p2x*p3x*xa*xb*ya*yb - 
        2*p1x*p1y*p2y*p3x*xa*xb*ya*yb - 2*p1y*p2x*p2y*p3x*xa*xb*ya*yb + 
        2*p1x*pow(p2y,2)*p3x*xa*xb*ya*yb - 
        2*pow(p1y,2)*pow(p3x,2)*xa*xb*ya*yb + 
        2*p1y*p2y*pow(p3x,2)*xa*xb*ya*yb - 2*p1x*p1y*p2x*p3y*xa*xb*ya*yb + 
        2*p1y*pow(p2x,2)*p3y*xa*xb*ya*yb + 2*pow(p1x,2)*p2y*p3y*xa*xb*ya*yb - 
        2*p1x*p2x*p2y*p3y*xa*xb*ya*yb + 4*p1x*p1y*p3x*p3y*xa*xb*ya*yb - 
        2*p1y*p2x*p3x*p3y*xa*xb*ya*yb - 2*p1x*p2y*p3x*p3y*xa*xb*ya*yb - 
        2*pow(p1x,2)*pow(p3y,2)*xa*xb*ya*yb + 
        2*p1x*p2x*pow(p3y,2)*xa*xb*ya*yb + pow(p1y,2)*pow(p3x,2)*pow(yb,2) - 
        p1y*p2y*pow(p3x,2)*pow(yb,2) - 2*p1x*p1y*p3x*p3y*pow(yb,2) + 
        p1y*p2x*p3x*p3y*pow(yb,2) + p1x*p2y*p3x*p3y*pow(yb,2) + 
        pow(p1x,2)*pow(p3y,2)*pow(yb,2) - p1x*p2x*pow(p3y,2)*pow(yb,2) + 
        pow(p1y,2)*p2x*p3x*xa*pow(yb,2) - p1x*p1y*p2y*p3x*xa*pow(yb,2) - 
        p1y*p2x*p2y*p3x*xa*pow(yb,2) + p1x*pow(p2y,2)*p3x*xa*pow(yb,2) - 
        2*pow(p1y,2)*pow(p3x,2)*xa*pow(yb,2) + 
        2*p1y*p2y*pow(p3x,2)*xa*pow(yb,2) - p1x*p1y*p2x*p3y*xa*pow(yb,2) + 
        p1y*pow(p2x,2)*p3y*xa*pow(yb,2) + pow(p1x,2)*p2y*p3y*xa*pow(yb,2) - 
        p1x*p2x*p2y*p3y*xa*pow(yb,2) + 4*p1x*p1y*p3x*p3y*xa*pow(yb,2) - 
        2*p1y*p2x*p3x*p3y*xa*pow(yb,2) - 2*p1x*p2y*p3x*p3y*xa*pow(yb,2) - 
        2*pow(p1x,2)*pow(p3y,2)*xa*pow(yb,2) + 
        2*p1x*p2x*pow(p3y,2)*xa*pow(yb,2) - 
        pow(p1y,2)*p2x*p3x*pow(xa,2)*pow(yb,2) + 
        p1x*p1y*p2y*p3x*pow(xa,2)*pow(yb,2) + 
        p1y*p2x*p2y*p3x*pow(xa,2)*pow(yb,2) - 
        p1x*pow(p2y,2)*p3x*pow(xa,2)*pow(yb,2) + 
        pow(p1y,2)*pow(p3x,2)*pow(xa,2)*pow(yb,2) - 
        p1y*p2y*pow(p3x,2)*pow(xa,2)*pow(yb,2) + 
        p1x*p1y*p2x*p3y*pow(xa,2)*pow(yb,2) - 
        p1y*pow(p2x,2)*p3y*pow(xa,2)*pow(yb,2) - 
        pow(p1x,2)*p2y*p3y*pow(xa,2)*pow(yb,2) + 
        p1x*p2x*p2y*p3y*pow(xa,2)*pow(yb,2) - 
        2*p1x*p1y*p3x*p3y*pow(xa,2)*pow(yb,2) + 
        p1y*p2x*p3x*p3y*pow(xa,2)*pow(yb,2) + 
        p1x*p2y*p3x*p3y*pow(xa,2)*pow(yb,2) + 
        pow(p1x,2)*pow(p3y,2)*pow(xa,2)*pow(yb,2) - 
        p1x*p2x*pow(p3y,2)*pow(xa,2)*pow(yb,2) + 
        pow(p0y,2)*(p1x - p2x)*(xb + ya - xb*ya + xa*(-1 + yb) - yb)*
         (-(p2x*xb*ya) + p1x*(ya - yb) + p2x*xa*yb + p3x*((-1 + xb)*ya + yb - xa*yb))
          + pow(p0x,2)*(p1y - p2y)*(xb + ya - xb*ya + xa*(-1 + yb) - yb)*
         (-(p2y*xb*ya) + p1y*(ya - yb) + p2y*xa*yb + p3y*((-1 + xb)*ya + yb - xa*yb))
          - p0y*((-(p2x*xb*ya) + p2x*xa*yb + p3x*((-1 + xb)*ya + yb - xa*yb))*
            (-((p2y*p3x - p2x*p3y)*((-1 + xb)*ya + yb - xa*yb)) + 
              p1y*(p2x*(xb*(-2 + ya) - ya - xa*(-2 + yb) + yb) + 
                 p3x*((-1 + xb)*ya + yb - xa*yb))) + 
           p0x*(xb + ya - xb*ya + xa*(-1 + yb) - yb)*
            (p2y*p3x*ya + p2x*p3y*ya + 2*p2x*p2y*xb*ya - p2y*p3x*xb*ya - 
              p2x*p3y*xb*ya - p2y*p3x*yb - p2x*p3y*yb - 2*p2x*p2y*xa*yb + 
              p2y*p3x*xa*yb + p2x*p3y*xa*yb + p1y*p3x*((-1 + xb)*ya + yb - xa*yb) + 
              p1y*p2x*(-((1 + xb)*ya) + (1 + xa)*yb) + 
              p1x*(2*p1y*(ya - yb) + p3y*((-1 + xb)*ya + yb - xa*yb) + 
                 p2y*(-((1 + xb)*ya) + (1 + xa)*yb))) + 
           pow(p1x,2)*(p3y*(pow(xb,2)*(-1 + ya)*ya - xb*(-1 + 3*ya)*(ya - yb) + 
                 2*pow(ya - yb,2) + xa*(1 + xb - 3*yb)*yb + 
                 pow(xa,2)*(-1 + yb)*yb + xa*ya*(-1 + xb + 3*yb - 2*xb*yb)) - 
              p2y*(pow(xa,2)*(-1 + yb)*yb + 
                 xb*((-1 + xb)*pow(ya,2) + yb + ya*(-1 - xb + yb)) + 
                 xa*((-1 + xb - yb)*yb + ya*(1 + xb + yb - 2*xb*yb)))) + 
           p1x*(p3y*(ya - xb*ya + (-1 + xa)*yb)*
               (p2x*(xa - xb - 2*ya + 2*yb) + p3x*((-1 + xb)*ya + yb - xa*yb)) + 
              p1y*(-(p3x*(pow(xb,2)*(-1 + ya)*ya - xb*(-1 + 3*ya)*(ya - yb) + 
                      2*pow(ya - yb,2) + xa*(1 + xb - 3*yb)*yb + 
                      pow(xa,2)*(-1 + yb)*yb + xa*ya*(-1 + xb + 3*yb - 2*xb*yb)))\
                  + p2x*(pow(xa,2)*(-1 + yb)*yb + 
                    xb*((-1 + xb)*pow(ya,2) + yb + ya*(-1 - xb + yb)) + 
                    xa*((-1 + xb - yb)*yb + ya*(1 + xb + yb - 2*xb*yb)))) + 
              p2y*(p3x*(pow(xb,2)*ya - xb*(1 + ya)*(ya - yb) + pow(ya - yb,2) + 
                    pow(xa,2)*yb - xa*(ya*(-1 + xb - yb) + yb*(1 + xb + yb))) + 
                 p2x*(pow(xa,2)*(-2 + yb)*yb + xb*ya*(xb*(-2 + ya) - ya + yb) + 
                    xa*((ya - yb)*yb + 2*xb*(ya + yb - ya*yb)))))) + 
        p0x*((-(p2y*xb*ya) + p2y*xa*yb + p3y*((-1 + xb)*ya + yb - xa*yb))*
            (-((p2y*p3x - p2x*p3y)*((-1 + xb)*ya + yb - xa*yb)) + 
              p1x*(p2y*(2*xb + ya - xb*ya + xa*(-2 + yb) - yb) + 
                 p3y*(ya - xb*ya + (-1 + xa)*yb))) + 
           pow(p1y,2)*(-(p3x*(pow(xb,2)*(-1 + ya)*ya - 
                   xb*(-1 + 3*ya)*(ya - yb) + 2*pow(ya - yb,2) + 
                   xa*(1 + xb - 3*yb)*yb + pow(xa,2)*(-1 + yb)*yb + 
                   xa*ya*(-1 + xb + 3*yb - 2*xb*yb))) + 
              p2x*(pow(xa,2)*(-1 + yb)*yb + 
                 xb*((-1 + xb)*pow(ya,2) + yb + ya*(-1 - xb + yb)) + 
                 xa*((-1 + xb - yb)*yb + ya*(1 + xb + yb - 2*xb*yb)))) + 
           p1y*(p3y*((-1 + xb)*ya + yb - xa*yb)*
               (p2x*(xa - xb + ya - yb) + p3x*((-1 + xb)*ya + yb - xa*yb)) + 
              p1x*(p3y*(pow(xb,2)*(-1 + ya)*ya - xb*(-1 + 3*ya)*(ya - yb) + 
                    2*pow(ya - yb,2) + xa*(1 + xb - 3*yb)*yb + 
                    pow(xa,2)*(-1 + yb)*yb + xa*ya*(-1 + xb + 3*yb - 2*xb*yb)) - 
                 p2y*(pow(xa,2)*(-1 + yb)*yb + 
                    xb*((-1 + xb)*pow(ya,2) + yb + ya*(-1 - xb + yb)) + 
                    xa*((-1 + xb - yb)*yb + ya*(1 + xb + yb - 2*xb*yb)))) - 
              p2y*(p3x*(pow(xb,2)*ya + xb*(-1 + 2*ya)*(ya - yb) - 
                    2*pow(ya - yb,2) + pow(xa,2)*yb - 
                    xa*((1 + xb - 2*yb)*yb + ya*(-1 + xb + 2*yb))) + 
                 p2x*(pow(xa,2)*(-2 + yb)*yb + xb*ya*(xb*(-2 + ya) - ya + yb) + 
                    xa*((ya - yb)*yb + 2*xb*(ya + yb - ya*yb)))))))*
      log(p1y*p3x - p1x*p3y + p1y*p2x*xa - p1x*p2y*xa - p1y*p3x*xa + p1x*p3y*xa + 
        p0x*(p3y - p3y*xa + p2y*(xa - ya) + p1y*(-1 + ya)) - p1y*p3x*ya + 
        p2y*p3x*ya + p1x*p3y*ya - p2x*p3y*ya + 
        p0y*(p1x + p3x*(-1 + xa) - p1x*ya + p2x*(-xa + ya))) + 
     2*(pow(p1y,2)*p2x*p3x*xa*ya - p1x*p1y*p2y*p3x*xa*ya - p1x*p1y*p2x*p3y*xa*ya + 
        pow(p1x,2)*p2y*p3y*xa*ya - pow(p1y,2)*p2x*p3x*xb*ya + 
        p1x*p1y*p2y*p3x*xb*ya + p1x*p1y*p2x*p3y*xb*ya - pow(p1x,2)*p2y*p3y*xb*ya + 
        pow(p1y,2)*pow(p2x,2)*xa*xb*ya - 2*p1x*p1y*p2x*p2y*xa*xb*ya + 
        pow(p1x,2)*pow(p2y,2)*xa*xb*ya - pow(p1y,2)*p2x*p3x*xa*xb*ya + 
        p1x*p1y*p2y*p3x*xa*xb*ya + p1x*p1y*p2x*p3y*xa*xb*ya - 
        pow(p1x,2)*p2y*p3y*xa*xb*ya - pow(p1y,2)*pow(p2x,2)*pow(xb,2)*ya + 
        2*p1x*p1y*p2x*p2y*pow(xb,2)*ya - 
        pow(p1x,2)*pow(p2y,2)*pow(xb,2)*ya + 
        pow(p1y,2)*p2x*p3x*pow(xb,2)*ya - p1x*p1y*p2y*p3x*pow(xb,2)*ya - 
        p1x*p1y*p2x*p3y*pow(xb,2)*ya + pow(p1x,2)*p2y*p3y*pow(xb,2)*ya - 
        pow(p1y,2)*pow(p3x,2)*pow(ya,2) + p1y*p2y*pow(p3x,2)*pow(ya,2) + 
        2*p1x*p1y*p3x*p3y*pow(ya,2) - p1y*p2x*p3x*p3y*pow(ya,2) - 
        p1x*p2y*p3x*p3y*pow(ya,2) - pow(p1x,2)*pow(p3y,2)*pow(ya,2) + 
        p1x*p2x*pow(p3y,2)*pow(ya,2) - pow(p1y,2)*p2x*p3x*xb*pow(ya,2) + 
        p1x*p1y*p2y*p3x*xb*pow(ya,2) + p1y*p2x*p2y*p3x*xb*pow(ya,2) - 
        p1x*pow(p2y,2)*p3x*xb*pow(ya,2) + 
        2*pow(p1y,2)*pow(p3x,2)*xb*pow(ya,2) - 
        2*p1y*p2y*pow(p3x,2)*xb*pow(ya,2) + p1x*p1y*p2x*p3y*xb*pow(ya,2) - 
        p1y*pow(p2x,2)*p3y*xb*pow(ya,2) - pow(p1x,2)*p2y*p3y*xb*pow(ya,2) + 
        p1x*p2x*p2y*p3y*xb*pow(ya,2) - 4*p1x*p1y*p3x*p3y*xb*pow(ya,2) + 
        2*p1y*p2x*p3x*p3y*xb*pow(ya,2) + 2*p1x*p2y*p3x*p3y*xb*pow(ya,2) + 
        2*pow(p1x,2)*pow(p3y,2)*xb*pow(ya,2) - 
        2*p1x*p2x*pow(p3y,2)*xb*pow(ya,2) + 
        pow(p1y,2)*p2x*p3x*pow(xb,2)*pow(ya,2) - 
        p1x*p1y*p2y*p3x*pow(xb,2)*pow(ya,2) - 
        p1y*p2x*p2y*p3x*pow(xb,2)*pow(ya,2) + 
        p1x*pow(p2y,2)*p3x*pow(xb,2)*pow(ya,2) - 
        pow(p1y,2)*pow(p3x,2)*pow(xb,2)*pow(ya,2) + 
        p1y*p2y*pow(p3x,2)*pow(xb,2)*pow(ya,2) - 
        p1x*p1y*p2x*p3y*pow(xb,2)*pow(ya,2) + 
        p1y*pow(p2x,2)*p3y*pow(xb,2)*pow(ya,2) + 
        pow(p1x,2)*p2y*p3y*pow(xb,2)*pow(ya,2) - 
        p1x*p2x*p2y*p3y*pow(xb,2)*pow(ya,2) + 
        2*p1x*p1y*p3x*p3y*pow(xb,2)*pow(ya,2) - 
        p1y*p2x*p3x*p3y*pow(xb,2)*pow(ya,2) - 
        p1x*p2y*p3x*p3y*pow(xb,2)*pow(ya,2) - 
        pow(p1x,2)*pow(p3y,2)*pow(xb,2)*pow(ya,2) + 
        p1x*p2x*pow(p3y,2)*pow(xb,2)*pow(ya,2) - pow(p1y,2)*p2x*p3x*xa*yb + 
        p1x*p1y*p2y*p3x*xa*yb + p1x*p1y*p2x*p3y*xa*yb - pow(p1x,2)*p2y*p3y*xa*yb - 
        pow(p1y,2)*pow(p2x,2)*pow(xa,2)*yb + 
        2*p1x*p1y*p2x*p2y*pow(xa,2)*yb - 
        pow(p1x,2)*pow(p2y,2)*pow(xa,2)*yb + 
        pow(p1y,2)*p2x*p3x*pow(xa,2)*yb - p1x*p1y*p2y*p3x*pow(xa,2)*yb - 
        p1x*p1y*p2x*p3y*pow(xa,2)*yb + pow(p1x,2)*p2y*p3y*pow(xa,2)*yb + 
        pow(p1y,2)*p2x*p3x*xb*yb - p1x*p1y*p2y*p3x*xb*yb - p1x*p1y*p2x*p3y*xb*yb + 
        pow(p1x,2)*p2y*p3y*xb*yb + pow(p1y,2)*pow(p2x,2)*xa*xb*yb - 
        2*p1x*p1y*p2x*p2y*xa*xb*yb + pow(p1x,2)*pow(p2y,2)*xa*xb*yb - 
        pow(p1y,2)*p2x*p3x*xa*xb*yb + p1x*p1y*p2y*p3x*xa*xb*yb + 
        p1x*p1y*p2x*p3y*xa*xb*yb - pow(p1x,2)*p2y*p3y*xa*xb*yb + 
        2*pow(p1y,2)*pow(p3x,2)*ya*yb - 2*p1y*p2y*pow(p3x,2)*ya*yb - 
        4*p1x*p1y*p3x*p3y*ya*yb + 2*p1y*p2x*p3x*p3y*ya*yb + 
        2*p1x*p2y*p3x*p3y*ya*yb + 2*pow(p1x,2)*pow(p3y,2)*ya*yb - 
        2*p1x*p2x*pow(p3y,2)*ya*yb + pow(p1y,2)*p2x*p3x*xa*ya*yb - 
        p1x*p1y*p2y*p3x*xa*ya*yb - p1y*p2x*p2y*p3x*xa*ya*yb + 
        p1x*pow(p2y,2)*p3x*xa*ya*yb - 2*pow(p1y,2)*pow(p3x,2)*xa*ya*yb + 
        2*p1y*p2y*pow(p3x,2)*xa*ya*yb - p1x*p1y*p2x*p3y*xa*ya*yb + 
        p1y*pow(p2x,2)*p3y*xa*ya*yb + pow(p1x,2)*p2y*p3y*xa*ya*yb - 
        p1x*p2x*p2y*p3y*xa*ya*yb + 4*p1x*p1y*p3x*p3y*xa*ya*yb - 
        2*p1y*p2x*p3x*p3y*xa*ya*yb - 2*p1x*p2y*p3x*p3y*xa*ya*yb - 
        2*pow(p1x,2)*pow(p3y,2)*xa*ya*yb + 2*p1x*p2x*pow(p3y,2)*xa*ya*yb + 
        pow(p1y,2)*p2x*p3x*xb*ya*yb - p1x*p1y*p2y*p3x*xb*ya*yb - 
        p1y*p2x*p2y*p3x*xb*ya*yb + p1x*pow(p2y,2)*p3x*xb*ya*yb - 
        2*pow(p1y,2)*pow(p3x,2)*xb*ya*yb + 2*p1y*p2y*pow(p3x,2)*xb*ya*yb - 
        p1x*p1y*p2x*p3y*xb*ya*yb + p1y*pow(p2x,2)*p3y*xb*ya*yb + 
        pow(p1x,2)*p2y*p3y*xb*ya*yb - p1x*p2x*p2y*p3y*xb*ya*yb + 
        4*p1x*p1y*p3x*p3y*xb*ya*yb - 2*p1y*p2x*p3x*p3y*xb*ya*yb - 
        2*p1x*p2y*p3x*p3y*xb*ya*yb - 2*pow(p1x,2)*pow(p3y,2)*xb*ya*yb + 
        2*p1x*p2x*pow(p3y,2)*xb*ya*yb - 2*pow(p1y,2)*p2x*p3x*xa*xb*ya*yb + 
        2*p1x*p1y*p2y*p3x*xa*xb*ya*yb + 2*p1y*p2x*p2y*p3x*xa*xb*ya*yb - 
        2*p1x*pow(p2y,2)*p3x*xa*xb*ya*yb + 
        2*pow(p1y,2)*pow(p3x,2)*xa*xb*ya*yb - 
        2*p1y*p2y*pow(p3x,2)*xa*xb*ya*yb + 2*p1x*p1y*p2x*p3y*xa*xb*ya*yb - 
        2*p1y*pow(p2x,2)*p3y*xa*xb*ya*yb - 2*pow(p1x,2)*p2y*p3y*xa*xb*ya*yb + 
        2*p1x*p2x*p2y*p3y*xa*xb*ya*yb - 4*p1x*p1y*p3x*p3y*xa*xb*ya*yb + 
        2*p1y*p2x*p3x*p3y*xa*xb*ya*yb + 2*p1x*p2y*p3x*p3y*xa*xb*ya*yb + 
        2*pow(p1x,2)*pow(p3y,2)*xa*xb*ya*yb - 
        2*p1x*p2x*pow(p3y,2)*xa*xb*ya*yb - pow(p1y,2)*pow(p3x,2)*pow(yb,2) + 
        p1y*p2y*pow(p3x,2)*pow(yb,2) + 2*p1x*p1y*p3x*p3y*pow(yb,2) - 
        p1y*p2x*p3x*p3y*pow(yb,2) - p1x*p2y*p3x*p3y*pow(yb,2) - 
        pow(p1x,2)*pow(p3y,2)*pow(yb,2) + p1x*p2x*pow(p3y,2)*pow(yb,2) - 
        pow(p1y,2)*p2x*p3x*xa*pow(yb,2) + p1x*p1y*p2y*p3x*xa*pow(yb,2) + 
        p1y*p2x*p2y*p3x*xa*pow(yb,2) - p1x*pow(p2y,2)*p3x*xa*pow(yb,2) + 
        2*pow(p1y,2)*pow(p3x,2)*xa*pow(yb,2) - 
        2*p1y*p2y*pow(p3x,2)*xa*pow(yb,2) + p1x*p1y*p2x*p3y*xa*pow(yb,2) - 
        p1y*pow(p2x,2)*p3y*xa*pow(yb,2) - pow(p1x,2)*p2y*p3y*xa*pow(yb,2) + 
        p1x*p2x*p2y*p3y*xa*pow(yb,2) - 4*p1x*p1y*p3x*p3y*xa*pow(yb,2) + 
        2*p1y*p2x*p3x*p3y*xa*pow(yb,2) + 2*p1x*p2y*p3x*p3y*xa*pow(yb,2) + 
        2*pow(p1x,2)*pow(p3y,2)*xa*pow(yb,2) - 
        2*p1x*p2x*pow(p3y,2)*xa*pow(yb,2) + 
        pow(p1y,2)*p2x*p3x*pow(xa,2)*pow(yb,2) - 
        p1x*p1y*p2y*p3x*pow(xa,2)*pow(yb,2) - 
        p1y*p2x*p2y*p3x*pow(xa,2)*pow(yb,2) + 
        p1x*pow(p2y,2)*p3x*pow(xa,2)*pow(yb,2) - 
        pow(p1y,2)*pow(p3x,2)*pow(xa,2)*pow(yb,2) + 
        p1y*p2y*pow(p3x,2)*pow(xa,2)*pow(yb,2) - 
        p1x*p1y*p2x*p3y*pow(xa,2)*pow(yb,2) + 
        p1y*pow(p2x,2)*p3y*pow(xa,2)*pow(yb,2) + 
        pow(p1x,2)*p2y*p3y*pow(xa,2)*pow(yb,2) - 
        p1x*p2x*p2y*p3y*pow(xa,2)*pow(yb,2) + 
        2*p1x*p1y*p3x*p3y*pow(xa,2)*pow(yb,2) - 
        p1y*p2x*p3x*p3y*pow(xa,2)*pow(yb,2) - 
        p1x*p2y*p3x*p3y*pow(xa,2)*pow(yb,2) - 
        pow(p1x,2)*pow(p3y,2)*pow(xa,2)*pow(yb,2) + 
        p1x*p2x*pow(p3y,2)*pow(xa,2)*pow(yb,2) + 
        pow(p0y,2)*(p1x - p2x)*(xa + xb*(-1 + ya) - ya + yb - xa*yb)*
         (-(p2x*xb*ya) + p1x*(ya - yb) + p2x*xa*yb + p3x*((-1 + xb)*ya + yb - xa*yb))
          + pow(p0x,2)*(p1y - p2y)*(xa + xb*(-1 + ya) - ya + yb - xa*yb)*
         (-(p2y*xb*ya) + p1y*(ya - yb) + p2y*xa*yb + p3y*((-1 + xb)*ya + yb - xa*yb))
          + p0y*((-(p2x*xb*ya) + p2x*xa*yb + p3x*((-1 + xb)*ya + yb - xa*yb))*
            (-((p2y*p3x - p2x*p3y)*((-1 + xb)*ya + yb - xa*yb)) + 
              p1y*(p2x*(xb*(-2 + ya) - ya - xa*(-2 + yb) + yb) + 
                 p3x*((-1 + xb)*ya + yb - xa*yb))) + 
           p0x*(xb + ya - xb*ya + xa*(-1 + yb) - yb)*
            (p2y*p3x*ya + p2x*p3y*ya + 2*p2x*p2y*xb*ya - p2y*p3x*xb*ya - 
              p2x*p3y*xb*ya - p2y*p3x*yb - p2x*p3y*yb - 2*p2x*p2y*xa*yb + 
              p2y*p3x*xa*yb + p2x*p3y*xa*yb + p1y*p3x*((-1 + xb)*ya + yb - xa*yb) + 
              p1y*p2x*(-((1 + xb)*ya) + (1 + xa)*yb) + 
              p1x*(2*p1y*(ya - yb) + p3y*((-1 + xb)*ya + yb - xa*yb) + 
                 p2y*(-((1 + xb)*ya) + (1 + xa)*yb))) + 
           pow(p1x,2)*(p3y*(pow(xb,2)*(-1 + ya)*ya - xb*(-1 + 3*ya)*(ya - yb) + 
                 2*pow(ya - yb,2) + xa*(1 + xb - 3*yb)*yb + 
                 pow(xa,2)*(-1 + yb)*yb + xa*ya*(-1 + xb + 3*yb - 2*xb*yb)) - 
              p2y*(pow(xa,2)*(-1 + yb)*yb + 
                 xb*((-1 + xb)*pow(ya,2) + yb + ya*(-1 - xb + yb)) + 
                 xa*((-1 + xb - yb)*yb + ya*(1 + xb + yb - 2*xb*yb)))) + 
           p1x*(p3y*(ya - xb*ya + (-1 + xa)*yb)*
               (p2x*(xa - xb - 2*ya + 2*yb) + p3x*((-1 + xb)*ya + yb - xa*yb)) + 
              p1y*(-(p3x*(pow(xb,2)*(-1 + ya)*ya - xb*(-1 + 3*ya)*(ya - yb) + 
                      2*pow(ya - yb,2) + xa*(1 + xb - 3*yb)*yb + 
                      pow(xa,2)*(-1 + yb)*yb + xa*ya*(-1 + xb + 3*yb - 2*xb*yb)))\
                  + p2x*(pow(xa,2)*(-1 + yb)*yb + 
                    xb*((-1 + xb)*pow(ya,2) + yb + ya*(-1 - xb + yb)) + 
                    xa*((-1 + xb - yb)*yb + ya*(1 + xb + yb - 2*xb*yb)))) + 
              p2y*(p3x*(pow(xb,2)*ya - xb*(1 + ya)*(ya - yb) + pow(ya - yb,2) + 
                    pow(xa,2)*yb - xa*(ya*(-1 + xb - yb) + yb*(1 + xb + yb))) + 
                 p2x*(pow(xa,2)*(-2 + yb)*yb + xb*ya*(xb*(-2 + ya) - ya + yb) + 
                    xa*((ya - yb)*yb + 2*xb*(ya + yb - ya*yb)))))) + 
        p0x*((-(p2y*xb*ya) + p2y*xa*yb + p3y*((-1 + xb)*ya + yb - xa*yb))*
            ((p2y*p3x - p2x*p3y)*((-1 + xb)*ya + yb - xa*yb) + 
              p1x*(p2y*(xb*(-2 + ya) - ya - xa*(-2 + yb) + yb) + 
                 p3y*((-1 + xb)*ya + yb - xa*yb))) + 
           pow(p1y,2)*(p3x*(pow(xb,2)*(-1 + ya)*ya - xb*(-1 + 3*ya)*(ya - yb) + 
                 2*pow(ya - yb,2) + xa*(1 + xb - 3*yb)*yb + 
                 pow(xa,2)*(-1 + yb)*yb + xa*ya*(-1 + xb + 3*yb - 2*xb*yb)) - 
              p2x*(pow(xa,2)*(-1 + yb)*yb + 
                 xb*((-1 + xb)*pow(ya,2) + yb + ya*(-1 - xb + yb)) + 
                 xa*((-1 + xb - yb)*yb + ya*(1 + xb + yb - 2*xb*yb)))) + 
           p1y*(p3y*(ya - xb*ya + (-1 + xa)*yb)*
               (p2x*(xa - xb + ya - yb) + p3x*((-1 + xb)*ya + yb - xa*yb)) + 
              p1x*(-(p3y*(pow(xb,2)*(-1 + ya)*ya - xb*(-1 + 3*ya)*(ya - yb) + 
                      2*pow(ya - yb,2) + xa*(1 + xb - 3*yb)*yb + 
                      pow(xa,2)*(-1 + yb)*yb + xa*ya*(-1 + xb + 3*yb - 2*xb*yb)))\
                  + p2y*(pow(xa,2)*(-1 + yb)*yb + 
                    xb*((-1 + xb)*pow(ya,2) + yb + ya*(-1 - xb + yb)) + 
                    xa*((-1 + xb - yb)*yb + ya*(1 + xb + yb - 2*xb*yb)))) + 
              p2y*(p3x*(pow(xb,2)*ya + xb*(-1 + 2*ya)*(ya - yb) - 
                    2*pow(ya - yb,2) + pow(xa,2)*yb - 
                    xa*((1 + xb - 2*yb)*yb + ya*(-1 + xb + 2*yb))) + 
                 p2x*(pow(xa,2)*(-2 + yb)*yb + xb*ya*(xb*(-2 + ya) - ya + yb) + 
                    xa*((ya - yb)*yb + 2*xb*(ya + yb - ya*yb)))))))*
      log(p1y*p3x - p1x*p3y + p1y*p2x*xb - p1x*p2y*xb - p1y*p3x*xb + p1x*p3y*xb + 
        p0x*(p3y - p3y*xb + p2y*(xb - yb) + p1y*(-1 + yb)) - p1y*p3x*yb + 
        p2y*p3x*yb + p1x*p3y*yb - p2x*p3y*yb + 
        p0y*(p1x + p3x*(-1 + xb) - p1x*yb + p2x*(-xb + yb))))/
   (2.*pow(-(p0x*p2y*xa) + p1x*p2y*xa + p0x*p3y*xa - p1x*p3y*xa + p0x*p2y*xb - 
       p1x*p2y*xb - p0x*p3y*xb + p1x*p3y*xb + p0x*p2y*ya - p2y*p3x*ya - p1x*p3y*ya + 
       p2x*p3y*ya - p0x*p2y*yb + p2y*p3x*yb + p1x*p3y*yb - p2x*p3y*yb + 
       p1y*(p2x*(-xa + xb) + p3x*(xa - xb + ya - yb) + p0x*(-ya + yb)) + 
       p0y*(p3x*(-xa + xb) + p1x*(ya - yb) + p2x*(xa - xb - ya + yb)),3));


    // flux along lower x edge
    this->projections[0] = ((xb - xa)*(p3x - p0x) + (yb - ya)*(p3y - p0y)) * integrals[0]
                         + ((xb - xa)*(p2x - p1x) + (yb - ya)*(p2y - p1y)) * integrals[1];

    // flux along upper x edge
    this->projections[1] = ((xb - xa)*(p3x - p0x) + (yb - ya)*(p3y - p0y)) * integrals[3]
                         + ((xb - xa)*(p2x - p1x) + (yb - ya)*(p2y - p1y)) * integrals[2];

    // flux along lower y edge
    this->projections[2] = ((xb - xa)*(p1x - p0x) + (yb - ya)*(p1y - p0y)) * integrals[0]
                         + ((xb - xa)*(p2x - p3x) + (yb - ya)*(p2y - p3y)) * integrals[3];

    // flux along upper y edge
    this->projections[2] = ((xb - xa)*(p1x - p0x) + (yb - ya)*(p1y - p0y)) * integrals[1]
                         + ((xb - xa)*(p2x - p3x) + (yb - ya)*(p2y - p3y)) * integrals[2];
  }


};
 
#ifdef __cplusplus
extern "C" {
#endif

int SgQuadLineFlows_new(SgQuadLineFlows_type** self);
                       
int SgQuadLineFlows_del(SgQuadLineFlows_type** self);

int SgQuadLineFlows_setQuadPoints(SgQuadLineFlows_type** self, 
                                      const double* quadPoints);

int SgQuadLineFlows_setLinePoints(SgQuadLineFlows_type** self, 
                                      const double* linePoints);

int SgQuadLineFlows_computeProjections(SgQuadLineFlows_type** self);

int SgQuadLineFlows_getProjection(SgQuadLineFlows_type** self, int edgeIndex, double* flux);

#ifdef __cplusplus
}
#endif

#endif // SG_QUAD_LINE_FLOWS_H
