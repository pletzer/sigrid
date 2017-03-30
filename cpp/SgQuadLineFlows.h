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
#include "SgFindPointInCell.h"

const int EDGE_LO_X = 0;
const int EDGE_HI_X = 1;
const int EDGE_LO_Y = 2;
const int EDGE_HI_Y = 3;

struct SgQuadLineFlows_type {

    // flat array (node, component) of the nodal coordinates
    double* quadCoords;
    double* lineCoords;

    double projections[4]; // 4 edges/basis functions

    // to find the indices given a target position
    SgFindPointInCell_type* indexFinder;

    // store the coordinates in struct grid format
    std::vector<double*> quadGridCoords;

    // the start/end positions of the line in index space
    std::vector<double> xia;
    std::vector<double> xib;

    /**
     * Constructor
     */
     SgQuadLineFlows_type() {

        this->quadCoords = NULL;
        this->lineCoords = NULL;
        this->quadGridCoords.resize(2);
        this->xia.resize(2);
        this->xib.resize(2);
        this->quadGridCoords[0] = new double[2]; // x
        this->quadGridCoords[1] = new double[2]; // y
        int nitermax = 10;
        double tolpos = 1.e-10;
        this->indexFinder = new SgFindPointInCell_type(nitermax, tolpos);
    }

    /**
     * Destructor
     */
     ~SgQuadLineFlows_type() {
        delete[] this->quadGridCoords[0];
        delete[] this->quadGridCoords[1];
        delete this->indexFinder;
    }

    /**
     * Set the quad's coordinates
     * @param quadCoords quad coordinates
     */
    void setQuadPoints(const double* quadCoords) {

      this->quadCoords = (double*)quadCoords;

      // set the grid
      this->quadGridCoords[0][0] = quadCoords[0];
      this->quadGridCoords[1][0] = quadCoords[1];

      this->quadGridCoords[0][1] = quadCoords[2];
      this->quadGridCoords[1][1] = quadCoords[3];

      this->quadGridCoords[0][2] = quadCoords[6];
      this->quadGridCoords[1][2] = quadCoords[7];

      this->quadGridCoords[0][3] = quadCoords[4];
      this->quadGridCoords[1][3] = quadCoords[5];

      const int ndims = 2;
      const int dims[] = {2, 2};
      const int periodicity[] = {0, 0};
      this->indexFinder->setGrid(ndims, dims, periodicity, (const double**) &this->quadGridCoords[0]);
    }

    /**
     * Set the line's coordinates
     * @param lineCoords line coordinates
     */
    void setLinePoints(const double* lineCoords) {

      int end;

      this->lineCoords = (double*) lineCoords;
      double initInds[] = {0.5, 0.5};

      // search xia
      this->indexFinder->reset(initInds, &lineCoords[0]);
      end = 0;
      while (end == 0) {
        end = this->indexFinder->next();
      }
      this->xia = this->indexFinder->getIndices();

      // search xib
      this->indexFinder->reset(initInds, &lineCoords[2]);
      end = 0;
      while (end == 0) {
        end = this->indexFinder->next();
      }
      this->xib = this->indexFinder->getIndices();

      // need to check error
      // HERE...
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

      double xa = this->xia[0];
      double ya = this->xia[1];
      double xb = this->xib[0];
      double yb = this->xib[1];

      double p0x = 0.;
      double p0y = 0.;

      double p1x = 1.;
      double p1y = 0.;

      double p2x = 1.;
      double p2y = 1.;

      double p3x = 0.;
      double p3y = 1.;

      double integrals[4];

      integrals[0] = (-3*(-2 + xb + ya + yb) + xa*(-3 + 2*ya + yb) + xb*(ya + 2*yb))/6.;

      integrals[1] = (-(xa*(-3 + 2*ya + yb)) - xb*(-3 + ya + 2*yb))/6.;

      integrals[2] = (xa*(2*ya + yb) + xb*(ya + 2*yb))/6.;

      integrals[3] = (-((-3 + 2*xa + xb)*ya) - (-3 + xa + 2*xb)*yb)/6.;

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
