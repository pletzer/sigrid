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

const int EDGE_LO_0 = 0;
const int EDGE_HI_0 = 1;
const int EDGE_LO_1 = 2;
const int EDGE_HI_1 = 3;

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
        this->quadGridCoords[0] = new double[4]; // x, 4 nodes
        this->quadGridCoords[1] = new double[4]; // y
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

      // set the grid. Note that we're getting the coordinates in counterclockwise numbering 
      // while the grid assumes structured node ordering, the first index varying faster
      for (size_t k = 0; k < NDIMS_2D_PHYS; ++k) {
        this->quadGridCoords[k][0] = quadCoords[0*NDIMS_2D_PHYS + k];
        this->quadGridCoords[k][1] = quadCoords[3*NDIMS_2D_PHYS + k];
        this->quadGridCoords[k][2] = quadCoords[1*NDIMS_2D_PHYS + k];
        this->quadGridCoords[k][3] = quadCoords[2*NDIMS_2D_PHYS + k];
      }
      const int dims[] = {2, 2};
      const int periodicity[] = {0, 0};
      this->indexFinder->setGrid(NDIMS_2D_TOPO, dims, periodicity,
                                (const double**) &this->quadGridCoords[0]);
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
    int computeProjections() {

      int res = this->computeStartEndLineIndices();

      double x0a = this->xia[0];
      double x1a = this->xia[1];
      double x0b = this->xib[0];
      double x1b = this->xib[1];

      double x0abMid = 0.5*(x0a + x0b);
      double x1abMid = 0.5*(x1a + x1b);
      double x0abDiff = x0b - x0a;
      double x1abDiff = x1b - x1a;

      this->projections[EDGE_LO_0] = x0abDiff*(1. - x1abMid);
      this->projections[EDGE_HI_0] = x0abDiff*x1abMid;
      this->projections[EDGE_LO_1] = x1abDiff*(1. - x0abMid);
      this->projections[EDGE_HI_1] = x1abDiff*x0abMid;

      return res;
  }

private:

    /**
     * Compute the start/end indices of the line
     * @return error (0 = OK, 1 lower point failed, 2 upper point failed, < 0 other error)
     */
     int computeStartEndLineIndices() {

      if (! this->lineCoords) {
        // must call setLinePoints before calling this method
        return -1;
      }

      int end, res;

      // initial guess indices
      res = 0;
      double initInds[] = {0.5, 0.5};

      // search index position of the start point
      this->indexFinder->reset(initInds, &lineCoords[0]);
      end = 0;
      while (end == 0) {
        end = this->indexFinder->next();
      }
      this->xia = this->indexFinder->getIndices();
      if (end != 1) {
        // not converged
        res = 1;
      }

      // search index position of the end point
      this->indexFinder->reset(initInds, &lineCoords[2]);
      end = 0;
      while (end == 0) {
        end = this->indexFinder->next();
      }
      this->xib = this->indexFinder->getIndices();
      if (end != 1) {
        // not converged
        res = 2;
      }

      return res;
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
