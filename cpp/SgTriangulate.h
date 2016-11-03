/**
 * A class that triangulates a set of points
 */
 
#ifndef SG_TRIANGULATE_H
#define SG_TRIANGULATE_H
 
#include <vector>
#include <cstdio> // size_t
#include <cmath>
#include <limits>

struct SgSortByDistanceSquareFunctor {
    std::vector<double> gravityCenter;
    const double** points;
    int numPoints;
    SgSortByDistanceSquareFunctor(int numPoints, const double** points) {
        this->numPoints = numPoints;
        this->points = points;
        this->gravityCenter.resize(2, 0);
        for (int i = 0; i < numPoints; ++i) {
            this->gravityCenter[0] += points[i][0];
            this->gravityCenter[1] += points[i][1];
        }
        this->gravityCenter[0] /= (double)(numPoints);
        this->gravityCenter[0] /= (double)(numPoints);
    }

    bool operator()(size_t i, size_t j) {
        double cg0 = this->gravityCenter[0];
        double cg1 = this->gravityCenter[1];
        double pi0 = this->points[i][0];
        double pi1 = this->points[i][0];
        double pj0 = this->points[j][0];
        double pj1 = this->points[j][1];
        double diSq = (pi0 - cg0)*(pi0 - cg0) + (pi1 - cg1)*(pi1 - cg1);
        double djSq = (pj0 - cg0)*(pj0 - cg0) + (pj1 - cg1)*(pj1 - cg1);
        // sort by increasing distance
        return diSq < djSq;
    }
};


 
struct SgTriangulate_type {

    // the vertices, as a flat array of 2D coordinates
 	std::vector<double> points;

    // sorted indices
    std::vector<size_t> sortedInds;

    // boundary edges as a flat array of two point indices
    std::vector<size_t> boundaryEdges;

    // flat array of triangle point indices
    std::vector<size_t> triIndices;

 	// a tolerance used to determine whether a triangle area is positive
    // or negative
 	double eps;

double inline getParallelipipedArea(size_t ip0, size_t ip1, size_t ip2) {
    double d1[] = {this->points[ip1*2 + 0] - this->points[ip0*2 + 0],
                   this->points[ip1*2 + 1] - this->points[ip1*2 + 1]};
    double d2[] = {this->points[ip2*2 + 0] - this->points[ip0*2 + 0],
                   this->points[ip2*2 + 1] - this->points[ip0*2 + 1]};
    return (d1[0]*d2[1] - d1[1]*d2[0]);
}

void inline makeCounterClockwise(int ips[]) {
    double area = this->getParallelipipedArea(ips[0], ips[1], ips[2]);
    if (area < -this->eps) {
      size_t ip1 = ips[1];
      size_t ip2 = ips[2];
      // swap
      ips[1] = ip2;
      ips[2] = ip1;
    }
}

bool inline isEdgeVisible(size_t ip, size_t ie0, size_t ie1) {
    double area = this->getParallelipipedArea(ip, ie0, ie1);
    if (area < this->eps) {
        return true;
    }
    return false;
}

};
 
#ifdef __cplusplus
extern "C" {
#endif

 int SgTriangulate_new(SgTriangulate_type** self,
                       int numPoints, const double** points);

 int SgTriangulate_getConvexHullArea(SgTriangulate_type** self, double* area);
                       
 int SgTriangulate_del(SgTriangulate_type** self);

#ifdef __cplusplus
}
#endif

#endif // SG_QUAD_INTERSECT_H