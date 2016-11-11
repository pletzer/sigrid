/**
 * A class that triangulates a set of points
 */
 
#ifndef SG_TRIANGULATE_H
#define SG_TRIANGULATE_H
 
#include <vector>
#include <set>
#include <cstdio> // size_t
#include <cmath>
#include <limits>

/**
 * Functor used to order points 
 */
struct SgSortByDistanceSquareFunctor {

    std::vector<double> gravityCenter;

    const double** points;

    int numPoints;

    /**
     * Constructor
     * @param numPoints number of points 
     * @param points points as numPoints sized array of 2 coordinates
     */
    SgSortByDistanceSquareFunctor(int numPoints, const double** points) {
        this->numPoints = numPoints;
        this->points = points;
        const size_t NDIMS = 2;
        this->gravityCenter.resize(NDIMS, 0);
        for (int i = 0; i < numPoints; ++i) {
            for (size_t j = 0; j < NDIMS; ++j)
                this->gravityCenter[j] += points[i][j] / (double)(numPoints);
        }
    }

    double inline distanceSquareFromCenter(size_t i) {
        double res = 0;
        for (size_t j = 0; j < this->gravityCenter.size(); ++j) {
            double sq = this->points[i][j] - this->gravityCenter[j];
            res += sq * sq;
        }
        return res;
    }

    bool operator()(size_t i, size_t j) {
        double diSq = this->distanceSquareFromCenter(i);
        double djSq = this->distanceSquareFromCenter(j);
        // sort by increasing distance
        return diSq < djSq;
    }
};
 
struct SgTriangulate_type {

    // the vertices, as a flat array of 2D coordinates
 	std::vector<double> points;

    // set of 2-tuples nodes indices
    std::set< std::vector<size_t> > boundaryEdges;

    // set of 3-tuples
    std::set< std::vector<size_t> > triIndices;

 	// a tolerance used to determine whether a triangle area is positive
    // or negative
 	double eps;

    size_t NDIMS;

double getDistanceSquare(size_t i0, size_t i1) {
    double res = 0;
    for (size_t j = 0; j < this->NDIMS; ++j) {
        double d = this->points[i1*this->NDIMS + j] - this->points[i0*this->NDIMS + j];
        res += d * d;
    }
    return res;
}

void removeDegenerateSegments() {

    size_t numPoints = this->points.size() / this->NDIMS;

    // start with the first point
    std::vector<double> pts;
    for (size_t j = 0; j < this->NDIMS; ++j) {
        pts.push_back(this->points[0*this->NDIMS + j]);
    }

    // iterate over the remaining points, making sure the distance square between this 
    // and the previous point is > 0
    size_t i0 = 0;
    for (size_t i1 = 1; i1 < numPoints; ++i1) {
        if (this->getDistanceSquare(i0, i1) > this->eps) {
            // not degenerate
            for (size_t j = 0; j < this->NDIMS; ++j) {
                pts.push_back(this->points[i1*this->NDIMS + j]);
            }
            // new baseline
            i0 = i1;
        }
    }
    // copy
    this->points = pts;
}

void removeDegenerateTriangles() {

    size_t numPoints = this->points.size() / this->NDIMS;

    // start with the first 2 points
    std::vector<double> pts;
    for(size_t j = 0; j < this->NDIMS; ++j) {
        pts.push_back(this->points[0*this->NDIMS + j]);
    }
    for (size_t j = 0; j < this->NDIMS; ++j) {
        pts.push_back(this->points[1*this->NDIMS + j]);
    }

    // iterate over the remaining points, making sure the area of the triangle
    // between this and the previous two points is != 0
    size_t i0 = 0;
    size_t i1 = 1;
    for (size_t i2 = 2; i2 < numPoints; ++i2) {
        double area = this->getParallelipipedArea(i0, i1, i2);
        if (std::fabs(area) > this->eps) {
            // not degenerate
            for (size_t j = 0; j < this->NDIMS; ++j) {
                pts.push_back(this->points[i2*this->NDIMS + j]);
            }
            i0 = i1;
            i1 = i2;
        }
    }
    // copy
    this->points = pts;
}

double inline getParallelipipedArea(size_t ip0, size_t ip1, size_t ip2) {
    double d1[] = {0, 0};
    double d2[] = {0, 0};
    for (size_t j = 0; j < this->NDIMS; ++j) {
        d1[j] = this->points[ip1*this->NDIMS + j] - this->points[ip0*this->NDIMS + j];
        d2[j] = this->points[ip2*this->NDIMS + j] - this->points[ip0*this->NDIMS + j];
    }
    return (d1[0]*d2[1] - d1[1]*d2[0]);
}

void inline makeCounterClockwise(size_t ips[]) {
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
    if (area < -this->eps) {
        return true;
    }
    return false;
}

void addPoint(size_t ip) {

    std::set< std::vector<size_t> > boundaryEdgesToRemove;
    std::set< std::vector<size_t> > boundaryEdgesToAdd;

    // iterate over boundary edges
    for (std::set< std::vector<size_t> >::const_iterator it = this->boundaryEdges.begin();
         it != this->boundaryEdges.end(); ++it) {
        size_t iea = (*it)[0];
        size_t ieb = (*it)[1];
        if (this->isEdgeVisible(ip, iea, ieb)) {

            // create new triangle
            size_t tri[] = {iea, ip, ieb};
            this->triIndices.insert(std::vector<size_t>(tri, tri + 3));

            // keep track of edges to remove
            std::vector<size_t> oldEdge(&(*it)[0], &(*it)[2]);
            boundaryEdgesToRemove.insert(oldEdge);

            // keep track of edges to add
            std::vector<size_t> newEdge(2);
            newEdge[0] = iea;
            newEdge[1] = ip;
            boundaryEdgesToAdd.insert(newEdge);
            newEdge[0] = ip;
            newEdge[1] = ieb;
            boundaryEdgesToAdd.insert(newEdge);
        }
    }

    // remove the old boundary edges
    for (std::set< std::vector<size_t> >::const_iterator it = boundaryEdgesToRemove.begin();
         it != boundaryEdgesToRemove.end(); ++it) {
        this->boundaryEdges.erase(*it);
    }

    // add the new boundary edges
    for (std::set< std::vector<size_t> >::const_iterator it = boundaryEdgesToAdd.begin();
         it != boundaryEdgesToAdd.end(); ++it) {
        this->boundaryEdges.insert(*it);
    }


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