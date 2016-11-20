/**
 * A class that triangulates a set of points
 */
 
#ifndef SG_TRIANGULATE_H
#define SG_TRIANGULATE_H

#include "SgNdims.h"
 
#include <vector>
#include <algorithm>
#include <set>
#include <cstdio> // size_t
#include <cmath>
#include <limits>

/**
 * Functor used to order points 
 */
struct SgSortByDistanceSquareFunctor {

    std::vector<double> gravityCenter;

    double* points;

    int numPoints;

    /**
     * Constructor
     * @param numPoints number of points 
     * @param points points as numPoints sized array of 2 coordinates
     */
    SgSortByDistanceSquareFunctor(int numPoints, const double points[]) {
        this->numPoints = numPoints;
        this->points = (double*) points;
        this->gravityCenter.resize(NDIMS_2D_PHYS, 0);
        for (int i = 0; i < numPoints; ++i) {
            for (size_t j = 0; j < NDIMS_2D_PHYS; ++j)
                this->gravityCenter[j] += points[NDIMS_2D_PHYS*i + j] / (double)(numPoints);
        }
    }

    double inline distanceSquareFromCenter(size_t i) {
        double res = 0;
        for (size_t j = 0; j < this->gravityCenter.size(); ++j) {
            double sq = this->points[NDIMS_2D_PHYS*i + j] - this->gravityCenter[j];
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


    SgTriangulate_type(int numPoints, const double points[]) {

        // tolerance for floating point comparisons
        this->eps = 1.e-12;

        this->points.resize(2 * numPoints);

        SgSortByDistanceSquareFunctor sortFunc(numPoints, points);
 
        // set the indices before the sort
        std::vector<size_t> sortedInds(numPoints);
        for (int i = 0; i < numPoints; ++i) {
            sortedInds[i] = i;
        }
        // sort the point indices by increasing distance from the centre of gravity
        std::sort(sortedInds.begin(), sortedInds.end(), sortFunc);

        // set the points in increasing distance from the centre of gravity
        this->points.resize(2 * numPoints);
        for (size_t i = 0; i < sortedInds.size(); ++i) {
            size_t indx = sortedInds[i];
            for (size_t j = 0; j < NDIMS_2D_PHYS; ++j) {
                this->points[NDIMS_2D_PHYS*i + j] = points[NDIMS_2D_PHYS*indx + j];
            }        
        }

        // remove degenerate points
        this->removeDegenerateSegments();

        // create the first, non-degenerate triangle
        size_t lastIndexOfFirstTriangle = this->makeFirstTriangle();

        // add the remaining points
        for (int ip = lastIndexOfFirstTriangle + 1; ip < numPoints; ++ip) {
            this->addPoint(ip);
        }
    }

    double getDistanceSquare(size_t i0, size_t i1) {
        double res = 0;
        for (size_t j = 0; j < NDIMS_2D_PHYS; ++j) {
            double d = this->points[i1*NDIMS_2D_PHYS + j];
            d -= this->points[i0*NDIMS_2D_PHYS + j];
            res += d * d;
        }
        return res;
    }

    void removeDegenerateSegments() {

        size_t numPoints = this->points.size() / NDIMS_2D_PHYS;

        if (numPoints < 1) {
            // nothing to do 
            return;
        }

        // start with the first point
        std::vector<double> pts;
        for (size_t j = 0; j < NDIMS_2D_PHYS; ++j) {
            pts.push_back(this->points[0*NDIMS_2D_PHYS + j]);
        }

        // iterate over the remaining points, making sure the distance square between this 
        // and the previous point is > 0
        size_t i0 = 0;
        for (size_t i1 = 1; i1 < numPoints; ++i1) {
            if (this->getDistanceSquare(i0, i1) > this->eps) {
                // not degenerate
                for (size_t j = 0; j < NDIMS_2D_PHYS; ++j) {
                    pts.push_back(this->points[i1*NDIMS_2D_PHYS + j]);
                }
                // new base point
                i0 = i1;
            }
        }
        // copy
        this->points = pts;
    }

    std::vector<size_t> getExtremaPointIndices(size_t i0, size_t i1, size_t i2) {
        // choose each point as the base point, compute the dot product of the 
        // other vectors minus this point. If the dot product is negative then 
        // the point is a middle point. Return the other two points
        double da[] = {0, 0};
        double db[] = {0, 0};
        std::vector<size_t> inds(3);
        inds[0] = i0; inds[1] = i1; inds[2] = i2;
        double dotProduct = 1;
        size_t ia, ib, iBase;
        size_t j = 0;
        while (dotProduct > this->eps && j < 3) {
            iBase = inds[j];
            ia = inds[(j + 1) % 3];
            ib = inds[(j + 2) % 3];
            dotProduct = 0;
            for (size_t k = 0; k < NDIMS_2D_PHYS; ++k) {
                double da = this->points[ia*NDIMS_2D_PHYS + k];
                da -= this->points[iBase*NDIMS_2D_PHYS + k];
                double db = this->points[ib*NDIMS_2D_PHYS + k];
                db -= this->points[iBase*NDIMS_2D_PHYS + k];
                dotProduct += da*db;
            }
            j++;
        }
        std::vector<size_t> res(2);
        res[0] = ia;
        res[1] = ib;
        return res;
    }

    size_t makeFirstTriangle() {

        size_t numPoints = this->points.size() / NDIMS_2D_PHYS;

        if (numPoints < 3) {
            return numPoints - 1;
        }

        // the first triangle cannot be degenerate
        size_t i0 = 0;
        size_t i1 = 1;
        size_t i2 = 2;
        double area = std::abs(this->getParallelipipedArea(i0, i1, i2));
        while (i2 < numPoints && area < this->eps) {
            // the points appear to be a line, get the two extremum point indices
            std::vector<size_t> extrema = this->getExtremaPointIndices(i0, i1, i2);
            i0 = extrema[0];
            i1 = extrema[1];
            // try creating a triangle using the next available point
            i2++;
            area = std::abs(this->getParallelipipedArea(i0, i1, i2));
        }
        size_t lastPointIndexOfFirstTriangle = i2;

        if (area < this->eps) {
            // all the points are on a line, cannot create triangle
            return lastPointIndexOfFirstTriangle;
        }

        // got a non-degenerate triangle
        // re-order the indices so the area is positive
        std::vector<size_t> inds(3);
        inds[0] = i0; inds[1] = i1; inds[2] = i2;
        this->makeCounterClockwise(&inds[0]);
        i0 = inds[0]; i1 = inds[1]; i2 = inds[2];

        // create the triangle
        this->triIndices.insert(inds);

        // create the boundary edges associated with this triangle
        size_t edge[] = {i0, i1};
        this->boundaryEdges.insert(std::vector<size_t>(edge, edge + 2));
        edge[0] = i1; edge[1] = i2;
        this->boundaryEdges.insert(std::vector<size_t>(edge, edge + 2));
        edge[0] = i2; edge[1] = i0;
        this->boundaryEdges.insert(std::vector<size_t>(edge, edge + 2));

        // return the first index that yields a non-degenerate triangle
        // so we can add the remaining points starting from this index
        return lastPointIndexOfFirstTriangle;
    }

    double inline getParallelipipedArea(size_t ip0, size_t ip1, size_t ip2) const {
        double d1[] = {0, 0};
        double d2[] = {0, 0};
        for (size_t j = 0; j < NDIMS_2D_PHYS; ++j) {
            const double& p0 = this->points[ip0*NDIMS_2D_PHYS + j];
            d1[j] = this->points[ip1*NDIMS_2D_PHYS + j] - p0;
            d2[j] = this->points[ip2*NDIMS_2D_PHYS + j] - p0;
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

    double getConvexHullArea() const {
        double area = 0;
        for (std::set< std::vector<size_t> >::const_iterator it = this->triIndices.begin();
             it != this->triIndices.end(); ++it) {
            area += this->getParallelipipedArea((*it)[0], (*it)[1], (*it)[2]);
        }
        area *= 0.5;
        return area;
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
        for (std::set< std::vector<size_t> >::const_iterator 
             it = boundaryEdgesToAdd.begin();
             it != boundaryEdgesToAdd.end(); ++it) {
            this->boundaryEdges.insert(*it);
        }
    }
};
 
#ifdef __cplusplus
extern "C" {
#endif

 int SgTriangulate_new(SgTriangulate_type** self,
                       int numPoints, const double points[]);

 int SgTriangulate_getConvexHullArea(SgTriangulate_type** self, double* area);
                       
 int SgTriangulate_del(SgTriangulate_type** self);

#ifdef __cplusplus
}
#endif

#endif // SG_QUAD_INTERSECT_H
