/**
 * Testing octree partitioning
 */

#include <vector>
#include <cassert>
#include <cstdio>
#include <iostream>
#include <cmath>
#include "SgOctreePoints.h"

bool testCart(size_t numLevels, size_t nx, size_t ny) {

    // number of dimensions
    const size_t ndims = 2;

    // grid
    std::vector<double> coords;
    coords.reserve(nx * ny * ndims);
    for (size_t j = 0; j < ny; ++j) {
        double y = 0.0 + j * 1.0/(double) (ny - 1);
        for (size_t i = 0; i < nx; ++i) {
            double x = 0.0 + i * 1.0/(double) (nx - 1);
            coords.push_back(x);
            coords.push_back(y);
        }
    }

    SgOctreePoints_type octree(numLevels, ndims, coords);


    return true;

}


bool testPolar(size_t numLevels, size_t nr, size_t nt) {

    // number of dimensions
    const size_t ndims = 2;

    // grid
    std::vector<double> coords;
    coords.reserve(nr * nt * ndims);
    for (size_t i = 0; i < nr; ++i) {
        double r = 0. + (1. - 0.) * i / double(nr - 1);
        for (size_t j = 0; j < nt; ++j) {
            double t = 0. + (2*M_PI - 0.)* j / double(nt - 1);
            double x = r * cos(t);
            coords.push_back(x);
            double y = r * sin(t);
            coords.push_back(y);
        }
    }

    SgOctreePoints_type octree(numLevels, ndims, coords);


    return true;

}

int main(int argc, char** argv) {

    if (!testCart(0, 11, 11)) return 1;
    if (!testPolar(0, 5, 9)) return 1;

    return 0;
}
