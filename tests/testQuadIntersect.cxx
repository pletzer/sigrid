/**
 * Testing computation of intersection points between two quads
 */

#include <vector>
#include <cassert>
#include <cstdio>
#include <iostream>
#include "SgQuadIntersect.h"

bool testNoOverlap() {

    double quad1Coords[] = {0., 0., 
                            1., 0.,
                            1., 1.,
                            0., 1.};

    double quad2Coords[] = {1.001, 0.,
                            2., 0.,
                            2., 1.,
                            1.001, 1.};

    SgQuadIntersect_type* qis = NULL;
    SgQuadIntersect_new(&qis);
    SgQuadIntersect_setQuadPoints(&qis, quad1Coords, quad2Coords);
    int numPoints;
    double* points = NULL;
    SgQuadIntersect_getIntersectPoints(&qis, &numPoints, &points);
    SgQuadIntersect_del(&qis);

    std::cout << "testNoOverlap: num intersection points = " << numPoints << '\n';

    if (numPoints != 0) {
        // error
        return false;
    }
    // OK
    return true;
}

bool testQuad2IsInsideQuad1() {

    double quad1Coords[] = {0., 0.,
                            1., 0.,
                            1., 1.,
                            0., 1.};

    double quad2Coords[] = {0.3, 0.2,
                            0.8, 0.2,
                            0.7, 0.6,
                            0.2, 0.5};

    SgQuadIntersect_type* qis = NULL;
    SgQuadIntersect_new(&qis);
    SgQuadIntersect_setQuadPoints(&qis, quad1Coords, quad2Coords);
    int numPoints;
    double* points = NULL;
    SgQuadIntersect_getIntersectPoints(&qis, &numPoints, &points);
    SgQuadIntersect_del(&qis);

    std::cout << "testQuad2IsInsideQuad1: num intersection points = " << numPoints << '\n';

    if (numPoints != 4) {
        // error
        return false;
    }
    // OK
    return true;
}

bool testQuad1IsInsideQuad2() {

    double quad2Coords[] = {0., 0.,
                            1., 0.,
                            1., 1.,
                            0., 1.};

    double quad1Coords[] = {0.3, 0.2,
                            0.8, 0.2,
                            0.7, 0.6,
                            0.2, 0.5};

    SgQuadIntersect_type* qis = NULL;
    SgQuadIntersect_new(&qis);
    SgQuadIntersect_setQuadPoints(&qis, quad1Coords, quad2Coords);
    int numPoints;
    double* points = NULL;
    SgQuadIntersect_getIntersectPoints(&qis, &numPoints, &points);
    SgQuadIntersect_del(&qis);

    std::cout << "testQuad1IsInsideQuad2: num intersection points = " << numPoints << '\n';

    if (numPoints != 4) {
        // error
        return false;
    }
    // OK
    return true;
}

bool testPartial3Points() {

    double quad2Coords[] = {0., 0.,
                            1., 0.,
                            1., 1.,
                            0., 1.};

    double quad1Coords[] = {1.1, 0.4,
                            2.0, 0.4,
                            1.3, 1.1,
                            0.6, 1.1};

    SgQuadIntersect_type* qis = NULL;
    SgQuadIntersect_new(&qis);
    SgQuadIntersect_setQuadPoints(&qis, quad1Coords, quad2Coords);
    int numPoints;
    double* points = NULL;
    SgQuadIntersect_getIntersectPoints(&qis, &numPoints, &points);
    SgQuadIntersect_del(&qis);

    std::cout << "testPartial3Points: num intersection points = " << numPoints << '\n';

    if (numPoints != 3) {
        // error
        return false;
    }
    // OK
    return true;
}

int main(int argc, char** argv) {

    if (!testNoOverlap()) return 1;
    if (!testQuad2IsInsideQuad1()) return 2;
    if (!testQuad1IsInsideQuad2()) return 3;
    if (!testPartial3Points()) return 4;

    std::cout << "SUCCESS\n";
    return 0;
}
