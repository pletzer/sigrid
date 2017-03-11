/**
 * Testing computation of intersection points between quad and line
 */

#include <vector>
#include <cassert>
#include <cstdio>
#include <iostream>
#include "SgQuadLineIntersect.h"

bool testNoOverlap() {

    double quadCoords[] = {0., 0., 
                           1., 0.,
                           1., 1.,
                           0., 1.};

    double lineCoords[] = {1.001, 0.,
                           2., 0.};

    SgQuadLineIntersect_type* qis = NULL;
    SgQuadLineIntersect_new(&qis);
    SgQuadLineIntersect_setQuadPoints(&qis, quadCoords);
    SgQuadLineIntersect_setLinePoints(&qis, lineCoords);
    int numPoints;
    double* points = NULL;
    SgQuadLineIntersect_getIntersectPoints(&qis, &numPoints, &points);
    SgQuadLineIntersect_del(&qis);

    std::cout << "testNoOverlap: num intersection points = " << numPoints << '\n';

    if (numPoints != 0) {
        // error
        return false;
    }
    // OK
    return true;
}

bool testLineIsInsideQuad() {

    double quadCoords[] = {0., 0.,
                           1., 0.,
                           1., 1.,
                           0., 1.};

    double lineCoords[] = {0.3, 0.2,
                           0.8, 0.2};

    SgQuadLineIntersect_type* qis = NULL;
    SgQuadLineIntersect_new(&qis);
    SgQuadLineIntersect_setQuadPoints(&qis, quadCoords);
    SgQuadLineIntersect_setLinePoints(&qis, lineCoords);
    int numPoints;
    double* points = NULL;
    SgQuadLineIntersect_getIntersectPoints(&qis, &numPoints, &points);
    SgQuadLineIntersect_del(&qis);

    std::cout << "test testLineIsInsideQuad: num intersection points = " << numPoints << '\n';

    if (numPoints != 2) {
        // error
        return false;
    }
    // OK
    return true;
}

bool testOneIntersection() {

    double quadCoords[] = {0., 0.,
                           1., 0.,
                           1., 1.,
                           0., 1.};

    double lineCoords[] = {-0.2, 0.1,
                           0.8, 0.5};

    SgQuadLineIntersect_type* qis = NULL;
    SgQuadLineIntersect_new(&qis);
    SgQuadLineIntersect_setQuadPoints(&qis, quadCoords);
    SgQuadLineIntersect_setLinePoints(&qis, lineCoords);
    int numPoints;
    double* points = NULL;
    SgQuadLineIntersect_getIntersectPoints(&qis, &numPoints, &points);
    SgQuadLineIntersect_del(&qis);

    std::cout << "test testOneIntersection: num intersection points = " << numPoints << '\n';

    if (numPoints != 1) {
        // error
        return false;
    }
    // OK
    return true;
}

bool testTwoIntersections() {

    double quadCoords[] = {0., 0.,
                           1., 0.,
                           1., 1.,
                           0., 1.};

    double lineCoords[] = {-0.2, 0.1,
                           0.8, 1.1};

    SgQuadLineIntersect_type* qis = NULL;
    SgQuadLineIntersect_new(&qis);
    SgQuadLineIntersect_setQuadPoints(&qis, quadCoords);
    SgQuadLineIntersect_setLinePoints(&qis, lineCoords);
    int numPoints;
    double* points = NULL;
    SgQuadLineIntersect_getIntersectPoints(&qis, &numPoints, &points);
    SgQuadLineIntersect_del(&qis);

    std::cout << "test testTwoIntersections: num intersection points = " << numPoints << '\n';

    if (numPoints != 2) {
        // error
        return false;
    }
    // OK
    return true;
}


int main(int argc, char** argv) {

    if (!testNoOverlap()) return 1;
    if (!testLineIsInsideQuad()) return 2;
    if (!testOneIntersection()) return 3;
    if (!testTwoIntersections()) return 4;

    std::cout << "SUCCESS\n";
    return 0;
}
