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
    SgQuadLineIntersect_collectIntersectPoints(&qis, &numPoints, &points);
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
    SgQuadLineIntersect_collectIntersectPoints(&qis, &numPoints, &points);
    SgQuadLineIntersect_del(&qis);

    std::cout << "testLineIsInsideQuad: num intersection points = " << numPoints << '\n';

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
    SgQuadLineIntersect_collectIntersectPoints(&qis, &numPoints, &points);
    SgQuadLineIntersect_del(&qis);

    std::cout << "testOneIntersection: num intersection points = " << numPoints << '\n';

    // two points, including the end point of the line that is inside the quad
    if (numPoints != 2) {
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
    SgQuadLineIntersect_collectIntersectPoints(&qis, &numPoints, &points);
    SgQuadLineIntersect_del(&qis);

    std::cout << "testTwoIntersections: num intersection points = " << numPoints << '\n';

    if (numPoints != 2) {
        // error
        return false;
    }
    // OK
    return true;
}

bool testParallelSegment() {

    double quadCoords[] = {0., 0.,
                           1., 0.,
                           1., 1.,
                           0., 1.};

    double lineCoords[] = {-0.2, 0.0,
                           -0.2, 1.1};

    SgQuadLineIntersect_type* qis = NULL;
    SgQuadLineIntersect_new(&qis);
    SgQuadLineIntersect_setQuadPoints(&qis, quadCoords);
    SgQuadLineIntersect_setLinePoints(&qis, lineCoords);
    int numPoints;
    double* points = NULL;
    SgQuadLineIntersect_collectIntersectPoints(&qis, &numPoints, &points);
    SgQuadLineIntersect_del(&qis);

    std::cout << "testParallelSegment: num intersection points = " << numPoints << '\n';

    if (numPoints != 0) {
        // error
        return false;
    }
    // OK
    return true;
}

bool testLineOnEdge() {

    double quadCoords[] = {0., 0.,
                           1., 0.,
                           1., 1.,
                           0., 1.};

    double lineCoords[] = {0.0, 0.0,
                           0.0, 1.0};

    SgQuadLineIntersect_type* qis = NULL;
    SgQuadLineIntersect_new(&qis);
    SgQuadLineIntersect_setQuadPoints(&qis, quadCoords);
    SgQuadLineIntersect_setLinePoints(&qis, lineCoords);
    int numPoints;
    double* points = NULL;
    SgQuadLineIntersect_collectIntersectPoints(&qis, &numPoints, &points);
    SgQuadLineIntersect_del(&qis);

    std::cout << "testLineOnEdge: num intersection points = " << numPoints << '\n';

    if (numPoints != 2) {
        // error
        return false;
    }

    for (int i = 0; i < numPoints; ++i) {
        std::cout << "point " << i << ": ";
        for (size_t j = 0; j < 2; ++j) {
            std::cout << points[2*i + j] << ", ";
        }
        std::cout << '\n';
    }
    // OK
    return true;
}

bool testSmallLineOnEdge() {

    double quadCoords[] = {0., 0.,
                           1., 0.,
                           1., 1.,
                           0., 1.};

    double lineCoords[] = {0.0, 0.4,
                           0.0, 0.6};

    SgQuadLineIntersect_type* qis = NULL;
    SgQuadLineIntersect_new(&qis);
    SgQuadLineIntersect_setQuadPoints(&qis, quadCoords);
    SgQuadLineIntersect_setLinePoints(&qis, lineCoords);
    int numPoints;
    double* points = NULL;
    SgQuadLineIntersect_collectIntersectPoints(&qis, &numPoints, &points);
    SgQuadLineIntersect_del(&qis);

    std::cout << "testSmallLineOnEdge: num intersection points = " << numPoints << '\n';

    if (numPoints != 2) {
        // error
        return false;
    }

    for (int i = 0; i < numPoints; ++i) {
        std::cout << "point " << i << ": ";
        for (size_t j = 0; j < 2; ++j) {
            std::cout << points[2*i + j] << ", ";
        }
        std::cout << '\n';
    }
    // OK
    return true;
}

bool testAcrossQuad() {

    double quadCoords[] = {0., 0.,
                           0.5, 0.,
                           0.5, 0.5,
                           0., 0.5};

    double lineCoords[] = {0.0, 0.0,
                           1.0, 1.0};

    SgQuadLineIntersect_type* qis = NULL;
    SgQuadLineIntersect_new(&qis);
    SgQuadLineIntersect_setQuadPoints(&qis, quadCoords);
    SgQuadLineIntersect_setLinePoints(&qis, lineCoords);
    int numPoints;
    double* points = NULL;
    SgQuadLineIntersect_collectIntersectPoints(&qis, &numPoints, &points);
    SgQuadLineIntersect_del(&qis);

    std::cout << "testAcrossQuad: num intersection points = " << numPoints << '\n';

    if (numPoints != 2) {
        // error
        return false;
    }

    for (int i = 0; i < numPoints; ++i) {
        std::cout << "point " << i << ": ";
        for (size_t j = 0; j < 2; ++j) {
            std::cout << points[2*i + j] << ", ";
        }
        std::cout << '\n';
    }
    // OK
    return true;
}

int main(int argc, char** argv) {

    if (!testNoOverlap()) return 1;
    if (!testLineIsInsideQuad()) return 2;
    if (!testOneIntersection()) return 3;
    if (!testTwoIntersections()) return 4;
    if (!testParallelSegment()) return 5;
    if (!testLineOnEdge()) return 6;
    if (!testSmallLineOnEdge()) return 7;
    if (!testAcrossQuad()) return 8;

    std::cout << "SUCCESS\n";
    return 0;
}
