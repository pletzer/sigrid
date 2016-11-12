/**
 * Testing computation of intersection points between two quads
 */

 #include <vector>
 #include <cassert>
 #include <cstdio>
 #include <iostream>
 #include "SgQuadIntersect.h"
 #include "SgTriangulate.h"

bool testNoOverlap() {

    double quad1Node0[] = {0., 0.};
    double quad1Node1[] = {1., 0.};
    double quad1Node2[] = {1., 1.};
    double quad1Node3[] = {0., 1.};

    double quad2Node0[] = {1.001, 0.};
    double quad2Node1[] = {2., 0.};
    double quad2Node2[] = {2., 1.};
    double quad2Node3[] = {1.001, 1.};

    const double* quad1Coords[] = {quad1Node0, quad1Node1, quad1Node2, quad1Node3};
    const double* quad2Coords[] = {quad2Node0, quad2Node1, quad2Node2, quad2Node3};

    SgQuadIntersect_type* qis = NULL;
    SgQuadIntersect_new(&qis);
    SgQuadIntersect_setQuadPoints(&qis, quad1Coords, quad2Coords);
    int numPoints;
    double** points;
    SgQuadIntersect_getIntersectPoints(&qis, &numPoints, points);

    SgTriangulate_type* tri = NULL;
    SgTriangulate_new(&tri, numPoints, (const double**) points);
    double area;
    SgTriangulate_getConvexHullArea(&tri, &area);

    SgTriangulate_del(&tri);
    SgQuadIntersect_del(&qis);

    std::cout << "testNoOverlap: area intersection = " << area << '\n';
    assert(fabs(area) < 1.e-10);

    return true;
}

bool testQuad2IsInsideQuad1() {

    double quad1Node0[] = {0., 0.};
    double quad1Node1[] = {1., 0.};
    double quad1Node2[] = {1., 1.};
    double quad1Node3[] = {0., 1.};

    double quad2Node0[] = {0.3, 0.2};
    double quad2Node1[] = {0.8, 0.2};
    double quad2Node2[] = {0.7, 0.6};
    double quad2Node3[] = {0.2, 0.5};

    const double* quad1Coords[] = {quad1Node0, quad1Node1, quad1Node2, quad1Node3};
    const double* quad2Coords[] = {quad2Node0, quad2Node1, quad2Node2, quad2Node3};

    SgQuadIntersect_type* qis = NULL;
    SgQuadIntersect_new(&qis);
    SgQuadIntersect_setQuadPoints(&qis, quad1Coords, quad2Coords);
    int numPoints;
    double* points = NULL;
    SgQuadIntersect_getIntersectPoints(&qis, &numPoints, &points);
    SgQuadIntersect_del(&qis);

    std::cout << "testQuad2IsInsideQuad1: num intersection points = " << numPoints << '\n';

    if (numPoints != 4) {
         /// error
        return false;
    }
    // OK
    return true;
}

bool testQuad1IsInsideQuad2() {

     double quad2Node0[] = {0., 0.};
     double quad2Node1[] = {1., 0.};
     double quad2Node2[] = {1., 1.};
     double quad2Node3[] = {0., 1.};

     double quad1Node0[] = {0.3, 0.2};
     double quad1Node1[] = {0.8, 0.2};
     double quad1Node2[] = {0.7, 0.6};
     double quad1Node3[] = {0.2, 0.5};

     const double* quad1Coords[] = {quad1Node0, quad1Node1, quad1Node2, quad1Node3};
     const double* quad2Coords[] = {quad2Node0, quad2Node1, quad2Node2, quad2Node3};

     SgQuadIntersect_type* qis = NULL;
     SgQuadIntersect_new(&qis);
     SgQuadIntersect_setQuadPoints(&qis, quad1Coords, quad2Coords);
     int numPoints;
     double* points = NULL;
     SgQuadIntersect_getIntersectPoints(&qis, &numPoints, &points);
     SgQuadIntersect_del(&qis);

     std::cout << "testQuad1IsInsideQuad2: num intersection points = " << numPoints << '\n';

     if (numPoints != 4) {
         /// error
         return false;
     }
     // OK
     return true;
}

bool testPartial3Points() {

    double quad2Node0[] = {0., 0.};
    double quad2Node1[] = {1., 0.};
    double quad2Node2[] = {1., 1.};
    double quad2Node3[] = {0., 1.};

    double quad1Node0[] = {1.1, 0.4};
    double quad1Node1[] = {2.0, 0.4};
    double quad1Node2[] = {1.3, 1.1};
    double quad1Node3[] = {0.6, 1.1};

    const double* quad1Coords[] = {quad1Node0, quad1Node1, quad1Node2, quad1Node3};
    const double* quad2Coords[] = {quad2Node0, quad2Node1, quad2Node2, quad2Node3};

    SgQuadIntersect_type* qis = NULL;
    SgQuadIntersect_new(&qis);
    SgQuadIntersect_setQuadPoints(&qis, quad1Coords, quad2Coords);
    int numPoints;
    double* points = NULL;
    SgQuadIntersect_getIntersectPoints(&qis, &numPoints, &points);
    SgQuadIntersect_del(&qis);

    std::cout << "testPartial3Points: num intersection points = " << numPoints << '\n';

    if (numPoints != 3) {
        /// error
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
