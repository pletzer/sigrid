/**
 * Testing cell search in 1D
 */

 #include <vector>
 #include <cassert>
 #include <cstdio>
 #include <iostream>
 #include "SgQuadIntersect.h"

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
 	SgQuadIntersect_new(&qis, quad1Coords, quad2Coords);
 	int numPoints;
 	double* points = NULL;
 	SgQuadIntersect_getPoints(&qis, &numPoints, &points);
 	SgQuadIntersect_del(&qis);

 	std::cout << "testNoOverlap: num intersection points = " << numPoints << '\n';

 	if (numPoints != 0) {
 		return false;
 	}
 	return true;
 }

int main(int argc, char** argv) {

	if (!testNoOverlap()) return 1;

	return 0;
}
