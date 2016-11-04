/** 
 * Test triangulation
 */

#include <SgTriangulate.h>
#include <iostream>
#include <cassert>

void test1Point() {
	const double p0[] = {0., 0.};
	const double* points[] = {p0};
	const int numPoints = 1;
	SgTriangulate_type* tri;
	SgTriangulate_new(&tri, numPoints, points);
    double area;
    SgTriangulate_getConvexHullArea(&tri, &area);
    std::cout << "test1Point: area = " << area << '\n';
	SgTriangulate_del(&tri);	
}

void test2Points() {
    const double p0[] = {0., 0.};
    const double p1[] = {1., 0.};
    const double* points[] = {p0, p1};
    const int numPoints = 1;
    SgTriangulate_type* tri;
    SgTriangulate_new(&tri, numPoints, points);
    double area;
    SgTriangulate_getConvexHullArea(&tri, &area);
    std::cout << "test2Points: area = " << area << '\n';
    SgTriangulate_del(&tri);    
}

void test3Points() {
    const double p0[] = {0., 0.};
    const double p1[] = {1., 0.};
    const double p2[] = {1., 1.};
    const double* points[] = {p0, p1, p2};
    const int numPoints = 1;
    SgTriangulate_type* tri;
    SgTriangulate_new(&tri, numPoints, points);
    double area;
    SgTriangulate_getConvexHullArea(&tri, &area);
    std::cout << "test3Points: area = " << area << '\n';
    SgTriangulate_del(&tri);    
}

int main(int argc, char** argv) {

	test1Point();
    test2Points();
    test3Points();

	return 0;
}
