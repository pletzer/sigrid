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
    assert(std::abs(area - 0.0) < 1.e-10);
	SgTriangulate_del(&tri);	
}

void test2Points() {
    const double p0[] = {0., 0.};
    const double p1[] = {1., 0.};
    const double* points[] = {p0, p1};
    const int numPoints = 2;
    SgTriangulate_type* tri;
    SgTriangulate_new(&tri, numPoints, points);
    double area;
    SgTriangulate_getConvexHullArea(&tri, &area);
    std::cout << "test2Points: area = " << area << '\n';
    assert(std::abs(area - 0.0) < 1.e-10);
    SgTriangulate_del(&tri);    
}

void test3Points() {
    const double p0[] = {0., 0.};
    const double p1[] = {1., 0.};
    const double p2[] = {1., 1.};
    const double* points[] = {p0, p1, p2};
    const int numPoints = 3;
    SgTriangulate_type* tri;
    SgTriangulate_new(&tri, numPoints, points);
    double area;
    SgTriangulate_getConvexHullArea(&tri, &area);
    std::cout << "test3Points: area = " << area << '\n';
    assert(std::abs(area - 0.5) < 1.e-10);
    SgTriangulate_del(&tri);    
}

void test4Points() {
    const double p0[] = {0., 0.};
    const double p1[] = {1., 0.};
    const double p2[] = {1., 1.};
    const double p3[] = {0., 1.};
    const double* points[] = {p0, p1, p2, p3};
    const int numPoints = 4;
    SgTriangulate_type* tri;
    SgTriangulate_new(&tri, numPoints, points);
    double area;
    SgTriangulate_getConvexHullArea(&tri, &area);
    std::cout << "test4Points: area = " << area << '\n';
    assert(std::abs(area - 1.0) < 1.e-10);
    SgTriangulate_del(&tri);    
}

void test6Points() {
    const double p0[] = {0., 0.};
    const double p1[] = {1., 0.};
    const double p2[] = {1., 1.};
    const double p3[] = {0., 1.};
    const double p4[] = {0.5, 1.};
    const double p5[] = {0.5, 0.};
    const double* points[] = {p0, p1, p2, p3, p4, p5};
    const int numPoints = 6;
    SgTriangulate_type* tri;
    SgTriangulate_new(&tri, numPoints, points);
    double area;
    SgTriangulate_getConvexHullArea(&tri, &area);
    std::cout << "test6Points: area = " << area << '\n';
    assert(std::abs(area - 1.0) < 1.e-10);
    SgTriangulate_del(&tri);    
}

int main(int argc, char** argv) {

	test1Point();
    test2Points();
    test3Points();
    test4Points();
    test6Points();

	return 0;
}
