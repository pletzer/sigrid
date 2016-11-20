/** 
 * Test triangulation
 */

#include <SgTriangulate.h>
#include <iostream>
#include <cassert>

void test1Point() {
    const double points[] = {0., 0.};
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
    const double points[] = {0., 0.,
                             1., 0.};
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
    const double points[] = {0., 0.,
                             1., 0.,
                             1., 1.};
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
    const double points[] = {0., 0.,
                             1., 0.,
                             1., 1.,
                             0., 1.};
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
    const double points[] = {0., 0.,
                             1., 0.,
                             1., 1.,
                             0., 1.,
                             0.5, 1.,
                             0.5, 0.};
    const int numPoints = 6;
    SgTriangulate_type* tri;
    SgTriangulate_new(&tri, numPoints, points);
    double area;
    SgTriangulate_getConvexHullArea(&tri, &area);
    std::cout << "test6Points: area = " << area << '\n';
    assert(std::abs(area - 1.0) < 1.e-10);
    SgTriangulate_del(&tri);    
}

void test7PointsDegenerate() {
    const double points[] = {0., 0.,
                             1., 0.,
                             1., 1.,
                             0., 1.,
                             0., 0.,
                             0.1, 0.,
                             0.2, 0.};
    const int numPoints = 7;
    SgTriangulate_type* tri;
    SgTriangulate_new(&tri, numPoints, points);
    double area;
    SgTriangulate_getConvexHullArea(&tri, &area);
    std::cout << "test7PointsDegenerate: area = " << area << '\n';
    assert(std::abs(area - 1.0) < 1.e-10);
    SgTriangulate_del(&tri);    
}

void test4PointsTriangle() {
    const double points[] = {6.12e-17, 1,
                             0, 0.5,
                             0, 1,
                             0.5, 0.49999999999};
    const int numPoints = 4;
    SgTriangulate_type* tri;
    SgTriangulate_new(&tri, numPoints, points);
    double area;
    SgTriangulate_getConvexHullArea(&tri, &area);
    std::cout << "test4PointsTriangle: area = " << area << '\n';
    assert(std::abs(area - 0.125) < 1.e-10);
    SgTriangulate_del(&tri);    
}

int main(int argc, char** argv) {

    test1Point();
    test2Points();
    test3Points();
    test4Points();
    test6Points();
    test7PointsDegenerate();
    test4PointsTriangle();

    return 0;
}
