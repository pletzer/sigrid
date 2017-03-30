/**
 * Testing computation of intersection points between quad and line
 */

#include <vector>
#include <cassert>
#include <cstdio>
#include <iostream>
#include "SgQuadLineFlows.h"

bool test0() {
    double quadCoords[] = {0., 0., 
                           1., 0.,
                           1., 1.,
                           0., 1.};
    const double lineCoords[] = {
        0.0, 0.0,
        1.0, 0.0
    };

    SgQuadLineFlows_type qlflows;
    qlflows.setQuadPoints(quadCoords);
    qlflows.setLinePoints(lineCoords);
    int err = qlflows.computeProjections();
    if (err != 0) return false;

    // check
    double flux;
    flux = qlflows.getProjection(EDGE_LO_0);
    std::cout << "test 0: flux on x lo: " << flux << '\n';
    if (fabs(flux - 1.0) > 1.e-10) return false;
    flux = qlflows.getProjection(EDGE_HI_0);
    std::cout << "test 0: flux on x hi: " << flux << '\n';
    if (fabs(flux - 0.0) > 1.e-10) return false;
    flux = qlflows.getProjection(EDGE_LO_1);
    std::cout << "test 0: flux on y lo: " << flux << '\n';
    if (fabs(flux - 0.0) > 1.e-10) return false;
    flux = qlflows.getProjection(EDGE_HI_1);
    std::cout << "test 0: flux on y hi: " << flux << '\n';
    if (fabs(flux - 0.0) > 1.e-10) return false;

    return true;
}

bool test1() {
    double quadCoords[] = {0., 0., 
                           1., 0.,
                           1., 1.,
                           0., 1.};

    const double lineCoords[] = {
        0.0, 1.0,
        1.0, 1.0
    };

    SgQuadLineFlows_type qlflows;
    qlflows.setQuadPoints(quadCoords);
    qlflows.setLinePoints(lineCoords);
    int err = qlflows.computeProjections();
    if (err != 0) return false;

    // check
    double flux;
    flux = qlflows.getProjection(EDGE_LO_0);
    if (fabs(flux - 0.0) > 1.e-10) return false;
    flux = qlflows.getProjection(EDGE_HI_0);
    if (fabs(flux - 1.0) > 1.e-10) return false;
    flux = qlflows.getProjection(EDGE_LO_1);
    if (fabs(flux - 0.0) > 1.e-10) return false;
    flux = qlflows.getProjection(EDGE_HI_1);
    if (fabs(flux - 0.0) > 1.e-10) return false;

    return true;
}

bool test2() {
    double quadCoords[] = {0., 0., 
                           1., 0.,
                           1., 1.,
                           0., 1.};

    const double lineCoords[] = {
        0.0, 0.0,
        0.0, 1.0
    };

    SgQuadLineFlows_type qlflows;
    qlflows.setQuadPoints(quadCoords);
    qlflows.setLinePoints(lineCoords);
    int err = qlflows.computeProjections();
    if (err != 0) return false;

    // check
    double flux;
    flux = qlflows.getProjection(EDGE_LO_0);
    if (fabs(flux - 0.0) > 1.e-10) return false;
    flux = qlflows.getProjection(EDGE_HI_0);
    if (fabs(flux - 0.0) > 1.e-10) return false;
    flux = qlflows.getProjection(EDGE_LO_1);
    if (fabs(flux - 1.0) > 1.e-10) return false;
    flux = qlflows.getProjection(EDGE_HI_1);
    if (fabs(flux - 0.0) > 1.e-10) return false;

    return true;
}

bool test3() {
    double quadCoords[] = {0., 0., 
                           1., 0.,
                           1., 1.,
                           0., 1.};

    const double lineCoords[] = {
        1.0, 0.0,
        1.0, 1.0
    };

    SgQuadLineFlows_type qlflows;
    qlflows.setQuadPoints(quadCoords);
    qlflows.setLinePoints(lineCoords);
    int err = qlflows.computeProjections();
    if (err != 0) return false;

    // check
    double flux;
    flux = qlflows.getProjection(EDGE_LO_0);
    if (fabs(flux - 0.0) > 1.e-10) return false;
    flux = qlflows.getProjection(EDGE_HI_0);
    if (fabs(flux - 0.0) > 1.e-10) return false;
    flux = qlflows.getProjection(EDGE_LO_1);
    if (fabs(flux - 0.0) > 1.e-10) return false;
    flux = qlflows.getProjection(EDGE_HI_1);
    if (fabs(flux - 1.0) > 1.e-10) return false;

    return true;
}

int main(int argc, char** argv) {

    if (!test0()) return 1;
    if (!test1()) return 1;
    if (!test2()) return 1;
    if (!test3()) return 1;

    std::cout << "SUCCESS\n";
    return 0;
}
