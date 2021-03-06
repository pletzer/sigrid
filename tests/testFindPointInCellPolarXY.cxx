/**
 * Testing cell search in 2D curvilinear grid
 */

#include <vector>
#include <cassert>
#include <cstdio>
#include <iostream>
#include <cmath>
#include "SgFindPointInCell.h"

int main(int argc, char** argv) {

    int ier;
    const int nitermax = 10;
    const double tolpos = 1.e-10;

    // number of dimensions
    const int ndims = 2;

    // initial guess
    double dIndices[] = {0.1, 0.2}; // cannot be exactly (0, 0) since it is a singular point

    // target position in x-y space
    const double targetPoint[] = {0.7, 0.3};
    std::cout << "target point is ";
    for (size_t i = 0; i < ndims; ++i) std::cout << targetPoint[i] << ' ';
    std::cout << '\n';

    // grid
    size_t nr = 11;
    size_t nt = 21;
    std::vector<double> x(nr * nt);
    std::vector<double> y(nr * nt);
    size_t k = 0;
    for (size_t i = 0; i < nr; ++i) {
        double r = 0. + (1. - 0.) * i / double(nr - 1);
        for (size_t j = 0; j < nt; ++j) {
            double t = 0. + (2*M_PI - 0.)* j / double(nt - 1);
            x[k] = r * cos(t);
            y[k] = r * sin(t);
            k++;
        }
    }

    const int dims[] = {nr, nt};
    std::vector< double* > coords(ndims);
    coords[0] = &x[0];
    coords[1] = &y[0];

    double pos[ndims];
    double oldPos[ndims];

    SgFindPointInCell_type* picf = NULL;
    ier = SgFindPointInCell_new(&picf, nitermax, tolpos);
    assert(ier == 0);

    // periodic in the second index
    const int periodicity[] = {0, 1};
    ier = SgFindPointInCell_setGrid(&picf, ndims, dims,
                                    periodicity, (const double**) &coords[0]);
    assert(ier == 0);

    ier = SgFindPointInCell_reset(&picf, dIndices, targetPoint);
    assert(ier == 0);

    int iterDone = 0;
    int icount = 0;
    while (iterDone == 0) {

        ier = SgFindPointInCell_getPosition(&picf, oldPos);
        assert(ier == 0);

        iterDone = SgFindPointInCell_next(&picf);
        if (iterDone < 0) {
            std::cout << "*** reached max number of iterations!\n";
        }

        ier = SgFindPointInCell_getPosition(&picf, pos);
        assert(ier == 0);

        std::cout << "iter " << icount << " position old = ";
        for (int j = 0; j < ndims; ++j) std::cout << oldPos[j] << " ";
        std::cout << " -> new = ";
        for (int j = 0; j < ndims; ++j) std::cout << pos[j] << " ";

        double error;
        ier = SgFindPointInCell_getError(&picf, &error);
        assert(ier == 0);
        std::cout << " (error = " << error << ")\n";

        icount++;
    }

    if (iterDone < 0) {
        return 1;
    }

    ier = SgFindPointInCell_del(&picf);
    assert(ier == 0);

    return 0;
}
