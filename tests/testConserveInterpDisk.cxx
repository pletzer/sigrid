/**
 * Testing cconservative interpolation in 2D
 */

#include <vector>
#include <cassert>
#include <cstdio>
#include <iostream>
#include "SgConserveInterp2D.h"
#include "createGrids2D.h"

bool test() {

    // source grid is disk
    const int srcDims[] = {5, 33};
    const int periodicity[] = {0, 1};
    int srcNumPoints = srcDims[0] * srcDims[1];
    double* srcCoords[] = {new double[srcNumPoints], new double[srcNumPoints]};
    const double radius = 1.0;
    createPolarGrid(srcDims, radius, srcCoords);

    // destination grid is rectangle encompassing src grid
    const int dstDims[] = {11, 11}; //{401, 401};
    const double dstXmins[] = {-1., -1.};
    const double dstXmaxs[] = {1., 1.};
    int dstNumPoints = dstDims[0] * dstDims[1];
    int dstNumCells = (dstDims[0] - 1) * (dstDims[1] - 1);
    double* dstCoords[] = {new double[dstNumPoints], new double[dstNumPoints]};
    createRectangularGrid(dstDims, dstXmins, dstXmaxs, dstCoords);    

    SgConserveInterp2D_type* interp = NULL;
    SgConserveInterp2D_new(&interp);
    SgConserveInterp2D_setSrcGrid(&interp, srcDims, periodicity, (const double**) srcCoords);
    SgConserveInterp2D_setDstGrid(&interp, dstDims, (const double**) dstCoords);
    SgConserveInterp2D_computeWeights(&interp);
    int srcIndex, dstIndex;
    double weight;
    double totalWeight = 0;
    int end = 0;
    SgConserveInterp2D_reset(&interp);
    while (!end) {
        SgConserveInterp2D_get(&interp, &srcIndex, &dstIndex, &weight);
        totalWeight += weight;
        end = SgConserveInterp2D_next(&interp);
    }
    SgConserveInterp2D_del(&interp);

    double ratioOfValidCells = (double) totalWeight / (double) dstNumCells;
    std::cout << "Total weight: " << totalWeight << " totalWeight/dstNumCells = " << 
        ratioOfValidCells << " ~ pi/4 = " << M_PI/4.0 << '\n';
    // expect partial cells
    if (fabs(ratioOfValidCells - 0.780361) > 1.e-5) {
        // error
        std::cout << "error = " << ratioOfValidCells - 0.780361 << '\n';
        return false;
    }

    // clean up
    for (size_t j = 0; j < 2; ++j) {
        delete[] srcCoords[j];
        delete[] dstCoords[j];
    }

    return true;
}


int main(int argc, char** argv) {

    if (!test()) return 1;

    std::cout << "SUCCESS\n";
    return 0;
}
