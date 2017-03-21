/**
 * Testing flow interpolation in 2D
 */

#include <vector>
#include <cassert>
#include <cstdio>
#include <iostream>
#include "SgFlowInterp2D.h"
#include "createGrids2D.h"

void createLineGrid(const int dims[], const double xmins[], const double xmaxs[], double* coords[]) {
    for (size_t i = 0; i < dims[0]; ++i) {
        for (size_t j = 0; j < 2; ++j) {
            coords[j][i] = xmins[j] + (xmaxs[j] - xmins[j]) * i / double(dims[0] - 1);
        }
    }
}


bool testSimple() {

    // source grid
    const int srcDims[] = {2, 2}; // number of nodes
    const double srcXmins[] = {0., 0.};
    const double srcXmaxs[] = {1., 1.};
    int srcNumPoints = srcDims[0] * srcDims[1];
    int srcNumEdgeX = srcDims[0] * (srcDims[1] - 1);
    int srcNumEdgeY = (srcDims[0] - 1) * srcDims[1];
    double* srcCoords[] = {new double[srcNumPoints], new double[srcNumPoints]};
    createRectangularGrid(srcDims, srcXmins, srcXmaxs, srcCoords);

    // src data, dimensioned number of nodes
    double* srcData[] = {new double[srcNumEdgeX], new double[srcNumEdgeY]};

    // src fluxes. Shown below are the edge values
    //  
    //      3   
    //   +-------+
    // 4 |       | 2
    //   |       |
    //   +-------+
    //       1

    // fluxes along x edges
    srcData[0][0] = 1.0; // flux along x edge
    srcData[0][1] = 3.0;

    // fluxes along y edges
    srcData[1][0] = 4.0;
    srcData[1][1] = 2.0;

    // destination grid, a segment
    const int dstDims[] = {2}; // number of nodes
    const double dstXmins[] = {0., 0.};
    const double dstXmaxs[] = {1., 1.};
    int dstNumPoints = dstDims[0];
    double* dstCoords[] = {new double[dstNumPoints], new double[dstNumPoints]};
    createLineGrid(dstDims, dstXmins, dstXmaxs, dstCoords);    

    // dst data, dimensioned number of nodes
    double* dstData[] = {new double[dstNumPoints], new double[dstNumPoints]};

    SgFlowInterp2D_type* interp = NULL;
    SgFlowInterp2D_new(&interp);
    SgFlowInterp2D_setSrcGrid(&interp, srcDims, (const double**) srcCoords);
    SgFlowInterp2D_setDstGrid(&interp, dstDims, (const double**) dstCoords);
    SgFlowInterp2D_computeWeights(&interp);
    SgFlowInterp2D_apply(&interp, 0, srcData[0], dstData[0]);
    SgFlowInterp2D_apply(&interp, 1, srcData[1], dstData[1]);
    SgFlowInterp2D_del(&interp);

    // check flux projected onto segment
    double deltaX = dstXmaxs[0] - dstXmins[0];
    double deltaY = dstXmaxs[1] - dstXmins[1];
    double avrgX = 0.5*(dstXmaxs[0] + dstXmins[0]);
    double avrgY = 0.5*(dstXmaxs[1] + dstXmins[1]);
    double xFlux = deltaX*avrgY*(1. + 3.);
    double yFlux = deltaY*avrgX*(2. + 4.);
    std::cout << "Flux in x: " << dstData[0][0] << 
                 " expected: " << xFlux << 
                 " error: " << dstData[0][0] - xFlux << '\n';
    assert(fabs(dstData[0][0] - xFlux) < 1.e-8);
    std::cout << "Flux in y: " << dstData[1][0] << 
                 " expected: " << yFlux << 
                 " error: " << dstData[1][0] - yFlux << '\n';
    assert(fabs(dstData[1][0] - yFlux) < 1.e-8);

    // clean up
    for (size_t j = 0; j < 2; ++j) {
        delete[] srcCoords[j];
        delete[] dstCoords[j];
        delete[] srcData[j];
        delete[] dstData[j];
    }

    return true;
}


int main(int argc, char** argv) {

    if (!testSimple()) return 1;

    std::cout << "SUCCESS\n";
    return 0;
}
