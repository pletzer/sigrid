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


bool testSimple(const double dstXmins[], const double dstXmaxs[]) {

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
    int dstNumPoints = dstDims[0];
    int dstNumCells = dstDims[0] - 1;
    double* dstCoords[] = {new double[dstNumPoints], new double[dstNumPoints]};
    createLineGrid(dstDims, dstXmins, dstXmaxs, dstCoords);    

    // dst data, dimensioned number of nodes
    double* dstData = new double[dstNumCells];

    SgFlowInterp2D_type* interp = NULL;
    SgFlowInterp2D_new(&interp);
    SgFlowInterp2D_setSrcGrid(&interp, srcDims, (const double**) srcCoords);
    SgFlowInterp2D_setDstGrid(&interp, dstDims, (const double**) dstCoords);
    SgFlowInterp2D_computeWeights(&interp);
    SgFlowInterp2D_apply(&interp, (const double**) srcData, dstData);
    SgFlowInterp2D_del(&interp);

    // check flux projected onto segment
    double hx = dstXmaxs[0] - dstXmins[0];
    double hy = dstXmaxs[1] - dstXmins[1];
    double avrgX = 0.5*(dstXmaxs[0] + dstXmins[0]);
    double avrgY = 0.5*(dstXmaxs[1] + dstXmins[1]);
    const double exactFlux = hx*( (1. - avrgY)*1. + avrgY*3. ) + hy*( (1. - avrgX)*4. + avrgX*2. );
    std::cout << "Flux: " << dstData[0] << 
                 " expected: " << exactFlux << 
                 " error: " << dstData[0] - exactFlux << '\n';

    // clean up
    for (size_t j = 0; j < 2; ++j) {
        delete[] srcCoords[j];
        delete[] dstCoords[j];
        delete[] srcData[j];
    }
    delete[] dstData;

    return true;
}

bool testLinear(const double dstXmins[], const double dstXmaxs[], const int srcNodeDims[]) {

    // source grid
    const double srcXmins[] = {0., 0.};
    const double srcXmaxs[] = {1., 1.};
    int srcNumPoints = srcNodeDims[0] * srcNodeDims[1];
    int srcNumEdgeX = srcNodeDims[0] * (srcNodeDims[1] - 1);
    int srcNumEdgeY = (srcNodeDims[0] - 1) * srcNodeDims[1];
    double* srcCoords[] = {new double[srcNumPoints], new double[srcNumPoints]};
    createRectangularGrid(srcNodeDims, srcXmins, srcXmaxs, srcCoords);

    // src data, dimensioned number of nodes
    double* srcData[] = {new double[srcNumEdgeX], new double[srcNumEdgeY]};

    double hx = (srcXmaxs[0] - srcXmins[0])/double(srcNodeDims[1] - 1);
    double hy = (srcXmaxs[1] - srcXmins[1])/double(srcNodeDims[0] - 1);
    size_t k;

    // fluxes along x
    k = 0;
    for (size_t j = 0; j < srcNodeDims[0]; ++j) {
        double y = srcXmins[1] + j*hy;
        for (size_t i = 0; i < srcNodeDims[1] - 1; ++i) {
            size_t index = j*(srcNodeDims[1] - 1) + i;
            double xLo = srcXmins[0] + i*hx;
            double xHi = xLo + hx;
            srcData[k][index] = 0.5*y*(xHi*xHi - xLo*xLo);
        }
    }

    // fluxes along y
    k = 1;
    for (size_t j = 0; j < srcNodeDims[0] - 1; ++j) {
        double yLo = srcXmins[1] + j*hy;
        double yHi = yLo + hy;
        for (size_t i = 0; i < srcNodeDims[1]; ++i) {
            size_t index = j*srcNodeDims[1] + i;
            double x = srcXmins[0] + i*hx;
            srcData[k][index] = x*hy + 0.5*(yHi*yHi - yLo*yLo);
        }
    }

    // destination grid, a segment
    const int dstDims[] = {2}; // number of nodes
    int dstNumPoints = dstDims[0];
    int dstNumCells = dstDims[0] - 1;
    double* dstCoords[] = {new double[dstNumPoints], new double[dstNumPoints]};
    createLineGrid(dstDims, dstXmins, dstXmaxs, dstCoords);    

    // dst data, dimensioned number of cells
    double* dstData = new double[dstNumCells];

    SgFlowInterp2D_type* interp = NULL;
    SgFlowInterp2D_new(&interp);
    SgFlowInterp2D_setSrcGrid(&interp, srcNodeDims, (const double**) srcCoords);
    SgFlowInterp2D_setDstGrid(&interp, dstDims, (const double**) dstCoords);
    SgFlowInterp2D_computeWeights(&interp);
    SgFlowInterp2D_apply(&interp, (const double**) srcData, dstData);
    SgFlowInterp2D_del(&interp);

    // check flux projected onto segment
    const double xa = dstXmins[0];
    const double ya = dstXmins[1];
    const double xb = dstXmaxs[0];
    const double yb = dstXmaxs[1];
    double exactFlux = hy*(xa + ya + 0.5*(hx + hy)) 
                     + hx*(xa*ya + 0.5*(xb*ya + xa*yb - 2*xa*ya) + (1./3.)*hx*hy);
    std::cout << "Flux: " << dstData[0] << 
                 " expected: " << exactFlux << 
                 " error: " << dstData[0] - exactFlux << '\n';
    assert(fabs(dstData[0] - exactFlux) < 1.e-8);

    // clean up
    for (size_t j = 0; j < 2; ++j) {
        delete[] srcCoords[j];
        delete[] dstCoords[j];
        delete[] srcData[j];
    }
    delete[] dstData;

    return true;
}

int main(int argc, char** argv) {

    double dstXmins[2];
    double dstXmaxs[2];

    dstXmins[0] = 0.; dstXmins[1] = 0.;
    dstXmaxs[0] = 1.; dstXmaxs[1] = 1.;
    if (!testSimple(dstXmins, dstXmaxs)) return 1;

    dstXmins[0] = 0.2; dstXmins[1] = 0.;
    dstXmaxs[0] = 0.2; dstXmaxs[1] = 1.;
    if (!testSimple(dstXmins, dstXmaxs)) return 1;

    dstXmins[0] = 0.2; dstXmins[1] = 1.;
    dstXmaxs[0] = 0.2; dstXmaxs[1] = 0.;
    if (!testSimple(dstXmins, dstXmaxs)) return 1;

    dstXmins[0] = 0.2; dstXmins[1] = 0.9;
    dstXmaxs[0] = 0.3; dstXmaxs[1] = 0.3;
    if (!testSimple(dstXmins, dstXmaxs)) return 1;

    std::cout << "SUCCESS\n";
    return 0;
}
