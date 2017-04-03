/**
 * Testing flow interpolation in 2D
 */

#include <vector>
#include <cassert>
#include <cstdio>
#include <iostream>
#include "SgFlowInterp2D.h"
#include "SgBoxIterator.h"
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
    SgFlowInterp2D_debug(&interp);
    SgFlowInterp2D_del(&interp);

    // check flux projected onto segment
    double u = dstXmaxs[0] - dstXmins[0];
    double v = dstXmaxs[1] - dstXmins[1];
    double xa = dstXmins[0];
    double ya = dstXmins[1];
    const double exactFlux = 1.*u*(1 - ya - v/2.) + 3.*u*(ya + v/2.) + 4.*v*(1. - xa - u/2.) + 2.*v*(xa + u/2.);
    std::cout << "tstSimple flux: " << dstData[0] << 
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

bool testSimple2x2(const double dstXmins[], const double dstXmaxs[], double expectedVal) {

    size_t k;
    int inds[2];
    const int zeros[] = {0, 0};

    // source grid
    const int srcNodeDims[] = {3, 3}; // number of nodes
    const double srcXmins[] = {0., 0.};
    const double srcXmaxs[] = {1., 1.};
    int srcNumPoints = srcNodeDims[0] * srcNodeDims[1];
    int srcNumEdgeX = srcNodeDims[0] * (srcNodeDims[1] - 1);
    int srcNumEdgeY = (srcNodeDims[0] - 1) * srcNodeDims[1];
    double* srcCoords[] = {new double[srcNumPoints], new double[srcNumPoints]};
    createRectangularGrid(srcNodeDims, srcXmins, srcXmaxs, srcCoords);

    // src data, dimensioned number of nodes
    double* srcData[] = {new double[srcNumEdgeX], new double[srcNumEdgeY]};

    double hx = (srcXmaxs[0] - srcXmins[0])/double(srcNodeDims[0] - 1);
    double hy = (srcXmaxs[1] - srcXmins[1])/double(srcNodeDims[1] - 1);

    k = 0;
    const int srcEdgeXDims[] = {srcNodeDims[0] - 1, srcNodeDims[1]}; // x, y
    SgBoxIterator_type edgeXIt(2, zeros, srcEdgeXDims);
    for (int index = 0; index < edgeXIt.getNumberOfElements(); ++index) {
        edgeXIt.getElement(index, inds);
        int i = inds[0]; 
        int j = inds[1];
        double y = srcXmins[1] + j*hy;
        double xLo = srcXmins[0] + i*hx;
        double xHi = xLo + hx;
        srcData[k][index] = double(index);
    }

    // fluxes along y
    k = 1;
    const int srcEdgeYDims[] = {srcNodeDims[0], srcNodeDims[1] - 1}; // x, y
    SgBoxIterator_type edgeYIt(2, zeros, srcEdgeYDims);
    for (int index = 0; index < edgeYIt.getNumberOfElements(); ++index) {
        edgeYIt.getElement(index, inds);
        int i = inds[0]; 
        int j = inds[1];
        double yLo = srcXmins[1] + j*hy;
        double yHi = yLo + hy;
        double x = srcXmins[0] + i*hx;
        srcData[k][index] = double(edgeXIt.getNumberOfElements() + index);
    }

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
    SgFlowInterp2D_setSrcGrid(&interp, srcNodeDims, (const double**) srcCoords);
    SgFlowInterp2D_setDstGrid(&interp, dstDims, (const double**) dstCoords);
    SgFlowInterp2D_computeWeights(&interp);
    SgFlowInterp2D_apply(&interp, (const double**) srcData, dstData);
    SgFlowInterp2D_debug(&interp);
    SgFlowInterp2D_del(&interp);

    // check flux projected onto segment
    double u = dstXmaxs[0] - dstXmins[0];
    double v = dstXmaxs[1] - dstXmins[1];
    double xa = dstXmins[0];
    double ya = dstXmins[1];
    std::cout << "tstSimple2x2 flux: " << dstData[0] << 
                 " expected: " << expectedVal << 
                 " error: " << dstData[0] - expectedVal << '\n';

    // clean up
    for (size_t j = 0; j < 2; ++j) {
        delete[] srcCoords[j];
        delete[] dstCoords[j];
        delete[] srcData[j];
    }
    delete[] dstData;

    return true;
}

bool testLinear(const double dstXmins[], const double dstXmaxs[], const int srcNodeDims[], double expectedVal) {

    // 2-form field is
    // x*y*(dz ^ dx) + (x + y)*(dy ^ dz)

    // source grid
    const double srcXmins[] = {0., 0.};
    const double srcXmaxs[] = {1., 1.};
    int srcNumPoints = srcNodeDims[0] * srcNodeDims[1];
    double* srcCoords[] = {new double[srcNumPoints], new double[srcNumPoints]};
    createRectangularGrid(srcNodeDims, srcXmins, srcXmaxs, srcCoords);

    // src data, dimensioned number of nodes
    int srcNumEdgeX = (srcNodeDims[0] - 1) * srcNodeDims[1];
    int srcNumEdgeY = srcNodeDims[0] * (srcNodeDims[1] - 1);
    double* srcData[] = {new double[srcNumEdgeX], new double[srcNumEdgeY]};

    double hx = (srcXmaxs[0] - srcXmins[0])/double(srcNodeDims[0] - 1);
    double hy = (srcXmaxs[1] - srcXmins[1])/double(srcNodeDims[1] - 1);
    size_t k;

    const int zeros[] = {0, 0};
    int inds[] = {-1, -1};

    // fluxes along x
    k = 0;
    const int srcEdgeXDims[] = {srcNodeDims[0], srcNodeDims[1] - 1}; // y, x
    SgBoxIterator_type edgeXIt(2, zeros, srcEdgeXDims);
    for (int index = 0; index < edgeXIt.getNumberOfElements(); ++index) {
        edgeXIt.getElement(index, inds);
        int i = inds[0]; 
        int j = inds[1];
        double y = srcXmins[1] + j*hy;
        double xLo = srcXmins[0] + i*hx;
        double xHi = xLo + hx;
        srcData[k][index] = 1 * y * 0.5*(xHi*xHi - xLo*xLo);
    }

    // fluxes along y
    k = 1;
    const int srcEdgeYDims[] = {srcNodeDims[0] - 1, srcNodeDims[1]}; // y, x
    SgBoxIterator_type edgeYIt(2, zeros, srcEdgeYDims);
    for (int index = 0; index < edgeYIt.getNumberOfElements(); ++index) {
        edgeYIt.getElement(index, inds);
        int i = inds[0]; 
        int j = inds[1];
        double yLo = srcXmins[1] + j*hy;
        double yHi = yLo + hy;
        double x = srcXmins[0] + i*hx;
        srcData[k][index] = 1 * ( x*hy + 0.5*(yHi*yHi - yLo*yLo) );
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
    SgFlowInterp2D_debug(&interp);
    SgFlowInterp2D_del(&interp);

    // check flux projected onto segment
    std::cout << "testLinear flux: " << dstData[0] << 
                 " expected: " << expectedVal << '\n';
    assert(fabs(dstData[0] - expectedVal) < 1.e-3);

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
    int srcNodeDims[2];

    // make sure the segment is inside the cell
    dstXmins[0] = 0.0000001; dstXmins[1] = 0.0000001;
    dstXmaxs[0] = 0.9999999; dstXmaxs[1] = 0.0000001;
    //if (!testSimple(dstXmins, dstXmaxs)) return 1;

    // going along y
    dstXmins[0] = 0.2; dstXmins[1] = 0.0000001;
    dstXmaxs[0] = 0.2; dstXmaxs[1] = 0.9999999;
    //if (!testSimple(dstXmins, dstXmaxs)) return 2;

    // going down along y
    dstXmins[0] = 0.2; dstXmins[1] = 0.9999999;
    dstXmaxs[0] = 0.2; dstXmaxs[1] = 0.0000001;
    //if (!testSimple(dstXmins, dstXmaxs)) return 3;

    // arbitrary segment
    dstXmins[0] = 0.2; dstXmins[1] = 0.9;
    dstXmaxs[0] = 0.3; dstXmaxs[1] = 0.3;
    //if (!testSimple(dstXmins, dstXmaxs)) return 4;

    // segment along y
    dstXmins[0] = 0.0000001; dstXmins[1] = 0.00000001;
    dstXmaxs[0] = 0.0000001; dstXmaxs[1] = 0.99999999;
    srcNodeDims[0] = 3; srcNodeDims[1] = 3;
    double xa = dstXmins[0];
    double xb = dstXmaxs[0];
    double ya = dstXmins[1];
    double yb = dstXmaxs[1];
    double u = xb - xa;
    double v = yb - ya;
    double expectedVal = 1*u*(xa*ya + (u + v)/2. + u*v/3.) + 1*v*(xa + ya + (u + v)/2.);
    //if (!testLinear(dstXmins, dstXmaxs, srcNodeDims, expectedVal)) return 5;

    // other segment along y
    dstXmins[0] = 0.99999999; dstXmins[1] = 0.00000001;
    dstXmaxs[0] = 0.99999999; dstXmaxs[1] = 0.99999999;
    srcNodeDims[0] = 3; srcNodeDims[1] = 3;
    xa = dstXmins[0];
    xb = dstXmaxs[0];
    ya = dstXmins[1];
    yb = dstXmaxs[1];
    u = xb - xa;
    v = yb - ya;
    expectedVal = 1*u*(xa*ya + (u + v)/2. + u*v/3.) + 1*v*(xa + ya + (u + v)/2.);
    //if (!testLinear(dstXmins, dstXmaxs, srcNodeDims, expectedVal)) return 6;

    // segment along x
    dstXmins[0] = 0.0000001; dstXmins[1] = 0.00000001;
    dstXmaxs[0] = 0.9999999; dstXmaxs[1] = 0.00000001;
    srcNodeDims[0] = 3; srcNodeDims[1] = 3;
    xa = dstXmins[0];
    xb = dstXmaxs[0];
    ya = dstXmins[1];
    yb = dstXmaxs[1];
    u = xb - xa;
    v = yb - ya;
    expectedVal = 1*u*(xa*ya + (u + v)/2. + u*v/3.) + 1*v*(xa + ya + (u + v)/2.);
    //if (!testLinear(dstXmins, dstXmaxs, srcNodeDims, expectedVal)) return 7;

    // 2x2 segment along x
    dstXmins[0] = 0.0000001; dstXmins[1] = 0.00000001;
    dstXmaxs[0] = 0.9999999; dstXmaxs[1] = 0.00000001;
    srcNodeDims[0] = 3; srcNodeDims[1] = 3;
    xa = dstXmins[0];
    xb = dstXmaxs[0];
    ya = dstXmins[1];
    yb = dstXmaxs[1];
    u = xb - xa;
    v = yb - ya;
    expectedVal = 0. + 3.;
    if (!testSimple2x2(dstXmins, dstXmaxs, expectedVal)) return 8;

    return 0;
}
