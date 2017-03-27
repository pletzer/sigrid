/**
 * Testing flow interpolation in 2D
 */

#include <vector>
#include <cassert>
#include <cstdio>
#include <cmath>
#include <iostream>
#include "SgFlowInterp2D.h"
#include "createGrids2D.h"
#include "CmdLineArgParser.h" 

void createLineGrid(const int dims[], const double xmins[], const double xmaxs[], double* coords[]) {
    for (size_t i = 0; i < dims[0]; ++i) {
        for (size_t j = 0; j < 2; ++j) {
            coords[j][i] = xmins[j] + (xmaxs[j] - xmins[j]) * i / double(dims[0] - 1);
        }
    }
}

double psi(const double pos[]) {
    double x = pos[0];
    double y = pos[1];
    return 0.5*(x*x + cos(2.*y));
}

int main(int argc, char** argv) {

    // parse command line arguments
    CmdLineArgParser prsr;
    prsr.set("--ni", 11, "Number of nodes in the x direction");
    prsr.set("--nj", 21, "Number of nodes in the y direction");
    prsr.set("--pa", "0.,0.", "Start position");
    prsr.set("--pb", "0.,1.", "End position");
    prsr.parse(argc, argv);

    // create src grid
    const double srcXmins[] = {-1.0, -M_PI};
    const double srcXmaxs[] = {+1.0, +M_PI};
    int srcNodeDims[] = {prsr.get<int>("--ni"), prsr.get<int>("--nj")};
    int srcNumPoints = srcNodeDims[0] * srcNodeDims[1];
    double* srcCoords[] = {new double[srcNumPoints], new double[srcNumPoints]};
    createRectangularGrid(srcNodeDims, srcXmins, srcXmaxs, srcCoords);

    // src field
    int srcNumEdgeX = (srcNodeDims[0] - 1) * srcNodeDims[1];
    int srcNumEdgeY = srcNodeDims[0] * (srcNodeDims[1] - 1);
    double* srcData[] = {new double[srcNumEdgeX], new double[srcNumEdgeY]};

    double hx = (srcXmaxs[0] - srcXmins[0])/double(srcNodeDims[0] - 1);
    double hy = (srcXmaxs[1] - srcXmins[1])/double(srcNodeDims[1] - 1);
    size_t k;

    // fluxes along x (y-field)
    k = 0;
    for (size_t i = 0; i < srcNodeDims[0] - 1; ++i) {
        double xLo = srcXmins[0] + i*hx;
        double xHi = xLo + hx;
        for (size_t j = 0; j < srcNodeDims[1] - 0; ++j) {
            double y = srcXmins[1] + j*hy;
            size_t index = i*(srcNodeDims[1] - 0) + j;
            srcData[k][index] = 0.5*(xHi*xHi - xLo*xLo);
        }
    }

    // fluxes along y (x field)
    k = 1;
    for (size_t i = 0; i < srcNodeDims[0] - 0; ++i) {
        double x = srcXmins[0] + i*hx;
        for (size_t j = 0; j < srcNodeDims[1] - 1; ++j) {
            double yLo = srcXmins[1] + j*hy;
            double yHi = yLo + hy;
            size_t index = i*(srcNodeDims[1] - 1) + j;
            srcData[k][index] = 0.5*(cos(2.*yHi) - cos(2.*yLo));
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

    // project src field onto dst grid
    SgFlowInterp2D_type* interp = NULL;
    SgFlowInterp2D_new(&interp);
    SgFlowInterp2D_setSrcGrid(&interp, srcNodeDims, (const double**) srcCoords);
    SgFlowInterp2D_setDstGrid(&interp, dstDims, (const double**) dstCoords);
    SgFlowInterp2D_computeWeights(&interp);
    SgFlowInterp2D_apply(&interp, (const double**) srcData, dstData);
    SgFlowInterp2D_debug(&interp);
    SgFlowInterp2D_del(&interp);

    // check
    double exactFlux = psi(dstXmaxs) - psi(dstXmins);
    std::cout << "Flux: " << dstData[0] << 
                 " exact: " << exactFlux <<
                 " error: " << dstData[0] - exactFlux << '\n';

    // cleanup
    for (size_t j = 0; j < 2; ++j) {
        delete[] srcCoords[j];
        delete[] dstCoords[j];
        delete[] srcData[j];
    }
    delete[] dstData;

    std::cout << "SUCCESS\n";
    return 0;
}
