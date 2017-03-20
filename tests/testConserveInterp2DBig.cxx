/**
 * Testing conservative interpolation in 2D
 */

#include <vector>
#include <cassert>
#include <cstdio>
#include <iostream>
#include "SgConserveInterp2D.h"
#include "createGrids2D.h"
#include "CmdLineArgParser.h"

bool testPolar(const int srcDims[], const int dstDims[]) {

    // source grid
    const double srcXmins[] = {-1., -1.};
    const double srcXmaxs[] = {1., 1.};
    const int periodicity[] = {0, 1};
    int srcNumPoints = srcDims[0] * srcDims[1];
    double* srcCoords[] = {new double[srcNumPoints], new double[srcNumPoints]};
    createRectangularGrid(srcDims, srcXmins, srcXmaxs, srcCoords);

    // destination grid
    const double center[] = {0., 0.};
    const double radius = 1.0;
    int dstNumPoints = dstDims[0] * dstDims[1];
    double* dstCoords[] = {new double[dstNumPoints], new double[dstNumPoints]};
    createPolarGrid(dstDims, center, radius, dstCoords);    

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

    // clean up
    for (size_t j = 0; j < 2; ++j) {
        delete[] srcCoords[j];
        delete[] dstCoords[j];
    }

    std::cout << "testPolar: total weight = " << totalWeight << '\n';
    int dstNumCells = (dstDims[0] - 1) * (dstDims[1] - 1);
    if (fabs(totalWeight -  1.0 * dstNumCells) > 1.e-10) {
        // sum of the weights should match number of dst cells
        // that fall within the src domain
        return false;
    }

    return true;
}

int main(int argc, char** argv) {

    CmdLineArgParser prsr;
    prsr.set("--src_nj", 101, "Number of source grid latitudes");
    prsr.set("--src_ni", 201, "Number of source grid longitudes");
    prsr.set("--dst_nj", 21, "Number of destination grid latitudes");
    prsr.set("--dst_ni", 21, "Number of destination grid longitudes");
    prsr.parse(argc, argv);

    if (prsr.get<bool>("-h") || prsr.get<bool>("--help")) {
        prsr.help();
        return 1;
    }

    int src_nj = prsr.get<int>("--src_nj");
    int src_ni = prsr.get<int>("--src_ni");
    int dst_nj = prsr.get<int>("--dst_nj");
    int dst_ni = prsr.get<int>("--dst_ni");

    int srcNodeDims[] = {src_nj, src_ni};
    int dstNodeDims[] = {dst_nj, dst_ni};

    if (!testPolar(srcNodeDims, dstNodeDims)) return 1;

    std::cout << "SUCCESS\n";
    return 0;
}
