/**
 * Testing conservative interpolation in 2D
 */

#include <vector>
#include <cassert>
#include <cstdio>
#include <iostream>
#include "SgConserveInterp2D.h"
#include "createGrids2D.h"


bool testPolar() {

    // source grid
    const int srcDims[] = {11, 41};
    const double center[] = {0., 0.};
    const double radius = sqrt(2.0);
    int srcNumPoints = srcDims[0] * srcDims[1];
    double* srcCoords[] = {new double[srcNumPoints], new double[srcNumPoints]};
    createPolarGrid(srcDims, center, radius, srcCoords);    

    // destination grid
    const int dstDims[] = {101, 201}; // number of nodes
    const double dstXmins[] = {-1., -1.};
    const double dstXmaxs[] = {1., 1.};
    const int periodicity[] = {0, 1};
    int dstNumPoints = dstDims[0] * dstDims[1];
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
    std::map<size_t, double> dstIndex2Weight;
    while (!end) {
        SgConserveInterp2D_get(&interp, &srcIndex, &dstIndex, &weight);
        totalWeight += weight;
        end = SgConserveInterp2D_next(&interp);

        // accummulate the weights for each dst cell
        std::map<size_t, double>::iterator it = dstIndex2Weight.find(dstIndex);
        if (it != dstIndex2Weight.end()) {
            it->second += weight;
        }
        else {
            dstIndex2Weight.insert(std::pair<size_t, double>(dstIndex, weight));
        }
    }
    SgConserveInterp2D_del(&interp);

    // check that sum of weights for each dst cell is 1
    for (size_t dstIndex = 0; dstIndex < dstIndex2Weight.size(); ++dstIndex) {
        double w = dstIndex2Weight[dstIndex];
        if (fabs(w - 1.0) > 1.e-10) {
            std::cout << "*** dst index : " << dstIndex << " has weight: " << w << '\n';
        }
    }

    // clean up
    for (size_t j = 0; j < 2; ++j) {
        delete[] srcCoords[j];
        delete[] dstCoords[j];
    }

    std::cout << "testPolar: total weight = " << totalWeight << '\n';
    int dstNumCells = (dstDims[0] - 1) * (dstDims[1] - 1);
    double error = totalWeight -  1.0 * dstNumCells;
    if (fabs(error) > 1.e-10) {
        // sum of the weights should match number of dst cells
        // that fall within the src domain
        std::cout << "error = " << error << '\n';
        return false;
    }

    return true;
}

int main(int argc, char** argv) {

    if (!testPolar()) return 1;

    std::cout << "SUCCESS\n";
    return 0;
}
