/**
 * Testing conservative interpolation in 2D
 */

#include <vector>
#include <cassert>
#include <cstdio>
#include <iostream>
#include "SgConserveInterp2D.h"
#include "createGrids2D.h"


bool testSimple() {

    // source grid
    const int srcDims[] = {2, 2};
    const double srcXmins[] = {0., 0.};
    const double srcXmaxs[] = {1., 1.};
    const int periodicity[] = {0, 0};
    int srcNumPoints = srcDims[0] * srcDims[1];
    double* srcCoords[] = {new double[srcNumPoints], new double[srcNumPoints]};
    createRectangularGrid(srcDims, srcXmins, srcXmaxs, srcCoords);

    // destination grid
    const int dstDims[] = {2, 2};
    const double dstXmins[] = {0., 0.};
    const double dstXmaxs[] = {1., 1.};
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
    int end = 0;
    SgConserveInterp2D_reset(&interp);
    while (!end) {
        SgConserveInterp2D_get(&interp, &srcIndex, &dstIndex, &weight);
        std::cout << "testSimple src index: " << srcIndex << " dst index: " 
                  << dstIndex << " weight: " << weight << '\n';
        end = SgConserveInterp2D_next(&interp);
    }
    SgConserveInterp2D_del(&interp);

    // clean up
    for (size_t j = 0; j < 2; ++j) {
        delete[] srcCoords[j];
        delete[] dstCoords[j];
    }

    return true;
}

bool testSrcInsideDst() {

    // source grid
    const int srcDims[] = {2, 2};
    const double srcXmins[] = {0.25, 0.25};
    const double srcXmaxs[] = {0.75, 0.75};
    const int periodicity[] = {0, 0};
    int srcNumPoints = srcDims[0] * srcDims[1];
    double* srcCoords[] = {new double[srcNumPoints], new double[srcNumPoints]};
    createRectangularGrid(srcDims, srcXmins, srcXmaxs, srcCoords);

    // destination grid
    const int dstDims[] = {2, 2};
    const double dstXmins[] = {0., 0.};
    const double dstXmaxs[] = {1., 1.};
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
    while (!end) {
        SgConserveInterp2D_get(&interp, &srcIndex, &dstIndex, &weight);
        std::cout << "testSrcInsideDst src index: " << srcIndex << " dst index: " 
                  << dstIndex << " weight: " << weight << '\n';
        totalWeight += weight;
        end = SgConserveInterp2D_next(&interp);
    }
    SgConserveInterp2D_del(&interp);

    // clean up
    for (size_t j = 0; j < 2; ++j) {
        delete[] srcCoords[j];
        delete[] dstCoords[j];
    }

    std::cout << "Total weight: " << totalWeight << '\n';
    if (fabs(totalWeight - 0.25) > 1.e-10) {
        // error
        return false;
    }

    return true;
}

bool testSrc10By10() {

    // source grid
    const int srcDims[] = {11, 11}; // number of nodes
    const double srcXmins[] = {0., 0.};
    const double srcXmaxs[] = {1., 1.};
    const int periodicity[] = {0, 0};
    int srcNumPoints = srcDims[0] * srcDims[1];
    double* srcCoords[] = {new double[srcNumPoints], new double[srcNumPoints]};
    createRectangularGrid(srcDims, srcXmins, srcXmaxs, srcCoords);

    // destination grid
    const int dstDims[] = {2, 2};
    const double dstXmins[] = {0., 0.};
    const double dstXmaxs[] = {1., 1.};
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

    std::cout << "testSrc10By10: total weight = " << totalWeight << '\n';
    if (fabs(totalWeight -  1.0) > 1.e-10) {
        // sum of the weights should match number of dst cells
        return false;
    }

    return true;
}

bool testSrc10By20Dst100By200() {

    // source grid
    const int srcDims[] = {11, 21}; // number of nodes
    const double srcXmins[] = {0., 0.};
    const double srcXmaxs[] = {1., 1.};
    const int periodicity[] = {0, 0};
    int srcNumPoints = srcDims[0] * srcDims[1];
    double* srcCoords[] = {new double[srcNumPoints], new double[srcNumPoints]};
    createRectangularGrid(srcDims, srcXmins, srcXmaxs, srcCoords);

    // destination grid
    const int dstDims[] = {101, 201};
    const double dstXmins[] = {0., 0.};
    const double dstXmaxs[] = {1., 1.};
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

    std::cout << "testSrc10By20Dst100By200: total weight = " << totalWeight << '\n';
    if (fabs(totalWeight -  200*100) > 1.e-10) {
        // sum of the weights should match number of dst cells
        return false;
    }

    return true;
}

bool testPolar() {

    // source grid
    const int srcDims[] = {2, 3}; // number of nodes
    const double srcXmins[] = {0., 0.};
    const double srcXmaxs[] = {1., 1.};
    const int periodicity[] = {0, 0};
    int srcNumPoints = srcDims[0] * srcDims[1];
    double* srcCoords[] = {new double[srcNumPoints], new double[srcNumPoints]};
    createRectangularGrid(srcDims, srcXmins, srcXmaxs, srcCoords);

    // destination grid
    const int dstDims[] = {2, 5};
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
        std::cout << "tstPolar src index: " << srcIndex << " dst index: " 
                  << dstIndex << " weight: " << weight << '\n';
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
    if (fabs(totalWeight -  0.25 * dstNumCells) > 1.e-10) {
        // sum of the weights should match number of dst cells
        // that fall within the src domain
        return false;
    }

    return true;
}

int main(int argc, char** argv) {

    if (!testSimple()) return 1;
    std::cout << "testSimple PASSED\n";
    if (!testSrcInsideDst()) return 5;
    std::cout << "testSrcInsideDst PASSED\n";
    if (!testSrc10By10()) return 2;
    std::cout << "testSrc10By10 PASSED\n";
    if (!testSrc10By20Dst100By200()) return 3;
    std::cout << "testSrc10By20Dst100By200 PASSED\n";
    if (!testPolar()) return 4;
    std::cout << "testPolar PASSED\n";

    std::cout << "SUCCESS\n";
    return 0;
}
