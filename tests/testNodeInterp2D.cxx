/**
 * Testing cconservative interpolation in 2D
 */

#include <vector>
#include <cassert>
#include <cstdio>
#include <iostream>
#include "SgNodeInterp2D.h"
#include "createGrids2D.h"

/**
 * Set a linear field
 * @param numNodes number of nodes
 * @param coords grid coordinates
 * @param data (output)
 */
void setLinearField(int numNodes, const double** coords, double data[]) {
    for (int i = 0; i < numNodes; ++i) {
        double x = coords[0][i];
        double y = coords[1][i];
        data[i] = 1.2 + 2.3*x + 3.4*y + 5.6*x*y;
    }
}

/**
 * Compute the interpolation error
 * @param numNodes number of nodes
 * @param data interpolated data
 * @param dataExact exact nodeal values
 * @return error
 */
double getInterpError(int numNodes, const double data[], const double dataExact[]) {
    double res = 0;
    for (int i = 0; i < numNodes; ++i) {
        res += fabs(data[i] - dataExact[i]);
    }
    return res;
}


bool testSimple() {

    // source grid
    const int srcDims[] = {2, 2};
    const double srcXmins[] = {0., 0.};
    const double srcXmaxs[] = {1., 1.};
    const int periodicity[] = {0, 0};
    int srcNumPoints = srcDims[0] * srcDims[1];
    double* srcCoords[] = {new double[srcNumPoints], new double[srcNumPoints]};
    createRectangularGrid(srcDims, srcXmins, srcXmaxs, srcCoords);
    std::vector<double> srcData(srcNumPoints);
    setLinearField(srcNumPoints, (const double**) srcCoords, &srcData[0]);

    // destination grid
    const int dstDims[] = {2, 2};
    const double dstXmins[] = {0., 0.};
    const double dstXmaxs[] = {1., 1.};
    int dstNumPoints = dstDims[0] * dstDims[1];
    double* dstCoords[] = {new double[dstNumPoints], new double[dstNumPoints]};
    createRectangularGrid(dstDims, dstXmins, dstXmaxs, dstCoords);
    std::vector<double> dstDataExact(dstNumPoints);
    std::vector<double> dstDataInterp(dstNumPoints);
    setLinearField(dstNumPoints, (const double**) dstCoords, &dstDataExact[0]);


    SgNodeInterp2D_type* interp = NULL;
    SgNodeInterp2D_new(&interp);
    SgNodeInterp2D_setSrcGrid(&interp, srcDims, periodicity, (const double**) srcCoords);
    SgNodeInterp2D_setDstGrid(&interp, dstDims, (const double**) dstCoords);
    SgNodeInterp2D_computeWeights(&interp);
    SgNodeInterp2D_apply(&interp, &srcData[0], &dstDataInterp[0]);
    SgNodeInterp2D_del(&interp);

    double absError = getInterpError(dstNumPoints, &dstDataInterp[0], &dstDataExact[0]);
    std::cout << "testSimple: abs interp error = " << absError << '\n';
    if (absError > 1.e-12) {
        return false;
    }

    // clean up
    for (size_t j = 0; j < 2; ++j) {
        delete[] srcCoords[j];
        delete[] dstCoords[j];
    }

    return true;
}

bool testRect2Rect() {

    // source grid
    const int srcDims[] = {101, 201};
    const double srcXmins[] = {0., 0.};
    const double srcXmaxs[] = {1., 1.};
    const int periodicity[] = {0, 0};
    int srcNumPoints = srcDims[0] * srcDims[1];
    double* srcCoords[] = {new double[srcNumPoints], new double[srcNumPoints]};
    createRectangularGrid(srcDims, srcXmins, srcXmaxs, srcCoords);
    std::vector<double> srcData(srcNumPoints);
    setLinearField(srcNumPoints, (const double**) srcCoords, &srcData[0]);

    // destination grid
    const int dstDims[] = {21, 31};
    const double dstXmins[] = {0., 0.};
    const double dstXmaxs[] = {1., 1.};
    int dstNumPoints = dstDims[0] * dstDims[1];
    double* dstCoords[] = {new double[dstNumPoints], new double[dstNumPoints]};
    createRectangularGrid(dstDims, dstXmins, dstXmaxs, dstCoords);
    std::vector<double> dstDataExact(dstNumPoints);
    std::vector<double> dstDataInterp(dstNumPoints);
    setLinearField(dstNumPoints, (const double**) dstCoords, &dstDataExact[0]);


    SgNodeInterp2D_type* interp = NULL;
    SgNodeInterp2D_new(&interp);
    SgNodeInterp2D_setSrcGrid(&interp, srcDims, periodicity, (const double**) srcCoords);
    SgNodeInterp2D_setDstGrid(&interp, dstDims, (const double**) dstCoords);
    SgNodeInterp2D_computeWeights(&interp);
    SgNodeInterp2D_apply(&interp, &srcData[0], &dstDataInterp[0]);
    SgNodeInterp2D_del(&interp);

    double absError = getInterpError(dstNumPoints, &dstDataInterp[0], &dstDataExact[0]);
    std::cout << "testRect2Rect: abs interp error = " << absError << '\n';
    if (absError > 1.e-8) {
        return false;
    }

    // clean up
    for (size_t j = 0; j < 2; ++j) {
        delete[] srcCoords[j];
        delete[] dstCoords[j];
    }

    return true;
}


bool testRect2Polar() {

    // source grid
    const int srcDims[] = {101, 201};
    const double srcXmins[] = {-1., -1.};
    const double srcXmaxs[] = {1., 1.};
    const int periodicity[] = {0, 0};
    int srcNumPoints = srcDims[0] * srcDims[1];
    double* srcCoords[] = {new double[srcNumPoints], new double[srcNumPoints]};
    createRectangularGrid(srcDims, srcXmins, srcXmaxs, srcCoords);
    std::vector<double> srcData(srcNumPoints);
    setLinearField(srcNumPoints, (const double**) srcCoords, &srcData[0]);

    // destination grid
    const int dstDims[] = {11, 41};
    int dstNumPoints = dstDims[0] * dstDims[1];
    double* dstCoords[] = {new double[dstNumPoints], new double[dstNumPoints]};
    const double center[] = {0., 0.};
    const double radius = 1.0;
    createPolarGrid(dstDims, center, radius, dstCoords);
    std::vector<double> dstDataExact(dstNumPoints);
    std::vector<double> dstDataInterp(dstNumPoints);
    setLinearField(dstNumPoints, (const double**) dstCoords, &dstDataExact[0]);


    SgNodeInterp2D_type* interp = NULL;
    SgNodeInterp2D_new(&interp);
    SgNodeInterp2D_setSrcGrid(&interp, srcDims, periodicity, (const double**) srcCoords);
    SgNodeInterp2D_setDstGrid(&interp, dstDims, (const double**) dstCoords);
    SgNodeInterp2D_computeWeights(&interp);
    SgNodeInterp2D_apply(&interp, &srcData[0], &dstDataInterp[0]);
    SgNodeInterp2D_del(&interp);

    double absError = getInterpError(dstNumPoints, &dstDataInterp[0], &dstDataExact[0]);
    std::cout << "testRect2Polar: abs interp error = " << absError << '\n';
    if (absError > 1.e-8) {
        return false;
    }

    // clean up
    for (size_t j = 0; j < 2; ++j) {
        delete[] srcCoords[j];
        delete[] dstCoords[j];
    }

    return true;
}

bool testPolar2Rect() {

    // source grid
    const int srcDims[] = {21, 65};
    const double center[] = {0., 0.};
    const double radius = 1.0;
    const int periodicity[] = {0, 1};
    int srcNumPoints = srcDims[0] * srcDims[1];
    double* srcCoords[] = {new double[srcNumPoints], new double[srcNumPoints]};
    createPolarGrid(srcDims, center, radius, srcCoords);
    std::vector<double> srcData(srcNumPoints);
    setLinearField(srcNumPoints, (const double**) srcCoords, &srcData[0]);

    // destination grid
    const int dstDims[] = {11, 11};
    int dstNumPoints = dstDims[0] * dstDims[1];
    double* dstCoords[] = {new double[dstNumPoints], new double[dstNumPoints]};
    // ensure the dst grid is insde the src grid
    const double dstXmins[] = {-1./sqrt(2.), -1./sqrt(2.)};
    const double dstXmaxs[] = {1./sqrt(2.), 1./sqrt(2.)};
    createRectangularGrid(dstDims, dstXmins, dstXmaxs, dstCoords);
    std::vector<double> dstDataExact(dstNumPoints);
    std::vector<double> dstDataInterp(dstNumPoints);
    setLinearField(dstNumPoints, (const double**) dstCoords, &dstDataExact[0]);


    SgNodeInterp2D_type* interp = NULL;
    SgNodeInterp2D_new(&interp);
    SgNodeInterp2D_setSrcGrid(&interp, srcDims, periodicity, (const double**) srcCoords);
    SgNodeInterp2D_setDstGrid(&interp, dstDims, (const double**) dstCoords);
    SgNodeInterp2D_computeWeights(&interp);
    SgNodeInterp2D_apply(&interp, &srcData[0], &dstDataInterp[0]);
    SgNodeInterp2D_del(&interp);

    double absError = getInterpError(dstNumPoints, &dstDataInterp[0], &dstDataExact[0]);
    std::cout << "testPolar2Rect: abs interp error = " << absError << '\n';
    if (absError > 1.e-8) {
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

    if (!testSimple()) return 1;
    if (!testRect2Rect()) return 2;
    if (!testRect2Polar()) return 3;
    if (!testPolar2Rect()) return 4;
    std::cout << "SUCCESS\n";
    return 0;
}
