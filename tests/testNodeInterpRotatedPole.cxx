/**
 * Testing nodal interpolation in 2D
 */

#include <vector>
#include <cassert>
#include <cstdio>
#include <iostream>
#include "SgNodeInterp2D.h"
#include "createGrids2D.h"
#include "CmdLineArgParser.h"

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

int main(int argc, char** argv) {

    CmdLineArgParser prsr;
    prsr.set("--src_nj", 41, "Number of source grid latitudes");
    prsr.set("--src_ni", 81, "Number of source grid longitudes");
    prsr.set("--dst_nj", 21, "Number of destination grid latitudes");
    prsr.set("--dst_ni", 41, "Number of destination grid longitudes");
    prsr.set("--delta_lat", 30.0, "Latitude pole displacement in deg.");
    prsr.set("--delta_lon", 20.0, "Longitude pole displacement in deg.");
    prsr.set("--nitermax", 100, "Max number of iterations");
    prsr.set("--tolpos", 1.e-10, "Tolerance in physical space");
    prsr.parse(argc, argv);

    if (prsr.get<bool>("-h") || prsr.get<bool>("--help")) {
        prsr.help();
        return 1;
    }

    // source grid
    const int srcDims[] = {prsr.get<int>("--src_nj"), prsr.get<int>("--src_ni")};
    int srcNumPoints = srcDims[0] * srcDims[1];
    double* srcCoords[] = {new double[srcNumPoints], new double[srcNumPoints]};
    double delta_lat = prsr.get<double>("--delta_lat");
    double delta_lon = prsr.get<double>("--delta_lon");
    createRotatedPoleGrid(srcDims, delta_lat, delta_lon, srcCoords);
    std::vector<double> srcData(srcNumPoints);
    setLinearField(srcNumPoints, (const double**) srcCoords, &srcData[0]);

    // destination grid
    const int dstDims[] = {prsr.get<int>("--dst_nj"), prsr.get<int>("--dst_ni")};
    int dstNumPoints = dstDims[0] * dstDims[1];
    double* dstCoords[] = {new double[dstNumPoints], new double[dstNumPoints]};
    // ensure that dst grid is inside the src grid
    const double dstXmins[] = {-90.0, -180.0};
    const double dstXmaxs[] = {90.0, 180.0};
    createRectangularGrid(dstDims, dstXmins, dstXmaxs, dstCoords);
    std::vector<double> dstDataExact(dstNumPoints);
    std::vector<double> dstDataInterp(dstNumPoints);
    setLinearField(dstNumPoints, (const double**) dstCoords, &dstDataExact[0]);

    const int nitermax = prsr.get<int>("--nitermax");
    const double tolpos = prsr.get<double>("--tolpos");

    SgNodeInterp2D_type* interp = NULL;
    SgNodeInterp2D_new(&interp, nitermax, tolpos);
    const int periodicity[] = {0, 1};
    SgNodeInterp2D_setSrcGrid(&interp, srcDims, periodicity, (const double**) srcCoords);
    SgNodeInterp2D_setDstGrid(&interp, dstDims, (const double**) dstCoords);
    SgNodeInterp2D_computeWeights(&interp);
    SgNodeInterp2D_apply(&interp, &srcData[0], &dstDataInterp[0]);
    SgNodeInterp2D_del(&interp);

    double absError = getInterpError(dstNumPoints, &dstDataInterp[0], &dstDataExact[0]);
    std::cout << "abs interp error = " << absError << '\n';
    if (absError > 1.e-4) {
        return 1;
    }

    // clean up
    for (size_t j = 0; j < 2; ++j) {
        delete[] srcCoords[j];
        delete[] dstCoords[j];
    }

    std::cout << "SUCCESS\n";
    return 0;
}
