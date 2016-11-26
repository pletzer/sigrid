/**
 * Testing cconservative interpolation in 2D
 */

 #include <vector>
 #include <cassert>
 #include <cstdio>
 #include <iostream>
 #include "SgNodeInterp2D.h"

/**
 * Create rectangular grid
 * @param nodeDims number of nodes in the two directions
 * @param xmins low corner point of the grid
 * @param xmaxs high corner point of the grid
 * @param coords coordinates (output)
 */
void createRectangularGrid(const int nodeDims[], 
                           const double xmins[],
                           const double xmaxs[], 
                           double** coords) {

    double deltas[] = {(xmaxs[0] - xmins[0])/double(nodeDims[0] - 1),
                       (xmaxs[1] - xmins[1])/double(nodeDims[1] - 1)};
    size_t index = 0;
    for (size_t i = 0; i < nodeDims[0]; ++i) {
        for (size_t j = 0; j < nodeDims[1]; ++j) {
            coords[0][index] = xmins[0] + deltas[0]*i;
            coords[1][index] = xmins[1] + deltas[1]*j;
            index++;
        }
    }
}

/**
 * Create polar grid
 * @param nodeDims number of nodes in the two directions
 * @param radius radius
 * @param coords coordinates (output)
 */
void createPolarGrid(const int nodeDims[], 
                     double radius,
                     double** coords) {

    double deltas[] = {radius / double(nodeDims[0] - 1),
                       2 * M_PI / double(nodeDims[1] - 1)};
    size_t index = 0;
    for (size_t i = 0; i < nodeDims[0]; ++i) {
        for (size_t j = 0; j < nodeDims[1]; ++j) {
            double rho = 0.0 + deltas[0]*i;
            double the = 0.0 + deltas[1]*j;
            coords[0][index] = rho*cos(the);
            coords[1][index] = rho*sin(the);
            index++;
        }
    }
}

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

int main(int argc, char** argv) {

    if (!testSimple())
    std::cout << "SUCCESS\n";
    return 0;
}
