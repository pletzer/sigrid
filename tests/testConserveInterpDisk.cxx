/**
 * Testing cconservative interpolation in 2D
 */

 #include <vector>
 #include <cassert>
 #include <cstdio>
 #include <iostream>
 #include "SgConserveInterp2D.h"

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


bool test() {

    // source grid is disk
    const int srcDims[] = {5, 33};
    int srcNumPoints = srcDims[0] * srcDims[1];
    double* srcCoords[] = {new double[srcNumPoints], new double[srcNumPoints]};
    const double radius = 1.0;
    createPolarGrid(srcDims, radius, srcCoords);

    // destination grid is rectangle encompassing src grid
    const int dstDims[] = {401, 401};
    const double dstXmins[] = {-1., -1.};
    const double dstXmaxs[] = {1., 1.};
    int dstNumPoints = dstDims[0] * dstDims[1];
    int dstNumCells = (dstDims[0] - 1) * (dstDims[1] - 1);
    double* dstCoords[] = {new double[dstNumPoints], new double[dstNumPoints]};
    createRectangularGrid(dstDims, dstXmins, dstXmaxs, dstCoords);    

    SgConserveInterp2D_type* interp = NULL;
    SgConserveInterp2D_new(&interp);
    SgConserveInterp2D_setSrcGrid(&interp, srcDims, (const double**) srcCoords);
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

    std::cout << "Total weight: " << totalWeight << " totalWeight/dstNumCells = " << 
        (double) totalWeight / (double) dstNumCells << " ~ pi/4 = " << M_PI/4.0 << '\n';
    // expect partial cells
    if (totalWeight >= dstNumCells - 1.e-8) {
        // error
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

    if (!test()) return 1;

    std::cout << "SUCCESS\n";
    return 0;
}
