/**
 * Testing cconservative interpolation in 2D
 */

 #include <vector>
 #include <cassert>
 #include <cstdio>
 #include <iostream>
 #include "SgFindOverlappingCells2D.h"

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

bool testRect() {

    SgFindOverlappingCells2D_type* foc = NULL;
    SgFindOverlappingCells2D_new(&foc);

    const int nodeDims[] = {2, 3};
    int numNodes = nodeDims[0] * nodeDims[1];
    int numCells = (nodeDims[0] - 1) * (nodeDims[1] - 1);
    double* srcGrdCoords[] = {new double[numNodes], new double[numNodes]};
    const double xmins[] = {-1.0, -1.0};
    const double xmaxs[] = {1.0, 1.0};
    const int periodicity[] = {0, 0};
    createRectangularGrid(nodeDims, xmins, xmaxs, srcGrdCoords);

    SgFindOverlappingCells2D_setSrcGrid(&foc, nodeDims, periodicity,
                                        (const double**) srcGrdCoords);

    double polyCoords[] = {-1., -1.,
                           -1., 1.,
                            1., 1.,
                            1., -1};
    int numPolyPoints = 4;

    SgFindOverlappingCells2D_setPolygonPoints(&foc, numPolyPoints, polyCoords);

    SgFindOverlappingCells2D_findSrcCellIndices(&foc);
    int numSrcCellsUnderPoly;
    SgFindOverlappingCells2D_getNumberOfSrcCellIndices(&foc, &numSrcCellsUnderPoly);
    int* srcInds = NULL;
    SgFindOverlappingCells2D_getSrcCellIndices(&foc, &srcInds);

    if (numSrcCellsUnderPoly != numCells) {
        return false;
    }

    return true;
}

bool testRectQuadrant() {

    SgFindOverlappingCells2D_type* foc = NULL;
    SgFindOverlappingCells2D_new(&foc);

    const int nodeDims[] = {101, 201};
    int numNodes = nodeDims[0] * nodeDims[1];
    int numCells = (nodeDims[0] - 1) * (nodeDims[1] - 1);
    double* srcGrdCoords[] = {new double[numNodes], new double[numNodes]};
    const double xmins[] = {-1.0, -1.0};
    const double xmaxs[] = {1.0, 1.0};
    const int periodicity[] = {0, 0};
    createRectangularGrid(nodeDims, xmins, xmaxs, srcGrdCoords);

    SgFindOverlappingCells2D_setSrcGrid(&foc, nodeDims, periodicity,
                                        (const double**) srcGrdCoords);

    double polyCoords[] = {0., 0.,
                           -0., 1.,
                            1., 1.,
                            1., 0};
    int numPolyPoints = 4;

    SgFindOverlappingCells2D_setPolygonPoints(&foc, numPolyPoints, polyCoords);

    SgFindOverlappingCells2D_findSrcCellIndices(&foc);
    int numSrcCellsUnderPoly;
    SgFindOverlappingCells2D_getNumberOfSrcCellIndices(&foc, &numSrcCellsUnderPoly);
    int* srcInds = NULL;
    SgFindOverlappingCells2D_getSrcCellIndices(&foc, &srcInds);

    // because of the alignment of the grid with the polygon,
    // we expect at least numCells / 4 overlapping cells
    if (numSrcCellsUnderPoly < numCells / 4) {
        return false;
    }

    return true;
}

bool testPolyIsOutsideGrid() {

    SgFindOverlappingCells2D_type* foc = NULL;
    SgFindOverlappingCells2D_new(&foc);

    const int nodeDims[] = {101, 201};
    int numNodes = nodeDims[0] * nodeDims[1];
    int numCells = (nodeDims[0] - 1) * (nodeDims[1] - 1);
    double* srcGrdCoords[] = {new double[numNodes], new double[numNodes]};
    const double xmins[] = {0.0, 0.0};
    const double xmaxs[] = {1.0, 1.0};
    const int periodicity[] = {0, 0};
    createRectangularGrid(nodeDims, xmins, xmaxs, srcGrdCoords);

    SgFindOverlappingCells2D_setSrcGrid(&foc, nodeDims, periodicity,
                                        (const double**) srcGrdCoords);

    double polyCoords[] = {-1., -1.,
                           -0.5, -0.5};
    int numPolyPoints = 2;

    SgFindOverlappingCells2D_setPolygonPoints(&foc, numPolyPoints, polyCoords);

    SgFindOverlappingCells2D_findSrcCellIndices(&foc);
    int numSrcCellsUnderPoly;
    SgFindOverlappingCells2D_getNumberOfSrcCellIndices(&foc, &numSrcCellsUnderPoly);
    std::cout << "testPolyIsOutsideGrid: number of overlap cells = " << numSrcCellsUnderPoly << '\n';
    int* srcInds = NULL;
    SgFindOverlappingCells2D_getSrcCellIndices(&foc, &srcInds);
    for (int i = 0; i < numSrcCellsUnderPoly; ++i) {
        std::cout << "src cell: " << srcInds[i] << '\n';
    }

    // no overlap
    if (numSrcCellsUnderPoly != 0) {
        // error
        return false;
    }

    return true;
}


int main(int argc, char** argv) {

    if (!testRect()) return 1;
    if (!testRectQuadrant()) return 1;
    if (!testPolyIsOutsideGrid()) return 2;

    std::cout << "SUCCESS\n";
    return 0;
}
