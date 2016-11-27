/**
 * Testing cconservative interpolation in 2D
 */

#include <vector>
#include <cassert>
#include <cstdio>
#include <iostream>
#include "SgFindOverlappingCells2D.h"
#include "createGrids2D.h"

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
