 
#include "SgFindOverlappingCells2D.h"
#include <cmath>
#include <iostream>

extern "C"
int SgFindOverlappingCells2D_new(SgFindOverlappingCells2D_type** self) {
     *self = new SgFindOverlappingCells2D_type();
     return 0;
}
      
extern "C"                   
int SgFindOverlappingCells2D_del(SgFindOverlappingCells2D_type** self) {
     delete *self;
     return 0;
}

extern "C"
int SgFindOverlappingCells2D_setSrcGrid(SgFindOverlappingCells2D_type** self,
                                        const int dims[], const int periodicity[],
                                        const double** coords) {
    (*self)->setSrcGrid(dims, periodicity, coords);
    return 0;
}

extern "C"
int SgFindOverlappingCells2D_setPolygonPoints(SgFindOverlappingCells2D_type** self, 
                                               int numPoints, const double coords[]) {
    (*self)->setPolygonPoints(numPoints, coords);
    return 0;
}

extern "C"
int SgFindOverlappingCells2D_findSrcCellIndices(SgFindOverlappingCells2D_type** self) {
    (*self)->findSrcCellIndices();
    return 0;
}

extern "C"
int SgFindOverlappingCells2D_getNumberOfSrcCellIndices(SgFindOverlappingCells2D_type** self,
                                                       int* numSrcCells) {
    *numSrcCells = (*self)->srcCellFlatInds.size();
    return 0;
}

extern "C"
int SgFindOverlappingCells2D_getSrcCellIndices(SgFindOverlappingCells2D_type** self,
                                               int** srcCellInds) {
    *srcCellInds = &(*self)->srcCellFlatInds.front();
    return 0;
}
