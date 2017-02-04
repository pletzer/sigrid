#include <cmath>
#include <cstddef>

#ifndef CREATE_GRIDS_H
#define CREATE_GRIDS_H

#ifdef __cplusplus
extern "C" {
#endif


void matrixDotMatrix(const double mat1[], const double mat2[], double matres[]);

void matrixDotVector(const double mat[], const double vec[], double vecres[]);


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
                           double** coords);

/**
 * Create polar grid
 * @param nodeDims number of nodes in the two directions
 * @param center center of grid
 * @param radius radius
 * @param coords coordinates (output)
 */
void createPolarGrid(const int nodeDims[], 
                     const double center[],
                     double radius,
                     double** coords);

/** 
 * Create rotated pole grid
 * @param nodeDims number of nodes in the two directions
 * @param delta_lat pole displacement in latitude (deg)
 * @param delta_lon pole displacement in longitude (deg)
 * @param coords coordinates (output)
 */
void createRotatedPoleGrid(const int nodeDims[], 
                           double delta_lat, double delta_lon,
                           double** coords);

#ifdef __cplusplus
}
#endif


#endif // CREATE_GRIDS_H
