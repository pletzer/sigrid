#include "createGrids2D.h"

extern "C"
void matrixDotMatrix(const double mat1[], const double mat2[], double matres[]) {
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            size_t k = i*3 + j;
            matres[k] = 0.;
            for (size_t el = 0; el < 3; ++el) {
                matres[k] += mat1[i*3 + el] * mat2[el*3 + j];
            }
        }
    }
}

extern "C"
void matrixDotVector(const double mat[], const double vec[], double vecres[]) {
    for (size_t i = 0; i < 3; ++i) {
        vecres[i] = 0.;
        for (size_t el = 0; el < 3; ++el) {
            vecres[i] += mat[i*3 + el] * vec[el];
        }
    }
}

extern "C"
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

extern "C"
void createPolarGrid(const int nodeDims[], 
                     const double center[],
                     double radius,
                     double** coords) {

    double deltas[] = {radius / double(nodeDims[0] - 1),
                       2 * M_PI / double(nodeDims[1] - 1)};
    size_t index = 0;
    for (size_t i = 0; i < nodeDims[0]; ++i) {
        for (size_t j = 0; j < nodeDims[1]; ++j) {
            double rho = center[0] + deltas[0]*i;
            double the = center[1] + deltas[1]*j;
            coords[0][index] = rho*cos(the);
            coords[1][index] = rho*sin(the);
            index++;
        }
    }
}

extern "C"
void createRotatedPoleGrid(const int nodeDims[], 
                           double delta_lat, double delta_lon,
                           double** coords) {
    int nj = nodeDims[0];
    int ni = nodeDims[1];
    double alpha = M_PI * delta_lat / 180.;
    double beta = M_PI * delta_lon / 180.;

    double cos_alp = cos(alpha);
    double sin_alp = sin(alpha);
    double cos_bet = cos(beta);
    double sin_bet = sin(beta);

    // rotate in latitude
    double rot_alp[] = {cos_alp,     0., sin_alp,
                             0.,     1.,      0.,
                       -sin_alp,     0., cos_alp};

    // rotate in longitude
    double rot_bet[] = {cos_bet, sin_bet, 0.,
                       -sin_bet, cos_bet, 0., 
                             0.,      0., 1.};

    // total rotation
    double transfMatrix[3*3];

    matrixDotMatrix(rot_bet, rot_alp, transfMatrix);

    // coordinates before rotation
    double xyzPrime[3];

    // coordinates after rotation
    double xyz[3];

    size_t k = 0;
    for (size_t j = 0; j < nj; ++j) {
        double lat = -90.0 + 180.0 * j / double(nj - 1);
        double the = M_PI * lat / 180.0;
        double cos_the = cos(the);
        double sin_the = sin(the);
        double rho = cos_the;
        for (size_t i = 0; i < ni; ++i) {
            double lon = -180.0 + 360.0 * i / double(ni - 1);
            double lam = M_PI * lon / 180.0;
            double cos_lam = cos(lam);
            double sin_lam = sin(lam);

            xyzPrime[0] = rho * cos_lam;
            xyzPrime[1] = rho * sin_lam;
            xyzPrime[2] = sin_the;

            matrixDotVector(transfMatrix, xyzPrime, xyz);

            coords[0][k] = 180. * asin(xyz[2]) / M_PI;
            coords[1][k] = 180. * atan2(xyz[1], xyz[0]) / M_PI;
            k++;
        }
    }
}

