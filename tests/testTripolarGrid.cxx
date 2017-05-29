/**
 * Testing cconservative interpolation in 2D
 */

#include <vector>
#include <cassert>
#include <fstream>
#include <iostream>
#include <string>
#include <cmath>
#include "SgTripolarGrid.h"

void saveGrid(const std::string& filename, const int dims[], double** coords) {
    int ntot = dims[0] * dims[1];
    std::ofstream f;
    f.open(filename.c_str());
    f << "# vtk DataFile Version 2.0\n";
    f << "testTripolarGrid\n";
    f << "ASCII\n";
    f << "DATASET STRUCTURED_GRID\n";
    f << "DIMENSIONS " << dims[0] << " " << dims[1] << " 1\n";
    f << "POINTS " << ntot << " float\n";
    size_t k = 0;
    for (int i = 0; i < ntot; ++i) {
        double lat = coords[0][k];
        double lon = coords[1][k];
        double x = cos(lat*M_PI/180.) * cos(lon*M_PI/180.);
        double y = cos(lat*M_PI/180.) * sin(lon*M_PI/180.);
        double z = sin(lat*M_PI/180.);
        f << x << " " << y << " " << z << '\n';
        k++;
    }
    f << "\n";

    f.close();
}

int main(int argc, char** argv) {

    const int dims[] = {101, 201};
    int capLatIndex = int(0.65 * dims[0]);
    SgTripolarGrid_type tri(dims, capLatIndex);
    int outdims[2];
    double** coords;
    tri.getGrid(outdims, &coords);
    size_t k = 0;
    for (int j = 0; j < outdims[0]; ++j) {
        for (int i = 0; i < outdims[1]; ++i) {
            std::cout << "j = " << j << " i = " << i << 
                      " lat = " << coords[0][k] << " lon = " << coords[1][k] << '\n';
            k++;
        }
    }
    saveGrid("testTripolarGrid.vtk", outdims, coords);

    std::cout << "SUCCESS\n";
    return 0;
}
