/**
 * Testing cconservative interpolation in 2D
 */

#include <vector>
#include <cassert>
#include <cstdio>
#include <iostream>
#include "SgTripolarGrid.h"

int main(int argc, char** argv) {

    const int dims[] = {11, 21};
    int capLatIndex = 6;
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

    std::cout << "SUCCESS\n";
    return 0;
}
