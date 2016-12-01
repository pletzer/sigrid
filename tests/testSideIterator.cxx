/**
 * Testing side iterator
 */

#include <iostream>
#include "SgSideIterator.h"

void test2D() {

    size_t ndims = 2;
    SgSideIterator_type sit(ndims);
    sit.reset();
    int stop = 0;
    while (!stop) {
        std::vector<int> side = sit.get();
        for (size_t i = 0; i < ndims; ++i) {
            std::cout << side[i] << ", ";
        }
        std::cout << '\n';
        stop = sit.next();
    }
}

void test3D() {

    size_t ndims = 3;
    SgSideIterator_type sit(ndims);
    sit.reset();
    int stop = 0;
    while (!stop) {
        std::vector<int> side = sit.get();
        for (size_t i = 0; i < ndims; ++i) {
            std::cout << side[i] << ", ";
        }
        std::cout << '\n';
        stop = sit.next();
    }
}

int main(int argc, char** argv) {

    test2D();
    test3D();

    std::cout << "SUCCESS\n";
    return 0;
}
