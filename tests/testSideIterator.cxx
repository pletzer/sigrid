/**
 * Testing side iterator
 */

#include <iostream>
#include "SgSideIterator.h"

void test2D() {

    size_t ndims = 2;
    SgSideIterator_type sit(ndims);
    sit.reset();
    size_t numElems = sit.getNumberOfElements();
    for (size_t i = 0; i < numElems; ++i) {
        std::vector<int> side = sit.get();
        for (size_t j = 0; j < ndims; ++j) {
            std::cout << side[j] << ", ";
        }
        std::cout << '\n';
        sit.next();
    }
}

void test3D() {

    size_t ndims = 3;
    SgSideIterator_type sit(ndims);
    sit.reset();
    size_t numElems = sit.getNumberOfElements();
    for (size_t i = 0; i < numElems; ++i) {
        std::vector<int> side = sit.get();
        for (size_t j = 0; j < ndims; ++j) {
            std::cout << side[j] << ", ";
        }
        std::cout << '\n';
        sit.next();
    }
}

int main(int argc, char** argv) {

    test2D();
    test3D();

    std::cout << "SUCCESS\n";
    return 0;
}
