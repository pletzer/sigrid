/**
 * A class that computes the interpolation weights for a nodal field
 */
 
#ifndef SG_SIDE_ITERATOR_H
#define SG_SIDE_ITERATOR_H
 
#include <iostream>
#include <vector>
#include <cmath>
 
struct SgSideIterator_type {

    size_t index;
    size_t numElems;
    std::vector<size_t> prodDims;

    /**
     * Constructor
     */
    SgSideIterator_type(size_t numDims) {
        this->prodDims.resize(numDims);
        this->prodDims[numDims - 1] = 1;
        for (int i = numDims - 2; i >= 0; --i) {
            this->prodDims[i] = this->prodDims[i + 1]*3;
        }
        this->index = 0;
        this->numElems = std::pow(3, numDims);
    }

    /**
     * Destructor (default)
     */
    ~SgSideIterator_type() {}

    /** 
     * Get the number of elements
     * return number
     */
    size_t getNumberOfElements() const {
        return this->numElems;
    }

    /**
     * Reset the iterator
     */
    void reset() {
        this->index = 0;
    }

    /** 
     * Advance the iterator
     * @return 0 if not finished, 1 if finished
     */
    int next() {
        this->index++;
        if (this->index < this->numElems) {
            return 0;
        }
        return 1;
    }

    /** 
     * Get the current side
     * @return array of -1 (= low side), 0 (= middle), and 1 (= high side)
     */
    std::vector<int> get() const {
        size_t ndims = this->prodDims.size();
        std::vector<int> res(ndims, 0);
        for (size_t i = 0; i < ndims; ++i) {
            res[i] = -1 + index / this->prodDims[i] % 3;
        }
        return res;
    }
};

#endif // SG_SIDE_ITERATOR_H
