/**
 * Box iterator class
 */
 
#ifndef SG_BOX_ITERATOR_H
#define SG_BOX_ITERATOR_H
 
#include <vector>
#include "SgNdims.h"
#include <iostream>
 
struct SgBoxIterator_type {

    // low end index set
    std::vector<int> lowEnd;

    // number of cells/elements along each axis
    std::vector<int> dims;

    // used to compute a flat index
 	std::vector<int> dimProd;

    // number of elements or cells
    int numElems;

    /**
     * Constructor
     * @param ndims number of axes
     * @param begInds begin index set 
     * @param endInds one past end index set
     */
    SgBoxIterator_type(int ndims, const int* begInds, const int* endInds) {

        this->lowEnd.resize(ndims);
        this->dims.resize(ndims);
        this->dimProd.resize(ndims);
        this->numElems = 1;
        for (int i = 0; i < ndims; ++i) {
            this->lowEnd[i] = begInds[i];
            this->dims[i] = endInds[i] - begInds[i];
            this->dimProd[i] = 1;
            this->numElems *= this->dims[i];
        }

        for (int i = ndims - 2; i >= 0; --i) {
            // last index varies fastest
            this->dimProd[i] =  this->dimProd[i + 1] * this->dims[i + 1];
        }
    }

    /**
     * Destructor
     */
    ~SgBoxIterator_type() {}

    /**
     * Get the number of elements
     * @return number
     */
    int getNumberOfElements() const {
        return this->numElems;
    }

    /**
     * Get the index set element
     * @param index flat index 
     * @param index set (output)
     */
    void getElement(int index, int inds[]) {
        for (size_t i = 0; i < this->dims.size(); ++i) {
            inds[i] = this->lowEnd[i] + 
                index / this->dimProd[i] % this->dims[i];
        }
    }
};
 
#ifdef __cplusplus
extern "C" {
#endif

int SgBoxIterator_new(SgBoxIterator_type** self,
                      int ndims,
                      const int* begInds, 
                      const int* endInds);
                       
int SgBoxIterator_del(SgBoxIterator_type** self);
 
int SgBoxIterator_getNumberOfElements(SgBoxIterator_type** self,
                                      int* num);
 
int SgBoxIterator_getElement(SgBoxIterator_type** self,
                             int index,
                             int inds[]);

#ifdef __cplusplus
}
#endif

#endif // SG_BOX_ITERATOR_H