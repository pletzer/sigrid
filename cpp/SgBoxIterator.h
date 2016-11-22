/**
 * Box iterator class
 */
 
#ifndef SG_BOX_ITERATOR_H
#define SG_BOX_ITERATOR_H
 
#include <vector>
#include "SgNdims.h"
#include <iostream>
 
struct SgBoxIterator_type {
    std::vector<int> lowEnd;
    std::vector<int> dims;
 	std::vector<int> dimProd;
    int numElems;

    SgBoxIterator_type(const int* begInds, const int* endInds) {
        this->lowEnd.resize(NDIMS_TOPO);
        this->dims.resize(NDIMS_TOPO);
        this->dimProd.resize(NDIMS_TOPO);
        this->numElems = 1;
        for (size_t i = 0; i < NDIMS_TOPO; ++i) {
            this->lowEnd[i] = begInds[i];
            this->dims[i] = endInds[i] - begInds[i];
            this->dimProd[i] = 1;
            this->numElems *= this->dims[i];
        }
        for (int i = (int) NDIMS_TOPO - 2; i >= 0; --i) {
            // last index varies fastest
            this->dimProd[i] =  this->dimProd[i + 1] * this->dims[i + 1];
        }
    }

    ~SgBoxIterator_type() {}

    int getNumberOfElements() const {
        return this->numElems;
    }

    void getElement(int index, int inds[]) {
        for (size_t i = 0; i < NDIMS_TOPO; ++i) {
            inds[i] = this->lowEnd[i] + 
                index / this->dimProd[i] % this->dims[i];
        }
    }
};
 
#ifdef __cplusplus
extern "C" {
#endif

int SgBoxIterator_new(SgBoxIterator_type** self,
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