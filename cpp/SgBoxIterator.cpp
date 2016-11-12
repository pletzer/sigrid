/**
 * Box iterator class
 */
 
 #include "SgNdims.h"
 #include "SgBoxIterator.h"
 #include <iostream>

 extern "C"
 int SgBoxIterator_new(SgBoxIterator_type** self,
                       const int* begInds, 
                       const int* endInds) {

  *self = new SgBoxIterator_type();

  (*self)->lowEnd.resize(NDIMS_TOPO);
  (*self)->dims.resize(NDIMS_TOPO);
  (*self)->dimProd.resize(NDIMS_TOPO);
  (*self)->numElems = 1;
  for (size_t i = 0; i < NDIMS_TOPO; ++i) {
    (*self)->lowEnd[i] = begInds[i];
    (*self)->dims[i] = endInds[i] - begInds[i];
    (*self)->dimProd[i] = 1;
    (*self)->numElems *= (*self)->dims[i];
  }
  for (int i = (int) NDIMS_TOPO - 2; i >= 0; --i) {
    // last index varies fastest
    (*self)->dimProd[i] =  (*self)->dimProd[i + 1] * (*self)->dims[i + 1];
  }
  return 0;
 }
                       
 extern "C"
 int SgBoxIterator_del(SgBoxIterator_type** self) {

  delete *self;
  return 0;
 }
 
 extern "C"
 int SgBoxIterator_getNumberOfElements(SgBoxIterator_type** self,
                                     int* num) {
 	*num = (*self)->numElems;
 	return 0;
 }
 
 extern "C"
 int SgBoxIterator_getElement(SgBoxIterator_type** self,
                              int index,
                              int inds[]) {
  for (size_t i = 0; i < NDIMS_TOPO; ++i) {
    inds[i] = (*self)->lowEnd[i] + 
        index / (*self)->dimProd[i] % (*self)->dims[i];
  }
 	return 0;
 }
 
