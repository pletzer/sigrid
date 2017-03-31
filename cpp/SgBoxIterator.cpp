/**
 * Box iterator class
 */
 
#include "SgBoxIterator.h"

extern "C"
int SgBoxIterator_new(SgBoxIterator_type** self,
                      int ndims,
                      const int* begInds, 
                      const int* endInds) {
  *self = new SgBoxIterator_type(ndims, begInds, endInds);
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
 	*num = (*self)->getNumberOfElements();
 	return 0;
}
 
extern "C"
int SgBoxIterator_getElement(SgBoxIterator_type** self,
                              int index,
                              int inds[]) {
  (*self)->getElement(index, inds);
 	return 0;
}
 
extern "C"
int SgBoxIterator_getIndex(SgBoxIterator_type** self,
                           const int inds[], int* index) {
  *index = (*self)->getIndex(inds);
  return 0;
}
