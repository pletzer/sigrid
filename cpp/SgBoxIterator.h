/**
 * Box iterator class
 */
 
 #ifndef SG_BOX_ITERATOR_H
 #define SG_BOX_ITERATOR_H
 
 #include <vector>
 
 struct SgBoxIterator_type {
 	std::vector<int> lowEnd;
 	std::vector<int> dims;
 	std::vector<int> dimProd;
  int numElems;
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