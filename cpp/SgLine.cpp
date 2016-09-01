/**
 * Line class
 */

#include <cstring> // for size_t 
#include "SgLine.h"

extern "C"
int SgLine_new(SgLine_type** self,
                      const double* p0, 
                      const double* p1) {
                       
    *self = new SgLine_type();
    
    (*self)->points.push_back(std::vector<double>(p0, p0 + 3));
    (*self)->points.push_back(std::vector<double>(p1, p1 + 3));
    
    return 0;
}

extern "C"                  
int SgLine_del(SgLine_type** self) {
 
 	delete *self;
 	
 	return 0;
}

extern "C"
int SgLine_getNumberOfPoints(SgLine_type** self,
                                     int* num) {
    *num = (int) (*self)->points.size();
    return 0;
}
 
extern "C"
int SgLine_getPoint(SgLine_type** self,
                            int i, double** p) {
    *p = &(*self)->points[(size_t) i].front();
    return 0;
}
 