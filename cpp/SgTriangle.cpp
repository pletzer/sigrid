/**
 * Triangle class
 */

#include <cstring> // for size_t 
#include "SgTriangle.h"

extern "C"
int SgTriangle_new(SgTriangle_type** self,
                      const double* p0, 
                      const double* p1, 
                      const double* p2) {
                       
    *self = new SgTriangle_type();
    
    (*self)->points.push_back(std::vector<double>(p0, p0 + 3));
    (*self)->points.push_back(std::vector<double>(p1, p1 + 3));
    (*self)->points.push_back(std::vector<double>(p2, p2 + 3));
    
    std::vector<int> e(2);
    e[0] = 0; e[1] = 1;
    (*self)->edges.push_back(e);
    e[0] = 1; e[1] = 2;
    (*self)->edges.push_back(e);
    e[0] = 2; e[1] = 0;
    (*self)->edges.push_back(e);
    
    return 0;
}

extern "C"                  
int SgTriangle_del(SgTriangle_type** self) {
 
 	delete *self;
 	
 	return 0;
}

extern "C"
int SgTriangle_getNumberOfPoints(SgTriangle_type** self,
                                     int* num) {
    *num = (int) (*self)->points.size();
    return 0;
}
 
extern "C"
int SgTriangle_getNumberOfEdges(SgTriangle_type** self,
                                    int* num) {
    *num = (int) (*self)->edges.size();
    return 0;
}
 
extern "C"
int SgTriangle_getPoint(SgTriangle_type** self,
                            int i, double** p) {
    *p = &(*self)->points[(size_t) i].front();
    return 0;
}
 
int SgTriangle_getEdge(SgTriangle_type** self,
                           int i, int** pids) {
    *pids = &(*self)->edges[(size_t) i].front();
    return 0;
}

 