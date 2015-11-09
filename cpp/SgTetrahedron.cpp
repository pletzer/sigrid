/**
 * Tetrahedron class
 */
 
#include "SgTetrahedron.h"

extern "C"
int SgTetrahedron_new(SgTetrahedron_type** self,
                      const double* p0, 
                      const double* p1, 
                      const double* p2,
                      const double* p3) {
                       
    *self = new SgTetrahedron_type();
    
    (*self)->points.push_back(std::vector<double>(p0, p0 + 3));
    (*self)->points.push_back(std::vector<double>(p1, p1 + 3));
    (*self)->points.push_back(std::vector<double>(p2, p2 + 3));
    (*self)->points.push_back(std::vector<double>(p3, p3 + 3));
    
    std::vector<int> e(2);
    e[0] = 0; e[1] = 1;
    (*self)->edges.push_back(e);
    e[0] = 1; e[1] = 2;
    (*self)->edges.push_back(e);
    e[0] = 2; e[1] = 0;
    (*self)->edges.push_back(e);
    e[0] = 0; e[1] = 3;
    (*self)->edges.push_back(e);
    e[0] = 1; e[1] = 3;
    (*self)->edges.push_back(e);
    e[0] = 2; e[1] = 3;
    (*self)->edges.push_back(e);
    
    std::vector<int> f(3);
    f[0] = 0; f[1] = 1; f[2] = 2;
    (*self)->faces.push_back(f);
    f[0] = 0; f[1] = 3; f[2] = 1;
    (*self)->faces.push_back(f);
    f[0] = 3; f[1] = 1; f[2] = 2;
    (*self)->faces.push_back(f);
    f[0] = 0; f[1] = 2; f[2] = 3;
    (*self)->faces.push_back(f);
    
    return 0;
}

extern "C"                  
int SgTetrahedron_del(SgTetrahedron_type** self) {
 
 	delete *self;
 	
 	return 0;
}

extern "C"
int SgTetrahedron_getNumberOfPoints(SgTetrahedron_type** self,
                                     int* num) {
    *num = (int) (*self)->points.size();
    return 0;
}
 
extern "C"
int SgTetrahedron_getNumberOfEdges(SgTetrahedron_type** self,
                                    int* num) {
    *num = (int) (*self)->edges.size();
    return 0;
}
 
extern "C"
int SgTetrahedron_getNumberOfFaces(SgTetrahedron_type** self,
                                    int* num) {
    *num = (int) (*self)->faces.size();
    return 0;
}
 
extern "C"
int SgTetrahedron_getPoint(SgTetrahedron_type** self,
                            int i, double** p) {
    *p = &(*self)->points[(size_t) i].front();
    return 0;
}
 
int SgTetrahedron_getEdge(SgTetrahedron_type** self,
                           int i, int** pids) {
    *pids = &(*self)->edges[(size_t) i].front();
    return 0;
}
 
int SgTetrahedron_getFace(SgTetrahedron_type** self,
                           int i, int** pids) {
    *pids = &(*self)->edges[(size_t) i].front();
    return 0;
}

 