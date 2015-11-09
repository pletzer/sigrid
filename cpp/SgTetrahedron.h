/**
 * Tetrahedron class
 */
 
 #ifndef SG_TETRAHEDRON_H
 #define SG_TETRAHEDRON_H
 
 #include <vector>
 
 struct SgTetrahedron_type {
 	std::vector< std::vector<double> > points;
 	std::vector< std::vector<int> > edges;
 	std::vector< std::vector<int> > faces;
 };
 
#ifdef __cplusplus
extern "C" {
#endif

 int SgTetrahedron_new(SgTetrahedron_type** self,
                       const double* p0, 
                       const double* p1, 
                       const double* p2,
                       const double* p3);
                       
 int SgTetrahedron_del(SgTetrahedron_type** self);
 
 int SgTetrahedron_getNumberOfPoints(SgTetrahedron_type** self,
                                     int* num);
 
 int SgTetrahedron_getNumberOfEdges(SgTetrahedron_type** self,
                                    int* num);
 
 int SgTetrahedron_getNumberOfFaces(SgTetrahedron_type** self,
                                    int* num);
 
 int SgTetrahedron_getPoint(SgTetrahedron_type** self,
                            int i, double** p);
 
 int SgTetrahedron_getEdge(SgTetrahedron_type** self,
                           int i, int** pids);
 
 int SgTetrahedron_getFace(SgTetrahedron_type** self,
                           int i, int** pids);
 
#ifdef __cplusplus
}
#endif

#endif // SG_TETRAHEDRON_H