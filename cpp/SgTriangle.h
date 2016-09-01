/**
 * Triangle class
 */
 
 #ifndef SG_TRIANGLE_H
 #define SG_TRIANGLE_H
 
 #include <vector>
 
 struct SgTriangle_type {
 	std::vector< std::vector<double> > points;
 	std::vector< std::vector<int> > edges;
 };
 
#ifdef __cplusplus
extern "C" {
#endif

 int SgTriangle_new(SgTriangle_type** self,
                       const double* p0, 
                       const double* p1, 
                       const double* p2);
                       
 int SgTriangle_del(SgTriangle_type** self);
 
 int SgTriangle_getNumberOfPoints(SgTriangle_type** self,
                                  int* num);
 
 int SgTriangle_getNumberOfEdges(SgTriangle_type** self,
                                 int* num);
 
 int SgTriangle_getPoint(SgTriangle_type** self,
                            int i, double** p);
 
 int SgTriangle_getEdge(SgTriangle_type** self,
                           int i, int** pids);
 
#ifdef __cplusplus
}
#endif

#endif // SG_TRIANGLE_H