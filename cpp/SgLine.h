/**
 * Line class
 */
 
 #ifndef SG_LINE_H
 #define SG_LINE_H
 
 #include <vector>
 
 struct SgLine_type {
 	std::vector< std::vector<double> > points;
 };
 
#ifdef __cplusplus
extern "C" {
#endif

 int SgLine_new(SgLine_type** self,
                       const double* p0, 
                       const double* p1);
                       
 int SgLine_del(SgLine_type** self);
 
 int SgLine_getNumberOfPoints(SgLine_type** self,
                                  int* num);
 int SgLine_getPoint(SgLine_type** self,
                            int i, double** p);

#ifdef __cplusplus
}
#endif

#endif // SG_LINE_H