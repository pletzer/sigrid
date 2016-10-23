/**
 * A class that finds the location of a point in index space
 */
 
 #ifndef SG_FIND_POINT_IN_CELL_H
 #define SG_FIND_POINT_IN_CELL_H
 
 #include <vector>
 #include <cstdio> // size_t
 #include "SgLinearSolve.h"
 
 struct SgFindPointInCell_type {
 	SgLinearSolve_type* slvr;
 	double tolpos;
 	std::vector<double> dIndices;
 	std::vector<double> targetPoint;
 	// the curvilinear coordinate points as ndims flat arrays
 	std::vector<std::vector<double> > coords;
 	std::vector<size_t> prodDims;
 	std::vector<size_t> dims;
 	int nitermax;
 	int iter;
 };
 
#ifdef __cplusplus
extern "C" {
#endif

 int SgFindPointInCell_new(SgFindPointInCell_type** self,
                       int nitermax, double tolpos);
                       
 int SgFindPointInCell_del(SgFindPointInCell_type** self);

 int SgFindPointInCell_setGrid(SgFindPointInCell_type** self, 
 	                           int ndims, const int dims[], 
 	                           const double** coords);

 int SgFindPointInCell_getPosition(SgFindPointInCell_type** self, 
 	                               double pos[]);

 int SgFindPointInCell_next(SgFindPointInCell_type** self);

 int SgFindPointInCell_setIndices(SgFindPointInCell_type** self,
 	                              double dIndices[]);

 int SgFindPointInCell_getIndices(SgFindPointInCell_type** self,
 	                              double dIndices[]);
 

#ifdef __cplusplus
}
#endif

#endif // SG_FIND_POINT_IN_CELL_H