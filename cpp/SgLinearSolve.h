/**
 * Linear solve class
 */
 
#ifndef SG_LINEAR_SOLVE_H
#define SG_LINEAR_SOLVE_H

#include <vector>

// Fortran name mangling
#define FC_GLOBAL(name,NAME) name##_
#define _GELS_ FC_GLOBAL(dgels,DGELS)

// least square solves
extern "C"
void _GELS_(char*,
            int*, int*, int*,
            double*, int*,
            double*, int*,
            double*, int*,
            int *);


struct SgLinearSolve_type {

  // the modified matrix
	std::vector<double> mat;

  // a copy of the matrix
	std::vector<double> matOri;

  // the right-hand side
	std::vector<double> b;

  // the solution vector
	std::vector<double> x;

  // additional work space
	std::vector<double> work;
	int lwork;
	int nrow;
	int ncol;
};

#ifdef __cplusplus
extern "C" {
#endif

int SgLinearSolve_new(SgLinearSolve_type** self,
                      int nrow, int ncol);
                      
int SgLinearSolve_del(SgLinearSolve_type** self);

int SgLinearSolve_setMatrix(SgLinearSolve_type** self,
                            const double mat[]);
                            
int SgLinearSolve_setRightHandSide(SgLinearSolve_type** self,
                            	     const double b[]);
                            	   
int SgLinearSolve_solve(SgLinearSolve_type** self);
                                     
int SgLinearSolve_getSolution(SgLinearSolve_type** self, 
                              double** sol);

int SgLinearSolve_getResidual(SgLinearSolve_type** self, 
                              double* res);

#ifdef __cplusplus
}
#endif

#endif // SG_LINEAR_SOLVE_H
