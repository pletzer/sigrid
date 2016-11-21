/**
 * Linear solve class
 */
 
#include "SgLinearSolve.h"

extern "C"
int SgLinearSolve_new(SgLinearSolve_type** self,
                      int nrow, int ncol) {
    *self = new SgLinearSolve_type(nrow, ncol);
    return 0;
}

extern "C"                    
int SgLinearSolve_del(SgLinearSolve_type** self) {
    delete *self;
    return 0;
}

extern "C"
int SgLinearSolve_setMatrix(SgLinearSolve_type** self,
                            const double mat[]) {
    (*self)->setMatrix(mat);
    return 0;
}

extern "C"                          
int SgLinearSolve_setRightHandSide(SgLinearSolve_type** self,
                                   const double b[]) {
    (*self)->setRightHandSide(b);
    return 0;   
}

extern "C"                                    
int SgLinearSolve_solve(SgLinearSolve_type** self) {
  return (*self)->solve();
}

extern "C"                                    
int SgLinearSolve_getSolution(SgLinearSolve_type** self, 
                              double** sol) {
    (*self)->getSolution(sol);
    return 0;
}

extern "C"
int SgLinearSolve_getResidual(SgLinearSolve_type** self,
                              double* res) {
    *res = (*self)->getResidual();
    return 0;
}

