/**
 * Linear solve class
 */
 
#include "SgLinearSolve.h"
#include <cstdlib> // std::max
#include <cmath> // fabs

extern "C"
int SgLinearSolve_new(SgLinearSolve_type** self,
                      int nrow, int ncol) {
    
    *self = new SgLinearSolve_type();
    (*self)->nrow = nrow;
    (*self)->ncol = ncol;
    int ldb = nrow > ncol? nrow: ncol;
    int mn = nrow < ncol? nrow: ncol;
    (*self)->mat.resize(nrow * ncol);
    (*self)->matOri.resize(nrow * ncol);
    (*self)->b.resize(nrow);
    (*self)->x.resize(ldb);
    int nb = 1; // optimal block size
    int nrhs = 1;
  	(*self)->lwork = std::max(1, mn + std::max(mn, nrhs)*nb);
  	(*self)->work.resize((size_t) (*self)->lwork);
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
                            
// row major storage
	for (size_t i = 0; i < (*self)->nrow; ++i) {
		for (size_t j = 0; j < (*self)->ncol; ++j) {
      // Fortran storage
			size_t k = i*(*self)->ncol + j;
			(*self)->mat[k] = mat[k];
			(*self)->matOri[k] = mat[k];
		}
	}
	return 0;
}

extern "C"                          
int SgLinearSolve_setRightHandSide(SgLinearSolve_type** self,
                            	   const double b[]) {
  size_t n = (*self)->b.size();
	for (size_t i = 0; i < (*self)->nrow; ++i) {
		(*self)->x[i] = b[i];
		(*self)->b[i] = b[i];
	} 
	return 0;   
}

extern "C"                                    
int SgLinearSolve_solve(SgLinearSolve_type** self) {

    int errCode = 0;
  	char t = 'T';
  	int nrhs = 1;
  	int mn = (int) (*self)->b.size();
  	int ldb = (int) (*self)->b.size();
  	_GELS_(&t, 
  	       &(*self)->nrow, 
  	       &(*self)->ncol, 
  	       &nrhs,
     	     &(*self)->mat.front(), 
     	     &(*self)->nrow,
     	     &(*self)->x.front(), &ldb,
     	     &(*self)->work.front(), &(*self)->lwork, 
     	     &errCode);

// if errCode == -i, the i-th argument had an illegal value
// if errCode == i, the i-th diagonal element of the
//                  triangular factor of A is zero, so that 
//                  A does not have full rank; the least 
//                  squares solution could not be computed. 
    return errCode;
}

extern "C"                                    
int SgLinearSolve_getSolution(SgLinearSolve_type** self, 
                              double** sol) {
    *sol = &(*self)->x.front();
    return 0;
}

extern "C"
int SgLinearSolve_getResidual(SgLinearSolve_type** self,
                              double* res) {
    *res = 0;
    for (size_t i = 0; i < (*self)->nrow; ++i) {
      double val = 0;
    	for (size_t j = 0; j < (*self)->ncol; ++j) {
    	  size_t k = i*(*self)->ncol + j;
    		val += (*self)->matOri[k]*(*self)->x[j];
    	}
    	val -= (*self)->b[i];
    	*res += fabs(val);
    }
    return 0;
}

