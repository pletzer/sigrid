/**
 * Linear solve class
 */
 
#ifndef SG_LINEAR_SOLVE_H
#define SG_LINEAR_SOLVE_H

#include <vector>
#include <cstdlib> // std::max
#include <cmath> // fabs

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

  SgLinearSolve_type(int nrow, int ncol) {
    this->nrow = nrow;
    this->ncol = ncol;
    int ldb = nrow > ncol? nrow: ncol;
    int mn = nrow < ncol? nrow: ncol;
    this->mat.resize(nrow * ncol);
    this->matOri.resize(nrow * ncol);
    this->b.resize(nrow);
    this->x.resize(ldb);
    int nb = 1; // optimal block size
    int nrhs = 1;
    this->lwork = std::max(1, mn + std::max(mn, nrhs)*nb);
    this->work.resize((size_t) this->lwork);
  }

  ~SgLinearSolve_type() {}

  void setMatrix(const double mat[]) {
                            
    // row major storage
    for (int i = 0; i < this->nrow; ++i) {
        for (int j = 0; j < this->ncol; ++j) {
        // C storage
            size_t k = i*this->ncol + j;
            this->mat[k] = mat[k];
            this->matOri[k] = mat[k];
        }
    }
  }

  void setRightHandSide(const double b[]) {
    for (int i = 0; i < this->nrow; ++i) {
        this->x[i] = b[i];
        this->b[i] = b[i];
    } 
  }

  int solve() {
    int errCode = 0;
    char t = 'T'; // C storage
    int nrhs = 1;
    int ldb = (int) this->b.size();
    _GELS_(&t, 
           &this->nrow, 
           &this->ncol, 
           &nrhs,
           &this->mat.front(), 
           &this->nrow,
           &this->x.front(), &ldb,
           &this->work.front(), &this->lwork, 
           &errCode);

    // if errCode == -i, the i-th argument had an illegal value
    // if errCode == i, the i-th diagonal element of the
    //                  triangular factor of A is zero, so that 
    //                  A does not have full rank; the least 
    //                  squares solution could not be computed. 
    return errCode;
  }

  void getSolution(double** sol) {
    *sol = &this->x.front();
  }

  double getResidual() const {
    double res = 0;
    for (int i = 0; i < this->nrow; ++i) {
      double val = 0;
        for (int j = 0; j < this->ncol; ++j) {
          int k = i*this->ncol + j;
          val += this->matOri[k]*this->x[j];
        }
        val -= this->b[i];
        res += fabs(val);
    }
    return res;
  }

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
