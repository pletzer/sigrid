/** 
 * Test linear solve
 */

#include <SgLinearSolve.h>
#include <iostream>
#include <cassert>

int main(int argc, char** argv) {
	const int nrow = 2;
	const int ncol = 2;
	const double mat[] = {1., 2., 3., 4.};
	const double b[] = {1., 2.};
	double* x;
	int ier;
	SgLinearSolve_type* slvr;
	ier = SgLinearSolve_new(&slvr, nrow, ncol);
	assert(ier == 0);
	ier = SgLinearSolve_setMatrix(&slvr, mat);
	assert(ier == 0);
	ier = SgLinearSolve_setRightHandSide(&slvr, b);
	assert(ier == 0);
	ier = SgLinearSolve_solve(&slvr);
	assert(ier == 0);
	ier = SgLinearSolve_getSolution(&slvr, &x);
	assert(ier == 0);
	std::cout << "solution = ";
	for (int i = 0; i < nrow; ++i) 
		std::cout << x[i] << ' ';
	std::cout << '\n';

	std::cout << "mat . x - b = "; 
	for (int i = 0; i < nrow; ++i) {
		double val = 0;
		for (int j = 0; j < ncol; ++j) {
			int k = i*ncol + j;
			val += mat[k] * x[j];
		}
		std::cout << val - b[i] << ' ';
	}
	std::cout << '\n';

	double residual;
	ier = SgLinearSolve_getResidual(&slvr, &residual);
	assert(ier == 0);
	std::cout << "residual = " << residual << '\n';
	ier = SgLinearSolve_del(&slvr);
	assert(ier == 0);

	return 0;
}
