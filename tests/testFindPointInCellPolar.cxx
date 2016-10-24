/**
 * Testing cell search in 1D
 */

 #include <vector>
 #include <cassert>
 #include <cstdio>
 #include <iostream>
 #include <cmath>
 #include "SgFindPointInCell.h"

int main(int argc, char** argv) {

	int ier;
	const int nitermax = 10;
	const double tolpos = 1.e-10;

	// number of dimensions
	const int ndims = 2;

	// initial guess
	double dIndices[] = {3.8, 4.9};

	// target position
	const double targetPoint[] = {0.7, M_PI/3.};

	// grid
	size_t nr = 11;
	size_t nt = 21;
	std::vector<double> rho(nr * nt);
	std::vector<double> the(nr * nt);
	for (size_t i = 0; i < nr; ++i) {
		double r = 0. + (1. - 0.) * i / double(nr - 1);
		for (size_t j = 0; j < nt; ++j) {
			double t = 0. + (2*M_PI - 0.)* j / double(nt - 1);
			size_t k = j * nt*i;
			rho[k] = r;
			the[k] = t;
		}
	}

	const int dims[] = {nr, nt};
	std::vector< double* > coords(ndims);
	coords[0] = &rho[0];
	coords[1] = &the[0];

	double pos[ndims];
	double oldPos[ndims];

	SgFindPointInCell_type* picf = NULL;
	ier = SgFindPointInCell_new(&picf, nitermax, tolpos);
	assert(ier == 0);

	ier = SgFindPointInCell_setGrid(&picf, ndims, dims, (const double**) &coords[0]);
	assert(ier == 0);

	ier = SgFindPointInCell_reset(&picf, dIndices, targetPoint);
	assert(ier == 0);

	bool iterFlag = true;
	size_t icount = 0;
	while (iterFlag) {

		ier = SgFindPointInCell_getPosition(&picf, oldPos);
		assert(ier == 0);

		ier = SgFindPointInCell_next(&picf);
		if (ier != 0) iterFlag = false;

		ier = SgFindPointInCell_getPosition(&picf, pos);
		assert(ier == 0);

		std::cout << "iter " << icount << " position old = ";
		for (int j = 0; j < ndims; ++j) std::cout << oldPos[j] << " ";
		std::cout << " -> new = ";
	    for (int j = 0; j < ndims; ++j) std::cout << pos[j] << " ";
	    std::cout << '\n';

		double error = 0;
		for (int j = 0; j < ndims; ++j) {
			double dp = targetPoint[j] - pos[j];
			error = dp * dp;
		}
		error = sqrt(error);
		std::cout << " error = " << error << "\n";

		icount++;
	}

	ier = SgFindPointInCell_del(&picf);
	assert(ier == 0);

	return 0;
}
