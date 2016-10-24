/**
 * Testing cell search in 1D
 */

 #include <vector>
 #include <cassert>
 #include <cstdio>
 #include <iostream>
 #include "SgFindPointInCell.h"

int main(int argc, char** argv) {

	int ier;
	const int nitermax = 1;
	const double tolpos = 1.e-10;

	// number of dimensions
	const int ndims = 1;

	// initial guess
	double dIndices[] = {3.8};

	// target position
	const double targetPoint[] = {0.3};

	// grid
	std::vector<double> xs(11);
	for (size_t i = 0; i < xs.size(); ++i) {
		xs[i] = 0. + (1. - 0.)*i/double(xs.size() - 1);
	}

	const int dims[] = {xs.size()};
	std::vector< double* > coords(ndims);
	coords[0] = &xs[0];

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

		std::cout << "iter " << icount << " old position = ";
		for (int j = 0; j < ndims; ++j) std::cout << oldPos[j] << " ";
		std::cout << " -> new position = ";
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
