/**
 * Testing cell search in 2D polar grid
 */

#include <vector>
#include <cassert>
#include <cstdio>
#include <iostream>
#include <cmath>
#include "SgFindPointInCell.h"
#include "createGrids2D.h"



int main(int argc, char** argv) {

	int ier;
	const int nitermax = 100;
	const double tolpos = 1.e-10;

	// number of dimensions
	const int ndims = 2;

	// initial guess
	double dIndices[] = {6.8, 3.4};

	// target position
	const double targetPoint[] = {0.7, M_PI/3.};
	std::cout << "target point is ";
	for (size_t i = 0; i < ndims; ++i) std::cout << targetPoint[i] << ' ';
	std::cout << '\n';

	// grid
	
	int nr = 11;
	int nt = 33;
	int dims[] = {nr, nt};
	double* coords[] = {new double[nr * nt], new double[nr * nt]};
	const double center[] = {0., 0.};
	const double radius = 1.0;
	createPolarGrid(dims, center, radius, coords);

	double pos[ndims];
	double oldPos[ndims];

	SgFindPointInCell_type* picf = NULL;
	ier = SgFindPointInCell_new(&picf, nitermax, tolpos);
	assert(ier == 0);

	// periodic in theta
	const int periodicity[] = {0, 1};
	ier = SgFindPointInCell_setGrid(&picf, ndims, dims, periodicity,
		                            (const double**) &coords[0]);
	assert(ier == 0);

	ier = SgFindPointInCell_reset(&picf, dIndices, targetPoint);
	assert(ier == 0);

	bool iterFlag = true;
	size_t icount = 0;
	while (iterFlag) {

		ier = SgFindPointInCell_getPosition(&picf, oldPos);
		assert(ier == 0);

		ier = SgFindPointInCell_next(&picf);
		if (ier != 0) {
			iterFlag = false;
			if (ier < 0) {
				std::cout << "*** reached max number of iterations!\n";
			}
		}

		ier = SgFindPointInCell_getPosition(&picf, pos);
		assert(ier == 0);

		std::cout << "iter " << icount << " position old = ";
		for (int j = 0; j < ndims; ++j) std::cout << oldPos[j] << " ";
		std::cout << " -> new = ";
	    for (int j = 0; j < ndims; ++j) std::cout << pos[j] << " ";

		double error;
		ier = SgFindPointInCell_getError(&picf, &error);
		assert(ier == 0);
		std::cout << " (error = " << error << ")\n";

		icount++;
	}

	ier = SgFindPointInCell_del(&picf);
	assert(ier == 0);

	return 0;
}
