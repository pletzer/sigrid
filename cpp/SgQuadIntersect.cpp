/**
 * A class that finds all intersection points between two quads
 */
 
#include "SgQuadIntersect.h"
#include <iostream>

extern "C"
int SgQuadIntersect_new(SgQuadIntersect_type** self) {

 	*self = new SgQuadIntersect_type();

 	// tolerance for floating point comparisons
 	(*self)->tol = 1.e-12;

 	(*self)->numIntersectionPoints = 0;

 	(*self)->slvr = NULL;
 	SgLinearSolve_new(&(*self)->slvr, 2, 2);

 	(*self)->quad1Coords = NULL;
 	(*self)->quad2Coords = NULL;

 	return 0;
}
      
extern "C"                   
int SgQuadIntersect_del(SgQuadIntersect_type** self) {

 	if ((*self)->slvr) SgLinearSolve_del(&(*self)->slvr);
 	delete *self;

 	return 0;
}

extern "C"
int SgQuadIntersect_setQuadPoints(SgQuadIntersect_type** self,
	                              const double** quad1Coords, const double** quad2Coords) {

 	(*self)->quad1Coords = (double**) quad1Coords;
 	(*self)->quad2Coords = (double**) quad2Coords;

}


extern "C"
int SgQuadIntersect_getIntersectPoints(SgQuadIntersect_type** self,
 	                                   int* numPoints, double** points) {

	// quickly check if there is any chance of overlap
	if (!(*self)->checkIfBoxesOverlap()) {
		// no chance to have an overlap
		return 0;
	}

	// seems like the quads are at least partially overlapping 
	(*self)->collectNodesInsideQuad((const double**)(*self)->quad1Coords, 
		                            (const double**)(*self)->quad2Coords);
	(*self)->collectNodesInsideQuad((const double**)(*self)->quad2Coords,
		                            (const double**)(*self)->quad1Coords);
	// iterate over edges
	for (size_t i = 0; i < 4; ++i) {
		// the edge of one of the first quad
		size_t quad1Indx0 = i;
		size_t quad1Indx1 = (i + 1) % 4;
		const double* quad1Coord0 = (*self)->quad1Coords[quad1Indx0];
		const double* quad1Coord1 = (*self)->quad1Coords[quad1Indx1];
		// iterate over the other quad's edges
		for (size_t j = 0; j < 4; ++j) {
			// the edges of the second quad
			size_t quad2Indx0 = j;
			size_t quad2Indx1 = (j + 1) % 4;
			const double* quad2Coord0 = (*self)->quad2Coords[quad2Indx0];
			const double* quad2Coord1 = (*self)->quad2Coords[quad2Indx1];
			(*self)->collectEdgeToEdgeIntersectionPoints(quad1Coord0, quad1Coord1,
				                                         quad2Coord0, quad2Coord1);
		}
	}

	*numPoints = (*self)->numIntersectionPoints;
	*points = (*self)->intersectionPoints;

 	return 0;
}
