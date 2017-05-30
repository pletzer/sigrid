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
#include "CmdLineArgParser.h"

std::vector<double> getDoubleVectorFromString(const std::string& s) {
    std::vector<double> res;
    size_t begPos = 0;
    bool cont = true;
    while (cont) {
        size_t endPos = s.find(',', begPos);
        size_t n = endPos - begPos;
        std::string sval = s.substr(begPos, n);
        double val = atof(sval.c_str());
        res.push_back(val);
        begPos = endPos + 1;
        if (endPos > s.size()) {
            cont = false;
        }
    }
    return res;
}

int main(int argc, char** argv) {

    CmdLineArgParser prsr;
    prsr.set("-p", std::string("0., 0."), "Target position");
    prsr.set("-i", std::string("0., 0."), "Initial index location (guess)");
    prsr.set("-nr", 11, "Number of grid radii");
    prsr.set("-nt", 33, "Number of grid poloidal points");
    prsr.set("-m", 100, "Max number of iterations");
    prsr.set("-t", 1.e-10, "Tolerance in physical space");
    prsr.parse(argc, argv);

    if (prsr.get<bool>("-h") || prsr.get<bool>("--help")) {
        prsr.help();\
        return 1;
    }


    int nitermax = prsr.get<int>("-m");
    double tolpos = prsr.get<double>("-t");
    int nt = prsr.get<int>("-nt");
    int nr = prsr.get<int>("-nr");

    std::vector<double> targetPoint = getDoubleVectorFromString(prsr.get<std::string>("-p"));
    std::vector<double> dIndices = getDoubleVectorFromString(prsr.get<std::string>("-i"));

    int ier;
    // number of dimensions
    const int ndims = 2;

    std::cout << "target point is ";
    for (size_t i = 0; i < ndims; ++i) std::cout << targetPoint[i] << ' ';
    std::cout << '\n';

    // grid
    
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

    ier = SgFindPointInCell_reset(&picf, &dIndices[0], &targetPoint[0]);
    assert(ier == 0);

    int end = 0;
    size_t icount = 0;
    double error = std::numeric_limits<double>::max();

    printf("iter    old_position             new_position               error\n");
    while (end == 0) {

        printf("%lu   ", icount);

        ier = SgFindPointInCell_getPosition(&picf, oldPos);
        assert(ier == 0);

        for (size_t i = 0; i < ndims; ++i) {
            printf("%10.4f ", oldPos[i]);
        }
        printf("    ");


        end = SgFindPointInCell_next(&picf);

        ier = SgFindPointInCell_getPosition(&picf, pos);
        assert(ier == 0);

        for (size_t i = 0; i < ndims; ++i) {
            printf("%10.5f ", pos[i]);
        }
        printf("    ");

        ier = SgFindPointInCell_getError(&picf, &error);
        assert(ier == 0);

        printf("    %e\n", error);

        icount++;
    }

    ier = SgFindPointInCell_del(&picf);
    assert(ier == 0);

    if (error > tolpos) {
        // did not converge
        std::cout << "ERROR " << end << '\n';
        return 1;
    }

    std::cout << "SUCCESS\n";
    return 0;
}
