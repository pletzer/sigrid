/**
 * Testing nodal interpolation in 2D
 */

#include <vector>
#include <cassert>
#include <cstdio>
#include <iostream>
#include "SgFindPointInCell.h"
#include "createGrids2D.h"
#include "CmdLineArgParser.h"
#include "SgNdims.h"

std::vector<double> getVecFromString(const std::string& expr) {
    std::vector<double> res;
    size_t posBeg = 0;
    size_t posEnd;
    bool go = true;
    while (go) {
        posEnd = expr.find(',', posBeg);
        std::string sn = expr.substr(posBeg, posEnd - posBeg);
        res.push_back(atof(sn.c_str()));
        if (posEnd == std::string::npos) {
            go = false;
        }
        posBeg = posEnd + 1;
    }
    return res;
}

int main(int argc, char** argv) {

    CmdLineArgParser prsr;
    prsr.set("--src_nj", 41, "Number of source grid latitudes");
    prsr.set("--src_ni", 81, "Number of source grid longitudes");
    prsr.set("--delta_lat", 30.0, "Latitude pole displacement in deg.");
    prsr.set("--delta_lon", 20.0, "Longitude pole displacement in deg.");
    prsr.set("--nitermax", 100, "Max number of iterations");
    prsr.set("--tolpos", 1.e-10, "Tolerance in physical space");
    prsr.set("--target", std::string("66,-36"), "Target point");
    prsr.parse(argc, argv);

    if (prsr.get<bool>("-h") || prsr.get<bool>("--help")) {
        prsr.help();
        return 1;
    }

    // source grid
    const int srcDims[] = {prsr.get<int>("--src_nj"), prsr.get<int>("--src_ni")};
    int srcNumPoints = srcDims[0] * srcDims[1];
    double* srcCoords[] = {new double[srcNumPoints], new double[srcNumPoints]};
    double delta_lat = prsr.get<double>("--delta_lat");
    double delta_lon = prsr.get<double>("--delta_lon");
    createRotatedPoleGrid(srcDims, delta_lat, delta_lon, srcCoords);

    const int nitermax = prsr.get<int>("--nitermax");
    const double tolpos = prsr.get<double>("--tolpos");

    // target point
    std::string tg = prsr.get<std::string>("--target");
    std::vector<double> targetPoint = getVecFromString(tg);

    // initial guess
    std::vector<double> dIndices(NDIMS_2D_TOPO);
    dIndices[0] = (srcDims[0] - 1)/2.14423;
    dIndices[1] = (srcDims[1] - 1)/1.93265;

    // initialize search
    SgFindPointInCell_type pointFinder(nitermax, tolpos);
    const int periodicity[] = {0, 0}; // no periodicity in lat-lon space
    pointFinder.setGrid(NDIMS_2D_PHYS, srcDims, periodicity, (const double**) srcCoords);

    // search
    pointFinder.reset(&dIndices[0], &targetPoint[0]);
    int stop = 0;
    while(!stop) {
        stop = pointFinder.next();
    }

    double error = pointFinder.getError();
    dIndices = pointFinder.getIndices();
    std::vector<double> position = pointFinder.getPosition();

    std::cout << "Result: indices              = " << dIndices[0] << ", " << dIndices[1] << '\n';
    std::cout << "        position             = " << position[0] << ", " << position[1] << '\n';
    std::cout << "        target               = " << targetPoint[0] << ", " << targetPoint[1] << '\n';
    std::cout << "        error in coord space = " << error << '\n';

    int niter;
    double* errorHistory;
    pointFinder.getErrorHistory(&niter, &errorHistory);
    for (int i = 0; i < niter; ++i) {
        std::cout << i << " error = " << errorHistory[i] << '\n';
    }

    // clean up
    for (size_t j = 0; j < 2; ++j) {
        delete[] srcCoords[j];
    }

    std::cout << "SUCCESS\n";
    return 0;
}
