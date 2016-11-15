/**
 * Testing cconservative interpolation in 2D
 */

 #include <vector>
 #include <cassert>
 #include <cstdio>
 #include <iostream>
 #include "SgConserveInterp2D.h"

void create2DGrid(const int nodeDims[], 
                  const double xmins[], const double xmaxs[], 
                  double** coords) {


    size_t index = 0;
    for (size_t i = 0; i < nodeDims[0]; ++i) {
        for (size_t j = 0; j < nodeDims[1]; ++j) {
            coords[0][index] = xmins[0] + (xmaxs[0] - xmins[0])*i/float(nodeDims[0] - 1);
            coords[1][index] = xmins[1] + (xmaxs[1] - xmins[1])*j/float(nodeDims[1] - 1);
            index++;
        }
    }

}

bool testSimple() {

    // source grid
    const int srcDims[] = {2, 2};
    const double srcXmins[] = {0., 0.};
    const double srcXmaxs[] = {1., 1.};
    int srcNumPoints = srcDims[0] * srcDims[1];
    double* srcCoords[] = {new double[srcNumPoints], new double[srcNumPoints]};
    create2DGrid(srcDims, srcXmins, srcXmaxs, srcCoords);

    // destination grid
    const int dstDims[] = {2, 2};
    const double dstXmins[] = {0., 0.};
    const double dstXmaxs[] = {1., 1.};
    int dstNumPoints = dstDims[0] * dstDims[1];
    double* dstCoords[] = {new double[dstNumPoints], new double[dstNumPoints]};
    create2DGrid(dstDims, dstXmins, dstXmaxs, dstCoords);    

    SgConserveInterp2D_type* interp = NULL;
    SgConserveInterp2D_new(&interp);
    SgConserveInterp2D_setSrcGrid(&interp, srcDims, (const double**) srcCoords);
    SgConserveInterp2D_setDstGrid(&interp, dstDims, (const double**) dstCoords);
    SgConserveInterp2D_computeWeights(&interp);
    SgConserveInterp2D_reset(&interp);
    int srcIndex, dstIndex;
    double weight;
    while (!SgConserveInterp2D_next(&interp)) {
        SgConserveInterp2D_get(&interp, &srcIndex, &dstIndex, &weight);
        std::cout << "src index: " << srcIndex << " dst index: " 
                  << dstIndex << " weight: " << weight << '\n';
    }
    SgConserveInterp2D_del(&interp);

    return true;
}

int main(int argc, char** argv) {

    if (!testSimple()) return 1;

    std::cout << "SUCCESS\n";
    return 0;
}
