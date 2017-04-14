/**
 * Testing flow interpolation in 2D
 */

#include <vector>
#include <cassert>
#include <cstdio>
#include <cmath>
#include <string>
#include <fstream>
#include <iostream>
#include "SgFlowInterp2D.h"
#include "SgBoxIterator.h"
#include "createGrids2D.h"
#include "CmdLineArgParser.h" 

void createLineGrid(const int dims[], const double xmins[], const double xmaxs[], double* coords[]) {
    double dx[] = {(xmaxs[0] - xmins[0]) / double(dims[0] - 1),  
                   (xmaxs[1] - xmins[1]) / double(dims[0] - 1)};
    for (size_t j = 0; j < 2; ++j) {
        for (size_t i = 0; i < dims[0]; ++i) {
            coords[j][i] = xmins[j] + i * dx[j];
        }
    }
}

/**
 * The stream function
 */
double psi(const std::vector<double>& pos) {
    double x = pos[0];
    double y = pos[1];
    return 0.5*(y*y + cos(2.*M_PI*x)/M_PI);
}

void saveStreamlinesVtk(const char* filename,
                        const int dims[], const double** coords, 
                        double (*psiFunc)(const std::vector<double>&)) {

    const int zeros[] = {0, 0};
    int inds[] = {-1, -1};

    std::fstream f;
    f.open(filename, std::ios_base::out);

    f << "# vtk DataFile Version 2.0\n";
    f << "streamFunctionInterp\n";
    f << "ASCII\n";
    f << "DATASET STRUCTURED_GRID\n";
    f << "DIMENSIONS 1 " << dims[1] << " " << dims[0] << "\n";
    int numPoints = dims[0] * dims[1];
    f << "POINTS " << numPoints << " float\n";
    for (size_t index = 0; index < numPoints; ++index) {
        f << coords[0][index] << " " << coords[1][index] << " 0.0\n";
    }
    f << "POINT_DATA " << numPoints << '\n';
    f << "SCALARS psi float\n";
    f << "LOOKUP_TABLE default\n";

     for (size_t index = 0; index < numPoints; ++index) {
        double x = coords[0][index];
        double y = coords[1][index];
        std::vector<double> pos(2);
        pos[0] = x; pos[1] = y;
        f << psi(pos) << '\n';
    }
    int numCells = (dims[0] - 1) * (dims[1] - 1);
    const int cellDims[] = {dims[0] - 1, dims[1] - 1};

    f << "CELL_DATA " << numCells << '\n';
    f << "SCALARS velocity float 3\n";
    f << "LOOKUP_TABLE default\n";

    SgBoxIterator_type nodeIt(2, zeros, dims);
    SgBoxIterator_type cellIt(2, zeros, cellDims);
    const int edgeDimsX[] = {dims[0] - 1, dims[1]};
    const int edgeDimsY[] = {dims[0], dims[1] - 1};
    SgBoxIterator_type edgeItX(2, zeros, edgeDimsX);
    SgBoxIterator_type edgeItY(2, zeros, edgeDimsY);

    for (size_t index = 0; index < cellIt.getNumberOfElements(); ++index) {

        cellIt.getElement(index, inds);

        size_t node00 = nodeIt.getIndex(inds);
        size_t edgeX0 = edgeItX.getIndex(inds);
        size_t edgeY0 = edgeItY.getIndex(inds);
        double x00 = coords[0][node00];
        double y00 = coords[1][node00];
        std::vector<double> pos00(2);
        pos00[0] = x00; pos00[1] = y00;

        inds[0] += 1;
        size_t node10 = nodeIt.getIndex(inds);
        size_t edgeY1 = edgeItY.getIndex(inds);
        double x10 = coords[0][node10];
        double y10 = coords[1][node10];
        std::vector<double> pos10(2);
        pos10[0] = x10; pos10[1] = y10;

        inds[1] += 1;
        size_t node11 = nodeIt.getIndex(inds);
        double x11 = coords[0][node11];
        double y11 = coords[1][node11];
        std::vector<double> pos11(2);
        pos11[0] = x11; pos11[1] = y11;

        inds[0] -= 1;
        size_t node01 = nodeIt.getIndex(inds);
        size_t edgeX1 = edgeItX.getIndex(inds);
        double x01 = coords[0][node01];
        double y01 = coords[1][node01];
        std::vector<double> pos01(2);
        pos01[0] = x01; pos01[1] = y01;
        
        double vx0 = -(psiFunc(pos01) - psiFunc(pos00))/(y01 - y00);
        double vx1 = -(psiFunc(pos11) - psiFunc(pos10))/(y11 - y10);
        double vy0 = +(psiFunc(pos10) - psiFunc(pos00))/(x10 - x00);
        double vy1 = +(psiFunc(pos11) - psiFunc(pos01))/(x11 - x01);

        double vx = 0.5*(vx0 + vx1);
        double vy = 0.5*(vy0 + vy1);

        f << vx << ' ' << vy << " 0.0\n";
    }

    f.close();
}

int main(int argc, char** argv) {

    // parse command line arguments
    CmdLineArgParser prsr;
    prsr.set("--dims", "21,11", "Number of nodes in the x and y");
    prsr.set("--pa", "0.,0.1", "Start position");
    prsr.set("--pb", "1.,0.1", "End position");
    prsr.set("--help", false, "Help");
    prsr.parse(argc, argv);

    if (prsr.get<bool>("--help")) {
        prsr.help();
        return 0;
    }

    std::vector<int> dims = prsr.get<std::vector<int> >("--dims");
    int srcNodeDims[] = {dims[0], dims[1]};
    std::vector<double> pa = prsr.get<std::vector<double> >("--pa");
    std::vector<double> pb = prsr.get<std::vector<double> >("--pb");

    std::cout << "Compute flow integral: " << srcNodeDims[0] << "x" << srcNodeDims[1] <<  " pa = ";
    for (size_t i = 0; i < pa.size(); ++i) {
        std::cout << pa[i] << ", ";
    }
    std::cout << " pb = ";
    for (size_t i = 0; i < pb.size(); ++i) {
        std::cout << pb[i] << ", ";
    }
    std::cout << '\n';


    // create src grid
    const double srcXmins[] = {-1.0, -1.0};
    const double srcXmaxs[] = {1.0, 1.0};
    int srcNumPoints = srcNodeDims[0] * srcNodeDims[1];
    double* srcCoords[] = {new double[srcNumPoints], new double[srcNumPoints]};
    createRectangularGrid(srcNodeDims, srcXmins, srcXmaxs, srcCoords);

    // src field
    int srcNumEdgeX = (srcNodeDims[0] - 1) * srcNodeDims[1];
    int srcNumEdgeY = srcNodeDims[0] * (srcNodeDims[1] - 1);
    double* srcData[] = {new double[srcNumEdgeX], new double[srcNumEdgeY]};

    double hx = (srcXmaxs[0] - srcXmins[0])/double(srcNodeDims[0] - 1);
    double hy = (srcXmaxs[1] - srcXmins[1])/double(srcNodeDims[1] - 1);
    size_t k;
    std::vector<double> posLo(2);
    std::vector<double> posHi(2);
    const int zeros[] = {0, 0};
    int inds[] = {-1, -1};

    // fluxes along x (y-field)
    k = 0;
    const int srcEdge0Dims[] = {srcNodeDims[0] - 1, srcNodeDims[1]}; // x, y
    SgBoxIterator_type edge0It(2, zeros, srcEdge0Dims);
    for (int index = 0; index < edge0It.getNumberOfElements(); ++index) {
        edge0It.getElement(index, inds);
        int i = inds[0]; 
        int j = inds[1];
        double xLo = srcXmins[0] + i*hx;
        double xHi = xLo + hx;
        double y = srcXmins[1] + j*hy;
        posLo[0] = xLo; posLo[1] = y;
        posHi[0] = xHi; posHi[1] = y;
        srcData[k][index] = psi(posHi) - psi(posLo);
    }

    // fluxes along y (x field)
    k = 1;
    const int srcEdge1Dims[] = {srcNodeDims[0], srcNodeDims[1] - 1}; // x, y
    SgBoxIterator_type edge1It(2, zeros, srcEdge1Dims);
    for (int index = 0; index < edge1It.getNumberOfElements(); ++index) {
        edge1It.getElement(index, inds);
        int i = inds[0]; 
        int j = inds[1];
        double x = srcXmins[1] + i*hx;
        double yLo = srcXmins[1] + j*hy;
        double yHi = yLo + hy;
        posLo[0] = x; posLo[1] = yLo;
        posHi[0] = x; posHi[1] = yHi;
        srcData[k][index] = psi(posHi) - psi(posLo);
    }

    saveStreamlinesVtk("streamFunctionInterp.vtk", srcNodeDims, (const double**) srcCoords, psi);

    // destination grid, a segment
    const int dstDims[] = {2}; // number of nodes
    int dstNumPoints = dstDims[0];
    int dstNumCells = dstDims[0] - 1;
    double* dstCoords[] = {new double[dstNumPoints], new double[dstNumPoints]};
    createLineGrid(dstDims, &pa[0], &pb[0], dstCoords);    

    // dst data, dimensioned number of cells
    double* dstData = new double[dstNumCells];

    // project src field onto dst grid
    SgFlowInterp2D_type* interp = NULL;
    SgFlowInterp2D_new(&interp);
    SgFlowInterp2D_setSrcGrid(&interp, srcNodeDims, (const double**) srcCoords);
    SgFlowInterp2D_setDstGrid(&interp, dstDims, (const double**) dstCoords);
    SgFlowInterp2D_computeWeights(&interp);
    SgFlowInterp2D_apply(&interp, (const double**) srcData, dstData);
    SgFlowInterp2D_debug(&interp);
    SgFlowInterp2D_del(&interp);

    // check
    double exactFlux = psi(pb) - psi(pa);
    std::cout << "Flux: " << dstData[0] << 
                 " exact: " << exactFlux <<
                 " error: " << dstData[0] - exactFlux << '\n';

    // cleanup
    for (size_t j = 0; j < 2; ++j) {
        delete[] srcCoords[j];
        delete[] dstCoords[j];
        delete[] srcData[j];
    }
    delete[] dstData;

    return 0;
}
