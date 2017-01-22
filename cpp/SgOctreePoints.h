/**
 * Classify location of points using an octree
 */
 
#ifndef SG_OCTREE_POINTS_H
#define SG_OCTREE_POINTS_H
 
#include <vector>
#include <cstdio> // size_t
#include <cmath>
#include <limits>
#include <map>
#include <vector>

// min/max coordinates of box
struct SgOctreeElement {
    std::vector<double> xmins;
    std::vector<double> xmaxs;
    SgOctreeElement(const std::vector<double>& xLo, const std::vector<double>& xHi) : xmins(xLo), xmaxs(xHi) {
    }  

};

struct SgOctreePoints_type {

    // maps a key to a box in space. A key is a list of ints, each
    // int being a flat index of a box
    std::map< std::vector<size_t>, SgOctreeElement > octree;

    // maps a node index to a key
    std::map< size_t, std::vector<size_t> > node2Key;

    // need this to convert flat box indices to index sets
    std::vector<size_t> prodDims;

    // domain sizes
    std::vector<double> xmins;
    std::vector<double> xmaxs;

    // number of space dimensions
    size_t ndims;

    // number of levels
    size_t nlevs; 

/**
 * Constructor
 * @param 
 * @param 
 */
SgOctreePoints_type(size_t numLevels, size_t ndims, const std::vector<double>& points) {

    // find the low/high corner points
    std::vector<double> xmins(ndims, std::numeric_limits<double>::max());
    std::vector<double> xmaxs(ndims, -std::numeric_limits<double>::max());
    for (size_t i = 0; i < points.size(); ++i) {
        size_t dim = i % ndims;
        double p = points[i];
        xmins[dim] = (p < xmins[dim]? p: xmins[dim]);
        xmaxs[dim] = (p > xmaxs[dim]? p: xmaxs[dim]);
    }

    // used to map the flat indices to index sets
    this->prodDims.resize(ndims, 1);
    for (int i = ndims - 2; i >= 0; --i) {
        this->prodDims[i] = prodDims[i + 1] * 2;
    }

    // build the octree
    std::vector<size_t> key(numLevels);
    for (size_t level = 0; level < numLevels; ++level) {

        // iterate over the quadrants
        for (size_t iBox = 0; iBox < std::pow(2, ndims); ++iBox) {
            key.push_back(iBox);

            // compute the box lo/hi corners in (0, 1) space
            std::vector<double> xLo(ndims, 0);
            std::vector<double> xHi(ndims, 0);
            for (size_t level2 = 0; level2 < numLevels; ++level2) {
                double dx = 1.0 / std::pow(2, level2);
                for (size_t j = 0; j < ndims; ++j) {
                    // index of the index set
                    size_t boxIndx = (key[level2] / this->prodDims[j] % 2);
                    xLo[j] += boxIndx * dx;
                    xHi[j] = xLo[j] + dx;
                }
            }
            // scale and translate
            for (size_t j = 0; j < ndims; ++j) {
                xLo[j] *= (xmaxs[j] - xmins[j]);
                xHi[j] *= (xmaxs[j] - xmins[j]);
                xLo[j] += xmins[j];
                xHi[j] += xmins[j];
            }

            // add entry to the tree
            SgOctreeElement elem(xLo, xHi);
            this->octree.insert(std::pair< std::vector<size_t>, SgOctreeElement>(key, elem));
        }
    }

    // attach each point to a sub-box
    size_t numPoints = points.size() / ndims;
    for (size_t i = 0; i < numPoints; ++i) {
        std::vector<size_t> key = this->getKey(&points[i*ndims], this->nlevs);
        this->node2Key.insert(std::pair<size_t, std::vector<size_t> >::pair(i, key));
    }

    this->ndims = ndims;
    this->nlevs = numLevels;
    this->xmins = xmins;
    this->xmaxs = xmaxs;
}

/**
 * Destructor
 */
~SgOctreePoints_type() {
}


const std::vector<double>& getLo(const std::vector<size_t>& key) const {
    std::map< std::vector<size_t>, SgOctreeElement >::const_iterator
        it = this->octree.find(key);
    if (it != this->octree.end()) {
        return it->second.xmins;
    }
    return this->xmins;
}

const std::vector<double>& getHi(const std::vector<size_t>& key) const {
    std::map< std::vector<size_t>, SgOctreeElement >::const_iterator
        it = this->octree.find(key);
    if (it != this->octree.end()) {
        return it->second.xmaxs;
    }
    return this->xmaxs;
}

std::vector<size_t> getKey(const double* pt, size_t lev) {
    // normalize
    std::vector<double> x(pt, pt + this->ndims);
    for (size_t j = 0; j < this->ndims; ++j) {
        x[j] = (pt[j] - this->xmins[j])/(this->xmaxs[j] - this->xmins[j]); // should not be zero
    }
    std::vector<size_t> key(lev, 0);
    for (size_t el = 0; el < lev; ++el) {
        double fact = std::pow(2, el);
        for (size_t j = 1; j < this->ndims; ++j) {
            size_t indx = (size_t) floor(x[j] * fact);
            key[j] += indx * this->prodDims[j] * indx;
        }   
    }
    return key;
}

};

 
#ifdef __cplusplus
extern "C" {
#endif

int SgOctreePoints_new(SgOctreePoints_type** self, 
                       int numLevels, int ndims, int npoints, const double* points);
                       
int SgOctreePoints_del(SgOctreePoints_type** self);


#ifdef __cplusplus
}
#endif

#endif // SG_OCTREE_POINTS_H
