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
 * @param numLevels number of levels
 * @param ndims number of dimensions
 * @param points flat array of points 
 */
SgOctreePoints_type(size_t numLevels, size_t ndims, const std::vector<double>& points) {

    this->ndims = ndims;
    this->nlevs = numLevels;

    // find the low/high corner points
    this->xmins.resize(ndims, std::numeric_limits<double>::max());
    this->xmaxs.resize(ndims, -std::numeric_limits<double>::max());
    for (size_t i = 0; i < points.size(); ++i) {
        size_t dim = i % ndims;
        double p = points[i];
        this->xmins[dim] = (p < this->xmins[dim]? p: this->xmins[dim]);
        this->xmaxs[dim] = (p > this->xmaxs[dim]? p: this->xmaxs[dim]);
    }

    // used to map the flat indices to index sets
    this->prodDims.resize(ndims, 1);
    for (int i = ndims - 2; i >= 0; --i) {
        this->prodDims[i] = prodDims[i + 1] * 2;
    }

    // recursively build the octree
    std::vector<size_t> key;
    this->refineOctree(key);

    // attach each point to a sub-box
    size_t numPoints = points.size() / ndims;
    for (size_t i = 0; i < numPoints; ++i) {
        std::vector<size_t> key = this->getKey(&points[i*ndims], this->nlevs);
        this->node2Key.insert(std::pair<size_t, std::vector<size_t> >(i, key));
    }
}

/**
 * Destructor
 */
~SgOctreePoints_type() {
}

/*
 * Get the low corner of the partition
 * @param key array of indices in the range 0.. 2^ndims - 1
 * @return the low coordinate point of the partition
 */
const std::vector<double>& getLo(const std::vector<size_t>& key) const {
    std::map< std::vector<size_t>, SgOctreeElement >::const_iterator
        it = this->octree.find(key);
    if (it != this->octree.end()) {
        return it->second.xmins;
    }
    return this->xmins;
}

/*
 * Get the high corner of the partition
 * @param key array of indices in the range 0...2^ndims - 1
 * @return the high coordinate point of the partition
 */
const std::vector<double>& getHi(const std::vector<size_t>& key) const {
    std::map< std::vector<size_t>, SgOctreeElement >::const_iterator
        it = this->octree.find(key);
    if (it != this->octree.end()) {
        return it->second.xmaxs;
    }
    return this->xmaxs;
}

/**
 * Get the key associated with a point
 * @param pt point
 * @param lev level
 * @return key, array of integers in the range 0...2^ndims - 1
 */
std::vector<size_t> getKey(const double* pt, size_t lev) {
    // normalize
    std::vector<double> x(pt, pt + this->ndims);
    for (size_t j = 0; j < this->ndims; ++j) {
        // xmaxs should be > xmins
        double len = this->xmaxs[j] - this->xmins[j];
        x[j] = (pt[j] - this->xmins[j])/len;
    }
    std::vector<size_t> key(lev, 0);
    for (size_t el = 0; el < lev; ++el) {
        double fact = std::pow(2, el + 1);
        for (size_t j = 0; j < this->ndims; ++j) {
            size_t indx = (size_t) floor(x[j] * fact);
            key[el] += this->prodDims[j] * indx;
        }   
    }
    return key;
}

/**
 * Refine an existing partition
 * @param key, array of integers in the range 0...2^ndims - 1
 */
void refineOctree(const std::vector<size_t>& key) {
    // iterate over the 2^ndims quadrants
    for (size_t iQ = 0; iQ < std::pow(2, this->ndims); ++iQ) {
        // copy the key
        std::vector<size_t> key2 = key;
        // append the quadrant index
	key2.push_back(iQ);
        // compute and fill in the low/high corner points
        std::vector<double> xLo, xHi;
        this->computeLoHi(key2, xLo, xHi);
        // add the quadrant
        this->octree.insert(std::pair< std::vector<size_t>, SgOctreeElement>(key2, SgOctreeElement(xLo, xHi)));
        if (key2.size() < this->nlevs) {
            // recursively call this method
            this->refineOctree(key2);
        }
        // done
    }
}

/**
 * Compute the box boundaries
 * @param key, array of integers in the range 0...2^ndims - 1
 * @param xLo low corner (output)
 * @param xHi high corner (output)
 */
void computeLoHi(const std::vector<size_t>& key, std::vector<double>& xLo, std::vector<double>& xHi) {
    xLo.resize(this->ndims, 0.0);
    xHi.resize(this->ndims, 1.0);
    double dx = 0.5;
    for (size_t i = 0; i < key.size(); ++i) {
        size_t iQuad = key[i];
        std::vector<size_t> inds = this->getIndices(iQuad);
        for (size_t j = 0; j < this->ndims; ++j) {
            xLo[j] += inds[j] * dx;
            xHi[j] = xLo[j] + dx;
        }
        dx /= 2.0;
    }
    // scale and translate
    for (size_t j = 0; j < this->ndims; ++j) {
        xLo[j] *= (this->xmaxs[j] - this->xmins[j]);
        xHi[j] *= (this->xmaxs[j] - this->xmins[j]);
        xLo[j] += this->xmins[j];
        xHi[j] += this->xmins[j];
    }
}

/**
 * Get the index set of the partition
 * @param indx falt index in the range 0...2^ndims - 1
 * @return index set, array of zeros and ones
 */
std::vector<size_t> getIndices(size_t indx) {
    std::vector<size_t> inds(this->ndims, 0);
    for (size_t j = 0; j < this->ndims; ++j) {
        inds[j] = (indx / this->prodDims[j]) % 2;
    }
    return inds;
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
