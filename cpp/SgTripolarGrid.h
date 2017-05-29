/**
 * A class that constructs a tripolar grid
 */
 
#ifndef SG_TRIPOLAR_GRID_H
#define SG_TRIPOLAR_GRID_H
 
#include <vector>
#include <cstdio> // size_t
#include <cmath>
 
struct SgTripolarGrid_type {

    // dimensions
    int dims[2];

    int capLatIndex;

    int ntot;

    // coordinates
    std::vector<double> lats;
    std::vector<double> lons;

    // pointer to lats and lons
    double* coords[2];

/**
 * Constructor
 * @param dims lats/lons dimensions
 * @param capLatIndex index where the latitudes transition to the tripolar grid
 */
SgTripolarGrid_type(const int dims[], int capLatIndex) {
    this->capLatIndex = capLatIndex;
    this->dims[0] = dims[0];
    this->dims[1] = dims[1];
    size_t ntot = dims[0] * dims[1];
    this->lats.resize(ntot);
    this->lons.resize(ntot);
    this->coords[0] = &this->lats.front();
    this->coords[1] = &this->lons.front();
}

/**
 * Destructor
 */
~SgTripolarGrid_type() {
}

/**
 * Get the grid 
 * @param dims number of nodes along each dimensions (output)
 * @param coords arrays of flat coordinates (component, coordinates)
 */
void getGrid(int dims[], double*** coords) {

  dims[1] = this->dims[0];
  dims[0] = this->dims[1];
  *coords = this->coords;

  int i,j;
  double L, P, Lc, Pc, Pc0, Lc0, latPerim;
  double L0 = 0.;
  double di = 2*M_PI / (dims[0]-1);
  double dj = M_PI / (dims[1]-1);

  Pc0 = -0.5*M_PI;
  Lc0 = -M_PI;

  for (j=0; j < this->capLatIndex; ++j) {
    for (i=0; i<dims[0]; ++i) {
      P = Pc0 + j*dj;
      L = Lc0 + i*di;
      this->lons[i+j*dims[0]] = L;
      this->lats[i+j*dims[0]] = P;

      this->lons[i+j*dims[0]] *= 180./M_PI;
      this->lats[i+j*dims[0]] *= 180./M_PI;
    }
  }

  latPerim = M_PI - this->capLatIndex*dj;
  //di = M_PI / (dims[0]-1);
  dj = (0.5*M_PI) / (dims[1] - this->capLatIndex-1);
  Pc0 = -(M_PI/2.);
  Lc0 = (M_PI/2.);
  for (j = this->capLatIndex; j<dims[1]; ++j) {
    for (i=0; i<dims[0]/2; ++i) {
      Pc = Pc0 + i*di;
      Lc = Lc0 - (j - this->capLatIndex)*dj;
      this->lons[i+j*dims[0]] = L0 - atan2(sin(Lc), tan(Pc));
      this->lats[i+j*dims[0]] = 0.5*M_PI  - 2.*atan(tan(latPerim/2.)*
          tan(0.5*acos((cos(Lc)*cos(Pc)))));

      this->lons[i+j*dims[0]] *= 180./M_PI;
      this->lats[i+j*dims[0]] *= 180./M_PI;
    }
  }

  Pc0 = -(M_PI/2.);
  Lc0 = -(M_PI/2.);
  for (j = this->capLatIndex; j<dims[1]; ++j) {
    for (i=dims[0]/2; i<dims[0]; ++i) {
      Pc = Pc0 + (dims[0]-i-1)*di;
      Lc = Lc0 + (1-1.e-10)*(j - this->capLatIndex)*dj;
      this->lons[i+j*dims[0]] = L0 - atan2(sin(Lc), tan(Pc));
      this->lats[i+j*dims[0]] = 0.5*M_PI  - 2.*atan(tan(latPerim/2.)*
          tan(0.5*acos((cos(Lc)*cos(Pc)))));

      this->lons[i+j*dims[0]] *= 180./M_PI;
      this->lats[i+j*dims[0]] *= 180./M_PI;
    }
  }

}


private:

};
 
#ifdef __cplusplus
extern "C" {
#endif

int SgTripolarGrid_new(SgTripolarGrid_type** self,
                       const int dims[], int capLatIndex);
                       
int SgTripolarGrid_del(SgTripolarGrid_type** self);

int SgTripolarGrid_getGrid(SgTripolarGrid_type** self,
                           int dims[], double*** coords);

#ifdef __cplusplus
}
#endif

#endif // SG_TRIPOLAR_GRID_H
