#ifndef _JECOLI_H_
#define _JECOLI_H_

//jMODS:
#include <iostream>
#include <fstream>

#include "jchipmunk_habitat.h"
#include "jchipmunk_ecoli.h"

#define DEFAULT_ECOLI_INIT_SIZE 1.57      // fl
#define DEFAULT_ECOLI_DIAMETER 1.0        // um
#define PI 3.1415


typedef enum {
    ECOLITYPE_01,
    ECOLITYPE_02,
    ECOLITYPE_SPEEDTEST,

    ECOLITYPE_GAP0,
    ECOLITYPE_GAP1,
    ECOLITYPE_GAP2,
    ECOLITYPE_GAP3,

    ECOLITYPE_STUCKCELL

} jEColi_t;

typedef struct
{
    cpVect f;
    cpVect v;
} cellFVdata_t;

// class jEColi : public Cell {
class jEColi  {

 public:
  inline void mark_for_death ( void ) { marked = true; }
  inline bool marked_for_death ( void ) { return marked; }
  inline cpFloat frand(void) { return (cpFloat)rand()/(cpFloat)RAND_MAX; }
  
  void set_id ( int i ) { id = i; }
  int get_id ( void ) { return id; }


  jEColi (jChipmunk_Habitat *habitat,
              float x, float y, float a, float v,
              jEColi_t);
  ~jEColi();

  int jupdate( void );
  jEColi * jdivide ( void );
  void compute_parameter_derivatives ( void );

  float jget_length ( void ) {
  //length in pixels
      return cpCell->length;
  }
  float jlengthFromVol ( float v ) {
  //volume in fL, length in pixels (*10)
      return 10*v/( 0.25 * PI * DEFAULT_ECOLI_DIAMETER * DEFAULT_ECOLI_DIAMETER );
  }
  float jget_volume ( void ) {
  //return vol. in fL, length in pixels
    return (1/10.0)*jget_length()*( 0.25 * PI * DEFAULT_ECOLI_DIAMETER * DEFAULT_ECOLI_DIAMETER );
  }

  void force_divide ( void ) { force_div = true; }

//10 MAR.2016:
    jEColi_t           cellType;
//13May.2016
    void setData(cpVect force, cpVect vel, int which);
//29Oct.2016:
    unsigned int itheta, izt, icomp, izc;
    double timeTare_thetaz, timeTare_comp;
    unsigned int thetaold, ztold, compold, zcold;

//02Feb.2017:
    unsigned int itd, iztd;
    double timeTare_tdz;
    unsigned int tdold, ztdold;

    double protein;

 private:
  int id;
  float volume, lambda, div_vol;
  int div_count;
  bool force_div;
  bool marked, divided, daughter, selected;

//jMODS:  Add a chipmunk structure for cell-local storage of its cp. parameters
  jChipmunk_Habitat     *habitat;
  jChipmunk_EColi       *cpCell;
//13May.2016
  cellFVdata_t              forceBuffer[2];// = {{cpvzero,cpvzero}, {cpvzero,cpvzero}};


  friend void jcpBodyUpdateVelocity(cpBody *body, cpVect gravity, cpFloat damping, cpFloat dt);
  friend int main(int argc, char **argv);
  friend void recordData(jEColi *thisCell, float t);
  friend int runSimulation(unsigned int height, unsigned int width, unsigned int numInitCells, double maxSimTime,
                    double &finalStrainRatio, double &elapsedTimeForSimulation, double &);
  friend bool setBinData(jEColi * c, int which);
  friend void recordDataMM(jEColi * thisCell, double &simTime);

};

#endif
