
#include "jEColi.h"


#define FMULT 0.125


jEColi::jEColi (jChipmunk_Habitat *h,
                float x, float y, float a, float v ,
                jEColi_t type)
{
  forceBuffer[0].f =forceBuffer[0].v = forceBuffer[1].f =forceBuffer[1].v = cpvzero;

  marked=false;    volume = v;    div_vol = 3.5;

  habitat = h;
  cellType = type;

  jChipmunk_EColiInit_t cell;
      cell.habitat = habitat;
      if(ECOLITYPE_STUCKCELL == cellType)
      {
          cell.mass = 1000000.0;  cell.moment = 1000000.0;
          cell.length = 80.0;  cell.width = 50.0;
      }
      else
      {
          cell.mass = 1.0;  cell.moment = 100.0;
         // cell.mass = 0.5;  cell.moment = 100.0;
//          cell.mass = 0.1;  cell.moment = 50.0;
//          cell.mass = 1.0;  cell.moment = 1.0;
          cell.length = jlengthFromVol(v);  cell.width = 10.0;
      }
      cell.x = x;  cell.y = y;  cell.angle = a;
      cell.vx = 0.0;  cell.vy = 0.0;
      cell.owner = this;
      //set the dRL variable:  = dt * 0.1 (um/minute) * 10  (pixels/um) = dt
      cell.dRL = habitat->dt;
  cpCell = new jChipmunk_EColi(&cell);

  div_count = 0;
  force_div = false;

  protein = 0.0;

}

jEColi::~jEColi()
{
    if ( cpCell != NULL ) {
      delete cpCell;
    }
}

int jEColi::jupdate( void ) {
     //call the chipmunk model--pass growth rate as a parameter:
     int status;
     status = cpCell->updateModel(1.0);

    if(-1 == status)
     {
        mark_for_death();
        return -1;
     }

     volume = jget_volume();
     if ( volume > div_vol )  div_count++;

  return 0;
}

jEColi * jEColi::jdivide ( void ) {

//jMODS:  not sure about the div count need:
//    if ( div_count >= 10 || force_div ) {
//    if (0 ) {
    if ( div_count >= 1 || force_div ) {

    int r = frand() > 0.5 ? 1 : -1;

    div_count = 0;
    force_div = false;

    float frac = 0.5 + 0.1 * ( frand() - 0.5 );
    float oldvol = volume;

    double oldProtien = protein;
    protein = frac*oldProtien;
 //============================================================
    //this should all move to the cp model:
    // like cpModel_divideCell(); --with perhaps the rand. parameters

    jChipmunk_EColiInit_t cell;

    volume = frac * oldvol;
    float a = cpCell->angle;
    //random delta-angle for each cell
//    float da = 0.5 * (frand()-0.5);
    // float da = 0.20 * (frand()-0.5);
   float da = 0.0;

    float oldsize = cpCell->length;
    cpVect oldpos = cpCell->center;

    cpVect newpos = cpCell->center + cpvmult ( cpv ( cos ( a - r*da ), sin ( a - r*da ) ),
                                               (-r)*0.5*oldsize*(1-frac) );
//    cpVect vel = cpvzero;
    //assign the old cell-half velocities to new center-of-mass velocities
    //NOTE:  cell will not be expanding at first subsequent time step
    cpVect vel = (1 == r) ?  cpCell->velB:cpCell->velA;

    cell.habitat = habitat;
    //eventually read these from the parent cell:
    cell.mass = cpCell->initialMass;
    cell.moment =cpCell->initialMoment;

    cell.x = newpos.x;  cell.y = newpos.y;  cell.angle = a - r*da;
    cell.length = jlengthFromVol(volume);  cell.width = 10.0;
//    cell.vx = 0.0;  cell.vy = 0.0;
    cell.vx = vel.x;  cell.vy = vel.y;
//    cell.vx = cpCell->velA.x;  cell.vy = cpCell->velA.y;
    cell.owner = this;
    cell.dRL = habitat->dt;

    jChipmunk_EColi *newCell = new jChipmunk_EColi(&cell);

    //set rest length: NEED TO MAKE SETTER WHICH CALLS cp FUNCTION!
    //need to use averaged values to avoid noise in daughter cells:

    float newRestLength = (cpCell->springRestLength - cpCell->length) + cell.length;

//    float newRestLength = 1.0*(cpCell->springRestLengthVar - cpCell->separationDistance);
//    float newRestLength = (cpCell->springRestLengthVar);
//    float newRestLength = (cpCell->shapeLength + cpCell->dRL);
//      float newRestLength = (cpCell->shapeLength * 2.0f);

    newCell->springRestLength = newRestLength;

    //determine whether it is the "right" or "left" cell;  r=1 is "right"=bodyB
    //need to use averaged valued!
    //need to split the left/right halves of each new cell also so inner half is like c.o.m. of old cell.
    cpVect newV1, newV2;
    if(1 == r)
    {
        newV1 = cpCell->velB; newV2 = cpCell->velA;
    }
    else
    {
        newV1 = cpCell->velA; newV2 = cpCell->velB;
    }
//    if( ( cpvlength(cpCell->velA) > 20.0)
//            || ( cpvlength(cpCell->velB) > 20.0))
//    {
//        newV1 = newV2 = cpvzero;
//    }
    newCell->setVelocity(newV1);

    //replace the old cell with the new one:
    delete cpCell;
    cpCell = newCell;

//============================================================
//  new daughter cell
//============================================================
    float dvol = (1-frac)*oldvol;

    //make cell #2 ("daughter") appear:
    // jEColi *daughter = new jEColi (habitat, world,
    jEColi *daughter = new jEColi (habitat,
                                   oldpos.x + r*0.5*oldsize*frac*cos ( a + r*da ),
                                   oldpos.y + r*0.5*oldsize*frac*sin ( a + r*da ),
                                   a+r*da, dvol,
                                   cellType);
    daughter->cpCell->setVelocity(newV2);
//    daughter->cpCell->setVelocity((-1 == r) ? cpCell->velB:cpCell->velA);

    //set rest length
    daughter->cpCell->springRestLength = newRestLength - cell.length + jlengthFromVol(dvol);

    daughter->protein = (1.0-frac)*oldProtien;


    // set_division_indicator(true);
    // daughter->set_division_indicator(true);
    // daughter->set_daughter_indicator(true);

    return daughter;

  } else return NULL;

}

void jEColi::setData(cpVect force, cpVect vel, int which)
{
    forceBuffer[which].f=force;
    forceBuffer[which].v=vel;
}
