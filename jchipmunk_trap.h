#ifndef JCHIPMUNK_TRAP_H
#define JCHIPMUNK_TRAP_H

#include "jchipmunk_habitat.h"

class jChipmunk_Habitat;
class jChipmunk_EColi;
class jEColi;

class jChipmunk_Trap
{
public:
    jChipmunk_Trap(jChipmunk_Habitat *habitat);
    jChipmunk_Trap(jChipmunk_Habitat *habitat, int which);

    jChipmunk_Trap(jChipmunk_Habitat *habitat,
                   int width, int height,
                   // float vecs[], unsigned int count,
                   float flowForceH, float flowForceL, bool removeOutside);
    ~jChipmunk_Trap();
    bool outsideTrap(jChipmunk_EColi *cell);
    bool updateModel(jChipmunk_EColi *cell);
    
    double      compAccumulator;

private:
    bool outsideTrap2(jChipmunk_EColi *cell);

    cpSpace     *space;
    cpBody      *staticBody;
    cpShape     *trap[128];  //prob. should set dynamically
    int         trapSegments;

    int         highChannel, lowChannel;
    int         rightBoundary, leftBoundary;
    int         height, width;

    cpVect      flowForceHigh, flowForceLow;

    bool        removeOutsideTrap;
    bool        noTrap;
    int         whichTrap;

    friend bool setBinData(jEColi * c, int which);

};

#endif // JCHIPMUNK_TRAP_H
