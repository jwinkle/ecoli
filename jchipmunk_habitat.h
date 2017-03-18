#ifndef JCHIPMUNK_HABITAT_H
#define JCHIPMUNK_HABITAT_H

// #include "../Chipmunk-7.0.1/include/chipmunk/chipmunk.h"
#include "../Chipmunk-7.0.1/include/chipmunk/chipmunk_private.h"
#include "../Chipmunk-7.0.1/include/chipmunk/chipmunk_unsafe.h"
extern "C" {
    #include "../Chipmunk-7.0.1/include/chipmunk/cpHastySpace.h"
}


class jChipmunk_Trap;
class jChipmunk_EColi;
class jEColi;

class jChipmunk_Habitat
{
public:
    jChipmunk_Habitat(float sim_dt, unsigned int numThreads, unsigned int numIterations);
    ~jChipmunk_Habitat();

    unsigned long getNumThreads(void);
    void stepSimulation(float dt);
    bool outsideTrap(jChipmunk_EColi *);
    bool updateCell(jChipmunk_EColi *);

private:
    cpSpace         *space;
    cpFloat         dt;
    unsigned long   numThreads;


    jChipmunk_Trap *trap;

//allow cells to access the space pointer:
    friend class jChipmunk_EColi;
    friend class jEColi;
    friend class jChipmunk_Trap;
    
    friend bool setBinData(jEColi * c, int which);

};

#endif // JCHIPMUNK_HABITAT_H
