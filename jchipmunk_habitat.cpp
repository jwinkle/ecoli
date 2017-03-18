#include "jchipmunk_habitat.h"
#include "jchipmunk_trap.h"

#define jITERATIONS 10
//#define jITERATIONS 20
// #define jITERATIONS 50

jChipmunk_Habitat::jChipmunk_Habitat(float sim_dt, unsigned int nthreads, unsigned int niterations)
{
    //parent = w;
    // space = cpSpaceNew();
    space = cpHastySpaceNew();
        cpHastySpaceSetThreads(space, (unsigned long)nthreads);
        numThreads = cpHastySpaceGetThreads(space);

    dt = sim_dt;

//    cpSpaceSetCollisionSlop(space, 0.01f);
    cpSpaceSetCollisionSlop(space, 0.1f);
    //    cpSpaceSetCollisionSlop(space, 0.2f);
//        cpSpaceSetCollisionSlop(space, 0.5f);

    cpSpaceSetCollisionBias(space, 0.1);
//    cpSpaceSetCollisionBias(space, 0.5);
//    cpSpaceSetCollisionPersistence(space, 6);//cpTimestamp value)

    cpSpaceUseSpatialHash(space, 10.0f, 20000);
//    cpSpaceUseSpatialHash(space, 15.0f, 10000);
    // cpSpaceSetIterations(space, jITERATIONS);
    cpSpaceSetIterations(space, niterations);

//a lower number for damping is more damping:
    //  1.0 ==> no damping.
    // epsilon->0.0 ==> full damping (must be > 0).
    // cpSpaceSetDamping(space, 0.8);
//    cpSpaceSetDamping(space, 1.0e-6f);
   cpSpaceSetDamping(space, 1.0);


//  trap is initialized if/when a trap is present (via friend class)
    trap = NULL;

}

jChipmunk_Habitat::~jChipmunk_Habitat()
{
    // cpSpaceFree(space);
    cpHastySpaceFree(space);
}
unsigned long jChipmunk_Habitat::getNumThreads(void)
{
    return numThreads;
}
bool jChipmunk_Habitat::outsideTrap(jChipmunk_EColi *cell)
{
    return (NULL != trap) ? trap->outsideTrap(cell) : false;
}

bool jChipmunk_Habitat::updateCell(jChipmunk_EColi *cell)
{
    return (NULL != trap) ? trap->updateModel(cell) : false;
}

void jChipmunk_Habitat::stepSimulation(float dt)
{
    // cpSpaceStep(space, cpFloat(dt));
    cpHastySpaceStep(space, cpFloat(dt));
}
