#ifndef JCHIPMUNK_ECOLI_H
#define JCHIPMUNK_ECOLI_H

#include "jchipmunk_habitat.h"
class jEColi;

//1/2 the size of the trap width:
#define jTRAP_SIZE 1000
// #define jTRAP_SIZE 500

typedef struct{
    int which;
    jEColi *parent;
} body_t;

typedef struct jChipmunk_EColiInit_t
{
    jChipmunk_Habitat       *habitat;
    jEColi                  *owner;
    float                   mass, moment;
    float                   x,y;
    float                   angle, length, width;
    float                   vx, vy, av;
    float                   dRL;

} jChipmunkEColiInit_t;

class jChipmunk_EColi
{
public:
    jChipmunk_EColi(jChipmunk_EColiInit_t *ecoli);
    ~jChipmunk_EColi(void);

    int     updateModel(float growthRate);
    void    setVelocity(cpVect vel);
    void    setBodyVelocities(cpVect vA, cpVect vB);
    void    applyForce(cpVect);
    float   getScaledExpansionForce(float scale);
    int     updateExpansionForce(void);

private:
    void    applyGrowthForce(cpFloat growthRate, bool expanding);


    cpVect              vertsA[4];
    cpVect              vertsB[4];
    cpFloat             springRestLength;
    cpFloat             dRL;
    cpFloat             compression;

    cpVect               posA, posB;
    cpVect               velA, velB;
    cpFloat              separationDistance, separationAngle, separationVelocityAngle;
    cpFloat              separationSpeed, projectionSpeed;
    bool                 cellExpanding;
    cpVect               center;
    cpFloat              angleA,angleB,angle;
    cpFloat              length;
    cpVect               velocity, separationVelocity;
    cpFloat              angularVelocity;
    cpVect              forceVE;
    body_t              pA,pB;

    cpFloat             radius, offset;
    cpFloat             nextRatchet;
    cpFloat             shapeLength;
    cpFloat             shapeHeight;
    cpFloat             initialLength;
    cpFloat             initialMass, initialMoment;

    jEColi              *parent;
    jChipmunk_Habitat   *habitat;
    cpSpace             *space;
    cpBody              *bodyA, *bodyB;

        cpShape           *shapes[32];
        unsigned int      shapeCount;
    cpShape             *capA, *capB;
    cpShape             *capA2, *capB2;
    cpShape             *boxA, *boxB;
        cpConstraint        *constraints[32];
        unsigned int      constraintCount;
    cpConstraint        *grooveJointA, *grooveJointB;
    cpConstraint        *grooveJointA2, *grooveJointB2;
    cpConstraint        *gearJoint;
    cpConstraint        *springJoint1, *springJoint2;

    bool        firstIterationFlag;
    cpFloat     vc, v0;

    friend class jChipmunk_Trap;
    friend class jEColi;
    friend class World;

    friend void jcpBodyUpdateVelocity(cpBody *body, cpVect gravity, cpFloat damping, cpFloat dt);
    friend int main(int argc, char **argv);
    friend void recordData(jEColi *thisCell, float t);
    friend int runSimulation(unsigned int height, unsigned int width, unsigned int numInitCells, double maxSimTime,
                    double &finalStrainRatio, double &elapsedTimeForSimulation, double &);
    friend bool setBinData(jEColi * c, int which);
    friend void recordDataMM(jEColi * thisCell, double &simTime);

};


#endif // JCHIPMUNK_ECOLI_H
