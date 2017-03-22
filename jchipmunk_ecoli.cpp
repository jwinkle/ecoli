#include "jchipmunk_ecoli.h"
#include "jchipmunk_habitat.h"
#include "jEColi.h"

//too low elast. results in overlap?:
#define SHAPE_ELASTICITY     0.1f
//#define SHAPE_ELASTICITY     1.0f
#define SHAPE_FRICTION     0.0f
//#define SHAPE_FRICTION     0.8f

#define GROOVE_JOINT_MAXFORCE   1.0e6
//#define GROOVE_JOINT_MAXFORCE   5000.0
//#define GROOVE_JOINT_MAXFORCE   2500.0
//#define GROOVE_JOINT_MAXFORCE   1024.0
#define GROOVE_JOINT_ERRORBIAS   0.5
#define GROOVE_JOINT_MAXBIAS    40.0


//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#define KSPRING_SECS    1.0         //normalize to 1.0  
#define KSPRING        (KSPRING_SECS*3600.0)    //convert to min^-2


// #define SCALE_DT    1.0
// #define SCALE_DT    2.0
// #define SCALE_DT    5.0
#define SCALE_DT    10.0
// #define SCALE_DT    100.0

// #define GAMMA_FLUID         (SCALE_DT * 0.25e0 * KSPRING_SECS*60.0) //dt=0.001
#define GAMMA_FLUID         (SCALE_DT  * KSPRING_SECS*60.0) //dt=0.001
//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||


#define RATCHET_QUANTUM    0.5
#define COMPRESSION_GAP    1.0
//#define RATCHET_QUANTUM    0.1
//#define COMPRESSION_GAP    0.1

extern double       g_lambdaFitness;
extern double       g_compressionLimit;
extern unsigned int g_indexValue;
extern unsigned int g_overComp;

//==============================================================================
//==============================================================================
void jcpBodyUpdateVelocity(cpBody *body, cpVect gravity, cpFloat damping, cpFloat dt)
{
    // Skip kinematic bodies.
    if(cpBodyGetType(body) == CP_BODY_TYPE_KINEMATIC) return;

    body_t *pBody =  (body_t *)cpBodyGetUserData(body);

    cpAssertSoft(body->m > 0.0f && body->i > 0.0f, 
        "Body's mass and moment must be positive to simulate. (Mass: %f Moment: %f)", 
        body->m, body->i);
    
    // cpVect Fbym = body->f * body->m_inv;
    cpVect Fbym     = body->f;
    body->v         = Fbym * (1.0/GAMMA_FLUID);
    body->w         = body->w*damping + body->t*body->i_inv*dt;
    
    //transfer the force for this step to be recorded:
        pBody->parent->setData(body->f, body->v, pBody->which);
    // Reset forces.
    body->f = cpvzero;
    body->t = 0.0f;
//    cpAssertSaneBody(body);
}
//==============================================================================
//==============================================================================

//==============================================================================
//==============================================================================
jChipmunk_EColi::jChipmunk_EColi(jChipmunk_EColiInit_t *ecoli)
{
/*
     * Constructor for the chipmunk model.
     * Input data is through Init structure pointer
     * Takes 'length' parameter as actual cell length and forms the
     * model from this and a 'width' -- fixed at 10 pixels = 1um.
     *
*/
    parent              = ecoli->owner;
    habitat             = ecoli->habitat;
    space               = ecoli->habitat->space;

    initialMass         = ecoli->mass;
    initialMoment       = ecoli->moment;
    shapeLength         = ecoli->length - ecoli->width;
    shapeHeight         = ecoli->width;
    length              = ecoli->length;
    initialLength       = ecoli->length;
    center              = cpv(ecoli->x, ecoli->y);
    angle               = ecoli->angle;
    velocity            = cpv(ecoli->vx, ecoli->vy);
    angularVelocity     = ecoli->av;

    posA                = center;
    posB                = center;
    separationDistance  = 0.0;
    separationVelocity  = cpvzero;

    compression         = 0.0;
    springRestLength    = length;

    dRL = ecoli->dRL;

    radius = shapeHeight * 0.5;
    offset = shapeLength * 0.5;
//    nextRatchet = 0.0;
    nextRatchet = RATCHET_QUANTUM;

//GROOVE JOINT:
    cpVect upperLeft = cpv(-offset, radius);
    cpVect upperRight = cpv(offset, radius);
    cpVect lowerRight = cpv(offset, -radius);
    cpVect lowerLeft = cpv(-offset, -radius);

//==============================================================================
//      CELL BODY HALVES:  mass,moment,center,angle,vel.,angular_vel
//==============================================================================
    bodyA = cpSpaceAddBody(
                space, cpBodyNew(
                    initialMass, initialMoment));
    bodyB = cpSpaceAddBody(
                space, cpBodyNew(
                    initialMass, initialMoment));
    cpBodySetPosition(bodyA, center);
        cpBodySetAngle(bodyA, angle);
            cpBodySetVelocity(bodyA, velocity);
    cpBodySetPosition(bodyB, center);
        cpBodySetAngle(bodyB, angle);
            cpBodySetVelocity(bodyB, velocity);

    cpBodySetVelocityUpdateFunc(bodyA,(cpBodyVelocityFunc) jcpBodyUpdateVelocity );
    cpBodySetVelocityUpdateFunc(bodyB,(cpBodyVelocityFunc) jcpBodyUpdateVelocity );

    pA.parent      = ecoli->owner; pA.which=0;
    pB.parent      = ecoli->owner; pB.which=1;
    cpBodySetUserData(bodyA,  (cpDataPointer *) &pA);
    cpBodySetUserData(bodyB,  (cpDataPointer *) &pB);

    velA = velocity;
    velB = velocity;
//==============================================================================
//      CELL SHAPES:  box(length,width), pole(radius,offset)
//==============================================================================
//vertices of initial box, for updates to growth:
    vertsA[0]=upperLeft;
    vertsA[1]=upperRight;
    vertsA[2]=lowerRight;
    vertsA[3]=lowerLeft;
        vertsB[0]=upperLeft;
        vertsB[1]=upperRight;
        vertsB[2]=lowerRight;
        vertsB[3]=lowerLeft;
//BOX CENTER FOR CELL:
    //create and add the box and circle shapes on the body:
    shapeCount = 0;
    boxA = cpSpaceAddShape(
                space, cpBoxShapeNew(
                  bodyA, shapeLength, shapeHeight, 0.0f));
    shapes[shapeCount++] = boxA;
    boxB = cpSpaceAddShape(
                space, cpBoxShapeNew(
                  bodyB, shapeLength, shapeHeight, 0.0f));
    shapes[shapeCount++] = boxB;
//CIRC ENDS FOR CELL:
//"primary" poles for the cell halves (expansion-side):
    capA = cpSpaceAddShape(
                space, cpCircleShapeNew(
                    bodyA, radius, cpv(-offset, 0.0)));
    shapes[shapeCount++] = capA;
    capB = cpSpaceAddShape(
                space, cpCircleShapeNew(
                    bodyB, radius, cpv(offset, 0.0)));
    shapes[shapeCount++] = capB;
//"secondary" poles (back-filled side):
    capA2 = cpSpaceAddShape(
                space, cpCircleShapeNew(
                    bodyA, radius, cpv(offset, 0.0)));
    shapes[shapeCount++] = capA2;
    capB2 = cpSpaceAddShape(
                space, cpCircleShapeNew(
                    bodyB, radius, cpv(-offset, 0.0)));
    shapes[shapeCount++] = capB2;
    for(unsigned int i=0;i<shapeCount;i++)
    {
        cpShapeSetElasticity(shapes[i], SHAPE_ELASTICITY);
        cpShapeSetFriction(shapes[i], SHAPE_FRICTION);
    }
//============================================================
//      CELL GROOVE JOINTS:
//============================================================
    constraintCount = 0;
//GROOVE JOINTS
    //lock the two bodies together using symmetric groove joints (and set to not collide the two bodies):
    grooveJointA = cpSpaceAddConstraint(
                space, cpGrooveJointNew(
                    bodyA,bodyB, upperLeft, upperRight, upperLeft));
    constraints[constraintCount++] = grooveJointA;
    grooveJointB = cpSpaceAddConstraint(
                space, cpGrooveJointNew(
                    bodyB,bodyA, lowerLeft, lowerRight, lowerRight));
    constraints[constraintCount++] = grooveJointB;

    grooveJointA2 = cpSpaceAddConstraint(
                space, cpGrooveJointNew(
                    bodyA,bodyB, lowerLeft, lowerRight, lowerLeft));
    constraints[constraintCount++] = grooveJointA2;
    grooveJointB2 = cpSpaceAddConstraint(
                space, cpGrooveJointNew(
                    bodyB,bodyA, upperLeft, upperRight, upperRight));
    constraints[constraintCount++] = grooveJointB2;

    for(unsigned int i=0;i<constraintCount;i++)
    {
        cpConstraintSetCollideBodies(constraints[i], cpFalse);
        cpConstraintSetMaxForce(constraints[i], GROOVE_JOINT_MAXFORCE);
        cpConstraintSetErrorBias(constraints[i], GROOVE_JOINT_ERRORBIAS);
        cpConstraintSetMaxBias(constraints[i], GROOVE_JOINT_MAXBIAS);
    }
}//end jChipmunk_EColi::jChipmunk_EColi()
//==============================================================================
//==============================================================================

//==============================================================================
//==============================================================================
jChipmunk_EColi::~jChipmunk_EColi(void)
{
    for(unsigned int i=0;i<constraintCount;i++)
    {
        cpSpaceRemoveConstraint(space, constraints[i]);
        cpConstraintFree(constraints[i]);
    }
    for(unsigned int i=0;i<shapeCount;i++)
    {
        cpSpaceRemoveShape(space, shapes[i]);
        cpShapeFree(shapes[i]);
    }
    cpSpaceRemoveBody(space, bodyA);
    cpSpaceRemoveBody(space, bodyB);
    cpBodyFree(bodyA);
    cpBodyFree(bodyB);
}


void jChipmunk_EColi::setVelocity(cpVect vel)
{
    cpBodySetVelocity(bodyA, vel); velA = vel;
    cpBodySetVelocity(bodyB, vel); velB = vel;
}
void jChipmunk_EColi::setBodyVelocities(cpVect vA, cpVect vB)
{
    cpBodySetVelocity(bodyA, vA); velA =vA;
    cpBodySetVelocity(bodyB, vB); velB=vB;
}

void jChipmunk_EColi::applyForce(cpVect force)
{//adds to the current force on the object
    cpBodyApplyForceAtWorldPoint(bodyA,  force, cpBodyLocalToWorld(bodyA, cpvzero));
    cpBodyApplyForceAtWorldPoint(bodyB,  force, cpBodyLocalToWorld(bodyB, cpvzero));
}
float jChipmunk_EColi::getScaledExpansionForce(float scale)
{
    return dRL*KSPRING*scale;
}

//==============================================================================
//==============================================================================
int jChipmunk_EColi::updateModel(float growthRate)
{
    //read the post-step data from the chipm. model:
    //position, angle, velocity, angularVelocity

    posA = cpBodyGetPosition(bodyA);
    posB = cpBodyGetPosition(bodyB);
        center = (posA + posB) * 0.5;                   //center of the cell wrt. the body positions
        separationDistance = cpvdist(posA, posB);       //distance between centers of mass
        length = separationDistance + initialLength;    //length of the whole cell, incl. poles
    angleA = (cpBodyGetAngle(bodyA));
    angleB = (cpBodyGetAngle(bodyB));
        angle = (angleA + angleB) * 0.5;        //compass heading wrt. body angles (average)
    velA = cpBodyGetVelocity(bodyA);
    velB = cpBodyGetVelocity(bodyB);
        velocity = (velA + velB) * 0.5;                 //vel. of center of mass of cell
    angularVelocity = (cpBodyGetAngularVelocity(bodyA)
                       + cpBodyGetAngularVelocity(bodyB)) * 0.5;

    separationAngle         = cpvtoangle(posB - posA);          //compass heading wrt. centers of mass
    separationVelocity      = (velB - velA); //vel. wrt. cell frame (expansion velocity); uses bodyB as + dir.
    separationVelocityAngle = cpvtoangle(velB - velA);//compass heading wrt sepraration vel.
    separationSpeed         = cpvlength(separationVelocity);
    projectionSpeed         = cpvdot(separationVelocity, cpvnormalize(velB));//signed speed of separation
    cellExpanding           = (projectionSpeed >= 0.0) ? true:false;
//    if(false == cellExpanding) separationSpeed = -separationSpeed;

//============================================================
//  RATCHET ALGORITHM:
//============================================================
    if(separationDistance > (nextRatchet + COMPRESSION_GAP))
//        if(separationDistance > (nextRatchet + gaps[icell]))
    {
        cpFloat newOffset = offset + nextRatchet;
        nextRatchet = nextRatchet + RATCHET_QUANTUM;
//update the rectangle for each cell half (only back-filled half):
    //        vertsA[0]=upperLeft;
            vertsA[1]=cpv(newOffset, radius);
            vertsA[2]=cpv(newOffset, -radius);
    //        vertsA[3]=lowerLeft;
            vertsB[0]=cpv(-newOffset, radius);;
    //        vertsB[1]=upperRight;
    //        vertsB[2]=lowerRight;
            vertsB[3]=cpv(-newOffset, -radius);;

// #include "chipmunk_unsafe.h"
//redraw the shapes. NB: these should not generate new collisions (on non-expanding side)
    cpPolyShapeSetVertsRaw(boxA, 4, vertsA);
    cpPolyShapeSetVertsRaw(boxB, 4, vertsB);
    cpCircleShapeSetOffset(capA2, cpv(newOffset, 0.0));
    cpCircleShapeSetOffset(capB2, cpv(-newOffset, 0.0));
    //expand the grooves and anchors to match new lengths for shapes:
    cpGrooveJointSetGrooveB(grooveJointA, vertsA[1]);
        cpGrooveJointSetAnchorB(grooveJointA, vertsB[0]);
    cpGrooveJointSetGrooveB(grooveJointA2, vertsA[2]);
        cpGrooveJointSetAnchorB(grooveJointA2, vertsB[3]);
    cpGrooveJointSetGrooveA(grooveJointB, vertsB[3]);
        cpGrooveJointSetAnchorB(grooveJointB, vertsA[2]);
    cpGrooveJointSetGrooveA(grooveJointB2, vertsB[0]);
        cpGrooveJointSetAnchorB(grooveJointB2, vertsA[1]);
    }

    compression = springRestLength - length;

    return 0;
}
//============================================================
//  SPRING EXPANSION FORCE ALGORITHM:
//============================================================
int jChipmunk_EColi::updateExpansionForce(void)
{//cell update will have been called before this (updates compression)
//at dt=0.001, dRL increment is 0.0001 um => 20k x dRL = 2um (doubling)
    cpVect springForce;
    // cpFloat climit = (1.0/SCALE_DT)*500.0*(4.0)*dRL;
    // cpFloat climit = 100.0; //=5x2um at dt=0.001
    cpFloat climit = g_compressionLimit; 

        //REMOVES SCLALING:
        // if(true)
        if (compression < climit)
        {
            springForce = cpv(KSPRING * (compression+dRL), 0.0f);
            springRestLength += dRL;
        }
        else
        {
            g_overComp++;
            cpFloat c = (compression-climit)/climit;
            if (c > 1.0) c=1.0;
            cpFloat scaled_dRL = (1.0-c)*dRL;
            springForce = cpv(KSPRING * (compression + scaled_dRL), 0.0f);
            springRestLength += scaled_dRL;
        }

    cpBodyApplyForceAtLocalPoint(bodyA,  cpvneg(springForce), cpvzero);
    cpBodyApplyForceAtLocalPoint(bodyB,  springForce, cpvzero);
    // cpBodyApplyForceAtLocalPoint(bodyA,  cpvneg(springForce)* 0.5, cpvzero);
    // cpBodyApplyForceAtLocalPoint(bodyB,  springForce* 0.5, cpvzero);

    return 0;
}
