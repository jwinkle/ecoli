#include "jchipmunk_trap.h"
#include "jchipmunk_ecoli.h"

//#define TRAP_ELASTICITY     1.0f
#define TRAP_ELASTICITY     0.0f
#define TRAP_FRICTION       0.0f

// #define TRAP_BORDER_WIDTH  2.0f
#define TRAP_BORDER_WIDTH  10.0f

static  inline cpFloat frand(void) { return (cpFloat)rand()/(cpFloat)RAND_MAX; }

jChipmunk_Trap::jChipmunk_Trap(jChipmunk_Habitat *habitat, int which)
{
    noTrap=(0==which) ? true:false;
    space = habitat->space;
    habitat->trap = this;

    trapSegments = 0;
    cpShape *shape;


    whichTrap = which;
//no trap condition
    if(0==which) return;
//else, add jamming trap skeleton:

    compAccumulator = 0.0;
//JAMMING TRAP:
    staticBody = cpSpaceGetStaticBody(space);
#define JAMHEIGHT 125
#define JAMWIDTH  150
#define EDGEOUTLETX 35//50
#define EDGE1X  55//75//150
#define EDGE1Y  (JAMHEIGHT+50) //350
#define EDGE2Y  (EDGE1Y+JAMHEIGHT+50)
//THREE WALLS TO TRAP:   
        shape = cpSpaceAddShape(
                    space, cpSegmentShapeNew(
                        staticBody, cpv(-JAMWIDTH/2,0), cpv(-JAMWIDTH/2,-JAMHEIGHT), TRAP_BORDER_WIDTH));
            trap[trapSegments++] = shape;
        shape = cpSpaceAddShape(
                    space, cpSegmentShapeNew(
                        staticBody, cpv(-JAMWIDTH/2,-JAMHEIGHT), cpv(JAMWIDTH/2,-JAMHEIGHT), TRAP_BORDER_WIDTH));
            trap[trapSegments++] = shape;
        shape = cpSpaceAddShape(
                    space, cpSegmentShapeNew(
                        staticBody, cpv(JAMWIDTH/2,-JAMHEIGHT), cpv(JAMWIDTH/2,0), TRAP_BORDER_WIDTH));
            trap[trapSegments++] = shape;
//lower edges to outlet:
        shape = cpSpaceAddShape(
                    space, cpSegmentShapeNew(
                        staticBody, cpv(-JAMWIDTH/2,0), cpv(-EDGEOUTLETX,0), TRAP_BORDER_WIDTH));
            trap[trapSegments++] = shape;
        shape = cpSpaceAddShape(
                    space, cpSegmentShapeNew(
                        staticBody, cpv(JAMWIDTH/2,0), cpv(EDGEOUTLETX,0), TRAP_BORDER_WIDTH));
            trap[trapSegments++] = shape;
//initial path from outlet:  0->250
        shape = cpSpaceAddShape(
                    space, cpSegmentShapeNew(
                        staticBody, cpv(-EDGEOUTLETX,0), cpv(-EDGEOUTLETX,JAMHEIGHT), TRAP_BORDER_WIDTH));
            trap[trapSegments++] = shape;
        shape = cpSpaceAddShape(
                    space, cpSegmentShapeNew(
                        staticBody, cpv(EDGEOUTLETX,0), cpv(EDGEOUTLETX,JAMHEIGHT), TRAP_BORDER_WIDTH));
            trap[trapSegments++] = shape;
//final outlet: 10 um  500->600
        shape = cpSpaceAddShape(
                    space, cpSegmentShapeNew(
                        staticBody, cpv(-EDGEOUTLETX,EDGE2Y), cpv(-EDGEOUTLETX,EDGE2Y+100), TRAP_BORDER_WIDTH));
            trap[trapSegments++] = shape;
        shape = cpSpaceAddShape(
                    space, cpSegmentShapeNew(
                        staticBody, cpv(EDGEOUTLETX,EDGE2Y), cpv(EDGEOUTLETX,EDGE2Y+100), TRAP_BORDER_WIDTH));
            trap[trapSegments++] = shape;

    if(1==which)
    {
//edges from outlet:  no valve:  250->500
        shape = cpSpaceAddShape(
                    space, cpSegmentShapeNew(
                        staticBody, cpv(-EDGEOUTLETX,250), cpv(-EDGEOUTLETX,500), TRAP_BORDER_WIDTH));
            trap[trapSegments++] = shape;
        shape = cpSpaceAddShape(
                    space, cpSegmentShapeNew(
                        staticBody, cpv(EDGEOUTLETX,250), cpv(EDGEOUTLETX,500), TRAP_BORDER_WIDTH));
            trap[trapSegments++] = shape;
    }
    else if(2==which)
    {
            //left side
                shape = cpSpaceAddShape(
                            space, cpSegmentShapeNew(
                                staticBody, cpv(-EDGEOUTLETX,JAMHEIGHT), cpv(-EDGE1X,EDGE1Y), TRAP_BORDER_WIDTH));
                    trap[trapSegments++] = shape;
                shape = cpSpaceAddShape(
                            space, cpSegmentShapeNew(
                                staticBody, cpv(-EDGE1X,EDGE1Y), cpv(-EDGE1X,EDGE2Y), TRAP_BORDER_WIDTH));
                    trap[trapSegments++] = shape;
                shape = cpSpaceAddShape(
                            space, cpSegmentShapeNew(
                                staticBody, cpv(-EDGE1X,EDGE2Y), cpv(-EDGEOUTLETX,EDGE2Y), TRAP_BORDER_WIDTH));
                    trap[trapSegments++] = shape;
            //right side:
                shape = cpSpaceAddShape(
                            space, cpSegmentShapeNew(
                                staticBody, cpv(EDGEOUTLETX,JAMHEIGHT), cpv(EDGE1X,EDGE1Y), TRAP_BORDER_WIDTH));
                    trap[trapSegments++] = shape;
                shape = cpSpaceAddShape(
                            space, cpSegmentShapeNew(
                                staticBody, cpv(EDGE1X,EDGE1Y), cpv(EDGE1X,EDGE2Y), TRAP_BORDER_WIDTH));
                    trap[trapSegments++] = shape;
                shape = cpSpaceAddShape(
                            space, cpSegmentShapeNew(
                                staticBody, cpv(EDGE1X,EDGE2Y), cpv(EDGEOUTLETX,EDGE2Y), TRAP_BORDER_WIDTH));
                    trap[trapSegments++] = shape;
    }
    for(int i=0; i<trapSegments; i++)
    {
        cpShapeSetElasticity(trap[i], TRAP_ELASTICITY);
        cpShapeSetFriction(trap[i], TRAP_FRICTION);
    }    
}
jChipmunk_Trap::jChipmunk_Trap(jChipmunk_Habitat *habitat,
                               int w, int h,
                               // float vecs[], unsigned int count,
                               float flowForceH, float flowForceL, bool removeOutside)
{
    noTrap=false;
    whichTrap = -1;

    space = habitat->space;
    staticBody = cpSpaceGetStaticBody(space);

    //set only in +/- y direction for now:
    flowForceHigh = cpv(0.0, flowForceH);
    flowForceLow = cpv(0.0, -flowForceL);
//    flowForceHigh = cpv(flowForceH, 0.0);
//    flowForceLow = cpv(flowForceL, 0.0);
    removeOutsideTrap = removeOutside;

    width = w;
    height = h;

    trapSegments = 0;

//NO WALLS CASE:
#if 0
    cpShape *shape;

    // shape = cpSpaceAddShape(
    //             space, cpSegmentShapeNew(
    //                 staticBody, cpv(-jTRAP_SIZE,h), cpv(-w,h), TRAP_BORDER_WIDTH));
    //     trap[trapSegments++] = shape;
    shape = cpSpaceAddShape(
                space, cpSegmentShapeNew(
                    staticBody, cpv(-w,h), cpv(-w,-h), TRAP_BORDER_WIDTH));
        trap[trapSegments++] = shape;


//close the top of the trap:
    shape = cpSpaceAddShape(
                space, cpSegmentShapeNew(
                    staticBody, cpv(-w,-h), cpv(w,-h), TRAP_BORDER_WIDTH));
        trap[trapSegments++] = shape;
//#endif

//CODE TO ADD "BUMPS" TO THE TOP OF THE TRAP:
  // #define TRAPBUMPS
#if 0
    #define NUMBUMPS 30
   int tw = 2*width;
   for(int i=0;i<tw/NUMBUMPS;i++)
   {
       shape = cpSpaceAddShape(
                   space, cpCircleShapeNew(
                       staticBody, TRAP_BORDER_WIDTH/2.0, cpv(-w + NUMBUMPS*i, -h+TRAP_BORDER_WIDTH)));
           trap[trapSegments++] = shape;
   }
#endif


    shape = cpSpaceAddShape(
                space, cpSegmentShapeNew(
                    staticBody, cpv(w,-h), cpv(w,h), TRAP_BORDER_WIDTH));
        trap[trapSegments++] = shape;

    // shape = cpSpaceAddShape(
    //             space, cpSegmentShapeNew(
    //                 staticBody, cpv(w,h), cpv(jTRAP_SIZE,h), TRAP_BORDER_WIDTH));
    //     trap[trapSegments++] = shape;
#endif//NO-WALLS CASE

    for(int i=0; i<trapSegments; i++)
    {
        cpShapeSetElasticity(trap[i], TRAP_ELASTICITY);
        cpShapeSetFriction(trap[i], TRAP_FRICTION);
    }

    habitat->trap = this;
    highChannel = height;
    lowChannel = -height;

    rightBoundary = width+10;
    leftBoundary = -(width+10);
}

jChipmunk_Trap::~jChipmunk_Trap()
{
    if(false == noTrap)
    {
        for(int i=0; i<trapSegments; i++)
        {
            cpShapeFree(trap[i]);
        }
    }
}

bool jChipmunk_Trap::outsideTrap(jChipmunk_EColi *cell)
{
    if(noTrap) return false;

    float y = (float)(cell->center.y);
    float x = (float)(cell->center.x);


    if(1 == whichTrap)
    {
        if (y>500) return true;
    }
    else if(2 == whichTrap)
    {
        if (y>EDGE2Y+100) return true;
    }

    // if ( (y>highChannel+100) || (y<lowChannel-100))
    // else if ( (y>highChannel) || (y<lowChannel))
    else if ( (y>highChannel) || (y<lowChannel)
            || (x>rightBoundary) || (x<leftBoundary))
    {
        return true;
    }
    return false;
}

// bool jChipmunk_Trap::outsideTrap2(jChipmunk_EColi *cell)
// {
//     if(noTrap) return false;

//     float y = (float)(cell->center.y);

// //    if ( (y>highChannel-75) || (y<lowChannel+75))
//     if ( (y>highChannel) || (y<lowChannel))
//     {
//         return true;
//     }
//     else return false;
// }

bool jChipmunk_Trap::updateModel(jChipmunk_EColi *cell)
{
    if(noTrap) return false;

    // if(outsideTrap2(cell))
    if(outsideTrap(cell))
    {
        if(removeOutsideTrap)
        {
            ;
        }
        //else:  TODO:
        cpVect flow;
        float y = (float)(cell->center.y);

        if (y>highChannel)
        {
            flow = flowForceHigh;
        }
        else if (y<lowChannel)
        {
            flow = flowForceLow;
        }
//        return true;
        cell->applyForce(flow);
        return true;
    }
    else
    {
        if(2==whichTrap)
        {
            if (cell->center.y < 0)
            {
                compAccumulator += cell->compression;
            }
        }

//space-filling parabolas: y = 2a(1 - x^2) - 1
//
        //normalize to [-1,1] X [-1,1]
        // float x = cell->center.x/width;
        float y = cell->center.y/height;
        // float yr = (2.0*frand()-1.0);

        // float alpha = (-y + 1.0)/(2.0*(1.0 - x*x));

        //y = 2a(1-x^2/a^2) - 1     y' = -4x/a
        // float alpha2 = 0.25*((-y + 1) + sqrt((-y + 1)*(-y+1) + 16.0*x*x));


        //only add force if point inside or on a=1
//        if(alpha2 <= 1.0)
//        if( (alpha2 <= 0.75) && (alpha2 > 0.05))
        if( (frand() > 0.5) && (-y < 0.5) )
        {
#define FORCE_FACTOR 50.0
    //        cpVect force = cpvnormalize(cpv(xr, yr*4.0*alpha*x)) * cell->getScaledExpansionForce(10.5);
//            cpVect force = cpvnormalize(cpv(1.0, 4.0*x/alpha2)) * cell->getScaledExpansionForce(FORCE_FACTOR);
            // cpVect force = cpv(xr, y+1.0+frand()) * cell->getScaledExpansionForce(50.0);

//FORCE OPTION (LATEST)
        // float xr = (2.0*frand()-1.0);
            // cpVect force = cpv(xr, y+1.0+frand()) * cell->getScaledExpansionForce(50.0);
            // cell->applyForce(force);



//            cell->applyForce(force*frand());
            //turn ON/OFF flow:
        }
        return false;
    }
}

