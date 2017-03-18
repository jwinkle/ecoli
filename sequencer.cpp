#include "main.h"

#include "jchipmunk_habitat.h"
#include "jchipmunk_trap.h"
#include "jEColi.h"

#include <omp.h>
#include <algorithm>
#include <math.h>

//use some globals for now, until code factored:
static jChipmunk_Habitat        *habitat;
static jChipmunk_Trap           *trap;
static std::list<jEColi *>      jcells;
static int                      jnext_id=0;

extern "C" {
    extern unsigned int getSynchError(void);
}

//#define TRAP_FORCE      1.0e15
#define TRAP_FORCE      1.0e14


// #define SIM_DT          0.001
// #define SIM_DT          0.002
// #define SIM_DT          0.004
// #define SIM_DT          0.005
#define SIM_DT          0.01
// #define SIM_DT          0.02
// #define SIM_DT          0.05

void initialize(float dt);//for non-trap simulations
void initialize(float, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int);
void destruct(void);
void jadd_cell (jEColi * c);
void jremove_cell ( std::list<jEColi *>::iterator &cellList, double simTime  );
void markBinTime(jEColi * c, double simTime, int which);
bool setBinData(jEColi * c, int which);
void recordDataMM(jEColi * thisCell, double &simTime);
void initDataMM(void);

extern unsigned int g_indexValue;
extern double       g_compressionLimit;
unsigned int        g_overComp;

#define MAX_SIM_THREADS     20
#define MAX_SIM_TIME        2000


// ################################
//      MAIN PROGRAM DEFINES:
// ################################   
//#define INIT_ON_FULLTRAP
// #define EXTINCTION_SIMULATION
// #define NOTRAP_SIMULATION
#define RECORD_BINDATA
// #define JAMMING_SIMULATION

// #define MOTHERMACHINE_SIM
    // #define MMCELLS  1
    // #define MMCELLS  2
    #define MMCELLS  4
    // #define MMCELLS  10


// #define TZ_BINSIZE 16
#define TZ_BINSIZE 32
// #define TZ_BINSIZE 64
double thetaZbin[TZ_BINSIZE][TZ_BINSIZE];
double compBin[TZ_BINSIZE][TZ_BINSIZE];
double compData[TZ_BINSIZE][TZ_BINSIZE];



#ifdef MOTHERMACHINE_SIM
    #define SPEED_SIZE  MMCELLS
    // bool notDoneSpeedRecord = true;
    std::ofstream speedFiles[SPEED_SIZE];
#endif

#ifdef JAMMING_SIMULATION
    std::ofstream compFile;
    unsigned int timeCounter = 0;
#endif

//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

//==============================================================================
//          int runSimulation():
//==============================================================================
int runSimulation(unsigned int height, unsigned int width, unsigned int numInitCells, double maxSimTime,
                    double &finalStrainRatio, double &elapsedSimTimeForSimulation, double &elapsedWallTimeForSimulation)
{

//===========================================
//  0.1  INITIALIZE :
//===========================================
#ifdef RECORD_BINDATA
    #define ST_RECORD_INTERVAL  20
    // #define ST_RECORD_INTERVAL  40
    // #define ST_RECORD_INTERVAL  10
    double nextSpacetime = ST_RECORD_INTERVAL;
    double timeStamps = maxSimTime/ST_RECORD_INTERVAL;
    int iStamps(ceil(timeStamps));
    double compSpaceTime[iStamps][TZ_BINSIZE + 1];
        for(int i=0;i<iStamps;i++) 
        {
            for(int j=0;j<TZ_BINSIZE+1;j++) 
            {
                compSpaceTime[i][j] = 0.0;
            }
        }
    int compIndex = 0;

#endif

    // g_compressionLimit = 10.0;// => 1.0um; RL-dot decays to 0 at 2um

#if defined(NOTRAP_SIMULATION) || defined(JAMMING_SIMULATION)  || defined(MOTHERMACHINE_SIM)
    initialize(SIM_DT);
#else   
    unsigned int numThreads = 1;
    unsigned int iterations = 10;
	initialize(SIM_DT, numThreads, iterations, 
                numInitCells, height, width);
#endif

    unsigned int    cellsDeleted = 0;
    double          cellEnergy, cEmax, cEmin;

    //establish time=0 for the simulation
    double          stepTime = 0.0;
    // double timeTare = omp_get_wtime();
    // double timeNow = 0.0;
    // double timeInterval = 0.0;
    // double minRecordTimeTare = 0.0;
    double          stepTimeTare = 0.0;
    // startTime = timeTare;

    double          simTime = 0.0;

    double          angleAccumulator = 0.0;
    unsigned int    cellCount = 0;
    unsigned int    countA = 0;
    unsigned int    countB = 0;

    double          pAlpha = 0.1;
    double          pBeta = 0.005;


//06OCT.2016:
//determine #cells for full-trap: (0.8 factor for partial density)
#ifdef INIT_ON_FULLTRAP        
        double fullTrapCells = 0.8*((width-20)*height)/(100.0*3.0);    
        bool fullTrapThresh = false;
        bool fullTrapInit = false;    
#endif        
//29OCT.2016:
        for(int i=0;i<TZ_BINSIZE;i++) 
        {
            for(int j=0;j<TZ_BINSIZE;j++) 
            {
                thetaZbin[i][j] = 0.0;
                compBin[i][j] = 0.0;                
                compData[i][j] = 0.0;                
            }
        }

#ifdef MOTHERMACHINE_SIM
    initDataMM();
#endif

    g_overComp = 0;
    std::list<jEColi *>::iterator jj;    

    #define DATA_STAMP_INTERVAL     100
    auto nextDataStamp = DATA_STAMP_INTERVAL;

//==============================================================================
//      MAIN SIMULATION LOOP:
//==============================================================================
	while (simTime <= maxSimTime)
	{
        #ifdef RECORD_BINDATA
           if(simTime > nextSpacetime)
            {
                nextSpacetime += ST_RECORD_INTERVAL;
                //average data()
                for ( jj=jcells.begin(); jj!=jcells.end(); jj++ )
                {
                    jEColi *thisCell = (*jj);
                    markBinTime(thisCell, simTime, 0);
                }
                //output matrix: first entry [0] is tims stamp, then N data for depths
                //average over "rows" for each depth bin
                for(int depth=0;depth<TZ_BINSIZE;depth++) 
                {
                    double rowSum = 0.0;
                    for(int j=0;j<TZ_BINSIZE;j++) 
                    {
                        //weighted average using 0-->1 as value of compression:
                        //depth dimension in binned data is second
                        //use depth+1 in output since time stamp lives in entry [0]
                        compSpaceTime[compIndex][depth+1] += j*compData[j][depth]/TZ_BINSIZE;
                        rowSum += compData[j][depth];
                        compData[j][depth] = 0.0;               
                    }
                    if(rowSum > 0.0)
                        compSpaceTime[compIndex][depth+1] /= rowSum;
                    else
                        compSpaceTime[compIndex][depth+1] = 0.0;
                }
                compSpaceTime[compIndex][0] = simTime;
                compIndex++;
                for ( jj=jcells.begin(); jj!=jcells.end(); jj++ )
                {
                    jEColi *thisCell = (*jj);
                    setBinData(thisCell, 0);
                }
            }
        #endif
        //some printing to stdout to check progress:
        if(simTime > nextDataStamp)
        {
            nextDataStamp += DATA_STAMP_INTERVAL;

            auto sr = -1.0;  
            if(0 != cellCount) sr = countA/cellCount;

            std::cout << std::fixed << std::setprecision(3)
                    <<"time = "<< simTime  << std::endl
                    <<"Ratio, CellCount: " << sr << ", " << cellCount << ";" 
                    <<" WallTime: " << stepTime/60.0 << std::endl;
        }
//===========================================
//  1.  INITIALIZE PER-ITERATION DATA:
//===========================================
        angleAccumulator        = 0.0;
        countA                  = 0;  
        countB                  = 0;
        cellCount               = 0;        
        cellEnergy              = 0.0;
        cEmax                   = 0.0;
        cEmin                   = 1.0e6;
//===========================================
//  2.  ITERATE CELL DATA STRUCTURE
//===========================================
        for ( jj=jcells.begin(); jj!=jcells.end(); jj++ )
        {
            jEColi *thisCell = (*jj);

    //CHECK NULL:
            if (NULL == thisCell)
            {
                continue;    
            } 
    //CELL UPDATE FUNCTION:
            // int status = thisCell->jupdate(); //I don't use status yet...
            thisCell->jupdate();

            thisCell->protein += SIM_DT*(pAlpha - pBeta*thisCell->protein); 

    //CELL DIVISION CHECK:
            jEColi * d = thisCell->jdivide();
            if ( d != NULL )
            {
                //should be put on a pending queue after the iterator is done?
                jadd_cell ( d );
                d->timeTare_thetaz = simTime;
                d->timeTare_comp = simTime;
                //std::cout << "divided cell added at: " << simTime << std::endl;
            }
    //DATA MARKING:
#ifdef NOTRAP_SIMULATION
            if(true == setBinData(thisCell, 3))
            {
                markBinTime(thisCell, simTime, 3);
            }
#endif
#ifdef RECORD_BINDATA
            if(true == setBinData(thisCell, 0))
            {
                markBinTime(thisCell, simTime, 0);
            }
            if(true == setBinData(thisCell, 1))
            {
                markBinTime(thisCell, simTime, 1);
            }
            if(true == setBinData(thisCell, 2))
            {
                markBinTime(thisCell, simTime, 2);
            }
#endif

    //ALLOW HABITAT TO ADD CELL FORCE (FLOW) OR REMOVE IF OUT OF BOUNDS
            // bool outsideTrap = habitat->updateCell(thisCell->cpCell);//I don't use the bool return value yet.
            habitat->updateCell(thisCell->cpCell);//I don't use the bool return value yet.
            //I use the return value of this one:
            if ( thisCell->marked_for_death() || habitat->outsideTrap(thisCell->cpCell) )
            {
                // std::cout << "Removed Cell:  " << thisCell->cpCell->center.x << ","<< thisCell->cpCell->center.y << std::endl;        
                //29OCT.2016:  added remove function
                jremove_cell(jj, simTime);
                // jj = jcells.erase ( jj );
                // delete thisCell;
                cellsDeleted++;
                //std::cout << "deleted cell("<<type<<")  at: " << simTime << ", ("<<x<<","<<y<<")"<<std::endl;
            }
            else
            {//cell not deleted, accumulate cell data for statistics:
                double comp = thisCell->cpCell->compression;  cellEnergy += comp * comp;
                    if (comp > cEmax) cEmax = comp;
                    if (comp < cEmin) cEmin = comp;


                float a = thisCell->cpCell->angle;  angleAccumulator += a;
#ifdef INIT_ON_FULLTRAP
                //assign cell strain when trap is full:
                if(fullTrapInit)
                {
                    thisCell->cellType = (0 == (rand() % 2)) ?  ECOLITYPE_01: ECOLITYPE_02;   
                }
#endif

                (thisCell->cellType == ECOLITYPE_01) ? countA++ : countB++;
                cellCount++;   

#ifdef MOTHERMACHINE_SIM
                recordDataMM(thisCell,  simTime);
#endif                 
            //UPDATE THE GROWTH FORCE:
                thisCell->cpCell->updateExpansionForce();

            }//end else (cell not deleted; per-cell data compute, record)
        }//end cell iterator loop
//===========================================
//  3.  COMPUTE PER TIME-STEP DATA, ACTIONS:
//===========================================
        if(0 != cellCount)
        {
            angleAccumulator /= float(cellCount);   

#ifdef JAMMING_SIMULATION       //RECORD DATA FOR JAMMING TRAP
    // #define COMPINTERVAL 1000            
    #define COMPINTERVAL 500            
            double compTotal = trap->compAccumulator;
            if((timeCounter++)%COMPINTERVAL == 0)
            {
                compTotal /= COMPINTERVAL;
                trap->compAccumulator = 0.0;
                        char str[1024];
                        // char fileStr[128];
                        sprintf(str, "%.3f,%.3f\n",
                                simTime, compTotal);
                        compFile.open("compFile.txt", std::ios::app);
                            compFile << std::string(str);
                        compFile.close();
            }
#endif
#ifdef MOTHERMACHINE_SIM        //END AFTER 2 DIVISIONS
            // if(cellCount > 2*MMCELLS) break;//break outside while loop
            if(simTime > 40.0) break;//break outside while loop
#endif
#ifdef NOTRAP_SIMULATION        //END AFTER N CELLS
            if(cellCount > 10000) break;//break outside while loop
#endif
#ifdef EXTINCTION_SIMULATION    //END AFTER EXTINCTION OF A STRAIN
            if( (0 == countA) || (0 == countB)) break;//break outside while loop
#endif
#ifdef INIT_ON_FULLTRAP        //INITIALIZE STRAINS ON FULL TRAP CONDITION:
            if(fullTrapInit)
            {
                fullTrapInit = false;
                fullTrapThresh = true;
                // std::cout<<". Done init. CellCount(thresh): "<<cellCount<<": "<<fullTrapCells<<std::endl;
            }
            //SIMULATION TERMINATION CHECK:
            if(fullTrapThresh)
                if( (0 == countA) || (0 == countB))
                    break;
            else
                if(cellCount > fullTrapCells)
                    fullTrapInit = true;//trigger one-time init sequence
#endif
        }//END if(0 != cellCount)


//===========================================
//  4.  UPDATE 2D PHYSICS ENGINE
//===========================================        
         stepTimeTare = omp_get_wtime();   
             habitat->stepSimulation(SIM_DT);
         stepTime += (omp_get_wtime() - stepTimeTare);
             
         simTime += SIM_DT;

	}// end    while (simTime <= maxSimTime)
//==============================================================================
//      END MAIN SIMULATION LOOP;  COMPUTE, RECORD SIMULATION DATA; TERMINATE:
//==============================================================================
    //WRITE OUTPUT DATA:
    finalStrainRatio = (0 != cellCount) ? (double)countA/cellCount : -1.0;    
    elapsedSimTimeForSimulation = simTime;
    elapsedWallTimeForSimulation = stepTime;
#ifdef RECORD_BINDATA
    for ( jj=jcells.begin(); jj!=jcells.end(); jj++ )
    {
        markBinTime((*jj),simTime,1);
        markBinTime((*jj),simTime,2);
    }

    std::cout << std::fixed;
    std::cout << std::setprecision(3); 
    for(int i=0;i<TZ_BINSIZE;i++) 
    {
        for(int j=0;j<TZ_BINSIZE;j++) 
            std::cout<<thetaZbin[i][j]<<", ";
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
    std::cout<<std::endl;
    for(int i=0;i<TZ_BINSIZE;i++) 
    {
        for(int j=0;j<TZ_BINSIZE;j++) 
            std::cout<<compBin[i][j]<<", ";
        std::cout<<std::endl;
    }
//THETA DATA IN SPACETIME
    std::string dateString = __TIME__;
    std::replace( dateString.begin(), dateString.end(), ':', '_'); // replace all 'x' to 'y'
    
    std::ofstream spaceTime;
    std::string fileStr = "./spaceTime/spaceTime_" 
                            + std::to_string(g_indexValue) + "_" 
                            + dateString + ".txt";
    spaceTime.open(fileStr, std::ios::trunc);
        for(int i=0;i<iStamps;i++) 
        {
            for(int j=0;j<TZ_BINSIZE+1;j++) 
            {
                spaceTime<<compSpaceTime[i][j]<<", ";
            }
            spaceTime<<std::endl;
        }
    spaceTime.close();

#endif    

    std::cout<<"__TIME__: "<<dateString<<std::endl;
    destruct();
    // return 0;
    return cellCount;
}
//==============================================================================
//      END    int runSimulation():
//==============================================================================

//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

//==============================================================================
//      SUPPORT FUNCTIONS:
//==============================================================================
#ifdef MOTHERMACHINE_SIM
void initDataMM(void)
{
    char fileStr[128];
    for (int i = 0; i < MMCELLS; ++i)
    {
            sprintf(fileStr, "./speedFolder/speedfile_%.4d.txt", i);
            speedFiles[i].open(fileStr, std::ios::trunc);
            speedFiles[i].close();
    }//DATA RECORD
}
void recordDataMM(jEColi * thisCell, double &simTime)
{

    double xsimTime = simTime;
    static bool notReset=true;

     if( (simTime >= 5.0) &&(simTime < 5.0+SIM_DT) &&(notReset))
     {
        static int count=0;
        count++; if(MMCELLS == count){ simTime = 0.0; notReset=false;}
            std::cout<<"SimTime=5.0: "<<count<<std::endl;
            xsimTime=0.0;
        thisCell->cpCell->springRestLength = thisCell->cpCell->length;
        thisCell->cpCell->compression = 0.0;
        thisCell->cpCell->setVelocity(cpvzero);
     }  
        char str[1024];
        char fileStr[128];
        int id = thisCell->get_id();
        if(id < SPEED_SIZE)
        {
            // unsigned int tint = (unsigned int)(t*1.0/get_sim_dt());
            // if(0 == tint%10)
            if(1)
            {
                sprintf(str, "%.3f,%.3f,%.3f,"
                             "%.3f,%.3f,%.3f,%.3f,"
                             "%.3f"
                             "\n",
                        // simTime,
                        xsimTime,
                        thisCell->cpCell->springRestLength,
                        thisCell->cpCell->length,
                        thisCell->cpCell->velA.x,//matlab(:,22)
                        thisCell->cpCell->velA.y,
                        thisCell->cpCell->velB.x,
                        thisCell->cpCell->velB.y,//matlab(:,25)
                        // thisCell->cpCell->angularVelocity
                        thisCell->cpCell->angle
                        );
                sprintf(fileStr, "./speedFolder/speedfile_%.4d.txt", id);
                speedFiles[id].open(fileStr, std::ios::app);
                    speedFiles[id] << std::string(str);
                speedFiles[id].close();
            }
        }//DATA RECORD
}
#endif

void markBinTime(jEColi * c, double simTime, int which)
{
    if(0 == which)
    {
        compData[c->tdold][c->ztdold] += (simTime - c->timeTare_tdz);
        c->timeTare_tdz = simTime;
        // std::cout<<"("<<c->get_id()<<") "<<c->thetaold<<","<<c->zold<<" : "<<simTime<<std::endl;
    }
    if(1 == which)
    {
        thetaZbin[c->thetaold][c->ztold] += (simTime - c->timeTare_thetaz);
        c->timeTare_thetaz = simTime;
        // std::cout<<"("<<c->get_id()<<") "<<c->thetaold<<","<<c->zold<<" : "<<simTime<<std::endl;
    }
    else if(2 == which)
    {
        compBin[c->compold][c->zcold] += (simTime - c->timeTare_comp);
        c->timeTare_comp = simTime;
    }
    else if(3 == which)
    {
        compBin[c->compold][c->zcold] += (simTime - c->timeTare_comp);
        c->timeTare_comp = simTime;
    }
}
bool setBinData(jEColi * c, int which)
{
    if(0 == which)
    {
        //this function interfaces the sim data matrix and the cell storage of theta/position
        double angle = c->cpCell->angle;        
        angle = fmod(angle, 3.141593);//2-fold symmetry (0-pi)
        if(angle < 0.0) angle = 3.141593 + angle;

        unsigned int binTheta = trunc((angle/3.141593) * TZ_BINSIZE);
            if(binTheta >= TZ_BINSIZE) binTheta = TZ_BINSIZE-1;
        // int traph = c->habitat->trap->height;
        // unsigned int binZ = trunc(((-c->cpCell->center.y + traph)/(2*traph)) * TZ_BINSIZE);
        //     if(binZ >= TZ_BINSIZE) binZ = TZ_BINSIZE-1;

    //oops, I recorded theta when I meant compression;  leaving all the vars for now...
        // unsigned int binComp = trunc((c->cpCell->compression /(2.0*g_compressionLimit)) * TZ_BINSIZE);
        //     if(binComp >= TZ_BINSIZE) binComp = TZ_BINSIZE-1;

        int traph = c->habitat->trap->height;
        unsigned int binZ = trunc(((-c->cpCell->center.y + traph)/(2*traph)) * TZ_BINSIZE);
            if(binZ >= TZ_BINSIZE) binZ = TZ_BINSIZE-1;

        c->tdold = c->itd;  c->ztdold  = c->iztd;
        // c->itd   = binComp;   c->iztd    = binZ;
        c->itd   = binTheta;   c->iztd    = binZ;
        //set return value to true if bins changed
        // if( (c->tdold == binComp) && (c->ztdold == binZ))
        if( (c->tdold == binTheta) && (c->ztdold == binZ))
            return false;
        else
            return true;
    }
    if(1 == which)
    {
        //this function interfaces the sim data matrix and the cell storage of theta/position
        double angle = c->cpCell->angle;        
        angle = fmod(angle, 3.141593);//2-fold symmetry (0-pi)
        if(angle < 0.0) angle = 3.141593 + angle;

        unsigned int binTheta = trunc((angle/3.141593) * TZ_BINSIZE);
            if(binTheta >= TZ_BINSIZE) binTheta = TZ_BINSIZE-1;
        int traph = c->habitat->trap->height;
        unsigned int binZ = trunc(((-c->cpCell->center.y + traph)/(2*traph)) * TZ_BINSIZE);
            if(binZ >= TZ_BINSIZE) binZ = TZ_BINSIZE-1;

        c->thetaold = c->itheta;  c->ztold  = c->izt;
        c->itheta   = binTheta;   c->izt    = binZ;
        //set return value to true if bins changed
        if( (c->thetaold == binTheta) && (c->ztold == binZ))
            return false;
        else
            return true;
    }
    else if(2 == which)
    {
        //this function interfaces the sim data matrix and the cell storage of theta/position
        // unsigned int binComp = trunc((c->cpCell->compression /(2000.0*SIM_DT)) * TZ_BINSIZE);
        unsigned int binComp = trunc((c->cpCell->compression /(2.0*g_compressionLimit)) * TZ_BINSIZE);
            if(binComp >= TZ_BINSIZE) binComp = TZ_BINSIZE-1;
        int traph = c->habitat->trap->height;
        unsigned int binZ = trunc(((-c->cpCell->center.y + traph)/(2*traph)) * TZ_BINSIZE);
            if(binZ >= TZ_BINSIZE) binZ = TZ_BINSIZE-1;

        c->compold  = c->icomp;  c->zcold   = c->izc;
        c->icomp    = binComp;   c->izc     = binZ;
        //set return value to true if bins changed
        if( (c->compold == binComp) && (c->zcold == binZ))
            return false;
        else
            return true;
    }
    else if(3 == which)
    {
        double x = c->cpCell->center.x;
        double y = c->cpCell->center.y;
        double radius = x*x + y*y;
        unsigned int binComp = trunc((c->cpCell->compression /(2.0*g_compressionLimit)) * TZ_BINSIZE);
            if(binComp >= TZ_BINSIZE) binComp = TZ_BINSIZE-1;

        // unsigned int binZ = trunc((radius/(500.0*500.0)) * TZ_BINSIZE);
        unsigned int binZ = trunc((radius/(800.0*800.0)) * TZ_BINSIZE);
            if(binZ >= TZ_BINSIZE) binZ = TZ_BINSIZE-1;

        c->compold=c->icomp;  c->zcold = c->izc;
        c->icomp = binComp;  c->izc = binZ;
        //set return value to true if bins changed
        if( (c->compold == binComp) && (c->zcold == binZ))
            return false;
        else
            return true;
    }
    return false;        
}
void jadd_cell ( jEColi * c )
{
    jcells.push_back ( c );
    c->set_id ( jnext_id++ );
#ifdef NOTRAP_SIMULATION
        setBinData(c,3);
#endif
#ifdef RECORD_BINDATA    
        setBinData(c,0);
        setBinData(c,1);
        setBinData(c,2);
#endif
}
void jremove_cell ( std::list<jEColi *>::iterator &pCell, double simTime )
{
    jEColi * cell = (*pCell); 
#ifdef NOTRAP_SIMULATION
        markBinTime(cell, simTime, 3);
#endif
#ifdef RECORD_BINDATA    
        markBinTime(cell, simTime, 0);
        markBinTime(cell, simTime, 1);
        markBinTime(cell, simTime, 2);
#endif
    pCell = jcells.erase ( pCell );
    delete cell;    
}

void initialize(float dt)
{    
    jcells.clear();
    jnext_id = 0;

    // habitat = new jChipmunk_Habitat(dt, nt, ni);
    habitat = new jChipmunk_Habitat(dt, 1, 10);
    // std::cout << (unsigned int) habitat->getNumThreads() << " thread(s) will be active." << std::endl;
#ifdef MOTHERMACHINE_SIM
    trap = new jChipmunk_Trap(habitat,
                              30/2, 240/2, //MM
                              // 300/2, 240/2,
                              // 30/2, 400/2,
                              // 30/2, 800/2,
                              TRAP_FORCE, TRAP_FORCE, false);
#endif
#ifdef JAMMING_SIMULATION        
    //init special trap for jamming:
    trap = new jChipmunk_Trap(habitat, 2);
#endif

#ifdef MOTHERMACHINE_SIM
    #define CL 20
    //each cell is 2um=20pixels long; +10 pixels for trap "edge"
    //  RIGHT EDGE:|10| +20 : +20 : +20 : +20  = 90 (< 1/2 = 120) 
    //  CENTER:    |10| +10 : +20 : +20 : +20  => mother center @ -100pixels 
    //4-CELL MOTHER MACHINE:
    // int cellys[MMCELLS] = {-5*CL,-4*CL,-3*CL,-2*CL};
    //4-CELL W/ 1UM GAP:
    // -90, -60, -30, 0 
    #define MMOFFSET 8    
    // #define MMOFFSET 15
    // int cellys[MMCELLS] = {-5*CL+MMOFFSET,-4*CL+MMOFFSET,-3*CL+MMOFFSET,-2*CL+MMOFFSET};
    // int cellys[MMCELLS] = {-5*CL+MMOFFSET,-4*CL+MMOFFSET};
    // int cellys[MMCELLS] = {-5*CL,-4*CL};
    // int cellys[MMCELLS] = {-5*CL+MMOFFSET};
    
    // int cellys[MMCELLS] = {-5*CL+MMOFFSET, 0};
    // int cellys[MMCELLS] = {-5*CL};
    // int cellys[MMCELLS] = {100};
    int cellys[MMCELLS] = {-5*CL,-4*CL,-3*CL,-2*CL};
    unsigned int cellCount = MMCELLS;
#else
    int rx=0, ry=0;
    unsigned int cellCount = 4;
    int rxv[cellCount] = {-15,15, 30,50};
#endif

    //TRAP INIT:
#ifdef JAMMING_SIMULATION        
    ry = -50;
#endif

    // int x2[2] = {0, 10};
    // int y2[2] = {0, -20};
    jEColi *firstCell;
    //INITIALIZE CELLS IN TRAP:
    for(unsigned int icell = 0;icell<cellCount;icell++)
    {
#ifdef MOTHERMACHINE_SIM
        firstCell = new jEColi(habitat,
            0, cellys[icell], 3.14159265/2.0,
            // x2[icell], y2[icell], double(icell)*(3.14159265/2.0),
            DEFAULT_ECOLI_INIT_SIZE,
            ECOLITYPE_01);
#else
        rx=rxv[icell];
        firstCell = new jEColi(habitat,
            rx, ry, 3.1415/2.0,
            // rx, ry, (rand()%31415)/10000.0,
            DEFAULT_ECOLI_INIT_SIZE,
            (icell%2) == 0 ? ECOLITYPE_01: ECOLITYPE_02);
#endif

        jadd_cell(firstCell);
     // std::cout << "Cell y,x:  " << ry <<", "<< rx << std::endl;
   }

}
void initialize(float dt, unsigned int nt, unsigned int ni, 
                    unsigned int cellCount, unsigned int trapHeight, unsigned int trapWidth)
{
	// std::cout 
 //        << "jinit -- Simulation dt = " << SIM_DT 
 //        << " Trap Size:" << TRAP_W << ", " << TRAP_H 
 //        << std::endl;

    jcells.clear();
    jnext_id = 0;

    habitat = new jChipmunk_Habitat(dt, nt, ni);
        // std::cout << (unsigned int) habitat->getNumThreads() << " thread(s) will be active." << std::endl;
    trap = new jChipmunk_Trap(habitat,
                              trapWidth/2, trapHeight/2,
                              TRAP_FORCE, TRAP_FORCE, false);
    int rx, ry, rw, rh;
    if(trapWidth > 28)     rw = trapWidth - 28; else rw = 0;
    if(trapHeight > 10)    rh = trapHeight - 10; else rh = 20;

    int cellArray[cellCount];

    jEColi *firstCell;
    // srand (time(NULL));

    //INITIALIZE CELLS IN TRAP:
    for(unsigned int icell = 0;icell<cellCount;icell++)
    {
        //This loop ensures cells don't overlap (< 1.5um vertical separation)
        int cellScan;
        do
        {
            rx = (rand() % rw) - rw/2;
            ry = (rand() % rh)- rh/2;

            cellScan = 1;
            for(unsigned int i=0;i<icell;i++)
            {
                if(abs(cellArray[i] - ry) < 15)
                {
                    cellScan = 0;
                    break;    //break the for loop, try new rx,ry
                } 
            }
        } while(0 == cellScan);
		//std::cout << rx << ", " << ry << std::endl;
        
        //good cell init location, record y-value for next scan:
        cellArray[icell] = ry;
        //create and add new cell at this position
        firstCell = new jEColi(habitat,
            // rx, ry, (rand()%31415)/10000.0,
            rx, ry, 3.1415/2.0,
               DEFAULT_ECOLI_INIT_SIZE,
               (icell%2) == 0 ? ECOLITYPE_01: ECOLITYPE_02);

        jadd_cell(firstCell);
    // std::cout << "Cell y,x:  " << ry <<", "<< rx << std::endl;
    }
	// std::cout << icell << " cells added to space." << std::endl;    
}
void destruct(void)
{
	// unsigned int numCells = 0;
    std::list<jEColi *>::iterator jj;
    for ( jj=jcells.begin(); jj!=jcells.end(); jj++ ) 
    {
        delete (*jj);
        // numCells++;
    }

	delete trap;
    delete habitat;
    // std::cout << "jdestruct: " << numCells << "cells destructed" << std::endl;
    // std::cout << "Elapsed time(hrs): " << (omp_get_wtime() - startTime)/3600.0 << std::endl;
}

