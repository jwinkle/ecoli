#include "main.h"

#include <omp.h>
#include <algorithm>


extern int runSimulation(unsigned int height, unsigned int width, unsigned int numInitCells, double maxSimTime,
                    double &finalStrainRatio, double &elapsedTimeForSimulation, double &elapsedWallTimeForSimulation);

unsigned int g_indexValue;
double g_lambdaFitness = 1.0;
double g_compressionLimit = 1.0;

// #define NUM_LAMBDAS 6
#define NUM_LAMBDAS 4
// double lambdas[NUM_LAMBDAS] = {1.0, 1.05, 1.1, 1.15, 1.2, 1.25};
//max compression cutoffs:
//600x600 initial set of data
// double lambdas[NUM_LAMBDAS] = {10.0, 20.0, 30.0, 40.0, 60.0, 80.0};
//800x800 data, new cut-offs:
// double lambdas[NUM_LAMBDAS] = {2.0, 4.0, 8.0, 10.0, 20.0, 40.0};
// double lambdas[NUM_LAMBDAS] = {1.0, 2.0, 4.0, 8.0};
double lambdas[NUM_LAMBDAS] = {2.0, 4.0, 6.0, 8.0};
// double lambdas[NUM_LAMBDAS] = {4.0, 8.0, 12.0, 16.0};
// double lambdas[NUM_LAMBDAS] = {2.0, 8.0, 16.0, 32.0};
// double lambdas[NUM_LAMBDAS] = {8.0, 12.0, 16.0, 20.0};
// double lambdas[NUM_LAMBDAS] = {20.0, 30.0, 40.0, 60.0};

// double lambdas[NUM_LAMBDAS] = {12.0, 14.0, 16.0, 32.0};

#define numArgs 5
const char *argErrorStrings[] = {
    "argv[0]",
    "BAD TRAP HEIGHT",
    "BAD TRAP WIDTH",
    "BAD INITIAL CELLS",
    "BAD SIM NUMBER",
    "BAD SIM TIME",
    "BAD LAMBDA INDEX",

    "END"
};


int main(int argc, char **argv)
{
// #threads, simTime, numIterations
 
    double strainRatio, elapsedtime, wallTime;
    unsigned int trapHeight = 50;
    unsigned int trapWidth = 50;
    unsigned int numInitCells = 2;
    unsigned int numSims = 1;
    double maxSimTime = 200.0;

    #define ERROR_BAD_PARAM -1
    int argCounter = 0;
    std::cout << "argc = "<<argc<<" argv[] = "<<std::endl;
    for(int i=0;i<argc;i++) {std::cout<<argv[i]<<" ";}
    std::cout << std::endl;

    if(argc > ++argCounter)
    {
        trapHeight = atoi(argv[argCounter]);
        if ((trapHeight >= 20) && (trapHeight <= 5000)) {}
        else{std::cout<< argv[argCounter] << argErrorStrings[argCounter] << std::endl;  return ERROR_BAD_PARAM;} 
    }
    if(argc > ++argCounter)
    {
        trapWidth = atoi(argv[argCounter]);
        if ((trapWidth >= 20) && (trapWidth <= 5000)) {}
        else{std::cout<< argv[argCounter] << argErrorStrings[argCounter] << std::endl;  return ERROR_BAD_PARAM;} 
    }
    if(argc > ++argCounter)
    {
        numInitCells = atoi(argv[argCounter]);
        if ((numInitCells >= 1) && (numInitCells <= 1000)) {}
        else{std::cout<< argv[argCounter] << argErrorStrings[argCounter] << std::endl;  return ERROR_BAD_PARAM;} 
    }
    if(argc > ++argCounter)
    {
        numSims = atoi(argv[argCounter]);
        if ((numSims >= 1) && (numSims <= 1000)) {}
        else{std::cout<< argv[argCounter] << argErrorStrings[argCounter] << std::endl;  return ERROR_BAD_PARAM;} 
    }
    if(argc > ++argCounter)
    {
        unsigned int maxTime = atoi(argv[argCounter]);
        if ((maxTime >= 100) && (maxTime <= 10000))
        {
            maxSimTime = (double)maxTime;
        }
        else{std::cout<< argv[argCounter] << argErrorStrings[argCounter] << std::endl;  return ERROR_BAD_PARAM;} 
    }
    if(argc > ++argCounter)
    {
        unsigned int lambdaIndex = atoi(argv[argCounter]);
        //NB:  index 0-5, for example (not 1-N)
        if ((lambdaIndex >= 0) && (lambdaIndex < NUM_LAMBDAS))
        {
            g_indexValue = lambdaIndex;
            g_lambdaFitness = lambdas[lambdaIndex];

            //compute the compression threshold: uses 0.5 = 500xdRL with dRL=dt=0.001
            //for index 0..3, will compute 1,2,4,8x = 0.5,1.0,2.0,4.0 as limits
            // g_compressionLimit = std::pow(2.0,(double)lambdaIndex) * 0.5;
            // g_compressionLimit = lambdaIndex * 0.5;

            g_compressionLimit = lambdas[lambdaIndex];

            // std::cout << "GROWTH RATE: "<<g_lambdaFitness << std::endl;
            std::cout << "g_compressionLimit: "<<g_compressionLimit << std::endl;
        }
        else{std::cout<< argv[argCounter] << argErrorStrings[argCounter] << std::endl;  return ERROR_BAD_PARAM;} 
    }

    srand (time(NULL));

    double avgSimTime=0.0;  double avgRatio=0.0;
    for (unsigned int i=1;i<=numSims;i++)
    {
        std::cout 
                <<"Running simulation: " << i << std::endl
                <<"trapHeight, trapWidth: " << trapHeight << ", " << trapWidth << std::endl 
                <<"numInitCells, maxSimTime: " << numInitCells << ", " << maxSimTime << std::endl 
        ;

        int cellCount = runSimulation(trapHeight, trapWidth, numInitCells, maxSimTime, 
            strainRatio, elapsedtime, wallTime);

        std::cout << std::fixed << std::setprecision(3)
                <<"Sim (cellCount): " << i << "("<<cellCount<<")    Ratio, time = " << strainRatio << "," << elapsedtime << "," 
                <<" WallTime: " << wallTime/60.0 << std::endl;

        avgRatio += strainRatio; 
        avgSimTime += elapsedtime;                
    }
    std::cout<<"SUMMARY DATA: avgSimTime, avgRatio: "<<avgSimTime/numSims<<", "<<avgRatio/numSims<<std::endl;

}

