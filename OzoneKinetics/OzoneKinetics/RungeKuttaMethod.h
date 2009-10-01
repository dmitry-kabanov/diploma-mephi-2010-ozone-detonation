#ifndef RUNGE_KUTTA_METHOD_H
#define RUNGE_KUTTA_METHOD_H

#include <iostream>
#include <fstream>
#include "RealType.h"
#include "Mixture.h"
#include "Reaction.h"
#include "Substance.h"

class RungeKuttaMethod
{
public:
    RungeKuttaMethod(RealType p0, RealType t0,
                     const char *fileOfSubstances,
                     const char *fileOfReactions,
                     const char *fileOfVolumeFractions);
    ~RungeKuttaMethod();
private:
    void performIntegration();

    RealType rightSideForO3(RealType t, 
                            RealType concOfO, 
                            RealType concOfO3,
                            RealType concOfO2);
    RealType rightSideForO(RealType t, 
                           RealType concOfO, 
                           RealType concOfO3, 
                           RealType concOfO2);
    RealType rightSideForO2(RealType t, 
                           RealType concOfO, 
                           RealType concOfO3, 
                           RealType concOfO2);
    RealType calculateRateForForwardReaction(int i);
    RealType calculateRateForBackReaction(int i, RealType kf);
    RealType h;
    Mixture *mixture;
    int fullTime;
    int timeStepForOutput;
    RealType time;
    std::ofstream outputFile;
    void printToFile();
    void printHeadingToFile();
};

#endif