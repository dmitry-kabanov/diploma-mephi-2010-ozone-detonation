#include <iostream>
#include <fstream>
#include "RungeKuttaMethod.h"
using namespace std;

int main()
{
    RealType initialPressure;
    RealType initialTemperature;
    RealType initialTimeStep;
    int fullTime;
    int timeStepForOutput;
    char fileOfSubstances[64];
    char fileOfReactions[64];
    char fileOfFractions[64];
    char paramName[64];
    ifstream configFile("Config.txt");

    configFile >> paramName >> initialPressure;
    initialPressure *= 101325.0; // Переводим из физ. атмосфер в паскали.
    configFile >> paramName >> initialTemperature;
    configFile >> paramName >> initialTimeStep;
    configFile >> paramName >> fullTime;
    configFile >> paramName >> timeStepForOutput;
    configFile >> paramName >> fileOfSubstances;
    configFile >> paramName >> fileOfReactions;
    configFile >> paramName >> fileOfFractions;

    RungeKuttaMethod rungeKuttaMethod(initialPressure,
                                      initialTemperature,
                                      initialTimeStep,
                                      fullTime,
                                      timeStepForOutput,
                                      fileOfSubstances,
                                      fileOfReactions,
                                      fileOfFractions
    );
}