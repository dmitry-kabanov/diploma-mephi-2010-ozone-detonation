#include "func_output.h"
#include <iostream>
#include <fstream>
#include "constants.h"
#include "RealType.h"

void outputAsCSVFile(std::ofstream &outFile,
                     const int *cellNumbers,
                     const RealType *xCenter,
                     const RealType *p,
                     const RealType *u,
                     const RealType *rho,
                     const RealType *fullEnergy,
                     const RealType *intEnergy,
                     RealType *volumeFractions[],
                     const bool *shockWaveVelocity)
{
    outFile << "Cells Numbers;xCenter;Pressure (Pa);Velocity (m/s);" <<
        "Density (kg m-3);Full energy (J kg-1);Internal Energy (J kg-1);" << 
        "Mole Fraction of O;Mole Fraction of O2;Mole Fraction of O3" << 
        std::endl;

    for (int i = 1; i < N; i++) {
        outFile << cellNumbers[i] << ";" <<
            xCenter[i]            << ";" << 
            p[i]                  << ";" << 
            u[i]                  << ";" << 
            rho[i]                << ";" <<
            fullEnergy[i]         << ";" << 
            intEnergy[i]          << ";" <<
            volumeFractions[i][0] << ";" <<
            volumeFractions[i][1] << ";" <<
            volumeFractions[i][2] << 
            std::endl;
        if (shockWaveVelocity[i] == true && shockWaveVelocity[i+1] == false) {
            break;
        }
    }
}