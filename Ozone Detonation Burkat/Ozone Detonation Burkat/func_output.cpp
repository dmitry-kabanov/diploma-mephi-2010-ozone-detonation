#include "func_output.h"
#include <iostream>
#include <fstream>
#include "constants.h"
#include "RealType.h"

void outputAsCSVFile(int timeStep,
                     const int *cellNumbers,
                     const RealType *x,
                     const RealType *xCenter,
                     const RealType *p,
                     const RealType *u,
                     const RealType *rho,
                     const RealType *fullEnergy,
                     const RealType *intEnergy,
                     RealType *volumeFractions[],
                     const bool *shockWaveVelocity)
{
    int num_digits;
    char time[32];
    char filename[] = "Output\\data_";
    char file_ext[5] = ".csv";
    char fullname[64];
    
    num_digits = sprintf_s(time, "%d", timeStep);
    strcpy_s(fullname, filename);
    strcat_s(fullname, time);
    strcat_s(fullname, file_ext);

    std::ofstream outFile(fullname);
    outFile.setf(std::ios::fixed, std::ios::floatfield);
    outFile.precision(9);
    
    outFile << "Cell No;Right Bound(m);Coordinate (m);Pressure (Pa);Velocity (m/s);" <<
        "Density (kg m-3);Full energy (J kg-1);Internal Energy (J kg-1);" << 
        "X(O);X(O2);X(O3)" << 
        std::endl;

    outFile << ";" << x[0] << std::endl;
    for (int i = 1; i < N; i++) {
        outFile << cellNumbers[i] << ";" <<
            x[i]                  << ";" <<
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

    outFile.close();
}