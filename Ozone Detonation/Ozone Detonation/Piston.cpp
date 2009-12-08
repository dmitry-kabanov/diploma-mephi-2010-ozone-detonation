#include <cmath>
#include <fstream>
#include "Piston.h"
#include "constants.h"
using namespace std;

Piston::Piston()
{
    readFileOfPiston("Piston.txt");
}

Piston::~Piston()
{
    delete [] pressures;
    delete [] densities;
    delete [] fractions;
}

RealType Piston::calculateVelocity(RealType f)
{
    RealType p;
    RealType rho;

    for (int i = 0; i < nRows; ++i) {
        if (f <= fractions[i]) {
            p = pressures[i];
            rho = densities[i];
            break;
        }
    }

    return sqrt((p - P2) * (1 / RHO2 - 1 / rho));
}

void Piston::readFileOfPiston(const char* filename)
{
    ifstream iFile(filename);
    
    iFile >> nRows;

    pressures = new RealType[nRows];
    densities = new RealType[nRows];
    fractions = new RealType[nRows];

    for (int i = 0; i < nRows; ++i) {
        iFile >> pressures[i] >> densities[i] >> fractions[i];
    }

    iFile.close();
}
