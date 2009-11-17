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
    delete [] fractions;
}

RealType Piston::calculateVelocity(RealType f)
{
    RealType p;
    //f = 100.1 - f;
    
    if (f > fractions[0]) {
        p = pressures[0];
    }
    else {
        for (int i = 1; i < nColumns; ++i) {
            if (f > fractions[i] && f < fractions[i-1]) {
                p = pressures[i];
                break;
            }
        }
    }

    return sqrt((p - P2) * (1 / RHO2 - 1 / RHO1));
}

void Piston::readFileOfPiston(const char* filename)
{
    ifstream iFile(filename);
    
    iFile >> nColumns;

    pressures = new RealType[nColumns];
    fractions = new RealType[nColumns];

    for (int i = 0; i < nColumns; ++i) {
        iFile >> pressures[i];
    }

    for (int i = 0; i < nColumns; ++i) {
        iFile >> fractions[i];
    }
}
