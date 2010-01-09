#include <cmath>
#include <fstream>
#include "Piston.h"
#include "constants.h"
using namespace std;

Piston::Piston(RealType pInitial, 
			   RealType rhoInitial, 
			   const char *fileOfPistonData)
:	pInitial_(pInitial),
	rhoInitial_(rhoInitial)
{
    readFileOfPiston(fileOfPistonData);
}

Piston::~Piston()
{
    delete [] pressures_;
    delete [] densities_;
    delete [] fractions_;
}

RealType Piston::calculateVelocity(RealType f)
{
    RealType p;
    RealType rho;

    for (int i = 0; i < nRows_; ++i) {
        if (f <= fractions_[i]) {
            p = pressures_[i];
            rho = densities_[i];
            break;
        }
    }

    return sqrt((p - P2) * (1 / RHO2 - 1 / rho));
}

void Piston::readFileOfPiston(const char* fileOfPistonData)
{
    ifstream iFile(fileOfPistonData);
    
    iFile >> nRows_;

    pressures_ = new RealType[nRows_];
    densities_ = new RealType[nRows_];
    fractions_ = new RealType[nRows_];

    for (int i = 0; i < nRows_; ++i) {
        iFile >> pressures_[i] >> densities_[i] >> fractions_[i];
    }

    iFile.close();
}
