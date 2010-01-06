#include "func_output.h"
#include <iostream>
#include <fstream>
#include "constants.h"
#include "RealType.h"
using namespace std;

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
	// TODO: добавить проверку на существование каталога.
    char filename[] = "Output\\data_";
    char file_ext[5] = ".csv";
    char fullname[64];
    
    num_digits = sprintf_s(time, "%d", timeStep);
    strcpy_s(fullname, filename);
    strcat_s(fullname, time);
    strcat_s(fullname, file_ext);
    ofstream outFile(fullname);

	if (!outFile) {
		cout << "Can't open file '" << fullname << "' for writing." << endl;
	}
    
    outFile.setf(ios::fixed, ios::floatfield);
    outFile.precision(9);
    outFile << "Cell No;Right Bound(m);Coordinate (m);Pressure (Pa);Velocity (m/s);" <<
        "Density (kg m-3);Full energy (J kg-1);Internal Energy (J kg-1);" << 
        "X(O);X(O2);X(O3)" << 
        endl;
    outFile << ";" << x[0] << endl;
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
            endl;
        if (shockWaveVelocity[i] == true && shockWaveVelocity[i+1] == false) {
            break;
        }
    }
    outFile.flush();
    return;
}