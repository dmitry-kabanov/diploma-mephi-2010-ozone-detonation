#include "ThermDataVar.h"
#include <cstdlib>

ThermDataVar::ThermDataVar()
{
	n    = 0;
    Tlow = NULL;
	Tup  = NULL;
	a    = NULL;
}

ThermDataVar::~ThermDataVar()
{
	if (n > 0) {
		for (int j = 0; j < n; j++) {
			delete [] a[j];
		}
		delete [] a;
		delete [] Tlow;
		delete [] Tup;
	}
}

int ThermDataVar::AllocateMemoryForTemperatureRange()
{
		Tlow = new double[n];
		Tup  = new double[n];
		a    = new double *[n];
        for (int j = 0; j < n; j++) {
			a[j] = new double[7];
        }
		return 0;			
}

