/**
 * @file
 *
 * @author  Кабанов Дмитрий <kabanovdmitry@gmail.com
 * @version %I%
 *
 * @section DESCRIPTION
 *
 * Реализация класса Reaction.
 */
#include "Reaction.h"

Reaction::Reaction()
{
	nEff = false;
	pEff = 0;
}

Reaction::~Reaction()
{
    delete [] reagents;
    delete [] products;

	delete [] temperatureLow;
	delete [] temperatureHigh;
	delete [] log10A;
	delete [] n;
	delete [] activationEnergy;
}

void Reaction::allocateMemoryForParameters()
{
	temperatureLow   = new RealType[nTemperatureRanges];
	temperatureHigh  = new RealType[nTemperatureRanges];
	log10A			 = new RealType[nTemperatureRanges];
	n				 = new RealType[nTemperatureRanges];
	activationEnergy = new RealType[nTemperatureRanges];
}

void Reaction::allocateMemoryForCollisionEfficiency(int n)
{
	pEff = new RealType[n];
	for (int i = 0; i < n; ++i) {
		pEff[i] = 1.0;
	}
}