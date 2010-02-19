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
#include <cmath>
#include "constants.h"
using namespace std;

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

RealType Reaction::calculateConstantRate(RealType t)
{
	RealType k;
	// TODO: сделать, чтоб возможно было выбрать несколько температурных диапазонов.
	int j = 0;

	k = exp(log(t) * n[j] - 
		activationEnergy[j] / (GAS_CONSTANT_KCAL_OVER_MOL_K * t) + 
		2.30258 * log10A[j]
	);

	return k;
}