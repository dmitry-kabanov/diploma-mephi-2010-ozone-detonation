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
	withThirdBody = false;
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

RealType Reaction::calculateReactionRate(double *Y, 
										 RealType t, 
										 RealType *gibbsEnergies, 
										 int nSpecies)
{
	kf = calculateConstantRate(t);
	multiplicationOfReagents = 1.0;
	multiplicationOfProducts = 1.0;
	nMoles = 0;
	q = 0;
	RealType sumConc;

	for (int j = 0; j < nReagents; ++j) {
		if (reagents[j] == -1) {
			continue;
		}
		multiplicationOfReagents *= Y[reagents[j]];
		q -= gibbsEnergies[reagents[j]];
		nMoles--;
	}
	for (int j = 0; j < nProducts; ++j) {
		if (products[j] == -1) {
			continue;
		}
		multiplicationOfProducts *= Y[products[j]];
		q += gibbsEnergies[products[j]];
		nMoles++;
	}

	// Вычисляем константу равновесия по концентрации.
	kc  = exp(-q / (R_J_OVER_KMOL_K * t));
	kc *= exp(-log(10 * K_BOLTZMANN * t) * nMoles);

	if (direction == 0) {
		// Реакция обратимая.
		//kr = kf / kc;
		reactionRate_ = kf * 
			(multiplicationOfReagents - multiplicationOfProducts / kc);
	}
	else if (direction == 1) {
		// Реакция необратимая.
		reactionRate_ = kf * multiplicationOfReagents;
	}
	else {
		cout << "Unknown direction of reaction '"
			<< nameOfReaction << "'." << endl;
		exit(-3);
	}

	sumConc = 0.0;
	if (withThirdBody) {
		if (nEff) {
			for (int j = 0; j < nSpecies; ++j) {
				sumConc += pEff[j] * Y[j];
			}
		}
		else {
			for (int j = 0; j < nSpecies; ++j) {
				sumConc += Y[j];
			}
		}
		reactionRate_ *= sumConc;
	}

	return reactionRate_;
}