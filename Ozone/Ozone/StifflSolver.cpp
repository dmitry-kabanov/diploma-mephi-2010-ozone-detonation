/**
* @file
*
* @author  Кабанов Дмитрий <kabanovdmitry@gmail.com>
* @version %I%
*
* @section DESCRIPTION
*
* Реализация класса RungeKuttaMethod.
*/
#include "StifflSolver.h"
#include <iomanip>
#include <iostream>
#include <cmath>
#include "constants.h"
using namespace std;

StifflSolver::StifflSolver(int NYDIM_PAR,
						   double *values,
						   double t_begin,
						   double t_end,
						   double t_step_begin)
	: Stiffl(NYDIM_PAR, values, t_begin, t_end, t_step_begin)
{
    h = t_step_begin;
    timeStepForOutput = 100;
    mixture = 0;
    //outputFile.open("output.txt");
    //outputFile.precision(6);
    //outputFile.setf(ios::scientific, ios::floatfield);
}

StifflSolver::~StifflSolver()
{
    //outputFile.close();
}

void StifflSolver::performIntegration(Mixture &mix, RealType afullTime)
{
    time = 0.0;
	int stifflCode;

	t = 0;
    h = 1.0e-13;
	tfin = afullTime;

	mixture = &mix;

	prepareStiffl();

    //printHeadingToFile();
    //printToFile();

	for (int i = 0; i < mixture->nSubstances; ++i) {
		(Y[0])[i] = mixture->concentrations[i];
	}
	(Y[0])[mixture->nSubstances] = mixture->temperature;

	// Производим интегрирование.
	stifflCode = STIFFL();

	for (int i = 0; i < mixture->nSubstances; ++i) {
		mixture->concentrations[i] = (Y[0])[i];
	}
	mixture->temperature = (Y[0])[mixture->nSubstances];

	mixture->molecularWeight = mixture->calculateMolecularWeight();
	mixture->pressure = mixture->calculatePressure();

}

RealType StifflSolver::rightSideForO(RealType concOfO, 
                                         RealType concOfO3, 
                                         RealType concOfO2)
{
    RealType m1;
    RealType m2;
    RealType m3;

    m1 = -k1f * concOfO * concOfO3 + k1r * concOfO2 * concOfO2;
    m2 = k2f * concOfM * concOfO3 - k2r * concOfO2 * concOfO * concOfM;
    m3 = 2 * k3f * concOfM * concOfO2 - 
        2 * k3r * concOfO * concOfO * concOfM;

    return m1 + m2 + m3;
}

RealType StifflSolver::rightSideForO3(RealType concOfO, 
                                          RealType concOfO3, 
                                          RealType concOfO2)
{
    RealType m1;
    RealType m2;

    m1 = -k1f * concOfO * concOfO3 + k1r * concOfO2 * concOfO2;
    m2 = -k2f * concOfM * concOfO3 + k2r * concOfO2 * concOfO * concOfM;

    return m1 + m2;
}

RealType StifflSolver::rightSideForO2(RealType concOfO, 
                                         RealType concOfO3, 
                                         RealType concOfO2)
{
    RealType m1;
    RealType m2;
    RealType m3;

    m1 = 2 * k1f * concOfO * concOfO3 - 2 * k1r * concOfO2 * concOfO2;
    m2 = k2f * concOfM * concOfO3 - k2r * concOfO2 * concOfO * concOfM;
    m3 = - k3f * concOfM * concOfO2 + k3r * concOfO * concOfO * concOfM;

    return m1 + m2 + m3;
}

RealType StifflSolver::calculateRateForForwardReaction(int i)
{
    RealType k;
    RealType t = mixture->temperature;
	int j = 0;

	//if (mixture->reactions[i].nTemperatureRanges > 1) {
	//	for (int kk = mixture->reactions[i].nTemperatureRanges; kk > 0; kk--) {
	//		if (t >= mixture->reactions[i].temperatureLow[kk-1]) {
	//			j = kk;
	//			break;
	//		}
	//	}
	//}

    // Универсальная газовая постоянная, ккал / (моль*К).
    const RealType r = 0.001985846;

    k = exp(log(t) * mixture->reactions[i].n[j] - 
        mixture->reactions[i].activationEnergy[j] / (r * t) + 
        2.30258 * mixture->reactions[i].log10A[j]);

    return k;
}

RealType StifflSolver::calculateRateForBackReaction(int i, RealType kf)
{
    RealType t = mixture->temperature;
    // Тепловой эффект реакции.
    RealType q = 0.0;
    int substanceNumber;
    int nMoles = 0;

    for (int j = 0; j < mixture->reactions[i].nProducts; j++) {
        substanceNumber = mixture->reactions[i].products[j];
        // Проверяем, что вещество не является веществом "М".
        if (substanceNumber == -1) {
            continue;
        }
        q += mixture->calculateSubstanceGibbsEnergy(substanceNumber, t);
        nMoles++;
    }

    for (int j = 0; j < mixture->reactions[i].nReagents; j++) {
        substanceNumber = mixture->reactions[i].reagents[j];
        // Проверяем, что вещество не является веществом "М".
        if (substanceNumber == -1) {
            continue;
        }
        q -= mixture->calculateSubstanceGibbsEnergy(substanceNumber, t);
        nMoles--;
    }

    // Константа равновесия.
    // Значение константы скорости обратной реакции получается с помощью 
    // константы скорости прямой реакции и константы равновесия.
    RealType kc;
    kc  = exp(-q / (mixture->R_J_OVER_KMOL_K * t));
    kc *= exp(-log(10 * mixture->K_BOLTZMANN * t) * nMoles);

    return kf / kc;
}

void StifflSolver::printToFile()
{
    RealType sumConc = 0.0;
    for (int i = 0; i < mixture->nSubstances; i++) {
        sumConc += mixture->concentrations[i];
    }
    outputFile << t << "\t" <<
        mixture->concentrations[2] / sumConc * 100 << "\t" <<
        mixture->concentrations[6] / sumConc * 100 << "\t" <<
        mixture->concentrations[1] / sumConc * 100 << "\t" <<
        sumConc << "\t" <<
        (mixture->fullEnergy + mixture->pressure * mixture->volume) * 1.0e-3 << "\t" <<
        mixture->fullEnergy * 1.0e-3 << "\t" <<
        mixture->temperature << "\t" <<
        mixture->pressure / ONE_ATM << "\t" <<
        mixture->molecularWeight << "\t\t" <<
        1.0 / mixture->volume * 1.0e-3 <<  "\t" <<
        mixture->volume << endl;
}

void StifflSolver::printHeadingToFile()
{
    outputFile << "t (s)\tH2O\tO2\tH2\tN (1/cm3)\tH (J/g)\tU (J/g)\tT (K)\tP (atm)\t";
    outputFile << "Mu (g/mole)\t\tRho (g/cm3)\tV(m3/kg)" << endl;
}

RealType StifflSolver::getPressure()
{
    return mixture->calculatePressure();
}

Mixture *StifflSolver::getMixture()
{
    return mixture;
}

void StifflSolver::PEDERV()
{
}

int StifflSolver::IFNSH()
{
	//mixture->molecularWeight = mixture->calculateMolecularWeight();
	//mixture->pressure = mixture->calculatePressure();
	//mixture->temperature = Y[0][mixture->nSubstances];

 //   for (int i = 0; i < mixture->nSubstances; ++i) {
 //       mixture->concentrations[i] = Y[0][i];
 //   }
	//printToFile();

	return 0;
}

int StifflSolver::DIFFUN(double **YY, double *F)
{
	RealType kf;
	RealType kr;
	RealType multiplicationOfReagents;
	RealType multiplicationOfProducts;
	RealType reactionRate;
	RealType sumConc;
	bool withThirdBody = false;
	// Порядковый номер вещества. Вспомогательная переменная для читабельности.
	int speciesNumber;
	// Тепловой эффект реакции.
	RealType q;
	int nMoles;
	RealType kc;

	// TODO: ввести в StifflSolver переменную-член 
	// с количеством уравнений nEquations.
	mixture->temperature = (*YY)[mixture->getNSpecies()];

	for (int i = 0; i < mixture->getNSpecies(); ++i) {
		F[i] = 0.0;
		mixture->gibbsEnergy[i] = mixture->calculateSubstanceGibbsEnergy(i, mixture->temperature);
	}
	F[mixture->getNSpecies()]   = 0.0;
	F[mixture->getNSpecies()+1] = 0.0;
	

	for (int i = 0; i < mixture->nReactions; ++i) {
		kf = calculateRateForForwardReaction(i);

		multiplicationOfReagents = 1.0;
		multiplicationOfProducts = 1.0;
		nMoles = 0;
		q = 0;

		for (int j = 0; j < mixture->reactions[i].nReagents; ++j) {
			speciesNumber = mixture->reactions[i].reagents[j];
			if (speciesNumber == -1) {
				withThirdBody = true;
				continue;
			}
			multiplicationOfReagents *= (*YY)[speciesNumber];
			q -= mixture->gibbsEnergy[speciesNumber];
			nMoles--;
		}
		for (int j = 0; j < mixture->reactions[i].nProducts; ++j) {
			speciesNumber = mixture->reactions[i].products[j];
			if (speciesNumber == -1) {
				withThirdBody = true;
				continue;
			}
			multiplicationOfProducts *= (*YY)[speciesNumber];
			q += mixture->gibbsEnergy[speciesNumber];
			nMoles++;
		}
		
		// Вычисляем константу равновесия по концентрации.
		kc  = exp(-q / (mixture->R_J_OVER_KMOL_K * mixture->temperature));
		kc *= exp(-log(10 * mixture->K_BOLTZMANN * mixture->temperature) * nMoles);

		if (mixture->reactions[i].direction == 0) {
			// Реакция обратимая.
			kr = kf / kc;
			reactionRate = kf * multiplicationOfReagents - kr * multiplicationOfProducts;
		}
		else if (mixture->reactions[i].direction == 1) {
			// Реакция необратимая.
			reactionRate = kf * multiplicationOfReagents;
		}
		else {
			cout << "Unknown direction of reaction '"
				 << mixture->reactions[i].nameOfReaction << "'." << endl;
			exit(-1);
		}

		sumConc = 0.0;
		if (withThirdBody) {
			if (mixture->reactions[i].nEff) {
				//cout << "Reaction " << i << endl;
				for (int j = 0; j < mixture->getNSpecies(); ++j) {
					sumConc += mixture->reactions[i].pEff[j] * (*YY)[j];
				}
			}
			else {
				for (int j = 0; j < mixture->getNSpecies(); ++j) {
					sumConc += (*YY)[j];
				}
			}
		}
		if (withThirdBody) {
			reactionRate *= sumConc;
		}
		withThirdBody = false;

		for (int j = 0; j < mixture->reactions[i].nReagents; ++j) {
			speciesNumber = mixture->reactions[i].reagents[j];
			if (speciesNumber == -1) {
				continue;
			}
			F[speciesNumber] -= reactionRate;
		}

		for (int j = 0; j < mixture->reactions[i].nProducts; ++j) {
			speciesNumber = mixture->reactions[i].products[j];
			if (speciesNumber == -1) {
				continue;
			}
			F[speciesNumber] += reactionRate;
		}
	}

	// Cчитаем правую часть уравнения для температуры.
	RealType u = 0;
    RealType v = 0;

	for (int i = 0; i < mixture->nSubstances; i++) {
		u += F[i] * (
            mixture->calculateSubstanceEnthalpy(i, mixture->temperature) * 
            mixture->substances[i]->molecularWeight - 
			mixture->R_J_OVER_KMOL_K * mixture->temperature);
		v += (mixture->R_J_OVER_KMOL_K - 
            mixture->substances[i]->molecularWeight *
            mixture->calculateSubstanceCp(i, mixture->temperature)) * 
            (*YY)[i];
	}

	F[mixture->nSubstances] = u / v;

	return 0;
}

