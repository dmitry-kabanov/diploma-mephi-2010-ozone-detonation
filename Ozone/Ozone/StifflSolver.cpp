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
	RealType reactionRate;
	// Порядковый номер вещества. Вспомогательная переменная для читабельности.
	int speciesNumber;

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
		reactionRate = mixture->reactions[i].calculateReactionRate(
			*YY,
			mixture->temperature,
			mixture->gibbsEnergy,
			mixture->getNSpecies());

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