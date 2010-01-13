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
#include <iostream>
#include <cmath>
#include "constants.h"
using namespace std;

StifflSolver::StifflSolver(int NYDIM_PAR, double *values, double t_begin, double t_end,
					 double t_step_begin,
					 const char *fileOfSubstances,
                     const char *fileOfReactions,
                     const char *fileOfMoleFractions)
	: Stiffl(NYDIM_PAR, values, t_begin, t_end, t_step_begin)
{
    h = t_step_begin;
    timeStepForOutput = 100;
    mixture = new Mixture(fileOfSubstances,
                          fileOfReactions,
                          fileOfMoleFractions);
    //outputFile.open("output.txt");
    //outputFile.precision(6);
    //outputFile.setf(ios::scientific, ios::floatfield);
}

StifflSolver::~StifflSolver()
{
    delete mixture;
    //outputFile.close();
}

void StifflSolver::performIntegration(RealType afullTime)
{
    time = 0.0;
	int stifflCode;

	t = 0;
    h = 1.0e-13;
	tfin = afullTime;

	prepareStiffl();

    //printHeadingToFile();
    //printToFile();

    // Производим интегрирование.
		
	(Y[0])[0] = mixture->concentrations[0];
	(Y[0])[1] = mixture->concentrations[1];
	(Y[0])[2] = mixture->concentrations[2];
	(Y[0])[3] = mixture->temperature;

	stifflCode = STIFFL();

	mixture->concentrations[0] = (Y[0])[0];
	mixture->concentrations[1] = (Y[0])[1];
	mixture->concentrations[2] = (Y[0])[2];
	mixture->temperature = (Y[0])[3];

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

    // Универсальная газовая постоянная, ккал / (моль*К).
    const RealType r = 0.001985846;

    k = exp(log(t) * mixture->reactions[i].n - 
        mixture->reactions[i].activationEnergy / (r * t) + 
        2.30258 * mixture->reactions[i].log10A);

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
        mixture->concentrations[1] / sumConc * 100 << "\t" <<
        mixture->concentrations[2] / sumConc * 100 << "\t" <<
        mixture->concentrations[0] / sumConc * 100 << "\t" <<
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
    outputFile << "t (s)\tO2\tO3\tO\tN (1/cm3)\tH (J/g)\tU (J/g)\tT (K)\tP (atm)\t";
    outputFile << "Mu (g/mole)\t\tRho (g/cm3)\tV(m3/kg)" << endl;
}

void StifflSolver::updateMoleFractions(RealType *vf)
{
    RealType sumConc = 0.0;

    for (int i = 0; i < mixture->nSubstances; i++) {
        sumConc += mixture->concentrations[i];
    }

    for (int i = 0; i < mixture->nSubstances; i++) {
        vf[i] = mixture->concentrations[i] / sumConc * 100;
    }
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
	//mixture->temperature = (Y[0])[3];
	//printToFile();

	return 0;
}

int StifflSolver::DIFFUN(double **YY, double *F)
{
	mixture->concentrations[0] = Y[0][0];
	mixture->concentrations[1] = Y[0][1];
	mixture->concentrations[2] = Y[0][2];
	mixture->temperature = Y[0][3];

	concOfM = 0;
	for (int i = 0; i < mixture->nSubstances; i++) {
		concOfM += mixture->concentrations[i];
	}
	k1f = calculateRateForForwardReaction(0);
	k2f = calculateRateForForwardReaction(1);
	k3f = calculateRateForForwardReaction(2);

	k1r = calculateRateForBackReaction(0, k1f);
	k2r = calculateRateForBackReaction(1, k2f);
	k3r = calculateRateForBackReaction(2, k3f);

    F[0] = rightSideForO(
        mixture->concentrations[0], 
        mixture->concentrations[2],
        mixture->concentrations[1]);
    F[1] = rightSideForO2(
        mixture->concentrations[0],
        mixture->concentrations[2],
        mixture->concentrations[1]);
    F[2] = rightSideForO3(
        mixture->concentrations[0],
        mixture->concentrations[2],
        mixture->concentrations[1]);

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
            mixture->concentrations[i];
	}

	F[3] = u / v;

	return 0;
}