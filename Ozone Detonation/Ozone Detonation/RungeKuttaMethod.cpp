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
#include <iostream>
#include <cmath>
#include "RungeKuttaMethod.h"
using namespace std;

RungeKuttaMethod::RungeKuttaMethod(const char* fileOfSubstances,
                                   const char* fileOfReactions)
{
    h = 1.0e-12;
    timeStepForOutput = 100;
    mixture = new Mixture(fileOfSubstances,
                          fileOfReactions);
    outputFile.open("output.txt");
    outputFile.precision(6);
    outputFile.setf(ios::scientific, ios::floatfield);
}

RungeKuttaMethod::~RungeKuttaMethod()
{
    delete mixture;
}

void RungeKuttaMethod::performIntegration(RealType afullTime)
{
    // Относительное изменение концентрации O за один временной шаг.
    RealType alpha = 0.01;
    RealType k1, k2, k3, k4;
    RealType q1, q2, q3, q4;
    RealType r1, r2, r3, r4;
    time = 0.0;

    printHeadingToFile();
    printToFile();

    int i;
    // Производим интегрирование.
    for (i = 0; time <= afullTime; i++) {
        k1 = rightSideForO(time,
            mixture->concentrations[0], 
            mixture->concentrations[2],
            mixture->concentrations[1]);
        q1 = rightSideForO3(time,
            mixture->concentrations[0],
            mixture->concentrations[2],
            mixture->concentrations[1]);
        r1 = rightSideForO2(time,
            mixture->concentrations[0],
            mixture->concentrations[2],
            mixture->concentrations[1]);

        k2 = rightSideForO(time + h / 2.0, 
            mixture->concentrations[0] + h * k1 / 2.0, 
            mixture->concentrations[2] + h * q1 / 2.0,
            mixture->concentrations[1] + h * r1 / 2.0);
        q2 = rightSideForO3(time + h / 2.0, 
            mixture->concentrations[0] + h * k1 / 2.0, 
            mixture->concentrations[2] + h * q1 / 2.0,
            mixture->concentrations[1] + h * r1 / 2.0);
        r2 = rightSideForO2(time + h / 2.0, 
            mixture->concentrations[0] + h * k1 / 2.0, 
            mixture->concentrations[2] + h * q1 / 2.0,
            mixture->concentrations[1] + h * r1 / 2.0);

        k3 = rightSideForO(time + h / 2.0, 
            mixture->concentrations[0] + h * k2 / 2.0, 
            mixture->concentrations[2] + h * q2 / 2.0,
            mixture->concentrations[1] + h * r2 / 2.0);
        q3 = rightSideForO3(time + h / 2.0, 
            mixture->concentrations[0] + h * k2 / 2.0, 
            mixture->concentrations[2] + h * q2 / 2.0,
            mixture->concentrations[1] + h * r2 / 2.0);
        r3 = rightSideForO2(time + h / 2.0, 
            mixture->concentrations[0] + h * k2 / 2.0, 
            mixture->concentrations[2] + h * q2 / 2.0,
            mixture->concentrations[1] + h * r2 / 2.0);

        k4 = rightSideForO(time + h, 
            mixture->concentrations[0] + h * k3, 
            mixture->concentrations[2] + h * q3, 
            mixture->concentrations[1] + h * r3);
        q4 = rightSideForO3(time + h, 
            mixture->concentrations[0] + h * k3, 
            mixture->concentrations[2] + h * q3, 
            mixture->concentrations[1] + h * r3);
        r4 = rightSideForO2(time + h, 
            mixture->concentrations[0] + h * k3, 
            mixture->concentrations[2] + h * q3, 
            mixture->concentrations[1] + h * r3);

        mixture->previousConcentrations[0] = mixture->concentrations[0];
        mixture->previousConcentrations[1] = mixture->concentrations[1];
        mixture->previousConcentrations[2] = mixture->concentrations[2];
        
        mixture->concentrations[0] += h * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;
        mixture->concentrations[2] += h * (q1 + 2 * q2 + 2 * q3 + q4) / 6.0;
        mixture->concentrations[1] += h * (r1 + 2 * r2 + 2 * r3 + r4) / 6.0;
        //mixture->concentrations[1] = (mixture->initialConcentrations[0] -
        //    mixture->concentrations[0] + 3 * mixture->initialConcentrations[2] -
        //    3 * mixture->concentrations[2] + 
        //    2 * mixture->initialConcentrations[1]) / 
        //    2.0;

        mixture->molecularWeight = mixture->calculateMolecularWeight();
        mixture->temperature     = mixture->calculateTemperature();

        time += h;

        //h = alpha * mixture->concentrations[0] / rightSideForO(time,
        //    mixture->concentrations[0],
        //    mixture->concentrations[2],
        //    mixture->concentrations[1]);

        //if ((time + h) > afullTime) {
        //    h = afullTime - time;
        //}

        if (i % timeStepForOutput == 0) {
            printToFile();
        }

        if (abs(
            mixture->concentrations[1] / mixture->previousConcentrations[1] - 
            1.0) <= 1.0e-6) {
                break;
        }
    }

    nIterations = i+1;
    cout << "Iterations: " << nIterations << endl;
}

RealType RungeKuttaMethod::rightSideForO(RealType t, 
                                         RealType concOfO, 
                                         RealType concOfO3, 
                                         RealType concOfO2)
{
    RealType concOfM = concOfO + concOfO2 + concOfO3;

    RealType m1;
    RealType m2;
    RealType m3;

    RealType k1f = calculateRateForForwardReaction(0);
    RealType k2f = calculateRateForForwardReaction(1);
    RealType k3f = calculateRateForForwardReaction(2);

    RealType k1r = calculateRateForBackReaction(0, k1f);
    RealType k2r = calculateRateForBackReaction(1, k2f);
    RealType k3r = calculateRateForBackReaction(2, k3f);

    RealType skor1 = k1f * concOfO * concOfO3;
    RealType skor2m = k2f * concOfM * concOfO3;
    RealType skor3m = k3f * concOfM;
    m1 = -k1f * concOfO * concOfO3 + k1r * concOfO2 * concOfO2;
    m2 = k2f * concOfM * concOfO3 - k2r * concOfO2 * concOfO * concOfM;
    m3 = 2 * k3f * concOfM * concOfO2 - 
        2 * k3r * concOfO * concOfO * concOfM;

    return m1 + m2 + m3;
}

RealType RungeKuttaMethod::rightSideForO3(RealType t, 
                                          RealType concOfO, 
                                          RealType concOfO3, 
                                          RealType concOfO2)
{
    RealType concOfM = concOfO + concOfO2 + concOfO3;

    RealType m1;
    RealType m2;

    RealType k1f = calculateRateForForwardReaction(0);
    RealType k2f = calculateRateForForwardReaction(1);

    RealType k1r = calculateRateForBackReaction(0, k1f);
    RealType k2r = calculateRateForBackReaction(1, k2f);

    m1 = -k1f * concOfO * concOfO3 + k1r * concOfO2 * concOfO2;
    m2 = -k2f * concOfM * concOfO3 + k2r * concOfO2 * concOfO * concOfM;

    return m1 + m2;
}

RealType RungeKuttaMethod::rightSideForO2(RealType t, 
                                         RealType concOfO, 
                                         RealType concOfO3, 
                                         RealType concOfO2)
{
    RealType concOfM = concOfO + concOfO2 + concOfO3;

    RealType m1;
    RealType m2;
    RealType m3;

    RealType k1f = calculateRateForForwardReaction(0);
    RealType k2f = calculateRateForForwardReaction(1);
    RealType k3f = calculateRateForForwardReaction(2);

    RealType k1r = calculateRateForBackReaction(0, k1f);
    RealType k2r = calculateRateForBackReaction(1, k2f);
    RealType k3r = calculateRateForBackReaction(2, k3f);

    m1 = 2 * k1f * concOfO * concOfO3 - 2 * k1r * concOfO2 * concOfO2;
    m2 = k2f * concOfM * concOfO3 - k2r * concOfO2 * concOfO * concOfM;
    m3 = - k3f * concOfM * concOfO2 + k3r * concOfO * concOfO * concOfM;

    return m1 + m2 + m3;
}

RealType RungeKuttaMethod::calculateRateForForwardReaction(int i)
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

RealType RungeKuttaMethod::calculateRateForBackReaction(int i, RealType kf)
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
        q += mixture->calculateGibbsEnergy(substanceNumber);
        nMoles++;
    }

    for (int j = 0; j < mixture->reactions[i].nReagents; j++) {
        substanceNumber = mixture->reactions[i].reagents[j];
        // Проверяем, что вещество не является веществом "М".
        if (substanceNumber == -1) {
            continue;
        }
        q -= mixture->calculateGibbsEnergy(substanceNumber);
        nMoles--;
    }

    // Константа равновесия.
    // Значение константы скорости обратной реакции получается с помощью 
    // константы скорости прямой реакции и константы равновесия.
    RealType kp;
    kp  = exp(-q / (mixture->R_J_OVER_MOL_K * t));
    kp *= exp(-log(10 * mixture->K_BOLTZMANN * t) * nMoles);

    return kf / kp;
}

void RungeKuttaMethod::printToFile()
{
    RealType sumConc = 0.0;
    for (int i = 0; i < mixture->nSubstances; i++) {
        sumConc += mixture->concentrations[i];
    }
    outputFile << time << "\t" <<
        mixture->concentrations[1] / sumConc * 100 << "\t" <<
        mixture->concentrations[2] / sumConc * 100 << "\t" <<
        mixture->concentrations[0] / sumConc * 100 << "\t" <<
        sumConc << "\t" <<
        mixture->fullEnergy * 1.0e-3 << "\t\t" <<
        mixture->temperature << "\t\t" <<
        mixture->molecularWeight << "\t\t" <<
        1.0 / mixture->volume * 1.0e-3 <<  "\t" <<
        mixture->volume << endl;
}

void RungeKuttaMethod::printHeadingToFile()
{
    outputFile << "t (s)\tO2\tO3\tO\tN (1/cm3)\tH (J/g)\t\tT (K)\t\t";
    outputFile << "Mu (g/mole)\t\tRho (g/cm3)\tV(m3/kg)" << endl;
}

RealType* RungeKuttaMethod::getVolumeFractions()
{
    RealType *volFracts = new RealType[mixture->nSubstances];
    RealType sumConc = 0.0;

    for (int i = 0; i < mixture->nSubstances; i++) {
        sumConc += mixture->concentrations[i];
    }

    for (int i = 0; i < mixture->nSubstances; i++) {
        volFracts[i] = mixture->concentrations[i] / sumConc * 100;
    }

    return volFracts;
}

RealType RungeKuttaMethod::getPressure()
{
    return mixture->calculatePressure();
}

Mixture *RungeKuttaMethod::getMixture()
{
    return mixture;
}