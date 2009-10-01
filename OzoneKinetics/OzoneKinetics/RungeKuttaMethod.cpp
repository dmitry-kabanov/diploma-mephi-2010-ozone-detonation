#include "RungeKuttaMethod.h"
#include <iostream>
#include <cmath>
using namespace std;

RungeKuttaMethod::RungeKuttaMethod(RealType p0, RealType t0, 
                                   const char *fileOfSubstances,
                                   const char *fileOfReactions,
                                   const char *fileOfVolumeFractions)
{
    h = 1.0e-12;
    time = 0.0;
    fullTime = 150000;
    timeStepForOutput = 500;
    mixture = new Mixture(p0, t0, 
                          fileOfSubstances,
                          fileOfReactions,
                          fileOfVolumeFractions
    );
    outputFile.open("output.txt");
    outputFile.precision(6);
    outputFile.setf(ios::scientific, ios::floatfield);

    performIntegration();
}

RungeKuttaMethod::~RungeKuttaMethod()
{
    delete mixture;
}

void RungeKuttaMethod::performIntegration()
{
    // ������������� ��������� ������������ O �� ���� ��������� ���.
    RealType alpha = 0.01;
    RealType k1, k2, k3, k4;
    RealType q1, q2, q3, q4;
    RealType r1, r2, r3, r4;

    mixture->calculateTemperature();
    printHeadingToFile();
    printToFile();

    // ���������� ��������������.
    for (int i = 0; i < fullTime; i++) {
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

        mixture->concentrations[0] = mixture->concentrations[0] + h * 
            (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;
        mixture->concentrations[2] = mixture->concentrations[2] + h * 
            (q1 + 2 * q2 + 2 * q3 + q4) / 6.0;
        mixture->concentrations[1] = mixture->concentrations[1] + h *
            (r1 + 2 * r2 + 2 * r3 + r4) / 6.0;
        mixture->concentrations[1] = (mixture->initialConcentrations[0] -
            mixture->concentrations[0] + 3 * mixture->initialConcentrations[2] -
            3 * mixture->concentrations[2] + 
            2 * mixture->initialConcentrations[1]) / 
            2.0;

        //mixture->concentrations[1] = mixture->concentrations[1] + 
        //    (-rightSideForO3(time, mixture->concentrations[0],
        //        mixture->concentrations[2], mixture->concentrations[1]) -
        //     rightSideForO(time, mixture->concentrations[0], 
        //        mixture->concentrations[2], mixture->concentrations[1])) *
        //     0.5 * h;

        //mixture->concentrations[1] = mixture->concentrations[1] + 
        //    (-k1 - q1) * 0.5 * h;
        
        mixture->assertConcentrationsArePositive();

        mixture->molecularWeight = mixture->calculateMolecularWeight();
        mixture->volume = mixture->calculateMixtureVolume();
        mixture->temperature = mixture->calculateTemperature();

        time += h;

        h = alpha * abs(mixture->concentrations[0]) / rightSideForO(time,
            mixture->concentrations[0],
            mixture->concentrations[2],
            mixture->concentrations[1]);

        if (i % timeStepForOutput == 0) {
            printToFile();
        }
    }
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

    // ������������� ������� ����������, ���� / (����*�).
    const RealType r = 0.001985846;

    k = exp(log(t) * mixture->reactions[i].n - 
        mixture->reactions[i].activationEnergy / (r * t) + 
        2.30258 * mixture->reactions[i].log10A);

    return k;
}

RealType RungeKuttaMethod::calculateRateForBackReaction(int i, RealType kf)
{
    RealType t = mixture->temperature;
    // �������� ������ �������.
    RealType q = 0.0;
    RealType multProducts = 1.0;
    RealType multReagents = 1.0;
    int substanceNumber;
    int nMoles = 0;

    RealType *gibbs_energies = new RealType[mixture->nSubstances];

    for (int j = 0; j < mixture->nSubstances; j++) {
        gibbs_energies[j] = mixture->GibbsCalc_ext(j, t);
    }

    delete [] gibbs_energies;

    for (int j = 0; j < mixture->reactions[i].nProducts; j++) {
        substanceNumber = mixture->reactions[i].products[j];
        // ���������, ��� �������� �� �������� ��������� "�".
        if (substanceNumber == -1) {
            continue;
        }
        multProducts *= mixture->concentrations[substanceNumber];
        q += mixture->GibbsCalc_ext(substanceNumber, t);
        nMoles++;
    }

    for (int j = 0; j < mixture->reactions[i].nReagents; j++) {
        substanceNumber = mixture->reactions[i].products[j];
        // ���������, ��� �������� �� �������� ��������� "�".
        if (substanceNumber == -1) {
            continue;
        }
        multReagents *= mixture->concentrations[substanceNumber];
        q -= mixture->GibbsCalc_ext(substanceNumber, t);
        nMoles--;
    }

    // ��������� ����������.
    RealType kp;
    kp  = exp(-q / (mixture->R_J_OVER_MOL_K * t));
    kp *= exp(-log(10 * mixture->K_BOLTZMANN * t) * nMoles);

    //if (multProducts != 0.0) {
    //    kp  = exp(-q / (mixture->R_J_OVER_MOL_K * t));
    //    kp *= exp(-log(10 * mixture->K_BOLTZMANN * t) * nMoles);
    //}
    //else {
    //    kp = 1.0;
    //}

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
        1.0 / mixture->volume * 1.0e-3 << endl;
}

void RungeKuttaMethod::printHeadingToFile()
{
    outputFile << "t (s)\tO2\tO3\tO\tN (cm3)\tH (J/g)\t\tT (K)\t\t";
    outputFile << "Mu (g/mole)\t\tRho (g/cm3)" << endl;
}