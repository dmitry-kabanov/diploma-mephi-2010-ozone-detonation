/**
* @file
*
* @author  Кабанов Дмитрий <kabanovdmitry@gmail.com>
* @version %I%
*
* @section DESCRIPTION
*
* Реализация класса Mixture.
*/
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include "Mixture.h"
#include "Substance.h"
#include "Reaction.h"
using namespace std;

const RealType Mixture::R_J_OVER_MOL_K  = 8.31441;
const RealType Mixture::R_J_OVER_KMOL_K = 8314.41;
const RealType Mixture::AVOGADRO_NUMBER = 6.022e23;
const RealType Mixture::K_BOLTZMANN     = 1.38e-23;

Mixture::Mixture(const char *fileOfSubstances, 
                 const char *fileOfReactions,
                 const char *fileOfMoleFractions)
{
    readFileOfSubstances(fileOfSubstances);
    readFileOfReactions(fileOfReactions);
    readFileOfVolumeFractions(fileOfMoleFractions);
}

Mixture::~Mixture()
{
    for (int i = 0; i < nSubstances; ++i) {
        delete substances[i];
    }
    delete [] substances;
    delete [] concentrations;
    delete [] initialConcentrations;
    delete [] previousConcentrations;
    delete [] reactions;
    delete [] volumeFractions;
}

void Mixture::readFileOfSubstances(const char *filename)
{
    char ch;
    std::string word;
    std::fstream iFile;

    iFile.open(filename, std::ios_base::in);    
    iFile >> nSubstances >> ch;

    allocateMemoryForSubstances();

    for (int i = 0; i < nSubstances; i++) {
        iFile >> substances[i]->nameOfSubstance;
        
        iFile >> ch; // Ноль
        iFile >> ch; // Ноль
        iFile >> substances[i]->enthalpyOfFormation;
        iFile >> substances[i]->molecularWeight;
        iFile >> substances[i]->nTemperatureRanges;

        substances[i]->allocateMemoryForTemperatureRanges(
            substances[i]->nTemperatureRanges
        );

        for (int j = 0; j < substances[i]->nTemperatureRanges; j++) {
            iFile >> substances[i]->temperatureLow[j];
            iFile >> substances[i]->temperatureHigh[j];

            for (int k = 0; k < 7; k++) {
                iFile >> substances[i]->a[j][k];
            }
        }

        iFile >> ch; // Ноль
        iFile >> ch; // Ноль
        iFile >> word; // User
    }
}

void Mixture::readFileOfReactions(const char *filename)
{
    char ch;
    ifstream iFile(filename);
    
    iFile >> nReactions >> ch;

    allocateMemoryForReactions();

    for (int i = 0; i < nReactions; i++) {
        iFile >> reactions[i].typeOfReaction;
        iFile >> ch; // Ноль
        iFile >> ch; // Количество температурных интервалов
        iFile >> ch; // Единица

        iFile >> reactions[i].temperatureLow;
        iFile >> reactions[i].temperatureHigh;
        iFile >> reactions[i].log10A;
        iFile >> reactions[i].n;
        iFile >> reactions[i].activationEnergy;
        iFile >> reactions[i].nameOfReaction;
    }
    fillReagents();
    fillProducts();

    iFile.close();
}

void Mixture::readFileOfVolumeFractions(const char *filename)
{
    volumeFractions = new RealType[nSubstances];

    ifstream iFile(filename);
    string substanceName;
    RealType fraction;

    for (int i = 0; i < nSubstances; i++) {
        iFile >> substanceName >> fraction;
        for (int j = 0; j < nSubstances; j++) {
            if (substanceName == substances[j]->nameOfSubstance) {
                volumeFractions[j] = fraction;
            }
        }
    }

    iFile.close();

    RealType sum = 0.0;
    for (int i = 0; i < nSubstances; i++) {
        sum += volumeFractions[i];
    }

    for (int i = 0; i < nSubstances; i++) {
        volumeFractions[i] /= sum;
    }

    allocateMemoryForConcentrations();
}

void Mixture::allocateMemoryForSubstances()
{
    substances = new Substance*[nSubstances];
    for (int i = 0; i < nSubstances; ++i) {
        substances[i] = new Substance();
    }
}

void Mixture::allocateMemoryForConcentrations()
{
    concentrations         = new RealType[nSubstances];
    initialConcentrations  = new RealType[nSubstances];
    previousConcentrations = new RealType[nSubstances];
}

void Mixture::allocateMemoryForReactions()
{
    reactions = new Reaction[nReactions];
}

RealType Mixture::calculateEnthalpy(int i, RealType t)
{
    int numberOfTemperatureRange;
    int j;

    for (j = 0; j < substances[i]->nTemperatureRanges; j++) {
        if (t <= substances[i]->temperatureHigh[j]) {
            numberOfTemperatureRange = j;
            break;
        }
        numberOfTemperatureRange = j;
    }
    

    RealType res;
    RealType *a = substances[i]->a[numberOfTemperatureRange];
    RealType x = t * 1E-4;
    res = -2*a[2]/x-a[3]+a[1]*x+a[4]*x*x+2*a[5]*x*x*x+3*a[6]*x*x*x*x;
    res *= 1E4;
    res *= 1E3;
    return res;
}

void Mixture::fillReagents()
{
    string s;
    string ch;

    for (int i = 0; i < nReactions; i++) {
        s = reactions[i].nameOfReaction;
        int lside = s.find("<");
        int k = 0;
        string subs[5];
        for (int j = 0; j < lside; j++) {
            ch = s.substr(j, 1);
            if (ch != "+") {
                subs[k] += ch;
            }
            else {
                k++;
            }
        }
        reactions[i].nReagents = k + 1;
        reactions[i].reagents = new int[reactions[i].nReagents];

        for (int j = 0; j < reactions[i].nReagents; j++) {
            if (subs[j] == "M") {
                reactions[i].reagents[j] = -1; // Признак вещества M
                continue;
            }
            for (int jj = 0; jj < nSubstances; jj++) {
                if (subs[j] == substances[jj]->nameOfSubstance) {
                    reactions[i].reagents[j] = jj;
                }
            }
        }
    }
}

void Mixture::fillProducts()
{
    string s;
    string ch;

    for (int i = 0; i < nReactions; i++) {
        s = reactions[i].nameOfReaction;
        int rside = s.find(">");
        int k = 0;
        string subs[5];
        for (unsigned int j = rside + 1; j < s.size(); j++) {
            ch = s.substr(j, 1);
            if (ch != "+") {
                subs[k] += ch;
            }
            else {
                k++;
            }
        }
        reactions[i].nProducts = k + 1;
        reactions[i].products = new int[reactions[i].nProducts];

        for (int j = 0; j < reactions[i].nProducts; j++) {
            if (subs[j] == "M") {
                reactions[i].products[j] = -1; // Признак вещества M
                continue;
            }
            for (int jj = 0; jj < nSubstances; jj++) {
                if (subs[j] == substances[jj]->nameOfSubstance) {
                    reactions[i].products[j] = jj;
                }
            }
        }
    }
}

RealType Mixture::calculateCp(RealType t)
{
    RealType res;
    sumPolynomialCoeffs(t);
    RealType *a = sumThermCoeffs;
    RealType x = t*1E-4;
    res=2*a[2]/x/x+a[1]+2*a[4]*x+6*a[5]*x*x+12*a[6]*x*x*x;	
    res/=AVOGADRO_NUMBER;
    return res;
}

RealType Mixture::calculateSubstanceCp(int i, RealType t)
{
    RealType res;
    int nInt;
    int j;

    for (j = 0; j < substances[i]->nTemperatureRanges; j++) {
        if (t <= substances[i]->temperatureHigh[j]) {
            nInt = j;
            break;
        }
        nInt = j;
    }

    RealType* a = substances[i]->a[nInt];
    RealType x  = t*1E-4;
    
    res  = 2*a[2]/x/x+a[1]+2*a[4]*x+6*a[5]*x*x+12*a[6]*x*x*x;	
    //res *= 1.0e4;
    res *= 1.0e3;
    //res /= AVOGADRO_NUMBER;

    return res;
}

void Mixture::sumPolynomialCoeffs(RealType t)
{
    int i,j,k, nInt ;

    for(k = 0; k < 7; k++)
    {
        sumThermCoeffs[k] = 0.0;

        for(i=0;i<nSubstances;i++)
        {
            for(j = 0; j < substances[i]->nTemperatureRanges; j++)
            {
                if(t<=substances[i]->temperatureHigh[j])
                {
                    nInt=j;
                    break;
                }
                nInt=j;
            }
            sumThermCoeffs[k]+=concentrations[i]*substances[i]->a[nInt][k];
        }
    }
}

RealType Mixture::calculateGibbsEnergy(int i)
{
    RealType result;
    RealType T = temperature;
    RealType *a;
    int nInt;
    int j;

    for (j = 0; j < substances[i]->nTemperatureRanges; j++) {
        if (T <= substances[i]->temperatureHigh[j]) {
            nInt = j;
            break;
        }
        nInt = j;
    }

    a = substances[i]->a[nInt];
    result  = -a[0]*T;
    result -= a[1]*T*(log(T)-4*log(10.0));
    result -= a[2]*1E8/T;
    result -= a[3]*1E4;
    result -= a[4]*T*T*1E-4;
    result -= a[5]*T*T*T*1E-8;
    result -= a[6]*T*T*T*T*1E-12;
    result  = substances[i]->enthalpyOfFormation * 1E3 + result;
    return result;
}

RealType Mixture::calculateMolecularWeight()
{
    RealType sumConc = 0.0;
    RealType weight = 0.0;

    for (int i = 0; i < nSubstances; i++) {
        weight  += substances[i]->molecularWeight * concentrations[i];
        sumConc += concentrations[i];
    }

    weight /= sumConc;

    return weight;
}

void Mixture::calculateInitialMolecularWeight()
{
    molecularWeight = 0.0;
    RealType sumOfFractions  = 0.0;

    for (int i = 0; i < nSubstances; i++) {
        molecularWeight += substances[i]->molecularWeight * volumeFractions[i];
        sumOfFractions  += volumeFractions[i];
    }

    molecularWeight /= sumOfFractions;
    initialMolecularWeight = molecularWeight;
}

RealType Mixture::calculateFullEnergy(RealType pressure, RealType volume)
{
    RealType energy = 0.0;
    RealType sumOfFractions = 0.0;

    for (int i = 0; i < nSubstances; i++) {
        // Коэффициент 1.0e6 для перевода из кДж / моль в Дж / кмоль.
        energy += (calculateEnthalpy(i, temperature) +
            substances[i]->enthalpyOfFormation * 1.0e6) *
            volumeFractions[i];
        sumOfFractions += volumeFractions[i];
    }

    energy = energy / (molecularWeight * sumOfFractions) -
        pressure * volume;

    return energy;
}

RealType Mixture::calculateTemperature()
{
    RealType sumConc = 0.0;
    RealType energy1 = 0.0;
    RealType energy2 = 0.0;
    RealType mixtureHeatCapacity;
    RealType mixtureHeatCapacityJOverKgK;
    RealType mixtureHeatCapacityJOverKmoleK;
    RealType mixEnthalpyOfFormation = 0.0;
    RealType mixEnthalpy = 0.0;

    for (int i = 0; i < nSubstances; i++) {
        sumConc += concentrations[i];
    }

    // Теплоёмкость в Дж / (см**3 * К)
    mixtureHeatCapacity = calculateCp(temperature);
    
    // Теплоёмкость в Дж / (кг * К)
    mixtureHeatCapacityJOverKgK = mixtureHeatCapacity * 1.0e6 * volume;
    
    // Теплоёмкость в Дж / (кмоль * К)
    mixtureHeatCapacityJOverKmoleK = mixtureHeatCapacityJOverKgK * 
        molecularWeight;

    for (int i = 0; i < nSubstances; i++) {
        // Коэффициент 1.0e6 для перевода из кДж/моль в Дж/кмоль.
        mixEnthalpyOfFormation += concentrations[i] * 
            substances[i]->enthalpyOfFormation * 1.0e6;
    }
    mixEnthalpyOfFormation /= sumConc;

    for (int i = 0; i < nSubstances; i++) {
        mixEnthalpy += concentrations[i] * calculateEnthalpy(i, temperature);
    }
    mixEnthalpy /= sumConc;
    mixEnthalpy += mixEnthalpyOfFormation;
    mixEnthalpy /= molecularWeight; // Переводим в Дж / кг.
    energy1 =  - mixtureHeatCapacityJOverKgK * temperature + mixEnthalpy;

    // Считаем выражение в знаменателе.
    for (int i = 0; i < nSubstances; i++) {
        energy2 += concentrations[i] * 
            (R_J_OVER_KMOL_K / substances[i]->molecularWeight);
    }
    energy2 /= sumConc;

    return (fullEnergy - energy1) / (mixtureHeatCapacityJOverKgK - energy2);
}

RealType Mixture::calculateInitialTemperature()
{
    // Номер итерации.
    int k = 0;
    // Максимальное количество итераций.
    const int maxIterations = 100;
    // Точность нахождения температуры.
    RealType precision = 1e-3;
    // Задаем начальное приближение.
    RealType t = 2000;
    // Уравнение, связывающее удельную внутреннюю энергию
    // с удельной энтальпией.
    RealType func = calculateMixtureEnthalpy(t) - R_J_OVER_KMOL_K * t / molecularWeight - fullEnergy;
    // Производная по температуре от функции func.
    RealType dfunc;

    while (abs(func) >= precision) {
        dfunc = (calculateMixtureCp(t) - R_J_OVER_KMOL_K / molecularWeight);
        t = t - func / dfunc;
        func = calculateMixtureEnthalpy(t) - R_J_OVER_KMOL_K * t / molecularWeight - fullEnergy;
        k++;

        if (k > maxIterations) {
            cout << "Initial temperature function: " << 
                "limit of iterations (" << maxIterations << 
                ") was reached." << endl;
            cout << "Internal energy = " << fullEnergy << " J/kg" << endl;
            cout << "Density         = " << 1.0 / volume << endl;
            cout << "X(O)            = "  << volumeFractions[0] << endl;
            cout << "X(O2)           = "  << volumeFractions[1] << endl;
            cout << "X(O3)           = "  << volumeFractions[2] << endl;
            cout << "T = " << t << endl;
            //exit(-1);
            return t;
        }
    }

    return t;
}

void Mixture::fillUpMixture(RealType internalEnergy, RealType density,
                                RealType *volFracts)
{
    volume = 1.0 / density;
    fullEnergy = internalEnergy; // Дж кг-1

    RealType sumFractions = 0.0;
    for (int i = 0; i < nSubstances; ++i) {
        sumFractions += volFracts[i];
    }

    for (int i = 0; i < nSubstances; ++i) {
        volumeFractions[i] = volFracts[i] / sumFractions;
    }

    calculateInitialMolecularWeight();

    temperature = calculateInitialTemperature();
    oldTemperature = temperature;

    pressure    = (R_J_OVER_KMOL_K * temperature) / (molecularWeight * volume);

    RealType densityKmoleOverM3 = pressure / 
        (temperature * R_J_OVER_KMOL_K); // кмоль м-3

    for (int i = 0; i < nSubstances; i++) {
        concentrations[i] = 1.0e-3 * densityKmoleOverM3 * 
            AVOGADRO_NUMBER * volumeFractions[i];
        initialConcentrations[i] = concentrations[i];
    }
}

RealType Mixture::calculatePressure()
{
    return (R_J_OVER_KMOL_K * temperature) / (molecularWeight * volume);
}

RealType Mixture::calculateMixtureCp(RealType t)
{
    RealType mixtureHeatCapacity = 0.0;
    RealType mixtureHeatCapacityJOverKgK;

    // Теплоёмкость в Дж / (кмоль * К)
    for (int i = 0; i < nSubstances; i++) {
        mixtureHeatCapacity += volumeFractions[i] * calculateSubstanceCp(i, t);
    }
    
    // Теплоёмкость в Дж / (кг * К)
    mixtureHeatCapacityJOverKgK = mixtureHeatCapacity / molecularWeight;
    
    // Теплоёмкость в Дж / (кмоль * К)
    // mixtureHeatCapacityJOverKmoleK = mixtureHeatCapacityJOverKgK * 
    //    molecularWeight;

    return mixtureHeatCapacityJOverKgK;
}

RealType Mixture::calculateMixtureEnthalpy(RealType t)
{
    RealType mixEnthalpyOfFormation = 0.0;
    RealType mixEnthalpy = 0.0;

    for (int i = 0; i < nSubstances; i++) {
        // Коэффициент 1.0e6 для перевода из кДж/моль в Дж/кмоль.
        mixEnthalpyOfFormation += volumeFractions[i] * 
            substances[i]->enthalpyOfFormation * 1.0e6;
    }

    for (int i = 0; i < nSubstances; i++) {
        mixEnthalpy += volumeFractions[i] * calculateEnthalpy(i, t);
    }
    
    mixEnthalpy += mixEnthalpyOfFormation;
    mixEnthalpy /= molecularWeight; // Переводим в Дж / кг.

    return mixEnthalpy;
}

RealType Mixture::getOldTemperature()
{
    return oldTemperature;
}

RealType Mixture::getTemperature()
{
    return temperature;
}

void Mixture::setConcentrations(RealType t, RealType density)
{
    volume = 1.0 / density;

    calculateInitialMolecularWeight();

    temperature = t;

    pressure = (R_J_OVER_KMOL_K * temperature) / (molecularWeight * volume);

    RealType densityKmoleOverM3 = pressure / 
        (temperature * R_J_OVER_KMOL_K); // кмоль м-3

    for (int i = 0; i < nSubstances; i++) {
        concentrations[i] = 1.0e-3 * densityKmoleOverM3 * 
            AVOGADRO_NUMBER * volumeFractions[i];
        initialConcentrations[i] = concentrations[i];
    }
}

void Mixture::setVolumeFractions(const RealType *fractions)
{
    RealType sumFractions = 0.0;
    for (int i = 0; i < nSubstances; ++i) {
        sumFractions += fractions[i];
    }

    for (int i = 0; i < nSubstances; ++i) {
        volumeFractions[i] = fractions[i] / sumFractions;
    }
}

RealType Mixture::getMolecularWeight()
{
    return molecularWeight;
}

RealType Mixture::calculateOldPressure()
{
    return (R_J_OVER_KMOL_K * oldTemperature) / (initialMolecularWeight * volume);
}

void Mixture::fillUpMixtureWithTAndP(RealType t, RealType p, RealType *vf)
{
    calculateInitialMolecularWeight();

    RealType sumFractions = 0.0;

    for (int i = 0; i < nSubstances; ++i) {
        volumeFractions[i] = vf[i];
        sumFractions += volumeFractions[i];
    }

    for (int i = 0; i < nSubstances; ++i) {
        volumeFractions[i] /= sumFractions;
    }

    temperature = t;
	oldTemperature = t;
    pressure    = p;

    volume = (R_J_OVER_KMOL_K * temperature) / (molecularWeight * pressure);

    RealType densityKmoleOverM3 = pressure / 
        (temperature * R_J_OVER_KMOL_K); // кмоль м-3

    for (int i = 0; i < nSubstances; i++) {
        concentrations[i] = 1.0e-3 * densityKmoleOverM3 * 
            AVOGADRO_NUMBER * volumeFractions[i];
        initialConcentrations[i] = concentrations[i];
    }

    fullEnergy = calculateMixtureEnthalpy(temperature) - pressure * volume;
}