/**
* @file
*
* @author  Кабанов Дмитрий <kabanovdmitry@gmail.com>
* @version %I%
*
* @section DESCRIPTION
*
* Класс смеси.
*/
#ifndef MIXTURE_H
#define MIXTURE_H

#include "RealType.h"

class Substance;
class Reaction;

class Mixture
{
public:
    /**
     * Конструктор класса.
     *
     * @param p0                    начальное давление, Па
     * @param t0                    начальная температура, К
     * @param fileOfSubstances      имя файла со свойствами веществ
     * @param fileOfReactions       имя файла с параметрами реакций
     * @param fileOfVolumeFractions имя файла с объемными долями газов
     */
    Mixture(RealType p0, RealType t0,
            const char *fileOfSubstances,
            const char *fileOfReactions,
            const char *fileOfVolumeFractions);
    /**
     * Деструктор класса.
     */
    ~Mixture();
    /**
     * Читает из файла filename вещества и их свойства.
     *
     * @param filename имя файла с веществами
     */
    void readFileOfSubstances(const char *filename);
    /**
     * Читает из файла filename реакции и их параметры.
     *
     * @param filename имя файла с реакциями
     */
    void readFileOfReactions(const char *filename);
    /**
     * Читает из файла filename имена веществ и их объемные доли.
     */
    void readFileOfVolumeFractions(const char *filename);
    /**
     * Объемные доли веществ.
     */
    RealType *volumeFractions;
    /**
     * Вычисляет изменение энтальпии от температуры 298.15 К 
     * до температуры t для i-го вещества.
     *
     * @param i порядковый номер вещества в массиве веществ
     * @param t температура, К
     * @return вычисленное приращение энтальпии, Дж/кмоль
     */
    RealType calculateEnthalpy(int i, RealType t);
    /**
     * Вычисляет энергию Гиббса для i-го вещества
     * при температуре t.
     *
     * @param i порядковый номер вещества в массиве веществ.
     * @return вычисленное значение энергии Гиббса, Дж / моль.
     */
    RealType Mixture::calculateGibbsEnergy(int i);
    /**
     * Вычисляет изобарную теплоемкость для смеси
     * при температуре t.
     *
     * @param t температура, К
     * @return вычисленное значение теплоемкости, Дж / (см**3 * К)
     */
    RealType calculateCp(RealType t);
    /**
     * Массив веществ.
     */
    Substance **substances;
    /**
     * Массив реакций.
     */
    Reaction *reactions;
    /**
     * Температура смеси, К.
     */
    RealType temperature;
    /**
     * Полная энергия системы.
     */
    RealType fullEnergy;
    /**
     * Количество веществ.
     */
    int nSubstances;
    /**
     * Количество реакций.
     */
    int nReactions;
    /**
     * Выделяет в куче память для массива веществ.
     */
    void allocateMemoryForSubstances();
    /**
     * Выделяет в куче память для массива реакций.
     */
    void allocateMemoryForReactions();
    /**
     * Заполняет массив реагентов.
     */
    void fillReagents();
    /**
     * Заполняет массив продуктов реакции.
     */
    void fillProducts();
    /**
     * Универсальная газовая постоянная, Дж / (моль * К).
     */
    static const RealType R_J_OVER_MOL_K;
    /**
     * Универсальная газовая постоянная, Дж / (кмоль * К).
     */
    static const RealType R_J_OVER_KMOL_K;
    /**
     * Константа Больцмана, Дж / К.
     */
    static const RealType K_BOLTZMANN;
    /**
    * Число Авогадро, 1 / моль.
    */
    static const RealType AVOGADRO_NUMBER;
    /**
     * Массив концентраций веществ.
     */
    RealType *concentrations;
    /**
     * Массив начальных концентраций.
     */
    RealType *initialConcentrations;
    /**
     * Выделяет в куче память для массива концентраций веществ в смеси.
     */
    void allocateMemoryForConcentrations();
    /**
     * Вычисляет молекулярный вес смеси в кг / кмоль.
     */
    RealType calculateMolecularWeight();
    /**
     *  Вычисляет начальный молекулярный вес смеси в кг / кмоль.
     */
    void calculateInitialMolecularWeight();
    /**
     * Молекулярный вес смеси, кг / кмоль.
     */
    RealType molecularWeight;
    /**
     * Начальный молекулярный вес, кг / кмоль.
     */
    RealType initialMolecularWeight;
    /**
     * Вычисляет полную энергию системы в Дж / кг.
     *
     * @param pressure давление в системе, Па
     * @param volume   удельный объем системы, м**3 / кг
     * @return внутренняя энергия системы, Дж / кг
     */
    RealType calculateFullEnergy(RealType pressure, RealType volume);
    /**
     * Вычисляет изобарную теплоемкость системы.
     */
    /**
     * Суммирует термодинамические коэффициенты при соответствующих степенях 
     * в полиномах для всех веществ в смеси.
     */
    void sumPolynomialCoeffs(RealType t);
    /**
     * Суммарные коэффициенты при соответствующих степенях полинома.
     */
    RealType sumThermCoeffs[8];
    /**
     * Давление смеси, Па.
     */
    RealType pressure;
    /**
     * Удельный объем смеси, м**3 / кг.
     */
    RealType volume;
    /**
     * Вычисляет температуру при заданном составе смеси.
     *
     * @return температура, К
     */
    RealType calculateTemperature();
};

#endif