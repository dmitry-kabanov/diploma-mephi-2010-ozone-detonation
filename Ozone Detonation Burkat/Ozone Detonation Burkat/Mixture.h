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
     * @param fileOfSubstances      имя файла со свойствами веществ
     * @param fileOfReactions       имя файла с параметрами реакций
     * @param fileOfMoleFractions   имя файла с значениями мольных долей
     */
    Mixture(const char *fileOfSubstances,
            const char *fileOfReactions,
            const char *fileOfMoleFractions);
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
     * Объемные доли компонентов смеси.
     */
    RealType *volumeFractions;
    /**
     * Массовые доли компонентов смеси.
     */
    RealType *y_;
    /**
     * Вычисляет изменение энтальпии от температуры 298.15 К 
     * до температуры t для i-го вещества.
     *
     * @param i порядковый номер вещества в массиве веществ
     * @param t температура, К
     * @return вычисленное приращение энтальпии, Дж/кмоль
     */
    RealType calculateSubstanceEnthalpy(int i, RealType t);
    /**
     * Вычисляет значение энтропии для i-го вещества 
     * при температуре t.
     *
     * @param i порядковый номер вещества в массиве веществ
     * @param t температура, К
     * @return вычисленное значение энтропии, Дж кг-1 К-1
     */
    RealType calculateSubstanceEntropy(int i, RealType t);
    /**
     * Вычисляет энергию Гиббса для i-го вещества
     * при температуре t.
     *
     * @param i порядковый номер вещества в массиве веществ.
     * @param t температура, К
     * @return вычисленное значение энергии Гиббса, Дж кг-1.
     */
    RealType calculateSubstanceGibbsEnergy(int i, RealType t);
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
     * Температура смеси в начале шага интегрирования, К. 
     * Вспомогательная переменная.
     */
    RealType oldTemperature;
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
    RealType *previousConcentrations;
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
    RealType getOldTemperature();
    RealType getTemperature();
    /**
     * Вычисляет температуру при заданном составе смеси.
     *
     * @return температура, К
     */
    RealType calculateTemperature();
    RealType calculateInitialTemperature();
    /**
    * Устанавливает термодинамическое состояние смеси 
    * по заданным внутренней энергии, плотности и мольным долям
    * компонентов.
     *
     * @param internalEnergy внутренняя энергия, Дж кг-1
     * @param density        плотность, кг м-3
     * @param x              массив мольных долей компонентов смеси
     */
    void setStateWithURhoX(RealType internalEnergy, 
        RealType density, 
        RealType* x);
    /**
     * Устанавливает термодинамическое состояние смеси 
     * по заданным температуре, давлению и мольным долям
     * компонентов.
     *
     * @param t  температура, К
     * @param p  давление, Па
     * @param vf массив мольных долей компонентов смеси
     */
    void setStateWithTPX(RealType t, RealType p, RealType *x);
    RealType calculatePressure();
    RealType calculateOldPressure();
    RealType calculateMixtureCp(RealType t);
    RealType calculateMixtureEnthalpy(RealType t);
    RealType calculateSubstanceCp(int i, RealType t);
    /**
     * Устанавливает концентрации компонентов смеси 
     * при заданной температуре и плотности.
     *
     * @param t       температура, К
     * @param density плотность, кг м-3
     */
    void setConcentrations(RealType t, RealType density);
    /**
     * Устанавливает объемные доли компонентов смеси.
     *
     * @param fractions объемные доли
     */
    void setVolumeFractions(const RealType *fractions);
    /**
     * Возвращает молекулярный вес смеси.
     */
    RealType getMolecularWeight();
};

#endif