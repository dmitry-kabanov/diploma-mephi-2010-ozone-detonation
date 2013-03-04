/**
 * @file
 *
 * @author  Кабанов Дмитрий <kabanovdmitry@gmail.com>
 * @version %I%
 *
 * @section DESCRIPTION
 *
 * Класс реакции для хранения параметров химической реакции.
 */
#ifndef REACTION_H
#define REACTION_H

#include <string>
#include "RealType.h"

class Reaction
{
public:
    /**
     * Конструктор класса.
     */
    Reaction();
    /**
     * Деструктор класса. Освобождает выделенную память.
     */
    ~Reaction();
	/**
	 * Выделяет память для массивов параметров реакции.
	 */
	void allocateMemoryForParameters();
	/**
	 * Выделяет память для массива коэффициентов эффективности столкновения.
	 *
	 * @param n количество веществ в смеси.
	 */
	void allocateMemoryForCollisionEfficiency(int n);
	/**
	 * Вычисляет константу скорости прямой реакции при заданной температуре.
	 *
	 * @param t температура, К
	 * @return вычисленное значение константы скорости прямой реакции
	 */
	RealType calculateConstantRate(RealType t);
	/**
	 * Вычисляет скорость реакции при заданной температуре.
	 *
	 * @param YY			массив значений концентраций компонентов смеси
	 * @param t				температура, К
	 * @param gibbsEnergies значения энергий Гиббса для всех компонентов смеси
	 * @param nSpecies		количество веществ в смеси
	 * @return вычисленное	значение скорости реакции
	 */
	RealType calculateReactionRate(double *Y, RealType t, RealType *gibbsEnergies, int nSpecies);
	RealType multiplicationOfReagents;
	RealType multiplicationOfProducts;
	RealType q;
	int nMoles;
	RealType kc;
	RealType kf;
	bool withThirdBody;
	RealType reactionRate_;
	/**
     * Тип реакции. Если равен 0, то реакция зависит от температуры,
     * если 1, то реакция зависит от давления.
     */
    int typeOfReaction;
	/**
	 * Количество температурных интервалов, для которых приводятся 
	 * параметры реакции.
	 */
	int nTemperatureRanges;
    /**
     * Нижняя температура для параметров реакции.
     */
    RealType *temperatureLow;
    /**
     * Верхняя температура для параметров реакции.
     */
    RealType *temperatureHigh;
    /**
     * Десятичный логарифм предэкспоненциального множителя.
     */
    RealType *log10A;
    /**
     * Показатель степени у температуры в модифицированном 
     * законе Аррениуса.
     */
    RealType *n;
    /**
     * Энергия активации для данной реакции, ккал моль-1.
     */
    RealType *activationEnergy;
    /**
     * Формула реакции.
     */
    std::string nameOfReaction;
    /**
     * Номера исходных веществ в данной реакции.
     */
    int *reagents;
    /**
     * Количество исходных веществ в данной реакции.
     */
    int nReagents;
    /**
     * Номера продуктов в данной реакции.
     */
    int *products;
    /**
     * Количество продуктов в данной реакции.
     */
    int nProducts;
	/**
	 * Направление реакции. Если 0, то реакция обратимая, если 1,
	 * то необратимая.
	 */
	int direction;
	/**
	 * Массив коэффициентов эффективности столкновения.
	 */
	RealType *pEff;
	/**
	 * Флаг, указывающий, есть ли коэффициенты эффективности столкновения.
	 */
	bool nEff;
};

#endif