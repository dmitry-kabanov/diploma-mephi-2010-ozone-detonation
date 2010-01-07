/**
 * @file
 *
 * @author  Кабанов Дмитрий <kabanovdmitry@gmail.com>
 * @version %I%
 *
 * @section DESCRIPTION
 *
 * Класс вещества для хранения его различных свойств.
 */
#ifndef SUBSTANCE_H
#define SUBSTANCE_H

#include <string>
#include "RealType.h"

class Substance
{
public:
    /**
     * Конструктор класса.
     */
    Substance();
    /**
     * Деструктор класса.
     */
    ~Substance();
    /**
     * Название вещества.
     */
    std::string nameOfSubstance;
    /**
     * Энтальпии образования веществ, кДж/моль.
     */
    RealType enthalpyOfFormation;
    /**
     * Молекулярные веса веществ, кг/кмоль.
     */
    RealType molecularWeight;
    /**
     * Количество аппроксимируемых температурных диапазонов.
     */
    int nTemperatureRanges;
    /**
     * Нижние границы температурных диапазонов, К.
     */
    RealType *temperatureLow;
    /**
     * Верхние границы температурных диапазонов, К.
     */
    RealType *temperatureHigh;
    /**
     * Коэффициенты для полиномиальных зависимостей термодинамических величин.
     */
    RealType **a;
    /**
     * Выделяет память для температурных диапазонов.
     * 
     * @param n количество температурных диапазонов.
     */
    void allocateMemoryForTemperatureRanges(int n);
};

#endif