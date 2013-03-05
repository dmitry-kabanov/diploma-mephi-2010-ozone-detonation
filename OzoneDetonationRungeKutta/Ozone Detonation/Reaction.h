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
     * Деструктор класса.
     */
    ~Reaction();
    /**
     * Тип реакции. Если равен 0, то реакция зависит от температуры,
     * если 1, то реакция зависит от давления.
     */
    int typeOfReaction;
    /**
     * Нижняя температура для параметров реакции.
     */
    RealType temperatureLow;
    /**
     * Верхняя температура для параметров реакции.
     */
    RealType temperatureHigh;
    /**
     * Десятичный логарифм предэкспоненциального множителя.
     */
    RealType log10A;
    /**
     * Показатель степени у температуры в модифицированном 
     * законе Аррениуса.
     */
    RealType n;
    /**
     * Энергия активации для данной реакции.
     */
    RealType activationEnergy;
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
};

#endif