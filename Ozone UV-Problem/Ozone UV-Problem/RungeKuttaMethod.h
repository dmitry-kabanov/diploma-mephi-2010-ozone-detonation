/**
* @file
*
* @author  Кабанов Дмитрий <kabanovdmitry@gmail.com>
* @version %I%
*
* @section DESCRIPTION
*
* Класс интегрирования системы ОДУ химической кинетики.
*/
#ifndef RUNGE_KUTTA_METHOD_H
#define RUNGE_KUTTA_METHOD_H

#include <iostream>
#include <fstream>
#include "RealType.h"
#include "Mixture.h"
#include "Reaction.h"
#include "Substance.h"

class RungeKuttaMethod
{
public:
    /**
     * Конструктор класса.
     *
     * @param initialPressure       начальное давление, Па
     * @param initialTemperature    начальная температура, К
     * @param initialTimeStep       начальный временной шаг, с
     * @param aFullTime             количество шагов интегрирования
     * @param aTimeStepForOutput    через сколько временных шагов выводить
     * @param fileOfSubstances      имя файла с веществами
     * @param fileOfReactions       имя файла с реакциями
     * @param fileOfVolumeFractions имя файла с объемными долями
     */
    RungeKuttaMethod(RealType    initialPressure,
                     RealType    initialTemperature,
                     RealType    initialTimeStep,
                     int         aFullTime,
                     int         aTimeStepForOutput,
                     const char *fileOfSubstances,
                     const char *fileOfReactions,
                     const char *fileOfVolumeFractions);
    /**
     * Деструктор класса.
     */
    ~RungeKuttaMethod();
private:
    /**
     * Производит интегрирование системы ОДУ.
     */
    void performIntegration();

    /**
     * Вычисляет значение скорости реакции по O3 при заданной температуре
     * и составе.
     *
     * @param t         температура, К
     * @param concOfO   концентрация O, молекул / см**3
     * @param concOfO3  концентрация O3, молекул / см**3
     * @param concOfO2  концентрация O2, молекул / см**3
     * @return          значение скорости реакции относительно O3,
     * молекул / (см**3 * с)
     */
    RealType rightSideForO3(RealType t, 
                            RealType concOfO, 
                            RealType concOfO3,
                            RealType concOfO2);
    /**
    * Вычисляет значение скорости реакции по O при заданной температуре
    * и составе.
    *
    * @param t         температура, К
    * @param concOfO   концентрация O, молекул / см**3
    * @param concOfO3  концентрация O3, молекул / см**3
    * @param concOfO2  концентрация O2, молекул / см**3
    * @return          значение скорости реакции относительно O,
    * молекул / (см**3 * с)
    */
    RealType rightSideForO(RealType t, 
                           RealType concOfO, 
                           RealType concOfO3, 
                           RealType concOfO2);
    /**
    * Вычисляет значение скорости реакции по O2 при заданной температуре
    * и составе.
    *
    * @param t         температура, К
    * @param concOfO   концентрация O, молекул / см**3
    * @param concOfO3  концентрация O3, молекул / см**3
    * @param concOfO2  концентрация O2, молекул / см**3
    * @return          значение скорости реакции относительно O2,
    * молекул / (см**3 * с)
    */
    RealType rightSideForO2(RealType t, 
                           RealType concOfO, 
                           RealType concOfO3, 
                           RealType concOfO2);
    /**
     * Вычисляет константу скорости прямой реакции.
     *
     * @param i  порядковый номер реакции в массиве реакций
     * @return   значение константы скорости прямой реакции.
     */
    RealType calculateRateForForwardReaction(int i);
    /**
    * Вычисляет константу скорости обратной реакции.
    *
    * @param i  порядковый номер реакции в массиве реакций
    * @param kf константа скорости прямой реакции
    * @return   значение константы скорости обратной реакции.
    */
    RealType calculateRateForBackReaction(int i, RealType kf);
    /**
     * Временной шаг интегрирования, с
     */
    RealType h;
    /**
     * Смесь газов.
     */
    Mixture *mixture;
    /**
     * Количество временных шагов интегрирования.
     */
    int fullTime;
    /**
     * Количество временных шагов, через которые происходит вывод в файл.
     */
    int timeStepForOutput;
    /**
     * Время, с
     */
    RealType time;
    /**
     * Выходной файл.
     */
    std::ofstream outputFile;
    /**
     * Выводит результаты расчёта в файл.
     */
    void printToFile();
    /**
     * Выводит заголовок в файл.
     */
    void printHeadingToFile();
};

#endif