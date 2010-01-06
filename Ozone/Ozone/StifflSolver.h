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
#include "Stiffl.h"

class StifflSolver : public Stiffl
{
public:
    /**
     * Конструктор класса.
     *
     * @param fileOfSubstances      имя файла с веществами
     * @param fileOfReactions       имя файла с реакциями
     */
	StifflSolver(int NYDIM_PAR, double *values, double t_begin, double t_end,
					 double t_step_begin,
					 const char *fileOfSubstances,
                     const char *fileOfReactions,
                     const char *fileOfMoleFractions);
    /**
     * Деструктор класса.
     */
    ~StifflSolver();
    /**
    * Производит интегрирование системы ОДУ.
    *
    * @param aFullTime временной интервал, на котором 
    * производится интегрирование.
    */
    void performIntegration(RealType aFullTime);
    /**
     * Обновляет значения мольных долей компонентов смеси.
     *
     * @param vf указатель на массив, в который пишутся
     * новые значения мольных долей.
     */
    void updateMoleFractions(RealType *vf);
    /**
     * Возвращает значение давления в смеси.
     *
     * @return давление в смеси, Па
     */
    RealType getPressure();
    /**
     * Смесь газов.
     */
    Mixture *mixture;
    /**
     * Возвращает объект смеси.
     *
     * @return указатель на объект класса Mixture.
     */
    Mixture *getMixture();

private:
    /**
     * Вычисляет значение скорости реакции по O3 при заданной температуре
     * и составе.
     *
     * @param concOfO   концентрация O, молекул / см**3
     * @param concOfO3  концентрация O3, молекул / см**3
     * @param concOfO2  концентрация O2, молекул / см**3
     * @return          значение скорости реакции относительно O3,
     * молекул / (см**3 * с)
     */
    RealType rightSideForO3(RealType concOfO, 
                            RealType concOfO3,
                            RealType concOfO2);
    /**
    * Вычисляет значение скорости реакции по O при заданной температуре
    * и составе.
    *
    * @param concOfO   концентрация O, молекул / см**3
    * @param concOfO3  концентрация O3, молекул / см**3
    * @param concOfO2  концентрация O2, молекул / см**3
    * @return          значение скорости реакции относительно O,
    * молекул / (см**3 * с)
    */
    RealType rightSideForO(RealType concOfO, 
                           RealType concOfO3, 
                           RealType concOfO2);
    /**
    * Вычисляет значение скорости реакции по O2 при заданной температуре
    * и составе.
    *
    * @param concOfO   концентрация O, молекул / см**3
    * @param concOfO3  концентрация O3, молекул / см**3
    * @param concOfO2  концентрация O2, молекул / см**3
    * @return          значение скорости реакции относительно O2,
    * молекул / (см**3 * с)
    */
    RealType rightSideForO2(RealType concOfO, 
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
    /**
     * Количество итераций в цикле интегрирования.
     */
    int nIterations;
	RealType k1f;
	RealType k2f;
	RealType k3f;
	RealType k1r;
	RealType k2r;
	RealType k3r;
	virtual int IFNSH();
	virtual int DIFFUN(double **YY, double *F);
	virtual void PEDERV();
};

#endif