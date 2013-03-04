/**
* @file
*
* @author  Кабанов Дмитрий <kabanovdmitry@gmail.com>
* @version $Id$
*
* @section DESCRIPTION
*
* Объявление класса StifflSolver.
*/
#ifndef STIFFL_SOLVER_H
#define STIFFL_SOLVER_H

#include <iostream>
#include <fstream>
#include "RealType.h"
#include "Mixture.h"
#include "Reaction.h"
#include "Substance.h"
#include "Stiffl.h"

/**
 * Класс интегрирования системы ОДУ химической кинетики. Для интегрирования 
 * использует метод STIFF, представляющий собой модификацию метода Гира.
 * См. Захаров А. Ю., Турчанинов В. П. STIFF - программа для решения жестких 
 * систем ОДУ. Препринт Института прикладной математики АН СССР. М., 1977 г.
 */
class StifflSolver : public Stiffl
{
public:
    /**
     * Конструктор класса.
     *
	 * @param NYDIM_PAR				количество ОДУ 1-го порядка
	 * @param values				начальные значения искомых переменных
	 * @param t_begin				начальное значение текущего времени, с
	 * @param t_end					конечное момент времени, до которого 
	 * нужно интегрировать систему ОДУ, с
	 * @param t_step_begin			начальный шаг интегрирования, с
     */
	StifflSolver(int NYDIM_PAR, 
		double *values,
		double t_begin,
		double t_end,
		double t_step_begin);
    /**
     * Деструктор класса.
     */
    ~StifflSolver();
    /**
    * Производит интегрирование системы ОДУ.
    *
	* @param mix       смесь, в которой происходят химические реакции
    * @param aFullTime временной интервал, на котором 
    * производится интегрирование
    */
    void performIntegration(Mixture &mix, RealType aFullTime);
    /**
     * Смесь газов.
     */
    Mixture *mixture;

private:
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
	RealType kf;
	RealType kr;
	/**
	 * Переопределяет IFNSH() из Stiffl. Служит для управления выводом 
	 * решения в процессе интегрирования. 
	 *
	 * @return 0, в случае, если нужно продолжить интегрирование, 
	 * 1, если по каким-либо причинам необходимо интегрирование прекратить.
	 */
	virtual int IFNSH();
	/**
	 * Переопределяет DIFFUN() из Stiffl. Вычисляет значения правых частей 
	 * систему ОДУ по заданным значениям YY[0], и записывает их в массив F.
	 *
	 * @param YY двумерный массив, YY[0] содержит текущие вычисленные значения
	 * искомых переменных
	 * @param F  массив, в который записываются значения правых частей 
	 * системы ОДУ
	 * @return 0
	 */
	virtual int DIFFUN(double **YY, double *F);
	/**
	 * Переопределяет PEDERV() из Stiffl. Вычисляет значения якобиана 
	 * интегрируемой системы ОДУ. Вычисление якобиана производится 
	 * по аналитическим формулам.
	 */
	virtual void PEDERV();
};

#endif // STIFFL_SOLVER_H