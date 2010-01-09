/**
* @file
*
* @author  Кабанов Дмитрий <kabanovdmitry@gmail.com>
* @version $Id$
*
* @section DESCRIPTION
*
* Объявление класса Piston.
*/
#ifndef PISTON_H
#define PISTON_H

#include "RealType.h"

/**
 * Класс поршня. Используется модель стационарного поршня.
 * См. Мейдер Ч. Численное моделирование детонации. Глава 1.
 */
class Piston
{
public:
	/**
	 * Конструктор класса Piston.
	 *
	 * @param pInitial			начальное давление, Па
	 * @param rhoInitial		начальная плотность, кг м-3
	 * @param fileOfPistonData	имя файла с табличными данными для поршня
	 */
	Piston(RealType pInitial, RealType rhoInitial, const char *fileOfPistonData);
	/**
	 * Деструктор класса Piston. Освобождает выделенную память.
	 */
    ~Piston();
    /**
     * Вычисляет скорость поршня в зависимости от заданной
	 * мольной доли исходного непрореагировавшего вещества.
     *
     * @param f мольная доля исходного вещества
     * @return новое значение скорости поршня
     */
    RealType calculateVelocity(RealType f);
private:
    /**
     * Читает из файла с именем fileOfPistonData табличные данные
	 * давления и плотности в зависимости от мольной доли
	 * исходного непрореагировавшего вещества.
     *
     * @param fileOfPistonData имя файла с табличными данными для поршня
     */
    void readFileOfPiston(const char *fileOfPistonData);
	/**
	 * Количество строк в файле с табличными данными для поршня.
	 */
    int nRows_;
	/**
	 * Массив давлений.
	 */
    RealType *pressures_;
	/**
	 * Массив плотностей.
	 */
    RealType *densities_;
	/**
	 * Массив мольных долей.
	 */
    RealType *fractions_;
	/**
	 * Начальное давление в смеси, Па.
	 */
	RealType pInitial_;
	/**
	 * Начальная плотность смеси, кг м-3.
	 */
	RealType rhoInitial_;
};

#endif // PISTON_H