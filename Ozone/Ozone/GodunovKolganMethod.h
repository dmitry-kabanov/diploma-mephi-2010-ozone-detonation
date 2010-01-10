/**
* @file
*
* @author  Кабанов Дмитрий <kabanovdmitry@gmail.com>
* @version $Id$
*
* @section DESCRIPTION
*
* Объявляет класс GodunovKolganMethod.
*/
#ifndef GODUNOVKOLGANMETHOD_H
#define GODUNOVKOLGANMETHOD_H

#include <vector>
#include "Config.h"
#include "Piston.h"
#include "RealType.h"

class StifflSolver;
class Output;

/**
 * Представляет собой решатель системы уравнений газовой динамики
 * методом Годунова-Колгана.
 */
class GodunovKolganMethod
{
public:
	/**
	 * Конструктор класса.
	 *
	 * @param config конфигурация программы
	 */
	GodunovKolganMethod(const Config &config);
	/**
	 * Деструктор класса. Освобождает память, выделенную с помощью
	 * оператора new.
	 */
	~GodunovKolganMethod();
	/**
	 * Интегрирует систему уравнений газовой динамики
	 * совместно с уравнениями химической кинетики.
	 */
	void run();

	int *getCellNumbers() { return &cells_numbers[0]; }
	RealType *getX() { return &x[0]; }
	RealType *getXCenter() { return &x_center[0]; }
	RealType *getP() { return &p[0]; }
	RealType *getU() { return &u[0]; }
	RealType *getRho() { return &rho[0]; }
	RealType *getFullEnergy() { return &e[0]; }
	RealType *getIntEnergy() { return &u_energy[0]; }
	RealType **getMoleFractions() { return volumeFractions; }
	bool *getShockWaveFront() { return shock_wave_front; }
	int getMeshSize() const { return config_->getMeshSize(); }

private:
	/**
	 * Инициализирует все газодинамические величины, 
	 * а также объект газовой смеси.
	 */
	void init_();

	/**
	 * Задает всем векторным величинам размер,
	 * соответствующий размеру разностной сетки.
	 */
	void resizeAllVectors();

	/**
	 * Вычисляет приращение величины f на границах ячейки.
	 * Так как в методе Годунова—Колгана используется линейная интерполяция
	 * газодинамических величин внутри ячейки, то необходимо найти приращение
	 * величин на границах ячейки по сравнение со значением в центре ячейки.
	 *
	 * @param f_left			значение f в ячейке слева от текущей
	 * @param f					значение f в текущей ячейке
	 * @param f_right			значение f в ячейке справа от текущей
	 * @param x_center_left		координата центра ячейки слева от текущей
	 * @param x_center			координата центра текущей ячейки
	 * @param x_center_right	координата центра ячейки справа от текущей
	 * @param x_bound_l			координата левой границы текущей ячейки
	 * @param x_bound_r			координата правой границы текущей ячейки
	 * @return					приращение величины f
	 */
	RealType calc_delta(RealType f_left,
		RealType f,
		RealType f_right,
		RealType x_center_left,
		RealType x_center,
		RealType x_center_right,
		RealType x_bound_l,
		RealType x_bound_r);

	/**
	 * Модифицирует ячейки слева и справа от фронта выделяемой ударной волны.
	 * Если ячейка слева от фронта становится много больше ячейки 
	 * справа от фронта, то ячейка слева делится пополам, а ячейка справа
	 * объединяется со следующей за ней ячейкой.
	 */
	void modifyShockWaveFront_();

	const Config *config_;
	// Массив номеров ячеек.
	std::vector<int> cells_numbers;

	RealType shock_wave_velocity;
	RealType delta_mass;
	Piston *piston;
	std::vector<RealType> x;
	std::vector<RealType> x_center;
	std::vector<RealType> m;
	std::vector<RealType> rho;
	std::vector<RealType> u;
	std::vector<RealType> p;
	RealType **volumeFractions;
	
	std::vector<RealType> internal_energy;
	std::vector<RealType> u_energy;
	std::vector<RealType> e;
	std::vector<RealType> deltaTemperature;
	std::vector<RealType> gamma;
	std::vector<RealType> p_contact;
	std::vector<RealType> u_contact;
	std::vector<RealType> delta_impulse;
	std::vector<RealType> delta_energy;
	
	std::vector<RealType> rho_bound_r;
	std::vector<RealType> rho_bound_l;
	std::vector<RealType> u_bound_r;
	std::vector<RealType> u_bound_l;
	std::vector<RealType> e_bound_r;
	std::vector<RealType> e_bound_l;
	std::vector<RealType> p_bound_r;
	std::vector<RealType> p_bound_l;
	std::vector<RealType> rho_u_bound_r;
	std::vector<RealType> rho_u_bound_l;
	std::vector<RealType> rho_e_bound_r;
	std::vector<RealType> rho_e_bound_l;

	std::vector<RealType> rho_delta;
	std::vector<RealType> rho_u_delta;
	std::vector<RealType> rho_e_delta;

	std::vector<RealType> rho_u;
	std::vector<RealType> rho_e;
	
	// Массив хранит данные о фронте ударной волны.
	// Если значение элемента равно true - на границе этой ячейки 
	// находится фронт ударной волны.
	// Если значение элемента равно false - то нет.
	bool *shock_wave_front;

	StifflSolver *kinetics;
	Output *plotter_;
	int frontCellNumber_;
	/**
	 * Шаг по времени, с.
	 */
	RealType dt;
	/**
	 * Размер разностной сетки.
	 */
	int meshSize_;

};

#endif // GODUNOVKOLGANMETHOD_H
