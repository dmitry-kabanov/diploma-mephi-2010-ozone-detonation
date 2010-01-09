/**
* @file
*
* @author  Кабанов Дмитрий <kabanovdmitry@gmail.com>
* @version $Id$
*
* @section DESCRIPTION
*
* Реализует функции класса GodunovKolganMethod.
*/
#include "GodunovKolganMethod.h"

#include <iostream>
#include "constants.h"
#include "func_riemann_solver.h"
#include "Output.h"
#include "StifflSolver.h"
using namespace std;

GodunovKolganMethod::GodunovKolganMethod(const Config &config)
{
	config_ = &config;

	int nSpecies = 3;
	double vf1[4];
	vf1[0] = 0.0;
	vf1[1] = 0.0;
	vf1[2] = 0.0;
	vf1[3] = 0.0;
	kinetics = new StifflSolver(
		nSpecies + 1,
		vf1,
		0.0,
		0.0,
		1.0e-13,
		"Substances.txt",
		"Reactions.txt",
		"MoleFractions.txt");

	piston = new Piston(config_->getPInitial(),
		config_->getRhoInitial(),
		"Piston.txt");
	plotter_ = new Output(*kinetics->getMixture(), "csv", "Output");

	volumeFractions  = new RealType*[config_->getMeshSize()+1];
	for (int i = 0; i <= config_->getMeshSize(); i++) {
		volumeFractions[i] = new RealType[kinetics->getMixture()->nSubstances];
	}

	shock_wave_front = new bool[config_->getMeshSize()+1];

	init_();
}

GodunovKolganMethod::~GodunovKolganMethod()
{
	delete kinetics;
	delete piston;
	delete plotter_;
	for (int i = 0; i <= config_->getMeshSize(); i++) {
		delete [] volumeFractions[i];
	}
	delete [] volumeFractions;
	delete [] shock_wave_front;
}

void GodunovKolganMethod::init_()
{
	int i;

	resizeAllVectors();

	for (i = 0; i <= config_->getMeshSize(); i++) {
		shock_wave_front[i] = false;
	}

	for (i = 0; i <= N; i++) {
		cells_numbers[i] = i;
	}

	x[0]        = 0;
	x_center[0] = 0;
	m[0]        = 0;
	rho[0]      = 0;
	u[0]        = 0;
	e[0]        = 0;
	p[0]        = 0;
	rho_u[0]    = 0;

	for (int i = 1; i <= N; i++) {
		if (i <= L) {
			x[i] = x[i-1] + DX1;
		}
		else {
			x[i] = x[i-1] + DX2;
		}

		if (i <= L) {
			p[i]   = P1;
			rho[i] = RHO1;
			u[i]   = U1;
		}
		else {
			p[i]   = P2;
			rho[i] = RHO2;
			u[i]   = U2;
		}

		m[i] = (x[i] - x[i-1]) * rho[i];

		x_center[i] = (x[i-1] + x[i]) / 2.0;
		rho_u[i] = rho[i] * u[i];
	}

	shock_wave_front[L]   = true;
	shock_wave_front[L+1] = true;

	for (int i = 0; i <= N; i++) {
		if (i <= L) {
			gamma[i] = GAMMA_BEHIND_FRONT;
		} else {
			gamma[i] = GAMMA_AHEAD_FRONT;
		}
	}

	for (int i = 0; i <= L; i++) {
		volumeFractions[i][0] = 0.0;
		volumeFractions[i][1] = 0.0;
		volumeFractions[i][2] = 100.0;
	}
	for (int i = L+1; i <= N; i++) {
		volumeFractions[i][0] = 0.0;
		volumeFractions[i][1] = 0.0;
		volumeFractions[i][2] = 100.0;
	}

	internal_energy[0] = 0.0;
	e[0]          = 0.0;
	u_energy[0]        = 0.0;

	// Температура в ячейке, К.
	RealType t;
	// Молекулярный вес смеси, кг кмоль-1.
	RealType mw;

	for (int i = 1; i <= N; i++) {
		if (i <= L) {
			internal_energy[i] = p[i] / ((GAMMA_BEHIND_FRONT - 1) * rho[i]);
		}
		else {
			internal_energy[i] = p[i] / ((GAMMA_AHEAD_FRONT - 1) * rho[i]);
		}

		kinetics->getMixture()->setVolumeFractions(volumeFractions[i]);
		kinetics->getMixture()->calculateInitialMolecularWeight();
		mw = kinetics->getMixture()->getMolecularWeight();
		t = p[i] * mw / (Mixture::R_J_OVER_KMOL_K * rho[i]);
		kinetics->getMixture()->setConcentrations(t, rho[i]);
		u_energy[i] = kinetics->getMixture()->calculateMixtureEnthalpy(t) -
			p[i] / rho[i];

		e[i] = u_energy[i] + u[i] * u[i] / 2.0;
	}

	rho_e[0] = 0.0;

	for (int i = 1; i <= N; i++) {
		rho_e[i] = rho[i] * internal_energy[i];
	}

	for (int i = 0; i <= N; i++) {
		p_contact[i] = p[i];
		u_contact[i] = u[i];
		delta_impulse[i] = 0;
		delta_energy[i] = 0;
	}

	for (int i = 0; i <= N; i++) {
		rho_bound_r[i] = 0;
		rho_bound_l[i] = 0;
		rho_u_bound_r[i] = 0;
		rho_u_bound_l[i] = 0;
		rho_e_bound_r[i] = 0;
		rho_e_bound_l[i] = 0;
		u_bound_r[i] = 0;
		u_bound_l[i] = 0;
		e_bound_r[i] = 0;
		e_bound_l[i] = 0;
		p_bound_r[i] = 0;
		p_bound_l[i] = 0;
	}

	for (int i = 0; i <= N; i++) {
		rho_tg_left[i]  = 0;
		rho_tg_right[i] = 0;
		rho_tg[i]       = 0;

		rho_u_tg_left[i]  = 0;
		rho_u_tg_right[i] = 0;
		rho_u_tg[i]       = 0;

		rho_e_tg_left[i]  = 0;
		rho_e_tg_right[i] = 0;
		rho_e_tg[i]       = 0;
	}

	// Пишем начальные значения в файл.
	if (! config_->getResume()) {
		plotter_->plotData(0, *this);
	}
}

void GodunovKolganMethod::resizeAllVectors()
{
	int n = config_->getMeshSize() + 1;

	cells_numbers.resize(n);
	x.resize(n);
	x_center.resize(n);
	m.resize(n);
	rho.resize(n);
	u.resize(n);
	p.resize(n);

	internal_energy.resize(n);
	u_energy.resize(n);
	e.resize(n);
	temperature.resize(n);
	gamma.resize(n);
	p_contact.resize(n);
	u_contact.resize(n);
	delta_impulse.resize(n);
	delta_energy.resize(n);

	rho_bound_r.resize(n);
	rho_bound_l.resize(n);
	u_bound_r.resize(n);
	u_bound_l.resize(n);
	e_bound_r.resize(n);
	e_bound_l.resize(n);
	p_bound_r.resize(n);
	p_bound_l.resize(n);
	rho_u_bound_r.resize(n);
	rho_u_bound_l.resize(n);
	rho_e_bound_r.resize(n);
	rho_e_bound_l.resize(n);

	rho_tg_left.resize(n);
	rho_tg_right.resize(n);
	rho_tg.resize(n);
	rho_u_tg_left.resize(n);
	rho_u_tg_right.resize(n);
	rho_u_tg.resize(n);
	rho_e_tg_left.resize(n);
	rho_e_tg_right.resize(n);
	rho_e_tg.resize(n);

	rho_u.resize(n);
	rho_e.resize(n);
}

void GodunovKolganMethod::run()
{
	int i;

	for (int j = config_->getStart() + 1; 
		j <= config_->getStart() + config_->getTimeSteps(); j++) {
		for (i = 2; i < N; i++) {
			// Считаем тангенс для rho на левой и правой границах ячейки.
			rho_tg_left[i] = (rho[i] - rho[i-1]) / (x_center[i] - x_center[i-1]);
			rho_tg_right[i] = (rho[i+1] - rho[i]) / (x_center[i+1] - x_center[i]);
			if (abs(rho_tg_left[i]) > abs(rho_tg_right[i])) {
				rho_tg[i] = rho_tg_right[i];
			}
			else {
				rho_tg[i] = rho_tg_left[i];
			}

			// Считаем тангенс для rho*u на левой и правой границах ячейки.
			rho_u_tg_left[i]  = (rho_u[i] - rho_u[i-1]) / 
				(x_center[i] - x_center[i-1]);
			rho_u_tg_right[i] = (rho_u[i+1] - rho_u[i]) / 
				(x_center[i+1] - x_center[i]);
			if (abs(rho_u_tg_left[i]) > abs(rho_u_tg_right[i])) {
				rho_u_tg[i] = rho_u_tg_right[i];
			}
			else {
				rho_u_tg[i] = rho_u_tg_left[i];
			}

			// Считаем тангенс для rho*e на левой и правой границах ячейки.
			rho_e_tg_left[i]  = (rho_e[i] - rho_e[i-1]) / (x_center[i] - x_center[i-1]);
			rho_e_tg_right[i] = (rho_e[i+1] - rho_e[i]) / (x_center[i+1] - x_center[i]);
			if (abs(rho_e_tg_left[i]) > abs(rho_e_tg_right[i])) {
				rho_e_tg[i] = rho_e_tg_right[i];
			}
			else {
				rho_e_tg[i] = rho_e_tg_left[i];
			}

			if (shock_wave_front[i] == true && shock_wave_front[i+1] == false)
			{
				break;
			}
		}

		rho_tg_left[1]    = rho_tg_left[2];
		rho_tg_right[1]   = rho_tg_right[2];
		rho_tg[1]         = rho_tg[2];

		rho_u_tg_left[1]  = rho_u_tg_left[2];
		rho_u_tg_right[1] = rho_u_tg_right[2];
		rho_u_tg[1]       = rho_u_tg[2];

		rho_e_tg_left[1]  = rho_e_tg_left[2];
		rho_e_tg_right[1] = rho_e_tg_right[2];
		rho_e_tg[1]       = rho_e_tg[2];

		for (i = 1; i < N; i++) {
			rho_bound_r[i] = rho[i] + (x[i] - x[i-1]) * rho_tg[i] / 2.0;
			rho_bound_l[i] = rho[i] - (x[i] - x[i-1]) * rho_tg[i] / 2.0;            

			rho_u_bound_r[i] = rho_u[i] + (x[i] - x[i-1]) * rho_u_tg[i] / 2.0;
			rho_u_bound_l[i] = rho_u[i] - (x[i] - x[i-1]) * rho_u_tg[i] / 2.0;

			rho_e_bound_r[i] = rho_e[i] + (x[i] - x[i-1]) * rho_e_tg[i] / 2.0;
			rho_e_bound_l[i] = rho_e[i] - (x[i] - x[i-1]) * rho_e_tg[i] / 2.0;

			// Считаем u на левой и правой границах ячейки.
			u_bound_r[i] = rho_u_bound_r[i] / rho_bound_r[i];
			u_bound_l[i] = rho_u_bound_l[i] / rho_bound_l[i];

			// Считаем e на левой и правой границах ячейки.
			e_bound_r[i] = rho_e_bound_r[i] / rho_bound_r[i];
			e_bound_l[i] = rho_e_bound_l[i] / rho_bound_l[i];

			// Считаем p на левой и правой границах ячейки.
			p_bound_r[i] = rho_bound_r[i] * (gamma[i] - 1) * e_bound_r[i];
			p_bound_l[i] = rho_bound_l[i] * (gamma[i] - 1) * e_bound_l[i];

			if (shock_wave_front[i] == true && shock_wave_front[i+1] == false)
			{
				break;
			}
		}

		// Решаем задачу Римана.
		for (i = 1; i < N; i++) {
			/*
			* Фронт ударной волны, которую мы выделяем, 
			* находится на правой границе i-й ячейки.
			* Задача Римана решается без линейной интерполяции
			* параметров в ячейках i и i+1.
			* После решения задачи Римана между этими ячейками
			* выходим из цикла, так как перед ударной волной ничего
			* не изменилось.
			*/
			if (shock_wave_front[i] == true && shock_wave_front[i+1] == true) {
				p_contact[i] = calc_p_contact(p[i], p[i+1],
					rho[i], rho[i+1],
					u[i], u[i+1], GAMMA_FRONT
					);
				u_contact[i] = calc_u_contact(p_contact[i],
					p[i], p[i+1],
					rho[i], rho[i+1],
					u[i], u[i+1], GAMMA_FRONT
					);
				shock_wave_velocity = calc_shock_wave_velocity(
					p_contact[i], u_contact[i],
					p[i], p[i+1],
					rho[i], rho[i+1],
					u[i], u[i+1], GAMMA_FRONT
					);
				if (shock_wave_velocity == 0) {
					cout << "Shock wave velocity is equals to 0." << endl;
					exit(1);
				}
				break;
			}
			/*
			* Фронт ударной волны, которую мы выделяем, 
			* находится на правой границе i+1-й ячейки.
			*/
			else if (shock_wave_front[i] == false && shock_wave_front[i+1] == true) {
				p_contact[i] = calc_p_contact(p_bound_r[i], p[i+1],
					rho_bound_r[i], rho[i+1],
					u_bound_r[i], u[i+1], GAMMA_BEHIND_FRONT
					);
				u_contact[i] = calc_u_contact(p_contact[i],
					p_bound_r[i], p[i+1],
					rho_bound_r[i], rho[i+1],
					u_bound_r[i], u[i+1], GAMMA_BEHIND_FRONT
					);
			}
			/*
			* Решение задачи Римана для соседних ячеек, ни одна из которых
			* не содержит ударную волну, которую мы выделяем.
			* Задача Римана решается через линейно интерполированные
			* параметры в ячейках (как и положено в методе Годунова-Колгана).
			*/
			else {
				p_contact[i] = calc_p_contact(p_bound_r[i], p_bound_l[i+1],
					rho_bound_r[i], rho_bound_l[i+1],
					u_bound_r[i], u_bound_l[i+1], GAMMA_BEHIND_FRONT);
				u_contact[i] = calc_u_contact(p_contact[i], p_bound_r[i], p_bound_l[i+1],
					rho_bound_r[i], rho_bound_l[i+1],
					u_bound_r[i], u_bound_l[i+1], GAMMA_BEHIND_FRONT);
			}
		}

		/*
		* На левой границе у нас находится поршень.
		* Считаем параметры для поршня.
		*/
		// Левая граница области расчета движется со скоростью поршня.
		p_contact[0] = p[1];
		u_contact[0] = piston->calculateVelocity(volumeFractions[1][2]);

		/**
		* Считаем импульс и энергию в ячейке.
		*/
		for (i = 1; i < N; i++) {
			/*
			* Фронт ударной волны, которую мы выделяем, 
			* находится на правой границе i-й ячейки.
			*/
			if (shock_wave_front[i] == true && shock_wave_front[i+1] == true) {
				delta_mass = (shock_wave_velocity - u[i+1]) * rho[i+1] * DT;

				delta_impulse[i] = p_contact[i-1] * DT + 
					delta_mass * u[i+1] - p[i+1] * DT;
				delta_energy[i]  = p_contact[i-1] * u_contact[i-1] * DT + 
					delta_mass * e[i+1] - p[i+1] * u[i+1] * DT;
			}
			/*
			* Фронт ударной волны, которую мы выделяем, 
			* находится на левой границе i-й ячейки.
			* После расчета delta_impulse[i] и delta_energy[i]
			* выходим из цикла, так как перед ударной волной 
			* ничего не изменилось.
			*/
			else if (shock_wave_front[i] == true 
				&& shock_wave_front[i+1] == false) {

					delta_impulse[i] = - delta_mass * u[i] + 
						p[i] * DT - 
						p_contact[i] * DT;
					delta_energy[i]  = - delta_mass * e[i] + 
						+ p[i] * u[i] * DT -
						p_contact[i] * u_contact[i] * DT;
					break;
			}
			/*
			* i-я ячейка не содержит ударную волну,
			* которую мы выделяем.
			* Считаем импульс и энергию обычным способом.
			*/
			else {
				delta_impulse[i] = p_contact[i-1] * DT -
					p_contact[i] * DT;
				delta_energy[i]  = p_contact[i-1] * u_contact[i-1] * DT - 
					p_contact[i] * u_contact[i] * DT;
			}
		}

		// Левая граница области расчета движется со скоростью поршня.
		x[0] = x[0] + u_contact[0] * DT;

		// Вычисляем основные величины.
		for (i = 1; i <= N; i++) {
			if (shock_wave_front[i] == true && shock_wave_front[i+1] == true) {
				// Фронт ударной волны, которую мы выделяем, 
				// находится на правой границе i-й ячейки.
				// Считаем изменение массы в i-й ячейке.
				// Правая граница i-й ячейки движется со скоростью ударной волны.
				x[i] = x[i] + shock_wave_velocity * DT;
				u[i] = (m[i] * u[i] + delta_impulse[i]) / (m[i] + delta_mass);
				e[i] = (m[i] * e[i] + delta_energy[i]) / (m[i] + delta_mass);
				m[i] = m[i] + delta_mass;
			}
			else if (shock_wave_front[i] == true && 
				shock_wave_front[i+1] == false) {
					// i-я ячейка содержит выделяемую 
					// ударную волну на левой границе.
					// Её правая граница движется со скоростью контактного разрыва.
					x[i] = x[i] + u_contact[i] * DT;
					x_center[i] = (x[i-1] + x[i]) / 2.0;
					u[i] = (m[i] * u[i] + delta_impulse[i]) / (m[i] - delta_mass);
					e[i] = (m[i] * e[i] + delta_energy[i]) / (m[i] - delta_mass);
					m[i] = m[i] - delta_mass;
					// cout << shock_wave_velocity << endl;
					break;
			}
			else {
				// Ячейка не содержит выделяемую ударную волну.
				x[i] = x[i] + u_contact[i] * DT;
				u[i] = (m[i] * u[i] + delta_impulse[i]) / m[i];
				e[i] = (m[i] * e[i] + delta_energy[i]) / m[i];
			}

			x_center[i] = (x[i-1] + x[i]) / 2.0;
			rho[i] = m[i] / (x[i] - x[i-1]);
			u_energy[i] = e[i] - u[i] * u[i] / 2.0;
			kinetics->getMixture()->setStateWithURhoX(u_energy[i], rho[i], volumeFractions[i]);
			kinetics->performIntegration(DT);
			kinetics->updateMoleFractions(volumeFractions[i]);
			p[i] = kinetics->getPressure();
			rho_u[i] = rho[i] * u[i];
			gamma[i] = p[i] / (rho[i] * u_energy[i]) + 1;
			rho_e[i] = rho[i] * p[i] / ((gamma[i] - 1) * rho[i]);
		}

		if (j % TIMEDIVISOR == 0) {
			cout << "j = " << j << endl;
			cout << "D = " << shock_wave_velocity << endl << endl;
			plotter_->plotData(j, *this);
		}

		// Модифицируем ячейки, связанные с фронтом ударной волны.
		for (i = 1; i < (N-1); i++) {
			if (shock_wave_front[i] == true && shock_wave_front[i+1] == true) {
				if ((x[i+1] - x[i]) <= K * (x[i] - x[i-1])) {
					RealType delta_x_right = x[i+2] - x[i+1];
					RealType delta_x_left  = x[i+1] - x[i];
					RealType delta_x_sum = delta_x_left + delta_x_right;
					RealType temp_rho = rho[i+2];
					RealType temp_x   = x[i];

					rho[i+2] = 1 / delta_x_sum * 
						(rho[i+2] * delta_x_right + 
						rho[i+1] * delta_x_left);

					m[i+2] = rho[i+2] * delta_x_sum;

					u[i+2] = 1 / (rho[i+2] * delta_x_sum) * 
						(temp_rho * delta_x_right * u[i+2] + 
						rho[i+1] * delta_x_left * u[i+1]);

					e[i+2] = 1 / (rho[i+2] * delta_x_sum) * 
						(temp_rho * delta_x_right * e[i+2] + 
						rho[i+1] * delta_x_left * e[i+1]);

					// Делим ячейку слева от фронта ударной волны пополам.
					p[i+1]   = p[i];
					rho[i+1] = rho[i];
					u[i+1]   = u[i];
					e[i+1]   = e[i];
					m[i+1]   = m[i] / 2.0;
					m[i]     = m[i] / 2.0;

					// Модифицируем координаты границ и центров ячеек.
					x[i]   = (x[i] + x[i-1]) / 2.0;
					x[i+1] = temp_x;

					//x_center[i]   = (x[i-1] + x[i])   / 2.0;
					//x_center[i+1] = (x[i]   + x[i+1]) / 2.0;
					//x_center[i+2] = (x[i+1] + x[i+2]) / 2.0;

					gamma[i+1] = gamma[i];

					shock_wave_front[i]   = 0;
					shock_wave_front[i+1] = true;
					shock_wave_front[i+2] = true;
				}
				break;
			}
		}
	}
}