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

#include <cmath>
#include <iomanip>
#include <iostream>
#include "constants.h"
#include "func_riemann_solver.h"
#include "Output.h"
#include "StifflSolver.h"
using namespace std;

GodunovKolganMethod::GodunovKolganMethod(const Config &config)
{
	config_ = &config;

	mixture = new Mixture(config_->getFileOfSubstances().c_str(),
		config_->getFileOfReactions().c_str(),
		config_->getFileOfFractions().c_str());

	double *vf1 = new double [mixture->getNSpecies()];
	for (int i = 0; i < mixture->getNSpecies(); ++i) {
		vf1[i] = 0.0;
	}
	kinetics = new StifflSolver(mixture->getNSpecies() + 1,
		vf1,
		0.0,
		0.0,
		1.0e-13);
	kinetics->Seteps(1e-2);

	piston = new Piston(config_->getPInitial(),
		config_->getRhoInitial(),
		config_->getFileOfPiston().c_str());
	plotter_ = new Output(*mixture, "csv", "Output");

	volumeFractions  = new RealType*[config_->getMeshSize()+1];
	for (int i = 0; i <= config_->getMeshSize(); i++) {
		volumeFractions[i] = new RealType[mixture->nSubstances];
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

	dt = config_->getDt();
	meshSize_ = config_->getMeshSize();

	resizeAllVectors();

	for (i = 0; i <= config_->getMeshSize(); i++) {
		shock_wave_front[i] = false;
	}

	for (i = 0; i <= meshSize_; i++) {
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

	for (int i = 1; i <= meshSize_; i++) {
		x[i] = x[i-1] + config_->getDx();

		if (i <= config_->getInitialShockWaveSize()) {
			p[i]   = config_->getPFront();
			rho[i] = config_->getRhoFront();
			u[i]   = config_->getUFront();
		}
		else {
			p[i]   = config_->getPInitial();
			rho[i] = config_->getRhoInitial();
			u[i]   = config_->getUInitial();
		}

		m[i] = (x[i] - x[i-1]) * rho[i];

		x_center[i] = (x[i-1] + x[i]) / 2.0;
		rho_u[i] = rho[i] * u[i];
	}

	shock_wave_front[config_->getInitialShockWaveSize()]   = true;
	shock_wave_front[config_->getInitialShockWaveSize()+1] = true;

	for (int i = 0; i <= meshSize_; i++) {
		if (i <= config_->getInitialShockWaveSize()) {
			gamma[i] = config_->getGammaBehindFront();
		} else {
			gamma[i] = config_->getGammaAheadFront();
		}
	}

	for (int i = 0; i <= meshSize_; i++) {
		mixture->returnMoleFractions(volumeFractions[i]);
	}

	internal_energy[0] = 0.0;
	e[0]          = 0.0;
	u_energy[0]        = 0.0;

	// Температура в ячейке, К.
	RealType t;
	// Молекулярный вес смеси, кг кмоль-1.
	RealType mw;

	for (int i = 1; i <= meshSize_; i++) {
		if (i <= config_->getInitialShockWaveSize()) {
			internal_energy[i] = p[i] / ((config_->getGammaBehindFront() - 1) * rho[i]);
		}
		else {
			internal_energy[i] = p[i] / ((config_->getGammaAheadFront() - 1) * rho[i]);
		}

		mixture->setVolumeFractions(volumeFractions[i]);
		mixture->calculateInitialMolecularWeight();
		mw = mixture->getMolecularWeight();
		t = p[i] * mw / (Mixture::R_J_OVER_KMOL_K * rho[i]);
		mixture->setConcentrations(t, rho[i]);
		u_energy[i] = mixture->calculateMixtureEnthalpy(t) -
			p[i] / rho[i];

		e[i] = u_energy[i] + u[i] * u[i] / 2.0;
	}

	rho_e[0] = 0.0;

	for (int i = 1; i <= meshSize_; i++) {
		rho_e[i] = rho[i] * internal_energy[i];
	}

	for (int i = 0; i <= meshSize_; i++) {
		p_contact[i] = p[i];
		u_contact[i] = u[i];
		delta_impulse[i] = 0;
		delta_energy[i] = 0;
	}

	for (int i = 0; i <= meshSize_; i++) {
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

	for (int i = 0; i <= meshSize_; i++) {
		rho_delta[i]   = 0;
		rho_u_delta[i] = 0;
		rho_e_delta[i] = 0;
	}

	if (config_->getStart() == 0)	{
		// Пишем начальные значения в файл.
		plotter_->plotData(0, *this);
		cout << "Starting from the very beginning." << endl;
	} 
	else
	{
		update_();
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
	deltaTemperature.resize(n);
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

	rho_delta.resize(n);
	rho_u_delta.resize(n);
	rho_e_delta.resize(n);

	rho_u.resize(n);
	rho_e.resize(n);
}

void GodunovKolganMethod::run()
{
	int i;

	ofstream outD("velocity.txt");
	if (!outD) {
		cout << "Can't create file for shock wave velocity" << endl;
	}

	for (int j = config_->getStart() + 1; j <= config_->getTimeSteps(); j++) {
		for (i = 2; i < meshSize_; i++) {
			// Считаем тангенс для rho на левой и правой границах ячейки.
			rho_delta[i] = calc_delta(rho[i-1], rho[i], rho[i+1], 
				x_center[i-1], x_center[i], x_center[i+1],
				x[i-1], x[i]);

			// Считаем тангенс для rho*u на левой и правой границах ячейки.
			rho_u_delta[i] = calc_delta(rho_u[i-1], rho_u[i], rho_u[i+1],
				x_center[i-1], x_center[i], x_center[i+1],
				x[i-1], x[i]);

			// Считаем тангенс для rho*e на левой и правой границах ячейки.
			rho_e_delta[i] = calc_delta(rho_e[i-1], rho_e[i], rho_e[i+1],
				x_center[i-1], x_center[i], x_center[i+1],
				x[i-1], x[i]);

			rho_bound_r[i] = rho[i] + rho_delta[i];
			rho_bound_l[i] = rho[i] - rho_delta[i];

			rho_u_bound_r[i] = rho_u[i] + rho_u_delta[i];
			rho_u_bound_l[i] = rho_u[i] - rho_u_delta[i];

			rho_e_bound_r[i] = rho_e[i] + rho_e_delta[i];
			rho_e_bound_l[i] = rho_e[i] - rho_e_delta[i];

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

		// В граничной левой ячейке вместо линейного распределения значений 
		// параметров по ячейке, как должно быть в методе Годунова-Колгана,
		// используем постоянное распределение, как в методе Годунова.
		rho_bound_r[1] = rho[1];
		rho_bound_l[1] = rho[1];

		rho_u_bound_r[1] = rho_u[1];
		rho_u_bound_l[1] = rho_u[1];

		rho_e_bound_r[1] = rho_e[1];
		rho_e_bound_l[1] = rho_e[1];

		// Считаем u на левой и правой границах ячейки.
		u_bound_r[1] = rho_u_bound_r[1] / rho_bound_r[1];
		u_bound_l[1] = rho_u_bound_l[1] / rho_bound_l[1];

		// Считаем e на левой и правой границах ячейки.
		e_bound_r[1] = rho_e_bound_r[1] / rho_bound_r[1];
		e_bound_l[1] = rho_e_bound_l[1] / rho_bound_l[1];

		// Считаем p на левой и правой границах ячейки.
		p_bound_r[1] = rho_bound_r[1] * (gamma[1] - 1) * e_bound_r[1];
		p_bound_l[1] = rho_bound_l[1] * (gamma[1] - 1) * e_bound_l[1];

		// Решаем задачу Римана.
		for (i = 1; i < meshSize_; i++) {
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
					u[i], u[i+1], config_->getGammaInsideFront()
					);
				u_contact[i] = calc_u_contact(p_contact[i],
					p[i], p[i+1],
					rho[i], rho[i+1],
					u[i], u[i+1], config_->getGammaInsideFront()
					);
				shock_wave_velocity = calc_shock_wave_velocity(
					p_contact[i], u_contact[i],
					p[i], p[i+1],
					rho[i], rho[i+1],
					u[i], u[i+1], config_->getGammaInsideFront()
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
				p_contact[i] = calc_p_contact(p_bound_r[i], p_bound_l[i+1],
					rho_bound_r[i], rho_bound_l[i+1],
					u_bound_r[i], u_bound_l[i+1], config_->getGammaBehindFront()
				);
				u_contact[i] = calc_u_contact(p_contact[i],
					p_bound_r[i], p_bound_l[i+1],
					rho_bound_r[i], rho_bound_l[i+1],
					u_bound_r[i], u_bound_l[i+1], config_->getGammaBehindFront()
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
					u_bound_r[i], u_bound_l[i+1], config_->getGammaBehindFront());
				u_contact[i] = calc_u_contact(p_contact[i], p_bound_r[i], p_bound_l[i+1],
					rho_bound_r[i], rho_bound_l[i+1],
					u_bound_r[i], u_bound_l[i+1], config_->getGammaBehindFront());
			}
		}

		/*
		* На левой границе у нас находится поршень.
		* Считаем параметры для поршня.
		*/
		// Левая граница области расчета движется со скоростью поршня.
		p_contact[0] = p[1];
		u_contact[0] = piston->calculateVelocity(volumeFractions[1]);

		/**
		* Считаем импульс и энергию в ячейке.
		*/
		for (i = 1; i < meshSize_; i++) {
			/*
			* Фронт ударной волны, которую мы выделяем, 
			* находится на правой границе i-й ячейки.
			*/
			if (shock_wave_front[i] == true && shock_wave_front[i+1] == true) {
				delta_mass = (shock_wave_velocity - u[i+1]) * rho[i+1] * dt;

				delta_impulse[i] = p_contact[i-1] * dt + 
					delta_mass * u[i+1] - p[i+1] * dt;
				delta_energy[i]  = p_contact[i-1] * u_contact[i-1] * dt + 
					delta_mass * e[i+1] - p[i+1] * u[i+1] * dt;
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
						p[i] * dt - 
						p_contact[i] * dt;
					delta_energy[i]  = - delta_mass * e[i] + 
						+ p[i] * u[i] * dt -
						p_contact[i] * u_contact[i] * dt;
					break;
			}
			/*
			* i-я ячейка не содержит ударную волну,
			* которую мы выделяем.
			* Считаем импульс и энергию обычным способом.
			*/
			else {
				delta_impulse[i] = p_contact[i-1] * dt -
					p_contact[i] * dt;
				delta_energy[i]  = p_contact[i-1] * u_contact[i-1] * dt - 
					p_contact[i] * u_contact[i] * dt;
			}
		}

		// Левая граница области расчета движется со скоростью поршня.
		x[0] = x[0] + u_contact[0] * dt;

		// Вычисляем основные величины.
		for (i = 1; i <= meshSize_; i++) {
			if (shock_wave_front[i] == true && shock_wave_front[i+1] == true) {
				// Фронт ударной волны, которую мы выделяем, 
				// находится на правой границе i-й ячейки.
				// Считаем изменение массы в i-й ячейке.
				// Правая граница i-й ячейки движется со скоростью ударной волны.
				x[i] = x[i] + shock_wave_velocity * dt;
				u[i] = (m[i] * u[i] + delta_impulse[i]) / (m[i] + delta_mass);
				e[i] = (m[i] * e[i] + delta_energy[i]) / (m[i] + delta_mass);
				m[i] = m[i] + delta_mass;
				frontCellNumber_ = i;
			}
			else if (shock_wave_front[i] == true && 
				shock_wave_front[i+1] == false) {
					// i-я ячейка содержит выделяемую 
					// ударную волну на левой границе.
					// Её правая граница движется со скоростью контактного разрыва.
					x[i] = x[i] + u_contact[i] * dt;
					x_center[i] = (x[i-1] + x[i]) / 2.0;
					u[i] = (m[i] * u[i] + delta_impulse[i]) / (m[i] - delta_mass);
					e[i] = (m[i] * e[i] + delta_energy[i]) / (m[i] - delta_mass);
					m[i] = m[i] - delta_mass;
					break;
			}
			else {
				// Ячейка не содержит выделяемую ударную волну.
				x[i] = x[i] + u_contact[i] * dt;
				u[i] = (m[i] * u[i] + delta_impulse[i]) / m[i];
				e[i] = (m[i] * e[i] + delta_energy[i]) / m[i];
			}

			x_center[i] = (x[i-1] + x[i]) / 2.0;
			rho[i] = m[i] / (x[i] - x[i-1]);
			u_energy[i] = e[i] - u[i] * u[i] / 2.0;
			mixture->setStateWithURhoX(u_energy[i], rho[i], volumeFractions[i]);
            //if (volumeFractions[i][2] > 0.01) {
			    kinetics->performIntegration(*mixture, dt);
            //}
			deltaTemperature[i] = mixture->getTemperature() -
				mixture->getOldTemperature();
			mixture->updateMoleFractions(volumeFractions[i]);
			p[i] = mixture->getPressure();
			rho_u[i] = rho[i] * u[i];
			rho_e[i] = p[i] / (gamma[i] - 1);
		}

		// Модифицируем ячейки, связанные с фронтом ударной волны.
		modifyShockWaveFront_();

		if (frontCellNumber_ == (meshSize_ - 1)) {
			cout << "Lead shock has reached right bound of the domain." << endl;
			exit(-128);
		}
		
		if (j % 10 == 0) {
			outD << j << " " << shock_wave_velocity << endl;
		}
		if (j % config_->getTimeStepForOutput() == 0) {
			cout << setw(6) << "j" << " = " << j << endl;
			cout << setw(6) << "D" << " = " << shock_wave_velocity << endl;
			cout << setw(6) << "xFront" << " = " << x[frontCellNumber_] << endl;

			RealType maxdT = deltaTemperature[1];
			for (i = 2; i < meshSize_; i++) {
				if (abs(maxdT) < abs(deltaTemperature[i])) {
					maxdT = deltaTemperature[i];
				}
				if (shock_wave_front[i] == true) {
					break;
				}
			}
			cout << setw(6) << "Max dT" << " = " << maxdT << endl;
			cout << endl;

            //plotter_->plotData(j, *this);
			modifyMesh();
			plotter_->plotData(j, *this);
		}
	}
}

RealType GodunovKolganMethod::calc_delta(RealType f_left, 
	RealType f, 
	RealType f_right, 
	RealType x_center_left, 
	RealType x_center, 
	RealType x_center_right,
	RealType x_bound_l,
	RealType x_bound_r)
{
	RealType f_tg_left, f_tg_right, f_tg;

	f_tg_left = (f - f_left) / (x_center - x_center_left);
	f_tg_right = (f_right - f) / (x_center_right - x_center);

	if (abs(f_tg_right) > abs(f_tg_left)) {
		f_tg = f_tg_left;
	}
	else {
		f_tg = f_tg_right;
	}

	return f_tg * (x_bound_r - x_bound_l) / 2.0;
}

void GodunovKolganMethod::modifyShockWaveFront_()
{
	int i = frontCellNumber_;

	if ((x[i+1] - x[i]) <= config_->getCellWidthCoeff() * (x[i] - x[i-1])) {
		//cout << "Modifying shock wave front." << endl;
		RealType delta_x_right = x[i+2] - x[i+1];
		RealType delta_x_left  = x[i+1] - x[i];
		RealType delta_x_sum = delta_x_left + delta_x_right;
		RealType m_left = rho[i+1] * delta_x_left;
		RealType m_right = rho[i+2] * delta_x_right;
		RealType new_rho;
		RealType temp_x   = x[i];

		new_rho = (m_left + m_right) / delta_x_sum;

		m[i+2] = new_rho * delta_x_sum;
		u[i+2] = (m_left * u[i+1] + m_right * u[i+2]) / m[i+2];
		e[i+2] = (m_left * e[i+1] + m_right * e[i+2]) / m[i+2];
		u_energy[i+2] = e[i+2] - u[i+2] * u[i+2] / 2.0;
		rho[i+2] = new_rho;

		mixture->setStateWithURhoX(u_energy[i+2], 
			rho[i+2], 
			volumeFractions[i+2]);
		p[i+2] = mixture->getPressure();

		// Делим ячейку слева от фронта ударной волны пополам.
		p[i+1]   = p[i];
		rho[i+1] = rho[i];
		u[i+1]   = u[i];
		e[i+1]   = e[i];
		rho_u[i+1] = rho_u[i];
		rho_e[i+1] = rho_e[i];

		// Модифицируем координаты границ и центров ячеек.
		x[i]   = (x[i] + x[i-1]) / 2.0;
		x[i+1] = temp_x;

		x_center[i]   = (x[i-1] + x[i])   / 2.0;
		x_center[i+1] = (x[i]   + x[i+1]) / 2.0;
		x_center[i+2] = (x[i+1] + x[i+2]) / 2.0;

		m[i]   = rho[i] * (x[i] - x[i-1]);
		m[i+1] = rho[i+1] * (x[i+1] - x[i]);

		gamma[i+1] = gamma[i];

		shock_wave_front[i]   = 0;
		shock_wave_front[i+1] = true;
		shock_wave_front[i+2] = true;
		frontCellNumber_ = i+1;
	}
}

void GodunovKolganMethod::modifyMesh()
{
	int nOfCellWithMinConc = 1;
	int nNoUnion	 = config_->getNCellsBehindLeadShock();
	int nUnitedCells = config_->getNCellsInReductionZone();
	int numberOfSpecies = config_->getWhatSpecies();
	RealType minConcentration = config_->getMinConcentration();
	int i;

	for (i = 1; i < frontCellNumber_; ++i) {
		if (volumeFractions[i][numberOfSpecies] >= minConcentration) {
			break;
		}
		nOfCellWithMinConc = i;
	}

	int reaction_start = min(nOfCellWithMinConc, frontCellNumber_ - nNoUnion);

	if (reaction_start < nUnitedCells) {
		//cout << "No modifying of the mesh." << endl;
		return;
	}

    //cout << "Objedineniye yacheek" << endl;

	if (reaction_start % 2 != 0) {
		reaction_start--;
	}

    RealType l_mass = 0.0;
    for (i = 1; i <= reaction_start; i++) {
        l_mass += m[i];
    }
    RealType delta_m0 = m[reaction_start+1];
    RealType a = (2 * l_mass - delta_m0 * reaction_start) / (RealType) reaction_start;
    RealType b = (2 * delta_m0 * reaction_start - 2 * l_mass) / 
        ((RealType) reaction_start * (reaction_start-1));

	RealType dx_left;
	RealType dx_right;
	RealType dx_new;
	RealType m_left;
	RealType m_right;
	RealType m_new;
	RealType rho_new;
	RealType *mf_left  = new RealType[mixture->getNSpecies()];
	RealType *mf_right = new RealType[mixture->getNSpecies()];
	RealType *mf_new   = new RealType[mixture->getNSpecies()];
	int offset;

	for (i = 1; i < reaction_start; i++) {
        if (m[i] < (a + b * (i-1))) {
		    dx_left  = x[i] - x[i-1];
		    dx_right = x[i+1] - x[i];
		    dx_new   = dx_left + dx_right;
		    m_left   = rho[i] * dx_left;
		    m_right  = rho[i+1] * dx_right;
    		
		    rho_new = (m_left + m_right) / dx_new;
		    m_new = rho_new * dx_new;
    		
		    mixture->setStateWithURhoX(u_energy[i], 
			    rho[i], 
			    volumeFractions[i]);
		    mixture->updateMassFractions(mf_left);
		    mixture->setStateWithURhoX(u_energy[i+1], 
			    rho[i+1], 
			    volumeFractions[i+1]);
		    mixture->updateMassFractions(mf_right);

		    // Величина, обратная молекулярному весу объединенной ячейки.
		    RealType invMWeightNew = 0;

		    for (int j = 0; j < mixture->getNSpecies(); ++j) {
			    mf_new[j] = (m_left * mf_left[j] + m_right * mf_right[j]) / m_new;
			    invMWeightNew += mf_new[j] / 
				    mixture->substances[j]->molecularWeight;
		    }
		    invMWeightNew /= 100.0;
		    for (int j = 0; j < mixture->getNSpecies(); ++j) {
			    volumeFractions[i][j] = mf_new[j] / invMWeightNew /
				    mixture->substances[j]->molecularWeight;
		    }

		    u[i]		= (m_left * u[i] + m_right * u[i+1]) / m_new;
		    e[i]		= (m_left * e[i] + m_right * e[i+1]) / m_new;
		    u_energy[i] = e[i] - u[i] * u[i] / 2.0;
		    rho[i]		= rho_new;
		    m[i]		= m_new;

		    mixture->setStateWithURhoX(u_energy[i], 
			    rho[i], 
			    volumeFractions[i]);
		    p[i] = mixture->getPressure();
		    rho_u[i] = rho[i] * u[i];
		    rho_e[i] = p[i] / (gamma[i] - 1);
		    x[i] = x[i+1];
		    x_center[i] = (x[i-1] + x[i+1]) / 2.0;
            for (int k = i+1; k < meshSize_; ++k) {
            	offset = k+1;
            	m[k] = m[offset];
            	u[k] = u[offset];
            	e[k] = e[offset];
            	rho[k] = rho[offset];
            	u_energy[k] = u_energy[offset];
            	p[k] = p[offset];
            	x[k] = x[offset];
            	x_center[k] = x_center[offset];
            	rho_u[k] = rho_u[offset];
            	rho_e[k] = rho_e[offset];
            	for (int j = 0; j < mixture->getNSpecies(); j++) {
            		volumeFractions[k][j] = volumeFractions[offset][j];
            	}
                gamma[k] = gamma[offset];
                shock_wave_front[k] = shock_wave_front[offset];
            }
            meshSize_--;
            frontCellNumber_--;
            reaction_start--;
        }
	}

	delete [] mf_left;
	delete [] mf_right;
	delete [] mf_new;

	for (i = frontCellNumber_ + 1; i <= meshSize_; ++i) {
		u_contact[i] = config_->getUInitial();
		p_contact[i] = config_->getPInitial();
	}
}

void GodunovKolganMethod::update_()
{
	int num_digits;
	char time[32];
	// TODO: добавить проверку на существование каталога.
	string fullname("Output\\data_");

	num_digits = sprintf_s(time, "%d", config_->getStart());
	fullname += time;
	fullname += ".txt";

	ifstream file(fullname.c_str());
	if (! file) {
		cout << "Cannot read file '" << fullname << "'." << endl;
		cout << "Calculations can not be done from time step" 
			 << config_->getStart() 
			 << endl;
		exit(-1);
	}
	else {
		cout << "Starting from time step " << config_->getStart() 
			 << "."
			 << endl;
	}

	string s;
	getline(file, s);
	int length = 0;

	file >> x[0];

	for (int i = 1; !file.eof(); ++i) {
		++length;
	    file >> s;
	    file >> x[i];
	    file >> x_center[i];
	    file >> p[i];
	    file >> u[i];
	    file >> rho[i];
		rho_u[i] = rho[i] * u[i];
		rho_e[i] = p[i] / (config_->getGammaBehindFront() - 1);
	    file >> e[i];
	    file >> u_energy[i];
	    internal_energy[i] = p[i] / 
			((config_->getGammaBehindFront() - 1) * rho[i]);
		for (int j = 0; j < mixture->getNSpecies(); ++j) {
			file >> volumeFractions[i][j];
		}
	    shock_wave_front[i] = false;
	    m[i] = (x[i] - x[i-1]) * rho[i];
	    gamma[i] = config_->getGammaBehindFront();
	}
	internal_energy[length] = p[length] / 
		((config_->getGammaAheadFront() - 1) * rho[length]);
	gamma[length] = config_->getGammaAheadFront();
	shock_wave_front[length-1] = true;
	shock_wave_front[length]   = true;

	for (int i = length + 1; i <= meshSize_; ++i) {
		x[i] = x[i-1] + config_->getDx();
		x_center[i] = (x[i-1] + x[i]) / 2.0;
	}
}