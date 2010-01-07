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

#include "constants.h"
#include "Output.h"
#include "StifflSolver.h"

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

}