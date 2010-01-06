#ifndef FUNC_INIT_H
#define FUNC_INIT_H

#include "RealType.h"
#include "StifflSolver.h"

void init_parameters(
    RealType *x, RealType *x_center,
    RealType *m, RealType *rho,
    RealType *u, RealType *e,
    RealType *p, RealType *rho_u,
    bool *shock_wave_front
);

/**
 * «адает начальное распределение термодинамических величин.
 * ‘ункци€ принимает также аргументы:
 * p (давление) - дл€ расчета внутренней энергии,
 * rho (плотность) - дл€ расчета внутренней энергии,
 * u (скорость) - дл€ расчета полной энергии.
 */
void init_thermodynamic_parameters(RealType *internal_energy,
                                   RealType *energy,
                                   RealType *p,
                                   RealType *rho,
                                   RealType *u,
                                   RealType *u_energy,
                                   RealType *volumeFractions[],
                                   StifflSolver &kinetics
);

/**
 * «адает начальное распределение rho * energy.
 * 
 * internal_energy Ч внутренн€€ теплова€ энерги€
 * rho             Ч плотность
 * rho_e           Ч внутренн€ энерги€ на единицу объема.
 */
void init_rho_e(const RealType *internal_energy,
                const RealType *rho,
                RealType *rho_e
);

void init_additional_parameters(
    RealType *p_contact, RealType *u_contact,
    RealType *impulse_flow, RealType *energy_flow,
    RealType *p, RealType *u
);

void init_boundary_parameters(
    RealType *rho_bound_r, RealType *rho_bound_l,
    RealType *rho_u_bound_r, RealType *rho_u_bound_l,
    RealType *rho_e_bound_r, RealType *rho_e_bound_l,
    RealType *u_bound_r, RealType *u_bound_l,
    RealType *e_bound_r, RealType *e_bound_l,
    RealType *p_bound_r, RealType *p_bound_l
);

void init_tangents(
    RealType *rho_tg_left, RealType *rho_tg_right, RealType *rho_tg,
    RealType *rho_u_tg_left, RealType *rho_u_tg_right, RealType *rho_u_tg,
    RealType *rho_e_tg_left, RealType *rho_e_tg_right, RealType *rho_e_tg
);

void init_gamma(RealType *gamma);

void init_volume_fractions(RealType *volumeFractions[]);

#endif