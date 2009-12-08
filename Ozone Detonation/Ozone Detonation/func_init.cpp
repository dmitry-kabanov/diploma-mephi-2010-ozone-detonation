#include <cmath>
#include "RealType.h"
#include "RungeKuttaMethod.h"
#include "constants.h"

void init_parameters(RealType *x, RealType *x_center,
                     RealType *m, RealType *rho,
                     RealType *u, RealType *e,
                     RealType *p, RealType *rho_u,
                     bool *shock_wave_front)
{
    for (int i = 0; i <= N; i++) {
        shock_wave_front[i] = false;
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
}

void init_thermodynamic_parameters(RealType *internal_energy,
    RealType *energy,
    RealType *p, 
    RealType *rho,
    RealType *u,
    RealType *u_energy,
    RealType *volumeFractions[],
    RungeKuttaMethod& kinetics)
{
    internal_energy[0] = 0.0;
    energy[0]          = 0.0;
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

        kinetics.getMixture()->setVolumeFractions(volumeFractions[i]);
        kinetics.getMixture()->calculateInitialMolecularWeight();
        mw = kinetics.getMixture()->getMolecularWeight();
        t = p[i] * mw / (Mixture::R_J_OVER_KMOL_K * rho[i]);
        kinetics.getMixture()->setConcentrations(t, rho[i]);
        u_energy[i] = kinetics.getMixture()->calculateMixtureEnthalpy(t) -
            p[i] / rho[i];

        energy[i] = u_energy[i] + u[i] * u[i] / 2.0;
    }
}

void init_rho_e(const RealType *internal_energy, 
                const RealType *rho,
                RealType *rho_e)
{
    rho_e[0] = 0.0;

    for (int i = 1; i <= N; i++) {
        rho_e[i] = rho[i] * internal_energy[i];
    }
}

void init_additional_parameters(
    RealType *p_contact, RealType *u_contact,
    RealType *impulse_flow, RealType *energy_flow,
    RealType *p, RealType *u)
{
    for (int i = 0; i <= N; i++) {
        p_contact[i] = p[i];
        u_contact[i] = u[i];
        impulse_flow[i] = 0;
        energy_flow[i] = 0;
	}
}

void init_boundary_parameters(
    RealType *rho_bound_r, RealType *rho_bound_l,
    RealType *rho_u_bound_r, RealType *rho_u_bound_l,
    RealType *rho_e_bound_r, RealType *rho_e_bound_l,
    RealType *u_bound_r, RealType *u_bound_l,
    RealType *e_bound_r, RealType *e_bound_l,
    RealType *p_bound_r, RealType *p_bound_l)
{
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
}

void init_tangents(
    RealType *rho_tg_left, RealType *rho_tg_right, RealType *rho_tg,
    RealType *rho_u_tg_left, RealType *rho_u_tg_right, RealType *rho_u_tg,
    RealType *rho_e_tg_left, RealType *rho_e_tg_right, RealType *rho_e_tg)
{
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
}

void init_gamma(RealType* gamma)
{
    for (int i = 0; i <= N; i++) {
        if (i <= L) {
            gamma[i] = GAMMA_BEHIND_FRONT;
        } else {
            gamma[i] = GAMMA_AHEAD_FRONT;
        }
    }
}

void init_volume_fractions(RealType* volumeFractions[])
{
    for (int i = 0; i <= L; i++) {
        //volumeFractions[i][0] = 7.16689;
        //volumeFractions[i][1] = 92.8312;
        //volumeFractions[i][2] = 0.0;
        volumeFractions[i][0] = 0.0;
        volumeFractions[i][1] = 0.0;
        volumeFractions[i][2] = 100.0;
    }
    for (int i = L+1; i <= N; i++) {
        volumeFractions[i][0] = 0.0;
        volumeFractions[i][1] = 0.0;
        volumeFractions[i][2] = 100.0;
    }
}