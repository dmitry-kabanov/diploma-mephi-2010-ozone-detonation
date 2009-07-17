#include <cmath>
#include "real_number_type.h"
#include "constants.h"

void init_parameters(
    long double *x, long double *x_center,
    long double *m, long double *rho,
    long double *u, long double *e,
    long double *p, long double *rho_u,
    bool *shock_wave_front)
{
    for (int i = 0; i <= N; i++) {
        shock_wave_front[i] = false;
    }

	for (int i = 0; i <= N; i++) {
        if (i == 0) {
            x[i] = 0;
            x_center[i] = 0;
            m[i] = 0;
            rho[i] = 0;
            u[i] = 0;
            e[i] = 0;
            p[i] = 0;
            rho_u[i] = 0;
        }
        else {
            if (x[i-1] + DX1 <= L) {
                x[i] = x[i-1] + DX1;
            }
            else {
                x[i] = x[i-1] + DX2;
            }

            shock_wave_front[L] = true;
            shock_wave_front[L+1] = true;
            
            if (i <= L) {
                p[i] = P1;
                rho[i] = RHO1;
                u[i]   = U1;
            }
            else {
		        p[i] = P2;
                rho[i] = RHO2;
                u[i] = U2;
            }

            m[i] = (x[i] - x[i-1]) * rho[i];

            x_center[i] = (x[i-1] + x[i]) / 2.0;
            rho_u[i] = rho[i] * u[i];
        }
	}
}

void init_thermodynamic_parameters(
    long double *internal_energy, long double *energy,
    long double *temperature, long double *p, 
    long double *rho, long double *u)
{
    internal_energy[0] = 0.0;
    energy[0]          = 0.0;
    temperature[0]     = 0.0;

    for (int i = 1; i <= N; i++) {
        internal_energy[i] = p[i] / ((GAMMA - 1) * rho[i]);
        energy[i]          = internal_energy[i] + pow(u[i], 2) / 2.0;
        temperature[i]     = internal_energy[i] / CV;
    }
}

void init_chemical_parameters(long double *w, long double *chemical_energy)
{
    w[0] = 0.0;
    chemical_energy[0] = 0.0;

    for (int i = 1; i <= N; i++) {
        w[i] = 1.0;
        chemical_energy[i] = (1 - w[i]) * Q;
    }
}

void init_rho_e(real_t *internal_energy, 
                real_t *chemical_energy, 
                real_t *rho,
                real_t *rho_e)
{
    rho_e[0] = 0.0;

    for (int i = 1; i <= N; i++) {
        rho_e[i] = rho[i] * (internal_energy[i] + chemical_energy[i]);
    }
}

void init_additional_parameters(
    long double *p_contact, long double *u_contact,
    long double *impulse_flow, long double *energy_flow,
    long double *p, long double *u)
{
    for (int i = 0; i <= N; i++) {
        p_contact[i] = p[i];
        u_contact[i] = u[i];
        impulse_flow[i] = 0;
        energy_flow[i] = 0;
	}
}

void init_boundary_parameters(
    long double *rho_bound_r, long double *rho_bound_l,
    long double *rho_u_bound_r, long double *rho_u_bound_l,
    long double *rho_e_bound_r, long double *rho_e_bound_l,
    long double *u_bound_r, long double *u_bound_l,
    long double *e_bound_r, long double *e_bound_l,
    long double *p_bound_r, long double *p_bound_l)
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
    long double *rho_tg_left, long double *rho_tg_right, long double *rho_tg,
    long double *rho_u_tg_left, long double *rho_u_tg_right, long double *rho_u_tg,
    long double *rho_e_tg_left, long double *rho_e_tg_right, long double *rho_e_tg)
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