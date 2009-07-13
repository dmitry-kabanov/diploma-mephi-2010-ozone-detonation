#ifndef FUNC_RIEMANN_SOLVER_H
#define FUNC_RIEMANN_SOLVER_H

long double calc_p_contact(long double p_left, long double p_right,
                             long double rho_left, long double rho_right,
                             long double u_left, long double u_right);

long double calc_u_contact(long double p, long double p_left, long double p_right,
                             long double rho_left, long double rho_right,
                             long double u_left, long double u_right);

long double calc_shock_wave_velocity(long double p, long double u,
                                     long double p_left, long double p_right,
                                     long double rho_left, long double rho_right,
                                     long double u_left, long double u_right);

#endif