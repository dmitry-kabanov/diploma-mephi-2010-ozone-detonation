#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "constants.h"
#include "func_riemann_solver.h"

// ‘лаг, определ€ющий, решаетс€ ли задача –имана точно (значение true)
// или же в акустическом приближении.
const bool EXACT_RIEMANN_SOLVER = true;

// “очность нахождени€ давлени€ на контактном разрыве.
const long double EPSILON = 1e-10;

long double p0 = 0;

long double calc_c(long double p, long double rho, RealType gamma)
{
    return sqrt(gamma * p / rho);
} 

long double calc_a(long double p, long double p_k, long double rho_k, RealType gamma)
{
    long double value;

    if (p >= p_k) {
        value = sqrt(rho_k * ((gamma + 1) * (p + p0) / 2.0 +
            (gamma - 1) * (p_k + p0) / 2.0));
    }
    else {
        long double c_k, pi_k, chisl, znam;
        c_k = calc_c(p_k, rho_k, gamma);
        pi_k = (p + p0) / (p_k + p0);
        pi_k = pow(pi_k, (gamma - 1.0) / (2.0 * gamma));

        chisl = ((gamma - 1) / (2.0 * gamma)) * rho_k * c_k * 
            (1.0 - (p + p0) / (p_k + p0));
        znam = 1.0 - pi_k;
        value = chisl / znam;
    }

    return value;
}

long double calc_u_contact(long double p, long double p_left, long double p_right,
              long double rho_left, long double rho_right,
              long double u_left, long double u_right, RealType gamma)
{
    long double u;

    if ((EXACT_RIEMANN_SOLVER == false)
        || (abs(p_left - p) <= 1e-3)
        || (abs(p_right - p) <= 1e-3)) {
        long double c_left, c_right;

        c_left = calc_c(p_left, rho_left, gamma);
        c_right = calc_c(p_right, rho_right, gamma);

        u = (p_left - p_right +
            rho_left * c_left * u_left + rho_right * c_right * u_right) /
            (rho_left * c_left + rho_right * c_right);
        return u;
    }

    long double a_left, a_right;

    a_left  = calc_a(p, p_left, rho_left, gamma);
    a_right = calc_a(p, p_right, rho_right, gamma);

    u = (a_left * u_left + a_right * u_right + p_left - p_right) /
        (a_left + a_right);

    return u;
}

long double calc_shock_wave_velocity(long double p, long double u,
                         long double p_left, long double p_right,
                         long double rho_left, long double rho_right,
                         long double u_left, long double u_right,
                         RealType gamma)
{
    long double a, d;

    // ѕрава€ волна €вл€етс€ ударной.
    if (p >= p_right) {
        a = calc_a(p, p_right, rho_right, gamma);
        d = u_right + a / rho_right;
    }
    // Ћева€ волна €вл€етс€ ударной.
    else if (p >= p_left) {
        a = calc_a(p, p_left, rho_left, gamma);
        d = u_left - a / rho_left;
    }
    else {
        d = 0.0;
    }

    return d;
}




long double f(long double p, long double p_k, long double rho_k, RealType gamma)
{
    long double value, c_k, pi_k;

    c_k  = calc_c(p_k, rho_k, gamma);
    pi_k = (p + p0) / (p_k + p0);

    if (p >= p_k) {
        value = (p - p_k) / (rho_k * c_k * 
            sqrt((gamma + 1) / (2.0 * gamma) * pi_k + (gamma - 1) / (2.0 * gamma)));
    }
    else {
        value = 2.0 * c_k * (pow(pi_k, (gamma - 1) / (2.0 * gamma)) - 1) / (gamma - 1);
    }
    
    return value;
}

long double df(long double p, long double p_k, long double rho_k, RealType gamma)
{
    long double c_k, pi_k, value;

    c_k  = calc_c(p_k + p0, rho_k, gamma);
    pi_k = (p + p0) / (p_k + p0);

    if (p >= p_k) {
        long double chisl, znam;

        chisl = (gamma + 1) * pi_k + (3 * gamma - 1);
        znam  = 4 * gamma * rho_k * c_k * sqrt(
            pow(
                (gamma + 1) / (2.0 * gamma) * pi_k + (gamma - 1) / (2.0 * gamma),
                3
            )
        );
        value = chisl / znam;
    }
    else {
        value = c_k * pow(pi_k, (gamma - 1) / (2.0 * gamma)) / (gamma * (p + p0));
    }

    return value;
}

long double calc_p_contact(long double p_left, long double p_right,
                      long double rho_left, long double rho_right,
                      long double u_left, long double u_right, RealType gamma)
{
    long double p;
	long double c_left, c_right;
    long double summa_f, delta_u;

    c_left  = calc_c(p_left, rho_left, gamma);
    c_right = calc_c(p_right, rho_right, gamma);

	p = (p_left * rho_right * c_right + p_right * rho_left * c_left + 
		(u_left - u_right) * rho_left * c_left * rho_right * c_right) / 
		(rho_left * c_left + rho_right * c_right);

    // ≈сли решаем в акустическом приближении.
    if (EXACT_RIEMANN_SOLVER == false) {
        long double u;
        u = (p_left - p_right + c_left * rho_left * u_left
            + c_right * rho_right *u_right)
            / (rho_left * c_left + rho_right * c_right);
        p = (u - u_right) * c_left * rho_left + p_right;

        return p;
    }

    summa_f = f(p, p_left, rho_left, gamma) + f(p, p_right, rho_right, gamma);
    delta_u = u_left - u_right;

    while (abs(summa_f - delta_u) > EPSILON) {
        p = p - (summa_f - delta_u) /
            (df(p, p_left, rho_left, gamma) + df(p, p_right, rho_right, gamma));

        summa_f = f(p, p_left, rho_left, gamma) + f(p, p_right, rho_right, gamma);
    }

    return p;
}