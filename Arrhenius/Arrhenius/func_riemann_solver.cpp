#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "constants.h"
#include "func_riemann_solver.h"

long double p0 = 0;

long double calc_c(long double p, long double rho)
{
    return sqrt(GAMMA * p / rho);
} 

long double calc_a(long double p, long double p_k, long double rho_k)
{
    long double value;

    if (p >= p_k) {
        value = sqrt(rho_k * ((GAMMA + 1) * (p + p0) / 2.0 +
            (GAMMA - 1) * (p_k + p0) / 2.0));
    }
    else {
        long double c_k, pi_k, chisl, znam;
        c_k = calc_c(p_k, rho_k);
        pi_k = (p + p0) / (p_k + p0);
        pi_k = pow(pi_k, (GAMMA - 1.0) / (2.0 * GAMMA));

        chisl = ((GAMMA - 1) / (2.0 * GAMMA)) * rho_k * c_k * 
            (1.0 - (p + p0) / (p_k + p0));
        znam = 1.0 - pi_k;
        value = chisl / znam;
    }

    return value;
}

long double calc_u_contact(long double p, long double p_left, long double p_right,
              long double rho_left, long double rho_right,
              long double u_left, long double u_right)
{
    long double u;

    if ((EXACT_RIEMANN_SOLVER == false)
        || (abs(p_left - p) <= 10e-10)
        || (abs(p_right - p) <= 10e-10)) {
        long double c_left, c_right;

        c_left = calc_c(p_left, rho_left);
        c_right = calc_c(p_right, rho_right);

        u = (p_left - p_right +
            rho_left * c_left * u_left + rho_right * c_right * u_right) /
            (rho_left * c_left + rho_right * c_right);
        return u;
    }

    long double a_left, a_right;

    a_left  = calc_a(p, p_left, rho_left);
    a_right = calc_a(p, p_right, rho_right);

    u = (a_left * u_left + a_right * u_right + p_left - p_right) /
        (a_left + a_right);

    return u;
}

long double calc_shock_wave_velocity(long double p, long double u,
                         long double p_left, long double p_right,
                         long double rho_left, long double rho_right,
                         long double u_left, long double u_right)
{
    long double a, d;

    // Правая волна является ударной.
    if (p >= p_right) {
        a = calc_a(p, p_right, rho_right);
        d = u_right + a / rho_right;
    }
    // Левая волна является ударной.
    else if (p >= p_left) {
        a = calc_a(p, p_left, rho_left);
        d = u_left - a / rho_left;
    }
    else {
        d = 0.0;
    }

    return d;
}




long double calc_u_ud(long double p_left, long double p_right, long double rho_left)
{
    long double chisl, znam;

    chisl = p_right - p_left;
    znam  = rho_left * ((GAMMA + 1) * (p_right + p0) / 2.0 +
        (GAMMA - 1) * (p_left + p0) / 2.0);
    
    return chisl / sqrt(znam);
}

long double calc_u_raz(long double p_left, long double p_right, long double rho_right)
{
    long double c_right, value;

    c_right = calc_c(p_right, rho_right);

    value = - 2.0 * c_right * (1 - pow(
        (p_left + p0) / (p_right + p0),
        (GAMMA - 1) / (2 * GAMMA)
        )) / (GAMMA - 1);

    return value;
}

long double calc_u_vac(long double p_left, long double p_right,
                  long double rho_left, long double rho_right)
{
    long double c_left, c_right;

    c_left  = calc_c(p_left, rho_left);
    c_right = calc_c(p_right, rho_right);

    return - 2.0 * (c_left + c_right) / (GAMMA - 1);
}

int calc_config(long double p_left, long double p_right,
                   long double rho_left, long double rho_right,
                   long double u_left, long double u_right)
{
    long double u_ud, u_raz, u_vac;

    if (p_left > p_right) {
        long double temp;

        temp = - u_left;
        u_left = - u_right;
        u_right = temp;

        temp = p_left;
        p_left = p_right;
        p_right = temp;

        temp = rho_left;
        rho_left = rho_right;
        rho_right = temp;
    }

    u_ud  = calc_u_ud(p_left, p_right, rho_left);
    u_raz = calc_u_raz(p_left, p_right, rho_right);
    u_vac = calc_u_vac(p_left, p_right, rho_left, rho_right);

    if ((u_left - u_right) > u_ud) {
        return 1;
    }
    else if (u_raz < (u_left - u_right) && (u_left - u_right) < u_ud) {
        return 2;
    }
    else if (u_vac < (u_left - u_right) && (u_left - u_right) < u_raz) {
        return 3;
    }
    else if ((u_left - u_right) < u_vac) {
        return 4;
    }
}

long double calc_rho_contact_right(int config,
                              long double p, long double p_left, long double p_right,
                              long double rho_right,
                              long double u_right, long double u)
{
    long double c_right, a_right, rho_contact_right;

    c_right = calc_c((p_right + p0), rho_right);
    a_right = calc_a(p, p_right, rho_right);

    switch (config) {
        case CONFIG_SHOCK_WAVES:
            // Вправо бежит ударная волна.
            rho_contact_right = rho_right * a_right /
                (a_right + rho_right * (u_right - u));
            break;
        case CONFIG_SHOCK_AND_RAREFACTION_WAVES:
            // Вправо бежит ударная волна.
            if (p_left > p_right) {
                rho_contact_right = rho_right * a_right /
                    (a_right + rho_right * (u_right - u));
            }
            // Вправо бежит волна разрежения.
            else {
                rho_contact_right = GAMMA * (p + p0) /
                    pow(c_right - (GAMMA - 1) * (u_right - u) / 2.0, 2);
            }
            break;
        case CONFIG_RAREFACTION_WAVES:
            // Вправо бежит волна разрежения.
            rho_contact_right = GAMMA * (p + p0) /
                    pow(c_right - (GAMMA - 1) * (u_right - u) / 2.0, 2);
            break;
        case CONFIG_VACUUM:
            // Истечение в вакуум.
            rho_contact_right = 0.0;
            break;
    }

    return rho_contact_right;
}

long double calc_rho_contact_left(int config,
                             long double p, long double p_left, long double p_right,
                             long double rho_left,
                             long double u_left, long double u)
{
    long double c_left, a_left, rho_contact_left;

    c_left = calc_c((p_left + p0), rho_left);
    a_left = calc_a(p, p_left, rho_left);

    switch (config) {
        case CONFIG_SHOCK_WAVES:
            // Влево бежит ударная волна.
            rho_contact_left = rho_left * a_left /
                (a_left - rho_left * (u_left - u));
            break;
        case CONFIG_SHOCK_AND_RAREFACTION_WAVES:
            // Влево бежит ударная волна
            if (p_left <= p_right) {
                rho_contact_left = rho_left * a_left /
                    (a_left - rho_left * (u_left - u));
            }
            // Влево бежит волна разрежения.
            else {
                rho_contact_left = GAMMA * (p + p0) /
                    pow(c_left + (GAMMA - 1) * (u_left - u) / 2.0, 2);
            }
            break;
        case CONFIG_RAREFACTION_WAVES:
            // Влево бежит волна разрежения.
            rho_contact_left = GAMMA * (p + p0) /
                    pow(c_left + (GAMMA - 1) * (u_left - u) / 2.0, 2);
            break;
        case CONFIG_VACUUM:
            // Истечение в вакуум.
            rho_contact_left = 0.0;
            break;
    }

    return rho_contact_left;
}


/**
 * Уравнение для ударного скачка.
 */
long double calc_shock_wave(long double p2, long double p1, long double rho1, long double u1, char direction)
{
    long double value, u2;

    value = (p2 - p1) * sqrt((2*(1/rho1)) / ((GAMMA + 1) * p2 + (GAMMA - 1) * p1));

    // Ударный скачок бежит вправо.
    if (direction == 'r') {
        u2 = u1 + value;
    }
    // Ударный скачок бежит влево.
    else {
        u2 = u1 - value;
    }

    return u2;
}

/**
 * Уравнение для волны разрежения.
 * p1        - начальное давление
 * rho1      - начальная плотность
 * u1        - начальная скорость
 * u         - скорость контактного разрыва
 * direction - направление движения волны
 */
long double calc_rarefaction_wave(long double p1, long double rho1, long double u1, long double u, char direction)
{
    long double c1, p;

    // Считаем скорость звука.
    c1 = calc_c(p1, rho1);

    // Волна разрежения бежит направо.
    if (direction == 'r') {
        p = p1 * pow(1 + (GAMMA - 1) * (u - u1) / (2 * c1), (2 * GAMMA) / (GAMMA - 1));
    }
    // Волная разрежения бежит налево.
    else {
        p = p1 * pow(1 - (GAMMA - 1) * (u - u1) / (2 * c1), (2 * GAMMA) / (GAMMA - 1));
    }

    return p;
}

void check_solution(long double p_left, long double p_right,
                    long double rho_left, long double rho_right,
                    long double u_left, long double u_right,
                    long double p, long double u,
                    int config)
{
    long double p_analit, u_analit;

    switch (config) {
        case CONFIG_SHOCK_WAVES:
            break;
        case CONFIG_SHOCK_AND_RAREFACTION_WAVES:
            if (p_left >= p_right) {
                p_analit = calc_rarefaction_wave(p_left, rho_left, u_left, u, 'l');
                u_analit = calc_shock_wave(p_analit, p_right, rho_right, u_right, 'r');
            }
            else {
                p_analit = calc_rarefaction_wave(p_right, rho_right, u_right, u, 'r');
                u_analit = calc_shock_wave(p_analit, p_left, rho_left, u_left, 'l');
            }
            printf("P analit = %f\nu analit = %f\n", p_analit, u_analit);
            break;
        case CONFIG_RAREFACTION_WAVES:
            break;
        case CONFIG_VACUUM:
            break;
    }
}

long double f(long double p, long double p_k, long double rho_k)
{
    long double value, c_k, pi_k;

    c_k  = calc_c(p_k, rho_k);
    pi_k = (p + p0) / (p_k + p0);

    if (p >= p_k) {
        value = (p - p_k) / (rho_k * c_k * 
            sqrt((GAMMA + 1) / (2.0 * GAMMA) * pi_k + (GAMMA - 1) / (2.0 * GAMMA)));
    }
    else {
        value = 2.0 * c_k * (pow(pi_k, (GAMMA - 1) / (2.0 * GAMMA)) - 1) / (GAMMA - 1);
    }
    
    return value;
}

long double df(long double p, long double p_k, long double rho_k)
{
    long double c_k, pi_k, value;

    c_k  = calc_c(p_k + p0, rho_k);
    pi_k = (p + p0) / (p_k + p0);

    if (p >= p_k) {
        long double chisl, znam;

        chisl = (GAMMA + 1) * pi_k + (3 * GAMMA - 1);
        znam  = 4 * GAMMA * rho_k * c_k * sqrt(
            pow(
                (GAMMA + 1) / (2.0 * GAMMA) * pi_k + (GAMMA - 1) / (2.0 * GAMMA),
                3
            )
        );
        value = chisl / znam;
    }
    else {
        value = c_k * pow(pi_k, (GAMMA - 1) / (2.0 * GAMMA)) / (GAMMA * (p + p0));
    }

    return value;
}

long double calc_p_contact(long double p_left, long double p_right,
                      long double rho_left, long double rho_right,
                      long double u_left, long double u_right)
{
    //int config;
    
    long double p;

	long double c_left, c_right;
    long double summa_f, delta_u;

    //config = calc_config(p_left, p_right, rho_left, rho_right, u_left, u_right);

    //switch (config) {
    //    case CONFIG_SHOCK_WAVES:
    //        printf("Вправо и влево распространяются ударные волны.\n");
    //        break;
    //    case CONFIG_SHOCK_AND_RAREFACTION_WAVES:
    //        if (p_left <= p_right) {
    //            printf("Влево - ударная волна, вправо - волна разрежения.\n");
    //        }
    //        else {
    //            printf("Влево - волна разрежения, вправо - ударная волна.\n");
    //        }
    //        break;
    //    case CONFIG_RAREFACTION_WAVES:
    //        printf("Вправо и влево распространяются волны разрежения.\n");
    //        break;
    //    case CONFIG_VACUUM:
    //        printf("Истечение в вакуум. Решения нет.\n");
    //        return 0;
    //        break;
    //}

    c_left  = calc_c(p_left, rho_left);
    c_right = calc_c(p_right, rho_right);

	p = (p_left * rho_right * c_right + p_right * rho_left * c_left + 
		(u_left - u_right) * rho_left * c_left * rho_right * c_right) / 
		(rho_left * c_left + rho_right * c_right);

    // Если решаем в акустическом приближении.
    if (EXACT_RIEMANN_SOLVER == false) {
        long double u;
        u = (p_left - p_right + c_left * rho_left * u_left
            + c_right * rho_right *u_right)
            / (rho_left * c_left + rho_right * c_right);
        p = (u - u_right) * c_left * rho_left + p_right;

        return p;
    }

    summa_f = f(p, p_left, rho_left) + f(p, p_right, rho_right);
    delta_u = u_left - u_right;

    while (abs(summa_f - delta_u) > EPSILON) {
        p = p - (summa_f - delta_u) /
            (df(p, p_left, rho_left) + df(p, p_right, rho_right));

        summa_f = f(p, p_left, rho_left) + f(p, p_right, rho_right);
    }

    //rho_contact_left  = calc_rho_contact_left(
    //    config, p, p_left, p_right, rho_left, u_left, u);
    //rho_contact_right = calc_rho_contact_right(
    //    config, p, p_left, p_right, rho_right, u_right, u);
    
    return p;
}