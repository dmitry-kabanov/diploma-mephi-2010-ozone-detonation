#include <math.h>
#include "constants.h"
#include "func_thermodynamics.h"
#include "func_piston.h"

long double calc_volume(long double w)
{
    long double a, b, c;
    long double d, x1, x2;
    a = (GAMMA + 1) * D_C_J * D_C_J;
    b = -2 * GAMMA * P2 * V2 - 2 * GAMMA * D_C_J * D_C_J;
    c = 2 * Q * w + 2 * GAMMA * I0 + GAMMA * D_C_J * D_C_J + 2 * GAMMA * Q - 
        2 * GAMMA * Q * w + 2 * GAMMA * P2 * V2 - D_C_J * D_C_J - 
        2 * P2 * V2 - 2 * Q - 2 * I0;

    d = b * b - 4 * a * c;
    //x1 = (-b + sqrt(d)) / (2 * a);
    x2 = (-b - sqrt(d)) / (2 * a) * V2;

    return x2;
}

long double internal_energy_func(long double v)
{
    return I0 + 0.5 * pow(RHO2 * D_C_J * (V2 - v), 2);
}

long double p_rayleigh(long double v)
{
    return P2 + pow(RHO2 * D_C_J, 2) * (V2 - v);
}


long double pressure_func(long double v, long double chemical_energy)
{
    return ((internal_energy_func(v) + chemical_energy) * 
        (GAMMA - 1)) / v - p_rayleigh(v);
}

long double d_pressure_func(long double v, long double chemical_energy)
{
    return -((internal_energy_func(v) + chemical_energy) * 
        (GAMMA - 1)) / pow(v, 2) - 
        pow(RHO2 * D_C_J, 2) * (V2 - v) / v + 
        pow(RHO2 * D_C_J, 2);
}

long double calc_piston_velocity(long double p, long double v)
{
    return sqrt((p-P2) * (V2 - v));
}

void calc_piston(long double *piston, long double w)
{
    long double v;
    long double chemical_energy = (1 - w) * Q;
    long double x, y, p_rayleigh, z, ph, uh, vh, c0, e, e2;
    long double mass, impulse, enthalpy;

    if (w < 10e-6) {
        w = 10e-6;
    }

    v = calc_volume(w);

    piston[0] = calc_pressure(internal_energy_func(v), (1 - w) * Q, 1.0 / v);
    piston[1] = calc_piston_velocity(piston[0], v);

    p_rayleigh = P2 + pow(RHO2 * D_C_J, 2) * (V2 - v);

    x = sqrt(GAMMA * piston[0] * v) + piston[1] - D_C_J;
    y = sqrt(GAMMA * piston[0] * v) / v;

    c0 = sqrt(GAMMA * P2 / RHO2);

    ph = P2 + RHO2 * D_C_J * D_C_J * (1 - (c0 * c0) / (D_C_J * D_C_J)) / (GAMMA + 1);
    uh = D_C_J * (1 - (c0 * c0) / (D_C_J * D_C_J)) / (GAMMA + 1);

    e = I0 + (piston[0] - P2) * (V2 - v) / 2.0 + Q;
    e2 = piston[0] * v / (GAMMA - 1);

    mass = RHO2 * D_C_J - (D_C_J - piston[1]) / v;
    impulse = P2 + RHO2 * D_C_J * D_C_J - piston[0] - pow(D_C_J - piston[1], 2) / v;
    enthalpy = I0 + D_C_J * D_C_J / 2.0 + P2 / RHO2 + Q - 
        v * piston[0] / (GAMMA - 1) - piston[0] * v - 
        pow(D_C_J - piston[1], 2) / 2.0;
}
