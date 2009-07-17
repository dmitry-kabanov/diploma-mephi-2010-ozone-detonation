#include <math.h>
#include "constants.h"
#include "func_thermodynamics.h"
#include "func_piston.h"

long double calc_volume(long double w)
{
    long double a, b, c;
    long double d, x2;
    a = (GAMMA + 1) * D_C_J * D_C_J;
    b = -2 * GAMMA * P2 * V2 - 2 * GAMMA * D_C_J * D_C_J;
    c = 2 * Q * w + 2 * GAMMA * I0 + GAMMA * D_C_J * D_C_J + 2 * GAMMA * Q - 
        2 * GAMMA * Q * w + 2 * GAMMA * P2 * V2 - D_C_J * D_C_J - 
        2 * P2 * V2 - 2 * Q - 2 * I0;

    d = b * b - 4 * a * c;
    x2 = (-b - sqrt(d)) / (2 * a) * V2;

    return x2;
}

long double internal_energy_func(long double v)
{
    return I0 + P2 * (V2 - v) + 0.5 * pow(RHO2 * D_C_J * (V2 - v), 2);
}

long double calc_piston_velocity(long double p, long double v)
{
    return sqrt((p-P2) * (V2 - v));
}

void calc_piston(long double *piston, long double w)
{
    long double v;

    if (w < 10e-6) {
        w = 10e-6;
    }

    v = calc_volume(w);

    piston[0] = calc_pressure(internal_energy_func(v), (1 - w) * Q, 1.0 / v);
    piston[1] = calc_piston_velocity(piston[0], v);
}
