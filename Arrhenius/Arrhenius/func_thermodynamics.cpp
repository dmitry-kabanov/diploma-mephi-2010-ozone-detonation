#include "func_thermodynamics.h"
#include "constants.h"

long double calc_temperature(
    long double internal_energy, 
    long double chemical_energy)
{
    return (internal_energy + chemical_energy) / CV;
}

long double calc_pressure(
    long double internal_energy,
    long double chemical_energy,
    long double density)
{
    return (internal_energy + chemical_energy) * (GAMMA - 1) * density;
}
