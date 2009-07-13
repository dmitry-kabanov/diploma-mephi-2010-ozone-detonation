#include <math.h>
#include "func_chemical_kinetics.h"
#include "constants.h"

long double calc_mass_fraction_of_initial_substance(
    long double w, long double t)
{
    return w - DT * Z * w * exp(-ACTIVATION_ENERGY / (R_GAS * t));              
} // long double calc_mass_fraction_of_initial_substance(long double, long double)
