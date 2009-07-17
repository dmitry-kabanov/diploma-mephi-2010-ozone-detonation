#include <cmath>
#include "real_number_type.h"
#include "func_chemical_kinetics.h"
#include "constants.h"

long double calc_mass_fraction_of_initial_substance(
    long double w, long double t)
{
    real_t delta_w;
    delta_w = DT * Z * w * exp(-ACTIVATION_ENERGY / (R_GAS * t));
    
    return (w - delta_w);
}
