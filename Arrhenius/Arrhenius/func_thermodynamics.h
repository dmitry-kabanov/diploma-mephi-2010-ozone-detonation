#ifndef FUNC_THERMODYNAMICS_H
#define FUNC_THERMODYNAMICS_H

/**
 * Возвращает значение температуры.
 */
long double calc_temperature(
    long double internal_energy,
    long double chemical_energy
);

long double calc_pressure(
    long double internal_energy,
    long double chemical_energy,
    long double density
);

#endif