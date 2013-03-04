#ifndef FUNC_CHEMICAL_KINETICS_H
#define FUNC_CHEMICAL_KINETICS_H

/**
 * Возвращает новое значение массовой доли непрореагировавшего вещества.
 * Для расчета используется закон Аррениуса.
 */
long double calc_mass_fraction_of_initial_substance(
    long double w, long double t
);

#endif