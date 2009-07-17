#ifndef FUNC_OUTPUT_H
#define FUNC_OUTPUT_H

#include <fstream>
#include "real_number_type.h"

void output_row(std::ofstream &output_file, long double *arr, char *row_name);

void output_row_int(std::ofstream &output_file, int *arr, char *row_name);

void output_row_bool(std::ofstream &output_file, bool *arr, char *row_name);

void output(std::ofstream &output_file, int *cells_numbers,
    real_t *x, real_t *p,
    real_t *u, real_t *w,
    real_t *rho, real_t *e,
    real_t *internal_energy, real_t *chemical_energy,
    real_t *p_bound_l, real_t *p_bound_r,
    real_t *rho_bound_l, real_t *rho_bound_r,
    real_t *u_bound_l, real_t *u_bound_r,
    bool *shock_wave_front
);

#endif