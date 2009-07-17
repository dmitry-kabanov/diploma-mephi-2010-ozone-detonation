#include <iostream>
#include <fstream>
#include "real_number_type.h"
#include "constants.h"
#include "func_output.h"

void output_row(std::ofstream &output_file, long double *arr, char *row_name)
{
    output_file << row_name << " ";
    for (int i = 0; i <= N; i++) {
        output_file << arr[i];
        if (i != N) {
            output_file << " ";
        }
        else {
            output_file << std::endl;
        }
    }
}

void output_row_int(std::ofstream &output_file, int *arr, char *row_name)
{
    output_file << row_name << " ";
    for (int i = 0; i <= N; i++) {
        output_file << arr[i];
        if (i != N) {
            output_file << " ";
        }
        else {
            output_file << std::endl;
        }
    }
}

void output_row_bool(std::ofstream &output_file, bool *arr, char *row_name)
{
    output_file << row_name << " ";
    for (int i = 0; i <= N; i++) {
        if (arr[i] == true) {
            output_file << "1";
        }
        else {
            output_file << "0";
        }

        if (i != N) {
            output_file << " ";
        }
        else {
            output_file << std::endl;
        }
    }
}


void output(std::ofstream &output_file, int *cells_numbers,
    real_t *x, real_t *p,
    real_t *u, real_t *w,
    real_t *rho, real_t *e,
    real_t *internal_energy, real_t *chemical_energy,
    real_t *p_bound_l, real_t *p_bound_r,
    real_t *rho_bound_l, real_t *rho_bound_r,
    real_t *u_bound_l, real_t *u_bound_r,
    bool *shock_wave_front)
{
    output_file << "Cells X Pressure Velocity MassFraction " <<
        "Density Energy IntEnergy ChemEnergy " << 
        "p_bound_l p_bound_r rho_bound_l rho_bound_r u_bound_l u_bound_r\n";
    for (int i = 1; i < N; i++) {
        output_file << cells_numbers[i] << " " <<
            x[i]               << " " << 
            p[i]               << " " << 
            u[i]               << " " << 
            w[i]               << " " <<
            rho[i]             << " " <<
            e[i]               << " " << 
            internal_energy[i] << " " <<
            chemical_energy[i] << " " <<
            p_bound_l[i]       << " " <<
            p_bound_r[i]       << " " <<
            rho_bound_l[i]     << " " <<
            rho_bound_r[i]     << " " <<
            u_bound_l[i]       << " " <<
            u_bound_r[i]       << " " <<
            std::endl;
        if (shock_wave_front[i] == true && shock_wave_front[i+1] == false) {
            break;
        }
    }
}