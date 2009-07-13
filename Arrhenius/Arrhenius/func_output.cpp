#include <iostream>
#include <fstream>
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


void output(std::ofstream &output_file, 
    int *cells_numbers,
    long double *x, 
    long double *p,
    long double *u,
    long double *w,
    bool *shock_wave_front)
{
    for (int i = 0; i < N; i++) {
        output_file << cells_numbers[i] << " " <<
            x[i] << " " << 
            p[i] << " " << 
            u[i] << " " << 
            w[i] << " " <<
            std::endl;
        if (shock_wave_front[i] == true && shock_wave_front[i+1] == false) {
            break;
        }
    }
}