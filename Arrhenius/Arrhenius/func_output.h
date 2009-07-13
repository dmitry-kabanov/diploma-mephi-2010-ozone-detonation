#ifndef FUNC_OUTPUT_H
#define FUNC_OUTPUT_H

#include <fstream>

void output_row(std::ofstream &output_file, long double *arr, char *row_name);

void output_row_int(std::ofstream &output_file, int *arr, char *row_name);

void output_row_bool(std::ofstream &output_file, bool *arr, char *row_name);

void output(std::ofstream &, int *cells_numbers, long double *, long double *, long double *, long double *, bool *);

#endif