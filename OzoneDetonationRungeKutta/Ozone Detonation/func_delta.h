#ifndef FUNC_DELTA_H
#define FUNC_DELTA_H

long double calc_delta(
                       long double f_left, long double f, long double f_right,
                       long double x_left, long double x, long double x_right,
                       long double x_bound_l, long double x_bound_r
);

#endif