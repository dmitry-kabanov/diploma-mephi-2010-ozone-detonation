#include <math.h>
#include "func_delta.h"

long double calc_delta(long double f_left, long double f, long double f_right,
                  long double x_center_left, long double x_center, long double x_center_right,
                  long double x_bound_left, long double x_bound_right)
{
    long double f_tg_left, f_tg_right, f_tg;

    f_tg_left = (f - f_left) / (x_center - x_center_left);
    f_tg_right = (f_right - f) / (x_center_right - x_center);

    if (abs(f_tg_right) > abs(f_tg_left)) {
        f_tg = f_tg_left;
    }
    else {
        f_tg = f_tg_right;
    }

    return f_tg * (x_bound_right - x_bound_left) / 2.0L;
}
