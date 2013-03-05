/**
 * @file
 * @author Кабанов Дмитрий <kabanovdmitry@gmail.com>
 * @version
 * 
 * @section DESCRIPTION
 * 
 * Реализация функции in_time_array(int j)
 */
#include "func_in_time_array.h"

bool in_time_array(int j)
{
    const static int SIZE = 7;

    const static int TIME_ARRAY[SIZE] = {
        729, 3925, 5069, 7741, 9955, 14253, 17501
    };

    for (int i = 0; i < SIZE; i++) {
        if (TIME_ARRAY[i] == j) {
            return true;
        }
    }

    return false;
}