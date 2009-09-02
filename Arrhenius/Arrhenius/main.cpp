#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include "real_number_type.h"
#include "constants.h"
#include "func_init.h"
#include "func_delta.h"
#include "func_riemann_solver.h"
#include "func_thermodynamics.h"
#include "func_chemical_kinetics.h"
#include "func_piston.h"
#include "func_output.h"
#include "func_in_time_array.h"

int main(int argc, char *argv[])
{
    int i, j;

    // Массив номеров ячеек.
    int cells_numbers[N+1];

    real_t shock_wave_velocity, delta_mass;
	
    // Объявление основных величин.
    real_t *x        = new real_t[N+1];
    real_t *x_center = new real_t[N+1];
    real_t *m        = new real_t[N+1];
    real_t *rho      = new real_t[N+1];
    real_t *u        = new real_t[N+1];
    real_t *p        = new real_t[N+1];

    // Объявление термодинамических величин.
    real_t *internal_energy = new real_t[N+1];
    real_t *e               = new real_t[N+1];
    real_t *temperature     = new real_t[N+1];

    // Объявление химических величин.
    real_t *w               = new real_t[N+1];
    real_t *chemical_energy = new real_t[N+1];

    // Объявление дополнительных величин.
    real_t *p_contact     = new real_t[N+1];
    real_t *u_contact     = new real_t[N+1];
    real_t *delta_impulse = new real_t[N+1];
    real_t *delta_energy  = new real_t[N+1];

    // Объявление величин для расчета на границах ячеек.    
    real_t *rho_bound_r = new real_t[N+1];
    real_t *rho_bound_l = new real_t[N+1];
    real_t *u_bound_r   = new real_t[N+1];
    real_t *u_bound_l   = new real_t[N+1];
    real_t *e_bound_r   = new real_t[N+1];
    real_t *e_bound_l   = new real_t[N+1];
    real_t *p_bound_r   = new real_t[N+1];
    real_t *p_bound_l   = new real_t[N+1];

    real_t *rho_tg_left     = new real_t[N+1];
    real_t *rho_tg_right    = new real_t[N+1];
    real_t *rho_tg          = new real_t[N+1];
    real_t *rho_u_tg_left   = new real_t[N+1];
    real_t *rho_u_tg_right  = new real_t[N+1];
    real_t *rho_u_tg        = new real_t[N+1];
    real_t *rho_e_tg_left   = new real_t[N+1];
    real_t *rho_e_tg_right  = new real_t[N+1];
    real_t *rho_e_tg        = new real_t[N+1];
    
    real_t *rho_u_bound_r   = new real_t[N+1];
    real_t *rho_u_bound_l   = new real_t[N+1];
    real_t *rho_e_bound_r   = new real_t[N+1];
    real_t *rho_e_bound_l   = new real_t[N+1];

    real_t *rho_u = new real_t[N+1];
    real_t *rho_e = new real_t[N+1];

    real_t piston[2];

    // Массив хранит данные о фронте ударной волны.
    // Если значение элемента равно true - на границе этой ячейки 
    // находится фронт ударной волны.
    // Если значение элемента равно false - то нет.
    bool shock_wave_front[N+1];

    int num_digits;
    char str_time[32];
    char filename[8] = "data_";
    char file_ext[5] = ".txt";
    char fullname[64];

    std::ofstream out_front("sw_front.txt");
    out_front.setf(std::ios::fixed, std::ios::floatfield);
    out_front.precision(6);

    // Инициализация массивов.
    for (i = 0; i <= N; i++) {
        cells_numbers[i] = i;
    }

    init_parameters(x, x_center, m, rho, u, e, p, rho_u, shock_wave_front);
    init_thermodynamic_parameters(
        internal_energy, e, 
        temperature, p, 
        rho, u
    );
    init_chemical_parameters(w, chemical_energy);
    init_rho_e(internal_energy, chemical_energy, rho, rho_e);
    init_additional_parameters(p_contact, u_contact, 
                               delta_impulse, delta_energy,
                               p, u
    );
    init_boundary_parameters(rho_bound_r, rho_bound_l,
                             rho_u_bound_r, rho_u_bound_l,
                             rho_e_bound_r, rho_e_bound_l,
                             u_bound_r, u_bound_l,
                             e_bound_r, e_bound_l,
                             p_bound_r, p_bound_l
    );
    init_tangents(
        rho_tg_left, rho_tg_right, rho_tg,
        rho_u_tg_left, rho_u_tg_right, rho_u_tg,
        rho_e_tg_left, rho_e_tg_right, rho_e_tg
    );

    /*
     * Основной цикл по времени
     */
    for (j = 1; j <= TIMESTEPS; j++) {
        for (i = 2; i < N; i++) {
            // Считаем тангенс для rho на левой и правой границах ячейки.
            rho_tg_left[i] = (rho[i] - rho[i-1]) / (x_center[i] - x_center[i-1]);
            rho_tg_right[i] = (rho[i+1] - rho[i]) / (x_center[i+1] - x_center[i]);
            if (abs(rho_tg_left[i]) > abs(rho_tg_right[i])) {
                rho_tg[i] = rho_tg_right[i];
            }
            else {
                rho_tg[i] = rho_tg_left[i];
            }
            
            // Считаем тангенс для rho*u на левой и правой границах ячейки.
            rho_u_tg_left[i]  = (rho_u[i] - rho_u[i-1]) / 
                (x_center[i] - x_center[i-1]);
            rho_u_tg_right[i] = (rho_u[i+1] - rho_u[i]) / 
                (x_center[i+1] - x_center[i]);
            if (abs(rho_u_tg_left[i]) > abs(rho_u_tg_right[i])) {
                rho_u_tg[i] = rho_u_tg_right[i];
            }
            else {
                rho_u_tg[i] = rho_u_tg_left[i];
            }

            // Считаем тангенс для rho*e на левой и правой границах ячейки.
            rho_e_tg_left[i]  = (rho_e[i] - rho_e[i-1]) / (x_center[i] - x_center[i-1]);
            rho_e_tg_right[i] = (rho_e[i+1] - rho_e[i]) / (x_center[i+1] - x_center[i]);
            if (abs(rho_e_tg_left[i]) > abs(rho_e_tg_right[i])) {
                rho_e_tg[i] = rho_e_tg_right[i];
            }
            else {
                rho_e_tg[i] = rho_e_tg_left[i];
            }

            if (shock_wave_front[i] == true && shock_wave_front[i+1] == false)
            {
                break;
            }
        }

        rho_tg_left[1]    = rho_tg_left[2];
        rho_tg_right[1]   = rho_tg_right[2];
        rho_tg[1]         = rho_tg[2];

        rho_u_tg_left[1]  = rho_u_tg_left[2];
        rho_u_tg_right[1] = rho_u_tg_right[2];
        rho_u_tg[1]       = rho_u_tg[2];

        rho_e_tg_left[1]  = rho_e_tg_left[2];
        rho_e_tg_right[1] = rho_e_tg_right[2];
        rho_e_tg[1]       = rho_e_tg[2];

        rho_tg_left[N]    = rho_tg_right[N-1];
        rho_tg_right[N]   = rho_tg_right[N-1];
        rho_tg[N]         = rho_tg[N-1];

        rho_u_tg_left[N]  = rho_u_tg_left[N-1];
        rho_u_tg_right[N] = rho_u_tg_right[N-1];
        rho_u_tg[N]       = rho_u_tg[N-1];

        rho_e_tg_left[N]  = rho_e_tg_left[N-1];
        rho_e_tg_right[N] = rho_e_tg_right[N-1];
        rho_e_tg[N]       = rho_e_tg[N-1];
        
        for (i = 1; i < N; i++) {
            rho_bound_r[i] = rho[i] + (x[i] - x[i-1]) * rho_tg[i] / 2.0;
            rho_bound_l[i] = rho[i] - (x[i] - x[i-1]) * rho_tg[i] / 2.0;            
            
            rho_u_bound_r[i] = rho_u[i] + (x[i] - x[i-1]) * rho_u_tg[i] / 2.0;
            rho_u_bound_l[i] = rho_u[i] - (x[i] - x[i-1]) * rho_u_tg[i] / 2.0;
            
            rho_e_bound_r[i] = rho_e[i] + (x[i] - x[i-1]) * rho_e_tg[i] / 2.0;
            rho_e_bound_l[i] = rho_e[i] - (x[i] - x[i-1]) * rho_e_tg[i] / 2.0;

            // Считаем u на левой и правой границах ячейки.
            u_bound_r[i] = rho_u_bound_r[i] / rho_bound_r[i];
            u_bound_l[i] = rho_u_bound_l[i] / rho_bound_l[i];

            // Считаем e на левой и правой границах ячейки.
            e_bound_r[i] = rho_e_bound_r[i] / rho_bound_r[i];
            e_bound_l[i] = rho_e_bound_l[i] / rho_bound_l[i];

            // Считаем p на левой и правой границах ячейки.
            p_bound_r[i] = rho_bound_r[i] * (GAMMA - 1) * e_bound_r[i];
            p_bound_l[i] = rho_bound_l[i] * (GAMMA - 1) * e_bound_l[i];
            
            if (shock_wave_front[i] == true && shock_wave_front[i+1] == false)
            {
                break;
            }
        }

        // Решаем задачу Римана.
        for (i = 1; i < N; i++) {
            /*
             * Фронт ударной волны, которую мы выделяем, 
             * находится на правой границе i-й ячейки.
             * Задача Римана решается без линейной интерполяции
             * параметров в ячейках i и i+1.
             * После решения задачи Римана между этими ячейками
             * выходим из цикла, так как перед ударной волной ничего
             * не изменилось.
             */
            if (shock_wave_front[i] == true && shock_wave_front[i+1] == true) {
                p_contact[i] = calc_p_contact(p[i], p[i+1],
                    rho[i], rho[i+1],
                    u[i], u[i+1]
                );
                u_contact[i] = calc_u_contact(p_contact[i],
                    p[i], p[i+1],
                    rho[i], rho[i+1],
                    u[i], u[i+1]
                );
                shock_wave_velocity = calc_shock_wave_velocity(
                    p_contact[i], u_contact[i],
                    p[i], p[i+1],
                    rho[i], rho[i+1],
                    u[i], u[i+1]
                );
                if (shock_wave_velocity == 0) {
                    std::cout << "Volny razrezheniya v obe storony" << std::endl;
                    std::cout << "j = " << j << ", i = " << i << std::endl;
                    std::cout << std::endl;
                }
                break;
            }
            /*
             * Фронт ударной волны, которую мы выделяем, 
             * находится на правой границе i+1-й ячейки.
             */
            else if (shock_wave_front[i] == false && shock_wave_front[i+1] == true) {
                p_contact[i] = calc_p_contact(p_bound_r[i], p[i+1],
                    rho_bound_r[i], rho[i+1],
                    u_bound_r[i], u[i+1]
                );
                u_contact[i] = calc_u_contact(p_contact[i],
                    p_bound_r[i], p[i+1],
                    rho_bound_r[i], rho[i+1],
                    u_bound_r[i], u[i+1]
                );
            }
            /*
             * Решение задачи Римана для соседних ячеек, ни одна из которых
             * не содержит ударную волну, которую мы выделяем.
             * Задача Римана решается через линейно интерполированные
             * параметры в ячейках (как и положено в методе Годунова-Колгана).
             */
            else {
                p_contact[i] = calc_p_contact(p_bound_r[i], p_bound_l[i+1],
                    rho_bound_r[i], rho_bound_l[i+1],
                    u_bound_r[i], u_bound_l[i+1]);
                u_contact[i] = calc_u_contact(p_contact[i], p_bound_r[i], p_bound_l[i+1],
                    rho_bound_r[i], rho_bound_l[i+1],
                    u_bound_r[i], u_bound_l[i+1]);
            }
        }

        p_contact[N] = p_contact[N-1];
        u_contact[N] = u_contact[N-1];

        /*
         * На левой границе у нас находится поршень.
         * Считаем параметры для поршня.
         */
        calc_piston(piston, w[1], sqrt(F) * D_C_J);
        p_contact[0] = p[1];
        // Левая граница области расчета движется со скоростью поршня.
        u_contact[0] = piston[1];

        /**
         * Считаем импульс и энергию в ячейке.
         */
        for (i = 1; i < N; i++) {
            /*
             * Фронт ударной волны, которую мы выделяем, 
             * находится на правой границе i-й ячейки.
             */
            if (shock_wave_front[i] == true && shock_wave_front[i+1] == true) {
                delta_mass = (shock_wave_velocity - u[i+1]) * rho[i+1] * DT;

                delta_impulse[i] = p_contact[i-1] * DT + 
                    delta_mass * u[i+1] - p[i+1] * DT;
                delta_energy[i]  = p_contact[i-1] * u_contact[i-1] * DT + 
                    delta_mass * e[i+1] - p[i+1] * u[i+1] * DT;
            }
            /*
             * Фронт ударной волны, которую мы выделяем, 
             * находится на левой границе i-й ячейки.
             * После расчета delta_impulse[i] и delta_energy[i]
             * выходим из цикла, так как перед ударной волной 
             * ничего не изменилось.
             */
            else if (shock_wave_front[i] == true 
                && shock_wave_front[i+1] == false) {
                
                delta_impulse[i] = - delta_mass * u[i] + 
                    p[i] * DT - 
                    p_contact[i] * DT;
                delta_energy[i]  = - delta_mass * e[i] + 
                    + p[i] * u[i] * DT -
                    p_contact[i] * u_contact[i] * DT;
                break;
            }
            /*
             * i-я ячейка не содержит ударную волну,
             * которую мы выделяем.
             * Считаем импульс и энергию обычным способом.
             */
            else {
                delta_impulse[i] = p_contact[i-1] * DT -
                    p_contact[i] * DT;
                delta_energy[i]  = p_contact[i-1] * u_contact[i-1] * DT - 
                    p_contact[i] * u_contact[i] * DT;
            }
        }

        delta_impulse[N] = delta_impulse[N-1];
        delta_energy[N]  = delta_energy[N-1];

        // Левая граница области расчета движется со скоростью поршня.
        x[0] = x[0] + u_contact[0] * DT;

        // Вычисляем основные величины.
        for (i = 1; i <= N; i++) {
            if (shock_wave_front[i] == true && shock_wave_front[i+1] == true) {
                // Фронт ударной волны, которую мы выделяем, 
                // находится на правой границе i-й ячейки.
                // Считаем изменение массы в i-й ячейке.
                // Правая граница i-й ячейки движется со скоростью ударной волны.
                x[i] = x[i] + shock_wave_velocity * DT;
                u[i] = (m[i] * u[i] + delta_impulse[i]) / (m[i] + delta_mass);
                e[i] = (m[i] * e[i] + delta_energy[i]) / (m[i] + delta_mass);
                m[i] = m[i] + delta_mass;
            }
            else if (shock_wave_front[i] == true && 
                shock_wave_front[i+1] == false) {
                // i-я ячейка содержит выделяемую 
                // ударную волну на левой границе.
                // Её правая граница движется со скоростью контактного разрыва.
                x[i] = x[i] + u_contact[i] * DT;
                u[i] = (m[i] * u[i] + delta_impulse[i]) / (m[i] - delta_mass);
                e[i] = (m[i] * e[i] + delta_energy[i]) / (m[i] - delta_mass);
                m[i] = m[i] - delta_mass;
                break;
            }
            else {
                // Ячейка не содержит выделяемую ударную волну.
                x[i] = x[i] + u_contact[i] * DT;
                u[i] = (m[i] * u[i] + delta_impulse[i]) / m[i];
                e[i] = (m[i] * e[i] + delta_energy[i]) / m[i];
            }
            x_center[i] = (x[i-1] + x[i]) / 2.0;
            rho[i] = m[i] / (x[i] - x[i-1]);
            internal_energy[i] = e[i] - pow(u[i], 2) / 2.0;
            w[i] = calc_mass_fraction_of_initial_substance(
                w[i], temperature[i]
            );
            chemical_energy[i] = (1.0 - w[i]) * Q;
            temperature[i] = calc_temperature(
                internal_energy[i], chemical_energy[i]
            );
            p[i] = rho[i] * (GAMMA - 1) * 
                (internal_energy[i] + chemical_energy[i]);
            rho_u[i] = rho[i] * u[i];
            rho_e[i] = rho[i] * (internal_energy[i] + chemical_energy[i]);
        }

        for (i = 1; i < (N-1); i++) {
            if (shock_wave_front[i] == true && shock_wave_front[i+1] == true) {
                if ((x[i+1] - x[i]) <= K * (x[i] - x[i-1])) {
                    long double delta_x_right = x[i+2] - x[i+1];
                    long double delta_x_left  = x[i+1] - x[i];
                    long double delta_x_sum = delta_x_left + delta_x_right;
                    long double temp_rho = rho[i+2];
                    long double temp_x   = x[i];
                    
                    rho[i+2] = 1 / delta_x_sum * 
                        (rho[i+2] * delta_x_right + 
                        rho[i+1] * delta_x_left);

                    m[i+2] = rho[i+2] * delta_x_sum;

                    u[i+2] = 1 / (rho[i+2] * delta_x_sum) * 
                        (temp_rho * delta_x_right * u[i+2] + 
                         rho[i+1] * delta_x_left * u[i+1]);

                    e[i+2] = 1 / (rho[i+2] * delta_x_sum) * 
                        (temp_rho * delta_x_right * e[i+2] + 
                         rho[i+1] * delta_x_left * e[i+1]);

                    p[i+2] = rho[i+2] * (GAMMA - 1) * 
                        (e[i+2] - pow(u[i+2], 2) / 2.0);

                    // Делим ячейку слева от фронта ударной волны пополам.
                    p[i+1]   = p[i];
                    rho[i+1] = rho[i];
                    u[i+1]   = u[i];
                    e[i+1]   = e[i];
                    m[i+1]   = m[i] / 2.0;
                    m[i]     = m[i] / 2.0;
                    x[i]     = (x[i] + x[i-1]) / 2.0;
                    x[i+1]   = temp_x;

                    shock_wave_front[i]   = 0;
                    shock_wave_front[i+1] = true;
                    shock_wave_front[i+2] = true;
                    break;
                }
            }
        }

        for (i = 1; i < N; i++) {
            /**
             * TODO: убрать магическое число 100.
             */
            if (i % 10 == 0) {
                if (shock_wave_front[i] == true 
                    && shock_wave_front[i+1] == true)
                {
                    out_front << (j * DT) << " " << p[i] << "\n";
                    break;
                }
            }
        }

        if (j % TIMEDIVISOR == 0 || in_time_array(j)) {
            num_digits = sprintf_s(str_time, "%d", j);
            strcpy_s(fullname, filename);
            strcat_s(fullname, str_time);
            strcat_s(fullname, file_ext);
            
            std::ofstream out(fullname);
            out.setf(std::ios::fixed, std::ios::floatfield);
            out.precision(6);
            
            output(
                out, cells_numbers, 
                x, p, 
                u, w, 
                rho, e,
                internal_energy, chemical_energy,
                p_bound_l, p_bound_r,
                rho_bound_l, rho_bound_r,
                u_bound_l, u_bound_r,
                shock_wave_front
            );
            
            out.close();
        }
    }

    /*
     * Освобождаем память.
     */
    delete [] x;
    delete [] x_center;
    delete [] m;
    delete [] rho;
    delete [] u;
    delete [] p;

    delete [] internal_energy;
    delete [] e;
    delete [] temperature;

    delete [] w;
    delete [] chemical_energy;

    delete [] p_contact;
    delete [] u_contact;
    delete [] delta_impulse;
    delete [] delta_energy;

    delete [] rho_bound_r;
    delete [] rho_bound_l;
    delete [] u_bound_r;
    delete [] u_bound_l;
    delete [] e_bound_r;
    delete [] e_bound_l;
    delete [] p_bound_r;
    delete [] p_bound_l;

    delete [] rho_tg_left;
    delete [] rho_tg_right;
    delete [] rho_tg;
    delete [] rho_u_tg_left;
    delete [] rho_u_tg_right;
    delete [] rho_u_tg;
    delete [] rho_e_tg_left;
    delete [] rho_e_tg_right;
    delete [] rho_e_tg;
    
    delete [] rho_u_bound_r;
    delete [] rho_u_bound_l;
    delete [] rho_e_bound_r;
    delete [] rho_e_bound_l;

    delete [] rho_u;
    delete [] rho_e;
}