/**
 * @file
 * @author Кабанов Дмитрий <kabanovdmitry@gmail.com
 * @version
 *
 * @section DESCRIPTION
 *
 * Главный файл проекта Ozone Detonation.
 */
#include "main.h"
#include <cstdio>
#include <cmath>
#include <ctime>
#include <iostream>
#include <fstream>
#include "constants.h"
#include "func_init.h"
#include "func_delta.h"
#include "func_riemann_solver.h"
#include "func_output.h"
#include "RealType.h"
#include "RungeKuttaMethod.h"
#include "Piston.h"
#include "Stiffl1.h"
using namespace std;

int main(/*int argc, char *argv[]*/)
{
    int i, j;

    // Массив номеров ячеек.
    int cells_numbers[N+1];

    RealType shock_wave_velocity, delta_mass;
	int nSpecies = 3;
	double vf1[4];
	vf1[0] = 0.0;
	vf1[1] = 0.0;
	vf1[2] = 0.0;
	vf1[3] = 0.0;
    RungeKuttaMethod kinetics(
		nSpecies + 1,
		vf1,
		0.0,
		0.0,
		1.0e-13,
		"Substances.txt",
        "Reactions.txt",
        "MoleFractions.txt");
    Piston aPiston;
	
    // Объявление основных величин.
    RealType *x        = new RealType[N+1];
    RealType *x_center = new RealType[N+1];
    RealType *m        = new RealType[N+1];
    RealType *rho      = new RealType[N+1];
    RealType *u        = new RealType[N+1];
    RealType *p        = new RealType[N+1];

    // Объемные доли компонентов смеси.
    RealType **volumeFractions  = new RealType*[N+1];
    for (i = 0; i <= N; i++) {
        volumeFractions[i] = new RealType[kinetics.getMixture()->nSubstances];
    }

    // Объявление термодинамических величин.
    RealType *internal_energy = new RealType[N+1];
    RealType *u_energy        = new RealType[N+1];
    RealType *e               = new RealType[N+1];
    RealType *temperature     = new RealType[N+1];
    RealType *gamma           = new RealType[N+1];

    // Объявление дополнительных величин.
    RealType *p_contact     = new RealType[N+1];
    RealType *u_contact     = new RealType[N+1];
    RealType *delta_impulse = new RealType[N+1];
    RealType *delta_energy  = new RealType[N+1];

    // Объявление величин для расчета на границах ячеек.    
    RealType *rho_bound_r = new RealType[N+1];
    RealType *rho_bound_l = new RealType[N+1];
    RealType *u_bound_r   = new RealType[N+1];
    RealType *u_bound_l   = new RealType[N+1];
    RealType *e_bound_r   = new RealType[N+1];
    RealType *e_bound_l   = new RealType[N+1];
    RealType *p_bound_r   = new RealType[N+1];
    RealType *p_bound_l   = new RealType[N+1];

    RealType *rho_tg_left     = new RealType[N+1];
    RealType *rho_tg_right    = new RealType[N+1];
    RealType *rho_tg          = new RealType[N+1];
    RealType *rho_u_tg_left   = new RealType[N+1];
    RealType *rho_u_tg_right  = new RealType[N+1];
    RealType *rho_u_tg        = new RealType[N+1];
    RealType *rho_e_tg_left   = new RealType[N+1];
    RealType *rho_e_tg_right  = new RealType[N+1];
    RealType *rho_e_tg        = new RealType[N+1];
    
    RealType *rho_u_bound_r   = new RealType[N+1];
    RealType *rho_u_bound_l   = new RealType[N+1];
    RealType *rho_e_bound_r   = new RealType[N+1];
    RealType *rho_e_bound_l   = new RealType[N+1];

    RealType *rho_u = new RealType[N+1];
    RealType *rho_e = new RealType[N+1];

    //RealType piston[2];

    // Массив хранит данные о фронте ударной волны.
    // Если значение элемента равно true - на границе этой ячейки 
    // находится фронт ударной волны.
    // Если значение элемента равно false - то нет.
    bool shock_wave_front[N+1];

    //std::ofstream out_front("sw_front.txt");
    //out_front.setf(std::ios::fixed, std::ios::floatfield);
    //out_front.precision(6);

    // Инициализация массивов.
    for (i = 0; i <= N; i++) {
        cells_numbers[i] = i;
    }

    init_parameters(x, x_center, m, rho, u, e, p, rho_u, shock_wave_front);
    init_gamma(gamma);
    init_volume_fractions(volumeFractions);
    init_thermodynamic_parameters(
        internal_energy, e, 
        p, rho, u,
        u_energy,
        volumeFractions,
        kinetics
    );
    init_rho_e(internal_energy, rho, rho_e);
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

    //RealType vf[3];
    //vf[0] = 0.0;
    //vf[1] = 0.0;
    //vf[2] = 100.0;

    //ofstream outFile("T_Cp.txt");
    //for (RealType t = 298; t < 3000; t += 1) {
    //    kinetics.getMixture()->fillUpMixtureWithTAndP(t, ONE_ATM, vf);
    //    outFile << t << " " << kinetics.getMixture()->calculateMixtureCp(t) << endl;
    //}
    //outFile.close();
    //exit(0);
    RealType vf[2];
    //vf[2] = 0.0;
    //vf[1] = 92.8312;
    //vf[0] = 7.16689;
    vf[2] = 100.0;
    vf[1] = 0.0;
    vf[0] = 0.0;
    clock_t start;
    clock_t finish;
    double workingTime;
 //   
 //////   kinetics.getMixture()->fillUpMixture(3297630, 3.58197, vf);
	start = clock();
    kinetics.getMixture()->setStateWithTPX(1914.18, 6543160, vf);
    cout << "oldT = " << kinetics.getMixture()->getOldTemperature() << endl;
	cout << "oldP = " << kinetics.getMixture()->calculateOldPressure() << endl;
    cout << "H(O3) = " << kinetics.getMixture()->calculateSubstanceEnthalpy(2, 1914.18) << endl;
    cout << "S(O3) = " << kinetics.getMixture()->calculateSubstanceEntropy(2, 1914.18) << endl;
    cout << "G(O3) = " << kinetics.getMixture()->calculateSubstanceGibbsEnergy(2, 1914.18) << endl;
	kinetics.performIntegration(1.0);
    kinetics.updateMoleFractions(vf);

    cout << "T = " << kinetics.getMixture()->getTemperature() << endl;
    cout << "P = " << kinetics.getMixture()->calculatePressure() << endl;
    cout << "X(O)  = " << vf[0]  << endl;
    cout << "X(O2) = " << vf[1] << endl;
    cout << "X(O3) = " << vf[2] << endl;
	finish = clock();

	workingTime = (double) (finish - start) / CLOCKS_PER_SEC;

	cout << "Calculations done in " << workingTime << " s." << endl;

	vf[2] = 100.0;
    vf[1] = 0.0;
    vf[0] = 0.0;
    
    kinetics.getMixture()->setStateWithURhoX(3297630, 3.58197, vf);
	start = clock();
    kinetics.getMixture()->setStateWithTPX(3445.82, 3326220, vf);
    cout << "oldT = " << kinetics.getMixture()->getOldTemperature() << endl;
	cout << "oldP = " << kinetics.getMixture()->calculateOldPressure() << endl;
	kinetics.performIntegration(1.0);
    kinetics.updateMoleFractions(vf);

    cout << "T = " << kinetics.getMixture()->getTemperature() << endl;
    cout << "P = " << kinetics.getMixture()->calculatePressure() << endl;
    cout << "X(O)  = " << vf[0]  << endl;
    cout << "X(O2) = " << vf[1] << endl;
    cout << "X(O3) = " << vf[2] << endl;
	finish = clock();

	workingTime = (double) (finish - start) / CLOCKS_PER_SEC;

	cout << "Calculations done in " << workingTime << " s." << endl;
    exit(0);

    // ********************************************************************
    // Пишем начальные значения в файл.
    outputAsCSVFile(0, cells_numbers, 
        x_center, p, 
        u, rho, 
        e, u_energy, 
        volumeFractions, shock_wave_front
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
            p_bound_r[i] = rho_bound_r[i] * (gamma[i] - 1) * e_bound_r[i];
            p_bound_l[i] = rho_bound_l[i] * (gamma[i] - 1) * e_bound_l[i];
            
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
                    u[i], u[i+1], GAMMA_FRONT
                );
                u_contact[i] = calc_u_contact(p_contact[i],
                    p[i], p[i+1],
                    rho[i], rho[i+1],
                    u[i], u[i+1], GAMMA_FRONT
                );
                shock_wave_velocity = calc_shock_wave_velocity(
                    p_contact[i], u_contact[i],
                    p[i], p[i+1],
                    rho[i], rho[i+1],
                    u[i], u[i+1], GAMMA_FRONT
                );
                if (shock_wave_velocity == 0) {
                    cout << "Shock wave velocity is equals to 0." << endl;
                    exit(1);
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
                    u_bound_r[i], u[i+1], GAMMA_BEHIND_FRONT
                );
                u_contact[i] = calc_u_contact(p_contact[i],
                    p_bound_r[i], p[i+1],
                    rho_bound_r[i], rho[i+1],
                    u_bound_r[i], u[i+1], GAMMA_BEHIND_FRONT
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
                    u_bound_r[i], u_bound_l[i+1], GAMMA_BEHIND_FRONT);
                u_contact[i] = calc_u_contact(p_contact[i], p_bound_r[i], p_bound_l[i+1],
                    rho_bound_r[i], rho_bound_l[i+1],
                    u_bound_r[i], u_bound_l[i+1], GAMMA_BEHIND_FRONT);
            }
        }

        /*
         * На левой границе у нас находится поршень.
         * Считаем параметры для поршня.
         */
        // Левая граница области расчета движется со скоростью поршня.
        p_contact[0] = p[1];
        u_contact[0] = aPiston.calculateVelocity(volumeFractions[1][2]);

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
                x_center[i] = (x[i-1] + x[i]) / 2.0;
                u[i] = (m[i] * u[i] + delta_impulse[i]) / (m[i] - delta_mass);
                e[i] = (m[i] * e[i] + delta_energy[i]) / (m[i] - delta_mass);
                m[i] = m[i] - delta_mass;
                // cout << shock_wave_velocity << endl;
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
            u_energy[i] = e[i] - u[i] * u[i] / 2.0;
            kinetics.getMixture()->setStateWithURhoX(u_energy[i], rho[i], volumeFractions[i]);
            kinetics.performIntegration(DT);
            kinetics.updateMoleFractions(volumeFractions[i]);
            p[i] = kinetics.getPressure();
            //cout << j << ": (" <<
            //    i << ") " <<
            //    volumeFractions[i][0] << " " << 
            //    volumeFractions[i][1] << " " << 
            //    volumeFractions[i][2] << " u=" << 
            //    u[i] << " rho=" << rho[i] << endl;
            //cout << "U=" << u_energy[i] << " " << " oldT=" <<
            //    kinetics.mixture->getOldTemperature() << " T=" <<
            //    kinetics.mixture->getTemperature() << " p=" << 
            //    p[i] << endl;
            rho_u[i] = rho[i] * u[i];
            gamma[i] = p[i] / (rho[i] * u_energy[i]) + 1;
            rho_e[i] = rho[i] * p[i] / ((gamma[i] - 1) * rho[i]);
        }

        // Записываем зависимость давления на фронте ударной волны
        // от времени в файл out_front.
        //for (i = 1; i < N; i++) {
        //    /**
        //     * TODO: убрать магическое число 100.
        //     */
        //    if (i % 500 == 0) {
        //        if (shock_wave_front[i] == true 
        //            && shock_wave_front[i+1] == true)
        //        {
        //            out_front << (j * DT) << " " << p[i] << "\n";
        //            break;
        //        }
        //    }
        //}

        if (j % TIMEDIVISOR == 0) {
            cout << "j = " << j << endl;
            cout << "D = " << shock_wave_velocity << endl << endl;
            outputAsCSVFile(j, cells_numbers, 
                            x_center, p, 
                            u, rho,  
                            e, u_energy, 
                            volumeFractions, shock_wave_front);
        }

        // Модифицируем ячейки, связанные с фронтом ударной волны.
        for (i = 1; i < (N-1); i++) {
            if (shock_wave_front[i] == true && shock_wave_front[i+1] == true) {
                if ((x[i+1] - x[i]) <= K * (x[i] - x[i-1])) {
                    RealType delta_x_right = x[i+2] - x[i+1];
                    RealType delta_x_left  = x[i+1] - x[i];
                    RealType delta_x_sum = delta_x_left + delta_x_right;
                    RealType temp_rho = rho[i+2];
                    RealType temp_x   = x[i];
                    
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
                    
                    //p[i+2] = rho[i+2] * (gamma[i+2] - 1) * 
                    //    (e[i+2] - pow(u[i+2], 2) / 2.0);

                    // Делим ячейку слева от фронта ударной волны пополам.
                    p[i+1]   = p[i];
                    rho[i+1] = rho[i];
                    u[i+1]   = u[i];
                    e[i+1]   = e[i];
                    m[i+1]   = m[i] / 2.0;
                    m[i]     = m[i] / 2.0;
                    
                    // Модифицируем координаты границ и центров ячеек.
                    x[i]   = (x[i] + x[i-1]) / 2.0;
                    x[i+1] = temp_x;

                    //x_center[i]   = (x[i-1] + x[i])   / 2.0;
                    //x_center[i+1] = (x[i]   + x[i+1]) / 2.0;
                    //x_center[i+2] = (x[i+1] + x[i+2]) / 2.0;

                    /*gamma[i]   = GAMMA_BEHIND_FRONT;*/
                    gamma[i+1] = gamma[i];

                    shock_wave_front[i]   = 0;
                    shock_wave_front[i+1] = true;
                    shock_wave_front[i+2] = true;
                }
                break;
            }
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

    for (int i = 0; i <= N; ++i) {
        delete [] volumeFractions[i];
    }
    delete [] volumeFractions;

    delete [] gamma;
    delete [] u_energy;

    delete [] p_contact;
    delete [] u_contact;
    delete [] delta_impulse;
    delete [] delta_energy;

    delete [] rho_u;
    delete [] rho_e;

    delete [] rho_bound_r;
    delete [] rho_bound_l;
    delete [] u_bound_r;
    delete [] u_bound_l;
    delete [] e_bound_r;
    delete [] e_bound_l;
    delete [] p_bound_r;
    delete [] p_bound_l;
    delete [] rho_u_bound_r;
    delete [] rho_u_bound_l;
    delete [] rho_e_bound_r;
    delete [] rho_e_bound_l;

    delete [] rho_tg_left;
    delete [] rho_tg_right;
    delete [] rho_tg;
    delete [] rho_u_tg_left;
    delete [] rho_u_tg_right;
    delete [] rho_u_tg;
    delete [] rho_e_tg_left;
    delete [] rho_e_tg_right;
    delete [] rho_e_tg;
}