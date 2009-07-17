#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <cmath>

const long double
    DX1 = 0.35,                 // ширина ячейки в КВД, м
    DX2 = 0.35,                 // ширина ячейки в КНД, м
    DT  = 0.0035,             // шаг по времени, с
    

    P1 = 67.345,
	P2 =  1.0,               // давление в КНД, Па

    RHO1 = 9.468,
    RHO2 = 1.0,               // плотность в КНД
    V2   = 1.0 / RHO2,        // начальный объем

    U1   = 7.703,
    U2   = 0.0,

    U_PISTON = 3.015,    // Скорость поршня.

    GAMMA = 1.2,               // показатель адиабаты

    R_GAS = 1.0,

    CV = R_GAS / (GAMMA - 1),      // Теплоемкость при постоянном объеме, Дж / К

    T0 = 1.0,                     // Начальная температура, К

    I0 = CV * T0,               // Начальная внутренняя энергия

    Q = 50 * R_GAS * T0,        // энергия, выделяющаяся в химической реакции, Дж

    D_C_J = sqrt((GAMMA * GAMMA - 1) * Q / 2.0 + GAMMA * P2 / RHO2) + 
        sqrt((GAMMA * GAMMA - 1) * Q / 2.0),              // Скорость в точке Чепмена-Жуге

 
    ACTIVATION_ENERGY = 50.0,

    Z = 206.0,
    K = 0.5;                   // используется для определения ширины ячейки,
                               // левая граница которой совпадает 
                              // с ударной волной

const int N = 3000,             // число ячеек
L = 3;

const bool EXACT_RIEMANN_SOLVER = true;

const int CONFIG_SHOCK_WAVES                  = 1;
const int CONFIG_SHOCK_AND_RAREFACTION_WAVES  = 2;
const int CONFIG_RAREFACTION_WAVES            = 3;
const int CONFIG_VACUUM                       = 4;

const long double EPSILON = 1e-10;

const int TIMESTEPS = 1000;
const int TIMEDIVISOR = 5000;

#endif