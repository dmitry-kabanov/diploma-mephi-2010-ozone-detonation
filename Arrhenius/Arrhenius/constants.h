#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <math.h>

const long double
    DX1 = 0.385,                 // ������ ������ � ���, �
    DX2 = 0.385,                 // ������ ������ � ���, �
    DT  = 0.00385,             // ��� �� �������, �
    

    P1 = 42.06,
	P2 =  1.0,               // �������� � ���, ��

    RHO1 = 1 / 0.11443,
    RHO2 = 1.0,               // ��������� � ���
    V2   = 1.0 / RHO2,        // ��������� �����

    U1   = 5.915,
    U2   = 0.0,

    U_PISTON = 3.015,    // �������� ������.

    GAMMA = 1.2,               // ���������� ��������

    R_GAS = 1.0,

    CV = R_GAS / (GAMMA - 1),      // ������������ ��� ���������� ������, �� / �

    T0 = 1.0,                     // ��������� �����������, �

    I0 = CV * T0,               // ��������� ���������� �������

    Q = 50 * R_GAS * T0,        // �������, ������������ � ���������� �������, ��

    D_C_J = sqrt((GAMMA * GAMMA - 1) * Q / 2.0 + GAMMA * P2 / RHO2) + 
        sqrt((GAMMA * GAMMA - 1) * Q / 2.0),              // �������� � ����� �������-����

 
    ACTIVATION_ENERGY = 50.0,

    Z = 82.5,
    K = 0.5;                   // ������������ ��� ����������� ������ ������,
                               // ����� ������� ������� ��������� 
                              // � ������� ������

const int N = 3000,             // ����� �����
L = 3;

const bool EXACT_RIEMANN_SOLVER = true;

const int CONFIG_SHOCK_WAVES                  = 1;
const int CONFIG_SHOCK_AND_RAREFACTION_WAVES  = 2;
const int CONFIG_RAREFACTION_WAVES            = 3;
const int CONFIG_VACUUM                       = 4;

const long double EPSILON = 1e-10;

const int TIMESTEPS = 20000;
const int TIMEDIVISOR = 5000;

#endif