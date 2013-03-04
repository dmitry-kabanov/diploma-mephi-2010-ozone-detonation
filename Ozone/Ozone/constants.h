#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <iostream>
#include "RealType.h"

const RealType
    // 1 физическая атмосфера, выраженная в паскалях.
    ONE_ATM = 101325.0,

    // 1 бар, выраженный в паскалях.
    ONE_BAR = 100000.0,
	
	// Комнатная температура, К
	NORMAL_TEMPERATURE = 298.15,

	// Универсальная газовая постоянная, ккал / (моль*К).
	GAS_CONSTANT_KCAL_OVER_MOL_K = 0.001985846,

	R_J_OVER_MOL_K  = 8.31441,
	R_J_OVER_KMOL_K = 8314.41,
	AVOGADRO_NUMBER = 6.022e23,
	K_BOLTZMANN     = 1.38e-23;

#endif