#ifndef VISCOSITY_H
#define VISCOSITY_H

#include "real_type.h"

class Viscosity
{
    REAL n1;
    REAL n2;
    REAL k;
public:
    Viscosity();
    REAL getN1();
    REAL getN2();
};

#endif