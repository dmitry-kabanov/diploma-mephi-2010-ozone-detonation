/**
* @file
*
* @author  Кабанов Дмитрий <kabanovdmitry@gmail.com>
* @version %I%
*
* @section DESCRIPTION
*
* Класс поршня.
*/
#ifndef Piston_h__
#define Piston_h__

#include "RealType.h"

class Piston
{
public:
    Piston();
    ~Piston();
    RealType calculateVelocity(RealType f);
private:
    void readFileOfPiston(const char* filename);
    int nColumns;
    RealType* pressures;
    RealType* fractions;
};

#endif // Piston_h__