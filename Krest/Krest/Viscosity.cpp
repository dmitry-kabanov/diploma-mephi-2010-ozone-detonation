#include "Viscosity.h"
#include <iostream>

Viscosity::Viscosity()
{
    std::cout << "Setting viscosity..." << std::endl;
    k = 0.5;
    n1 = 1;
    n2 = 1;
}