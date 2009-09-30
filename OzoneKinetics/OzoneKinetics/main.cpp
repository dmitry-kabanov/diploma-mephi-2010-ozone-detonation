#include <iostream>
#include "RungeKuttaMethod.h"

int main()
{
    // Давление в Па.
    RealType p0 = 64 * 101325;
    // Температура в К.
    RealType t0 = 1200.00;

    RungeKuttaMethod rungeKuttaMethod(p0, t0,
                                      "Substances.txt",
                                      "Reactions.txt",
                                      "VolumeFractions.txt"
    );
}