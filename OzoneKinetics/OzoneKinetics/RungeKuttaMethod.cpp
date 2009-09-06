#include "RungeKuttaMethod.h"
#include <iostream>

RungeKuttaMethod::RungeKuttaMethod()
{
    h = 0.01;
    t = 0.0;
    TIME = 100;
    mixture = new Mixture();
    mixture->readFileOfReactions("Reactions.txt");
    mixture->readFileOfSubstances("Substances.txt");
}

RungeKuttaMethod::~RungeKuttaMethod()
{
    delete mixture;
}

void RungeKuttaMethod::performIntegration()
{
    /**
     * Концентрация молекул O.
     */
    RealType concOfO = 0.0;
    /**
     * Концентрация молекул O2.
     */
    RealType concOfO2 = 0.0;
    /**
     * Концентрация молекул O3.
     */
    RealType concOfO3 = 1.0;
    /**
     * Шаг интегрирования по времени.
     */
    RealType h = 0.01;

    // Производим интегрирование в цикле.
    for (int i = 0; i < TIME; i++) {
        k1 = rightSideForO(t, concOfO, concOfO3, concOfO2);
        q1 = rightSideForO3(t, concOfO, concOfO3, concOfO2);

        k2 = rightSideForO(t + h / 2.0, 
            concOfO + h * k1 / 2.0, 
            concOfO3 + h * q1 / 2.0,
            concOfO2);
        q2 = rightSideForO3(t + h / 2.0, 
            concOfO + h * k1 / 2.0, 
            concOfO3 + h * q1 / 2.0,
            concOfO2);

        k3 = rightSideForO(t + h / 2.0, 
            concOfO + h * k2 / 2.0, 
            concOfO3 + h * q2 / 2.0,
            concOfO2);
        q3 = rightSideForO3(t + h / 2.0, 
            concOfO + h * k2 / 2.0, 
            concOfO3 + h * q2 / 2.0,
            concOfO2);

        k4 = rightSideForO(t + h, 
            concOfO + h * k3, 
            concOfO3 + h * q3, 
            concOfO2);
        q4 = rightSideForO3(t + h, 
            concOfO + h * k3, 
            concOfO3 + h * q3, 
            concOfO2);

        concOfO  = concOfO + h * 
            (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;
        concOfO3 = concOfO3 + h * 
            (q1 + 2 * q2 + 2 * q3 + q4) / 6.0;
        concOfO2 = concOfO2 + 
            (-3 * rightSideForO3(t, concOfO, concOfO3, concOfO2) -
            rightSideForO(t, concOfO, concOfO3, concOfO2)) * 0.5 * h;
        
        std::cout << concOfO << " " << concOfO2 << " " << concOfO3 << std::endl;
    }
}

RealType RungeKuttaMethod::rightSideForO(RealType t, 
                                         RealType concOfO, 
                                         RealType concOfO3, 
                                         RealType concOfO2)
{
    //RealType concOfM = concOfO + concOfO2 + concOfO3;

    //RealType m1;
    //RealType m2;
    //RealType m3;

    //m1 = -k1 * concOfO * concOfO3 + k1r * confOfO2 * confOfO2;
    //m2 = k2 * concOfM * concOfO3 - k2r * concOfO2 * concOfO * concOfM;
    //m3 = 2 * k3 * concOfM * concOfO2 - 
    //    2 * k3r * concOfO * concOfO * concOfM;

    //return m1 + m2 + m3;
    return 0;
}

RealType RungeKuttaMethod::rightSideForO3(RealType t, 
                                          RealType concOfO, 
                                          RealType concOfO3, 
                                          RealType concOfO2)
{
    //RealType concOfM = concOfO + concOfO2 + concOfO3;

    //RealType m1;
    //RealType m2;

    //m1 = -k1 * concOfO * concOfO3 + k1r * concOfO2 * concOfO2;
    //m2 = -k2 * concOfM * concOfO3 + k2r * concOfO2 * concOfO3 * concOfM;

    //return m1 + m2;
    return 0;
}