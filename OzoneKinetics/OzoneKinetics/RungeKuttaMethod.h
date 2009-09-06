#ifndef RUNGE_KUTTA_METHOD_H
#define RUNGE_KUTTA_METHOD_H

#include "RealType.h"
#include "Mixture.h"

class RungeKuttaMethod
{
public:
    RungeKuttaMethod();
    ~RungeKuttaMethod();
private:
    void performIntegration();

    RealType rightSideForO3(RealType t, 
                            RealType concOfO, 
                            RealType concOfO3,
                            RealType concOfO2);
    RealType rightSideForO(RealType t, 
                           RealType concOfO, 
                           RealType concOfO3, 
                           RealType concOfO2);
    RealType h;
    Mixture *mixture;
    RealType k1, k2, k3, k4;
    RealType q1, q2, q3, q4;
    int TIME;
    RealType t;
};

#endif