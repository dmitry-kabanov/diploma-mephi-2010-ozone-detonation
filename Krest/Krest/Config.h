#ifndef CONFIG_H
#define CONFIG_H

#include "real_type.h"
#include "Viscosity.h"

class Config
{
    REAL p1;
    REAL p2;
    REAL t1;
    REAL t2;
    REAL rho1;
    REAL rho2;
    REAL u1;
    REAL u2;
    REAL mu1;
    REAL mu2;
    REAL gamma;
    REAL l1;
    unsigned int n1;
    REAL l2;
    unsigned int n2;
    REAL dx;
    REAL dt;
    unsigned int timesteps;
    REAL r;
    Viscosity viscosity;
public:
    Config();
    REAL getP1()
    {
        return p1;
    }
    REAL getP2()
    {
        return p2;
    }
    REAL getT1()
    {
        return t1;
    }
    REAL getT2()
    {
        return t2;
    }
    REAL getRho1()
    {
        return rho1;
    }
    REAL getRho2()
    {
        return rho2;
    }
    REAL getU1()
    {
        return u1;
    }
    REAL getU2()
    {
        return u2;
    }
    REAL getMu1()
    {
        return mu1;
    }
    REAL getMu2()
    {
        return mu2;
    }
    REAL getGamma()
    {
        return gamma;
    }
    REAL getL1()
    {
        return l1;
    }
    unsigned int getN1()
    {
        return n1;
    }
    REAL getL2()
    {
        return l2;
    }
    unsigned int getN2()
    {
        return n2;
    }
    REAL getDx()
    {
        return dx;
    }
    REAL getDt()
    {
        return dt;
    }
    unsigned int getTimesteps()
    {
        return timesteps;
    }
    REAL getR()
    {
        return r;
    }
    //Viscosity getViscosity()
    //{
    //    return Viscosity;
    //}
};

#endif