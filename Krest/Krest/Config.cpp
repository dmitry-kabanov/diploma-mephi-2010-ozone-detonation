#include "Config.h"
#include <iostream>

Config::Config()
{
    std::cout << "Configuring..." << std::endl;
    p1 = 892035.3678;
    p2 = 1.0e5;
    t1 = 298;
    t2 = 298;
    rho1 = 4.27542699;
    rho2 = 1.17;
    u1   = 701.2135039;
    u2   = 0.0;
    mu1 = 28.9;
    mu2 = 28.9;
    gamma = 1.4;
    l1    = 0.5;
    n1    = 50;
    l2    = 1.0;
    n2    = 100;
    dx    = 1.0;
    dt    = 0.00005;
    timesteps = 1500;
    r     = 8314.0;
}