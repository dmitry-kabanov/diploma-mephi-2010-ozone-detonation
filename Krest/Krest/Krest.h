#ifndef KREST_H
#define KREST_H

#include "real_type.h"
#include "Config.h"
#include <vector>

class Krest
{
    std::vector<int> cellsNumbers;
    std::vector<REAL> mass;
    std::vector<REAL> x;
    std::vector<REAL> p;
    std::vector<REAL> rho;
    std::vector<REAL> u;
    std::vector<REAL> e;
    std::vector<REAL> iEnergy;
    std::vector<REAL> omega;
    std::vector<REAL> p_bound;
    std::vector<REAL> omega_bound;
    Config config;
    void init();
    void solve();
    void display();
public:
    Krest();
};

#endif