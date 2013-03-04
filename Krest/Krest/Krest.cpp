#include "Krest.h"
#include <iostream>
#include <fstream>
#include <cmath>

Krest::Krest()
{
    init();
    solve();
    display();
}

void Krest::init()
{
    unsigned int i, numCells;
    std::cout << "Init..." << std::endl;
    
    numCells = config.getN1() + config.getN2() + 2;

    cellsNumbers.resize(numCells);
    x.resize(numCells);
    p.resize(numCells);
    mass.resize(numCells);
    rho.resize(numCells);
    u.resize(numCells);
    e.resize(numCells);
    iEnergy.resize(numCells);
    omega.resize(numCells);
    p_bound.resize(numCells);
    omega_bound.resize(numCells);

    cellsNumbers[0] = 0;
    x[0]            = 0.0;
    p[0]            = config.getP1();
    rho[0]          = config.getRho1();
    u[0]            = 0.0;
    iEnergy[0]      = 0.0;
    e[0]            = 0.0;
    omega[0]        = 0.0;
    p_bound[0]      = p[0];
    omega_bound[0]  = omega[0];

    for (i = 1; i <= config.getN1(); i++) {
        cellsNumbers[i] = i;
        x[i]            = x[i-1] + config.getDx();
        p[i]            = config.getP1();
        rho[i]          = config.getRho1();
        u[i]            = config.getU1();
        iEnergy[i]      = p[i] / (rho[i] * (config.getGamma() - 1));
        e[i]            = iEnergy[i] + pow(u[i], 2) / 2.0;
        omega[i]        = 0.0;
        p_bound[i]      = p[i];
        omega_bound[i]  = omega[i];
    }

    for (i = config.getN1() + 1; i < numCells; i++) {
        cellsNumbers[i] = i;
        x[i]            = x[i-1] + config.getDx();
        p[i]            = config.getP2();
        rho[i]          = config.getRho2();
        u[i]            = config.getU2();
        iEnergy[i]      = p[i] / (rho[i] * (config.getGamma() - 1));
        e[i]            = iEnergy[i] + pow(u[i], 2);
        omega[i]        = 0.0;
        p_bound[i]      = p[i];
        omega_bound[i]  = omega[i];
    }

    for (i = 0; i < (numCells-1); i++) {
        mass[i] = (x[i+1] - x[i]) * rho[i];
    }
    mass[numCells-1] = mass[numCells-2];
}

void Krest::solve()
{
    std::cout << "Solving..." << std::endl;

    unsigned int numCells = cellsNumbers.size();
    unsigned int i;

    for (unsigned int j = 1; j <= config.getTimesteps(); j++) {
        for (i = 1; i < (numCells-1); i++) {
            u[i] = u[i] - config.getDt() * 
                (p[i] + omega[i] - p[i-1] - omega[i-1]) /
                (0.5 * (mass[i] + mass[i-1]));
            p_bound[i] = (p[i] * mass[i-1] + p[i-1] * mass[i]) / 
                (mass[i-1] + mass[i]);
            omega_bound[i] = (omega[i] * mass[i-1] + omega[i-1] * mass[i]) / 
                (mass[i-1] + mass[i]);
        }

        p_bound[0] = p_bound[1];
        omega_bound[0] = omega_bound[1];

        p_bound[numCells-1] = p_bound[numCells-2];
        omega_bound[numCells-1] = omega_bound[numCells-2];

        u[numCells-1] = u[numCells-2];

        for (i = 1; i < (numCells-1); i++) {
            // Считаем полную энергию.
            e[i] = e[i] - config.getDt() * 
                ((p_bound[i+1] + omega_bound[i+1]) * u[i+1] - 
                (p_bound[i] + omega_bound[i]) * u[i]) /
                mass[i];

            // Считаем внутренную энергию.
            iEnergy[i] = e[i] - pow(u[i] + u[i+1], 2) / 8.0;

            // Считаем новое положение левого узла текущей ячейки.
            x[i] = x[i] + u[i] * config.getDt();
        }

        for (i = 1; i < (numCells-1); i++) {
            // Считаем плотность.
            rho[i] = mass[i] / (x[i+1] - x[i]);

            // Считаем давление.
            p[i] = (config.getGamma() - 1) * iEnergy[i] * rho[i];

            // Считаем вязкость.
            if ((u[i] - u[i+1]) > 0) {
                omega[i] = abs(
                    0.5 * 1 * (u[i] + u[i+1]) * (u[i] - u[i+1]) * rho[i]
                );
            }
            else {
                omega[i] = 0.0;
            }
        }

        p[0] = p[1];
        omega[0] = omega[1];
    }
}

void Krest::display()
{
    std::cout << "Displaying..." << std::endl;

    std::ofstream out("data.txt");
    out.setf(std::ios::scientific, std::ios::floatfield);
    
    for (unsigned int i = 0; i < (cellsNumbers.size()-1); i++) {
        out << cellsNumbers[i] << " " << 
            0.5*(x[i] + x[i+1])  << " " <<
            p[i]                 << " " << 
            rho[i]               << " " <<
            u[i]                 << " " <<
            e[i]                 << " " <<
            iEnergy[i]           << " " <<
            omega[i]             << " " <<
            p_bound[i]           << " " <<
            omega_bound[i]       << " " <<
            std::endl;
    }

    out.close();
}