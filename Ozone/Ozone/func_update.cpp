#include "func_update.h"
#include <fstream>
#include <string>
#include "constants.h"
using namespace std;

void update(RealType *x, RealType *x_center,
            RealType *m, RealType *rho,
            RealType *u, RealType *e,
            RealType *u_energy, RealType *internal_energy,
            RealType *p, bool *shock_wave_front,
            RealType *mf[], RealType *gamma)
{
    ifstream file("data_650000.txt");
    string s;
    getline(file, s);
    int length = 23329;

    file >> x[0];

    for (int i = 1; i <= length; i++) {
        file >> s;
        file >> x[i];
        file >> x_center[i];
        file >> p[i];
        file >> u[i];
        file >> rho[i];
        file >> e[i];
        file >> u_energy[i];
        internal_energy[i] = p[i] / ((GAMMA_BEHIND_FRONT - 1) * rho[i]);
        file >> mf[i][0];
        file >> mf[i][1];
        file >> mf[i][2];
        shock_wave_front[i] = false;
        m[i] = (x[i] - x[i-1]) * rho[i];
        gamma[i] = GAMMA_BEHIND_FRONT;
    }
    internal_energy[length] = p[length] / ((GAMMA_AHEAD_FRONT - 1) * rho[length]);
    gamma[length] = GAMMA_AHEAD_FRONT;
    shock_wave_front[length-1] = true;
    shock_wave_front[length]   = true;
    
    file.close();
}

void update_rho_u_and_rho_e(RealType *rho, RealType *u,
                            RealType *internal_energy,
                            RealType *rho_u, RealType *rho_e)
{
    for (int i = 1; i <= N; i++) {
        rho_u[i] = rho[i] * u[i];
        rho_e[i] = rho[i] * internal_energy[i];
    }
}