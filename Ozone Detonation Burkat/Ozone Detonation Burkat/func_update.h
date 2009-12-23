#ifndef func_update_h__
#define func_update_h__

#include "RealType.h"

void update(RealType *x, RealType *x_center,
            RealType *m, RealType *rho,
            RealType *u, RealType *e,
            RealType *u_energy, RealType *internal_energy,
            RealType *p, bool *shock_wave_front,
            RealType *mf[], RealType *gamma);

void update_rho_u_and_rho_e(RealType *rho, RealType *u,
                            RealType *internal_energy,
                            RealType *rho_u,
                            RealType *rho_e);

#endif // func_update_h__