#include "connector.h"
#include "Integral.h"

#ifdef __cplusplus
extern "C" {
#endif

void compute_U0(double x, double y, double z, double a, double Z0){
    compute_V0(x, y, z, a, Z0);
}

double U_x0(void) {
    return get_Vx();
}

double U_y0(void) {
    return get_Vy();
}

double U_z0(void) {
    return get_Vz();
}

#ifdef __cplusplus
}
#endif
