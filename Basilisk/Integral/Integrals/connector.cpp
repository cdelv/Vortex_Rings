#include "connector.h"
#include "Integral.h"

#ifdef __cplusplus
extern "C" {
#endif

double U_x0(double x, double y, double z, double Z0, double a) {
    return Integral_x(x, y, z, Z0, a);
}

double U_y0(double x, double y, double z, double Z0, double a) {
    return Integral_y(x, y, z, Z0, a);
}

double U_z0(double x, double y, double z, double Z0, double a) {
    return Integral_z(x, y, z, Z0, a);
}

#ifdef __cplusplus
}
#endif
