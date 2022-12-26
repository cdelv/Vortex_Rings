#ifndef AAA_C_CONNECTOR_H
#define AAA_C_CONNECTOR_H

#ifdef __cplusplus
extern "C" {
#endif

double W_0(double rho, double z);
void compute_U0(double x, double y, double z, double a, double Z0);
double U_x0(void);
double U_y0(void);
double U_z0(void);

#ifdef __cplusplus
}
#endif


#endif
