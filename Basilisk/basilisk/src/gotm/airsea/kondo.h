// Generated automatically by gotm.awk from
// /home/popinet/local/src/code-5.2.1//src/airsea/kondo.F90

extern void kondo_ (
  realtype * sst,
  realtype * airt,
  realtype * u10,
  realtype * v10,
  realtype * precip,
  realtype * evap,
  realtype * taux,
  realtype * tauy,
  realtype * qe,
  realtype * qh
);
static inline void airsea_kondo (
  realtype * sst,
  realtype * airt,
  realtype * u10,
  realtype * v10,
  realtype * precip,
  realtype * evap,
  realtype * taux,
  realtype * tauy,
  realtype * qe,
  realtype * qh) {
  kondo_ (sst, airt, u10, v10, precip, evap, taux, tauy, qe, qh);
}
