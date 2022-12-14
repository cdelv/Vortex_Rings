// Compile with
// qcc -source -D_MPI=1 -O2 two_rings.c
// mpicc -Wall -O2 -std=c99 _two_rings.c -lm -lmpi -L$BASILISK/gl -I$BASILISK/gl -lglutils -lfb_osmesa -lGLU -lOSMesa
// mpirun -np 4 ./a.out

/*
  This ones are basilisk includes.
  see http://basilisk.fr/Front%20Page
*/
#include "grid/octree.h"
#include "poisson.h"
#include "navier-stokes/centered.h"
#include "navier-stokes/perfs.h"
#include "fractions.h"
#include "view.h"
#include "lambda2.h"

/*
  This include is for Paraview visualization. We took it from Sander Sandox.
  Thank you very much for your help Maximilian Sander. 
  See http://basilisk.fr/sandbox/sander/output_htg.h.

  This exports htg (Hyper Tree Grid) data from the simulation to Paraview.
  We use this because basilisk uses octrees for the AMR. 
  
  Known HTG Problems
  - HyperTreeGridToDualGrid stops working when advancing time step (Paraview related).
  - Contour Filter does not work on HTG with only one tree (like this exporter creates) (Paraview related).
  - x-z Axis are swapped (3D), x-y Axis are swapped (2D).
*/
#include "output_htg.h"

/*
  Time Variables:
  - tend: Finallization time.
*/
double tend = 20. + 0.1;

/*
  Fluid Parameters: CHANGE FOR VISCOSITY
  - Re: Reynolds number
*/
double Re = 1750.0;

/*
  Vortex Ring Parameters
  -Gamma:
  -aa:
  -RR:
  -LL:
  -ttol:
*/
double Gamma = 1.0;
double aa = 0.2;
double RR = 1.5;
double LL = -8.0;
double ttol = 0.000001;

/*
  AMR variables:
  - maxlevel: maximun refinementh depth in the tree.
  - np: 
  - ue: 
*/
int maxlevel = 8;
int np = 2e5;
double ue = 0.008;

/*
  For vorticity initial conditions.
*/
double W(double r, double z){
  double scale = Gamma/(M_PI*pow(aa,2));
  double T1 = pow(r-RR,2);
  double T2 = pow(z-LL,2);
  return scale*exp(-(T1+T2)/pow(aa,2));
}

double curl_w0_x(double x, double y, double z){
  double r = hypot(x,y)+1.0e-6;
  double T1 = 2.0*(z-LL)/(r*pow(aa,2));
  return W(r,z)*T1*x;
}

double curl_w0_y(double x, double y, double z){
  double r = hypot(x,y)+1.0e-6;
  double T1 = 2.0*(z-LL)/(r*pow(aa,2));
  return W(r,z)*T1*y;
}

double curl_w0_z(double x, double y, double z){
  double r = hypot(x,y)+1.0e-6;
  double T1 = 1.0/r;
  double T2 = z*(r-RR)/pow(aa,2);
  return W(r,z)*(T1-T2);
}

// Velocity boundary conditions
u.n[left]   = neumann (0.0);
u.n[right]  = neumann (0.0);
u.n[top]    = neumann (0.0);
u.n[bottom] = neumann (0.0);
u.n[front]  = neumann (0.0);
u.n[back]   = neumann (0.0);

// Pressure Boundary Conditions
p[top]      = dirichlet (0.0);
pf[top]     = dirichlet (0.0);
p[bottom]   = dirichlet (0.0);
pf[bottom]  = dirichlet (0.0);

// Laplace Boundary Conditions
scalar vx0[], vy0[], vz0[];
scalar bx[], by[], bz[]; // forcing terms

vx0[left] = dirichlet(0.0);
vy0[left] = dirichlet(0.0);
vz0[left] = dirichlet(0.0);

int main() {
  init_grid (128);
  size(32.0);
  X0 = Y0 = Z0 = -L0/2;
  const face vector muc[] = {1./Re, 1./Re, 1./Re};
  mu = muc;
  run();
}

/*
  Initial Condition:
*/
event init(t = 0.0) {
  foreach(){
    bx[] = curl_w0_x(x,y,z);
    by[] = curl_w0_y(x,y,z);
    bz[] = curl_w0_z(x,y,z);
  }

  poisson (vx0, bx, tolerance = ttol);
  poisson (vy0, by, tolerance = ttol);
  poisson (vz0, bz, tolerance = ttol);

  foreach(){
    u.x[] = vx0[];
    u.y[] = vy0[];
    u.z[] = vz0[];
  }
}

/*
  Perform the AMR.
  - We use the velocity as the criteria for refinement. 
*/
event adapt (i++){
  astats s = adapt_wavelet ((scalar*){u}, (double[]){1.6*ue, ue, ue}, maxlevel);
  fprintf (stderr, "# Time %3f step %d -> refined %d cells, coarsened %d cells\n", t, i, s.nf, s.nc);
}

/*
  Paraview Output.
*/
event snapshots (t += 0.1) {
  scalar l2[];
  lambda2 (u, l2); // vorticity.

  // Paraview
  char path[]="htg"; // no slash at the end!!
  char prefix[80];
  sprintf(prefix, "data_%03d_%06d", (int) t, i);
  output_htg((scalar *){l2},(vector *){u}, path, prefix, i, t);
}

event stop (t = tend);