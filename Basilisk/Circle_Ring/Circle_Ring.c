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
#include "lambda2.h"
#include "utils.h"

/*
  This include is for Paraview visualization. We took it from Sander Sandbox.
  Thank you very much for your help, Maximilian Sander.
  See http://basilisk.fr/sandbox/sander/output_htg.h.

  This exports htg (Hyper Tree Grid) data from the simulation to Paraview.
  We use this because Basilisk uses octrees for the AMR.

  Known HTG Problems
  - HyperTreeGridToDualGrid stops working when advancing the time step (Paraview related).
  - Contour Filter does not work on HTG with only one tree (like this exporter creates) (Paraview related).
  - x-z Axis is swapped (3D), x-y Axis is swapped (2D).
*/
#include "output_htg.h"

struct Config {
  // Finalization time.
  double tmax;

  // Box dimensions.
  double L;

  // Ring Parameters, radius, Initial position in the z-direction,
  //Sigma of the gaussian (thickness of the ring), Magnitude of the vortex, and
  // Tolerance for the Laplace solver.
  // R=1 and Gamma=1 always.
  double R, Z0, a, Gamma, tol;

  // Fluid parameters
  double Re, viscosity;

  // AMR parameters
  int max_level, initial_level;
  double threshold;

  // File saving parameters
  char *path;
} conf;

double save_dt = 0.1;

// Initialize the values of the Config struct
void init_values(int argc, char *argv[]);

/*
  For vorticity initial conditions.
  - W: Scale factor of the initial vorticity curl.
  - curl_w0_x: X component of the initial vorticity curl.
  - curl_w0_y: Y component of the initial vorticity curl.
  - curl_w0_z: Z component of the initial vorticity curl.
*/
double W(double r, double z);
double curl_w0_x(double x, double y, double z);
double curl_w0_y(double x, double y, double z);
double curl_w0_z(double x, double y, double z);

/*
  Velocity boundary conditions:
  The order of the boundaries in 3D is
  { right,  left,   top, bottom, front, back  }
  { X max, X min, Y max,  Y min, Z max, Z min }
*/
u.n[right]  = neumann(0.0);
u.n[left]   = neumann(0.0);
u.n[top]    = neumann(0.0);
u.n[bottom] = neumann(0.0);
u.n[front]  = neumann(0.0);
u.n[back]   = neumann(0.0);

/*
  Pressure Boundary Conditions
*/
p[front]  = neumann(0.0);
p[back]   =  neumann(0.0);
p[top]    = dirichlet(0.0);
p[bottom] = dirichlet(0.0);
p[right]  = dirichlet(0.0);
p[left]   = dirichlet(0.0);

/*
  Laplace Boundary Conditions
  The Dirichlet condition is necessary for convergence.
  It's the farthest boundary from the ring.
*/
scalar vx0[], vy0[], vz0[];
scalar bx[], by[], bz[]; // forcing terms

vx0[back] = dirichlet(0.0);
vy0[back] = dirichlet(0.0);
vz0[back] = dirichlet(0.0);

/*
  Function for computing the curl of a vector field.
*/
void curl(const vector v, vector curl);

int main(int argc, char *argv[]) {
  init_values(argc, argv);
  init_grid ((int)pow(2,conf.initial_level));
  size(conf.L);
  X0 = Y0 = Z0 = -L0 / 2;
  const face vector muc[] = {conf.viscosity, conf.viscosity, conf.viscosity};
  mu = muc;
  run();
}

/*
  Initial Condition: we initialize the velocity from a vorticity field.
*/
event init(t = 0.0) {
  // Restart simulation from previus run
  //if (!restore (file = "restart")){
  // Initialize the forcing terms of the Poisson equation
  foreach () {
    bx[] = curl_w0_x(x, y, z);
    by[] = curl_w0_y(x, y, z);
    bz[] = curl_w0_z(x, y, z);
  }

  // Solve Poisson equation
  // \nabla( \alpha \nabla(a)) - \lamba a = b
  poisson (vx0, bx, tolerance = conf.tol);
  poisson (vy0, by, tolerance = conf.tol);
  poisson (vz0, bz, tolerance = conf.tol);

  // Initialize velocity
  foreach () {
    u.x[] = vx0[];
    u.y[] = vy0[];
    u.z[] = vz0[];
  }
  //}
}

/*
  Perform the AMR.
  - We use the velocity as the criteria for refinement.
*/
event adapt (i++) {
  astats s = adapt_wavelet ((scalar*) {u}, (double[]) {conf.threshold , conf.threshold, conf.threshold}, conf.max_level);
  fprintf (stderr, "# Time %3f step %d -> refined %d cells, coarsened %d cells\n", t, i, s.nf, s.nc);
}

/*
  Paraview Output:
  - W_vec: Vectorial vorticity field (W)
  - W_mag: Magnitude of the vorticity
  - W_r: W*r where r = sqrt(x^2 + y^2)
  - W_r2: W*r^2 where r = sqrt(x^2 + y^2)
  - W_Z: W*(z-L) where L is the initial position of the ring
  - W_Z2: W*(z-L)^2 where L is the initial position of the ring
  - PID: Process ID of the CPU in charge of that cell
  - LEVEL: Level of refinement. The number of cells is 2^LEVEL
  CHECK FOR FILE EXISTANCE
*/
event snapshots (t += save_dt) {
  vector W_vec[], W_r[], W_r2[], W_Z[], W_Z2[];
  scalar W_mag[], PID[], LEVEL[];

  curl(u, W_vec);

  foreach () {
    double r = hypot(x, y);
    W_mag[] = norm(W_vec);
    PID[] = pid();
    LEVEL[] = level;
    foreach_dimension() {
      W_r.x[] = r * W_vec.x[];
      W_r2.x[] = pow(r , 2) * W_vec.x[];
      W_Z.x[] = (z - conf.Z0) * W_vec.x[];
      W_Z2.x[] = pow(z - conf.Z0, 2) * W_vec.x[];
    }
  }

  // Paraview output
  char path[] = "htg"; // no slash at the end!!
  char prefix[80];
  sprintf(prefix, "data_%03d_%06d", (int) t, i);
  output_htg((scalar *) {W_mag, PID, LEVEL}, (vector *) {uf, W_vec, W_r, W_r2, W_Z, W_Z2}, conf.path, prefix, i, t);
}

event stop (t = conf.tmax);

void init_values(int argc, char *argv[]) {
  char *ptr;
  conf.tmax = strtod(argv[1], &ptr);
  conf.L = strtod(argv[2], &ptr);
  conf.R = 1.0;
  conf.Z0 = strtod(argv[3], &ptr);
  conf.a = strtod(argv[4], &ptr);
  conf.Gamma = 1.0;
  conf.tol = strtod(argv[5], &ptr);
  conf.Re = strtod(argv[6], &ptr);
  conf.viscosity = 1.0 / conf.Re;
  conf.initial_level = strtod(argv[7], &ptr);
  conf.max_level = strtod(argv[8], &ptr);
  conf.threshold = strtod(argv[9], &ptr);
  save_dt = strtod(argv[10], &ptr);
  conf.path = argv[11];
  if (pid() == 0) {
    printf("tmax = %g \n", conf.tmax);
    printf("L = %g \n", conf.L);
    printf("Z0 = %g \n", conf.Z0);
    printf("a = %g \n", conf.a);
    printf("tol = %g \n", conf.tol);
    printf("Re = %g \n", conf.Re);
    printf("initial_level = %d \n", conf.initial_level);
    printf("max_level = %d \n", conf.max_level);
    printf("threshold = %g \n", conf.threshold);
    printf("save_dt = %g \n", save_dt);
    printf("path = %s \n", conf.path);
  }
}

double W(double r, double z) {
  double scale = conf.Gamma / (M_PI * pow(conf.a, 2));
  double T1 = pow(r - conf.R, 2);
  double T2 = pow(z - conf.Z0, 2);
  return scale * exp(-(T1 + T2) / pow(conf.a, 2));
}

double curl_w0_x(double x, double y, double z) {
  double r = hypot(x, y) + 1.0e-6;
  double T1 = 2.0 * (z - conf.Z0) / (r * pow(conf.a, 2));
  return W(r, z) * T1 * x;
}

double curl_w0_y(double x, double y, double z) {
  double r = hypot(x, y) + 1.0e-6;
  double T1 = 2.0 * (z - conf.Z0) / (r * pow(conf.a, 2));
  return W(r, z) * T1 * y;
}

double curl_w0_z(double x, double y, double z) {
  double r = hypot(x, y) + 1.0e-6;
  double T1 = 1.0 / r;
  double T2 = z * (r - conf.R) / pow(conf.a, 2);
  return W(r, z) * (T1 - T2);
}

void curl(const vector v, vector curl) {
  vector dvx[], dvy[], dvz[];
  scalar vx[], vy[], vz[];

  foreach () {
    vx[] = v.x[];
    vy[] = v.y[];
    vz[] = v.z[];
  }

  gradients ({vx}, {dvx});
  gradients ({vy}, {dvy});
  gradients ({vz}, {dvz});

  foreach () {
    curl.x[] = dvz.y[] - dvy.z[];
    curl.y[] = dvx.z[] - dvz.x[];
    curl.z[] = dvy.x[] - dvx.y[];
  }
}