// Compile with
// qcc -source -D_MPI=1 -O2 two_rings.c
// mpicc -Wall -O2 -std=c99 _two_rings.c -lm -lmpi -L$BASILISK/gl -I$BASILISK/gl -lglutils -lfb_osmesa -lGLU -lOSMesa
// mpirun -np 4 ./a.out

// Example modified from http://basilisk.fr/sandbox/Antoonvh/two_rings.c

/*
  This ones are basilisk includes.
  see http://basilisk.fr/Front%20Page
*/
#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "navier-stokes/perfs.h"
#include "fractions.h"
#include "lambda2.h"
#include "utils.h"

#define RADIUS (sqrt(sq(y) + sq(z)))

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

  // Injection time
  double ti;
  // File saving parameters
  char *path;
} conf;

double save_dt = 0.5;

// Initialize the values of the Config struct
void init_values(int argc, char *argv[]);

/*
  Time Variables:
  - ti: Injection time.
  - tend: Finalization time.
*/

scalar f[];

// Velocity boundary conditions
u.n[left]   = dirichlet( f[]  *(1.) * (t <= conf.ti));
u.n[right]  = neumann (0.0);
u.n[top]    = neumann (0.0);
u.n[bottom] = neumann (0.0);

// Pressure Boundary Conditions
p[top]      = dirichlet (0.0);
pf[top]     = dirichlet (0.0);
p[bottom]   = dirichlet (0.0);
pf[bottom]  = dirichlet (0.0);

/*
  Function for computing the curl of a vector field.
*/
void curl(const vector v, vector curl);


int main(int argc, char *argv[]) {
  init_values(argc, argv);
  init_grid ((int)pow(2,conf.initial_level));
  size (conf.L);
  X0 = Y0 = Z0 = -L0/2;
  const face vector muc[] = {conf.viscosity, conf.viscosity, conf.viscosity};
  mu = muc;
  run();
}

/*
  Initial Condition:
  - Everything starts as 0.
  - We refine the mesh near the inlets.
*/
event init(t = 0.0) {
  refine (RADIUS < 2.5 && x < -9.0*L0/20.   && level < (conf.max_level - 1));
  refine (RADIUS < 1.5 && x < -19.5*L0/40. && level < (conf.max_level));
  f.refine = f.prolongation = fraction_refine;
  fraction (f, 1. - RADIUS);
  boundary ({f});
}

/*
  When the inflow jet is active we recompute the inlet shape.
*/
event inject(i++; t <= conf.ti) {
  fraction (f, 1. - RADIUS);
  boundary ({f});
}

/*
  Perform the AMR.
  - We use the velocity as the criteria for refinement.
*/
event adapt (i++){
  astats s = adapt_wavelet ((scalar*){u}, (double[]){1.6*conf.threshold, conf.threshold, conf.threshold}, conf.max_level);
  fprintf (stderr, "# Time %3f step %d -> refined %d cells, coarsened %d cells\n", t, i, s.nf, s.nc);
}

/*
  Paraview Output.
*/
event snapshots (t += 0.5) {
  scalar l2[], W_mag[];
  vector W_vec[];
  lambda2 (u, l2);
  curl(u, W_vec);

  foreach ()
    W_mag[] = norm(W_vec);

  // Paraview
  char path[] = "htg"; // no slash at the end!!
  char prefix[80];
  sprintf(prefix, "data_%03d_%06d", (int) t, i);
  output_htg((scalar *){l2, W_mag}, (vector *){u, W_vec}, path, prefix, i, t);
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
  conf.ti = 4.;
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
