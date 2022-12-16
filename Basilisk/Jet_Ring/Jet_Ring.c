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

/*
  AMR variables:
  - maxlevel: maximum refinement depth in the tree.
  - ue: threshold for the refinement.
*/
int maxlevel = 9;
double ue = 0.008;

/*
  Time Variables:
  - ti: Injection time.
  - tend: Finalization time.
*/
double ti = 4.;
double tend = 50. + 0.1;

/*
  Fluid Parameters:
  - Re: Reynolds number.
*/
double Re = 1750.0;


scalar f[];

// Velocity boundary conditions
u.n[left]   = dirichlet( f[]  * (1.) * (t <= ti));
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


int main() {
  init_grid (64);
  size (32.0);
  X0 = Y0 = Z0 = -L0 / 2;
  const face vector muc[] = {1. / Re, 1. / Re, 1. / Re};
  mu = muc;
  run();
}

/*
  Initial Condition:
  - Everything starts as 0.
  - We refine the mesh near the inlets.
*/
event init(t = 0.0) {
  refine (RADIUS < 2.5 && x < -9.0 * L0 / 20.   && level < (maxlevel - 1));
  refine (RADIUS < 1.5 && x < -19.5 * L0 / 40. && level < (maxlevel));
  f.refine = f.prolongation = fraction_refine;
  fraction (f, 1. - RADIUS);
  boundary ({f});
}

/*
  When the inflow jet is active we recompute the inlet shape.
*/
event inject(i++; t <= ti) {
  fraction (f, 1. - RADIUS);
  boundary ({f});
}

/*
  Perform the AMR.
  - We use the velocity as the criteria for refinement.
*/
event adapt (i++) {
  astats s = adapt_wavelet ((scalar*) {u}, (double[]) {1.6 * ue, ue, ue}, maxlevel);
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
  output_htg((scalar *) {l2, W_mag}, (vector *) {u,W_vec}, path, prefix, i, t);
}

event stop (t = tend);

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