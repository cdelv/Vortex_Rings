// Compile with
// qcc -source -D_MPI=1 -O2 two_rings.c
// mpicc -Wall -O2 -std=c99 _two_rings.c -lm -lmpi -L$BASILISK/gl -I$BASILISK/gl -lglutils -lfb_osmesa -lGLU -lOSMesa
// mpirun -np 4 ./a.out

/*
  This ones are basilisk includes.
  see http://basilisk.fr/Front%20Page
*/
#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "navier-stokes/perfs.h"
#include "fractions.h"
#include "view.h"
#include "lambda2.h"

#define RADIUS (sqrt(sq(y) + sq(z)))

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
  AMR variables:
  - maxlevel: maximun refinementh depth in the tree.
  - np: 
  - ue: 
*/
int maxlevel = 9;
int np = 2e5;
double ue = 0.008;

/*
  Time Variables:
  - ti: Injection time.
  - tend: Finallization time.
*/
double ti = 4.;
double tend = 110. + 0.1;

/*
  Fluid Parameters:
  - 
*/
double Re = 1750.0;

/*
  What is this?
*/
scalar f[];

// Velocity boundary conditions
u.n[left]   = dirichlet( f[]  *(1.) * (t <= ti));
u.n[right]  = neumann (0.0);
u.n[top]    = neumann (0.0);
u.n[bottom] = neumann (0.0);

// Pressure Boundary Conditions
p[top]      = dirichlet (0.0);
pf[top]     = dirichlet (0.0);
p[bottom]   = dirichlet (0.0);
pf[bottom]  = dirichlet (0.0);

int main() {
  init_grid (64);
  size (32.0);
  X0 = Y0 = Z0 = -L0/2;
  const face vector muc[] = {1./Re, 1./Re, 1./Re};
  mu = muc;
  run();
}

/*
  Initial Condition:
  - Everithing starts as 0.
  - We refine the mesh near the inlets.  
*/
event init(t = 0.0) {
  refine (RADIUS < 2.5 && x < -9.0*L0/20.   && level < (maxlevel - 1));
  refine (RADIUS < 1.5 && x < -19.5*L0/40. && level < (maxlevel));
  f.refine = f.prolongation = fraction_refine;
  fraction (f, 1. - RADIUS);
  boundary ({f});
}

/*
  The inflow jet is active, we recompute the inlet shape.
*/
event inject(i++; t <= ti) {
  fraction (f, 1. - RADIUS);
  boundary ({f});
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
event snapshots (t += 0.5) {
  scalar l2[];
  lambda2 (u, l2); // vorticity.

  // Paraview
  char path[]="htg"; // no slash at the end!!
  char prefix[80];
  sprintf(prefix, "data_%03d_%06d", (int) t, i);
  output_htg((scalar *){l2},(vector *){u}, path, prefix, i, t);
}

event stop (t = tend);