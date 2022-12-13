#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "navier-stokes/perfs.h"
#include "fractions.h"
#include "view.h"
#include "lambda2.h"
#include "output_htg.h"

// qcc -source -D_MPI=1 -O2 two_rings.c
// mpicc -Wall -O2 -std=c99 _two_rings.c -lm -lmpi -L$BASILISK/gl -I$BASILISK/gl -lglutils -lfb_osmesa -lGLU -lOSMesa
// mpirun -np 4 ./a.out

#define RADIUS (sqrt(sq(y) + sq(z)))

scalar f[];
int maxlevel = 9;
double ti = 4., ue = 0.008;
double Re = 1750.;
double tend = 120. + 0.1;
int np = 2e5;
/**
## Boundary conditions

Initially, there are two opposing jets at the `left` an `right`
wall. The `top` and `bottom` boundary can leak fluid.
 */
u.n[left]   = dirichlet( f[]  *(1.) * (t <= ti));
u.n[right]  = neumann (0.); //dirichlet(-f[-1]*(1.) * (t <= ti)); quitar 1 anillo
u.n[top]    = neumann (0.);
p[top]      = dirichlet (0.);
pf[top]     = dirichlet (0.);
u.n[bottom] = neumann (0.);
p[bottom]   = dirichlet (0.);
pf[bottom]  = dirichlet (0.);

int main() {
  L0 = 32.;
  X0 = Y0 = Z0 = -L0/2;
  const face vector muc[] = {1./Re, 1./Re, 1./Re};
  mu = muc;
  run();
}
/**
The simulations starts when the jets are triggered. We have a guess at
an initial grid.
 */
event init (t = 0) {
  refine (RADIUS < 2.5 && fabs(x) > 9.*L0/20.   && level < (maxlevel - 1));
  //refine (RADIUS < 1.5 && fabs(x) > 19.5*L0/40. && level < (maxlevel));
  f.refine = f.prolongation = fraction_refine;
  fraction (f, 1. - RADIUS);
  boundary ({f});
}
/**
During the injection phase, the orifice shape is recomputed.
 */
event inject (i++; t <= ti) {
  fraction (f, 1. - RADIUS);
  boundary ({f});
}
/**
## Output

Two movies are generated, plotting an adaptive $\lambda_2$-isosurface
value, a slice of the the cells and the the particles.
 */
event snapshots (t += 0.1) {
  char str[99];
  sprintf (str,"Re = %g", Re);
  double val = -0.0001;
  scalar l2[];
  lambda2 (u, l2);
  foreach()
    l2[] = l2[] < val ? l2[] : nodata;
  stats l2s = statsf(l2);
  foreach()
    l2[] = l2[] < val ? l2[] : 0;
  boundary ({l2});

  char path[]="htg"; // no slash at the end!!
  char prefix[80];
  sprintf(prefix, "data_%03d_%06d", (int) t, i);
  output_htg((scalar *){l2},(vector *){uf}, path, prefix, i, t);
  

  /*int w = 1280;
  double tzoom = 100;
  double fov = t < tzoom ? 15 - t/80 : 15 - t/80 - (t - tzoom)/2.5;
  view (fov = 15 - t/80, width = w, height = 9*w/16, bg = {0.3, 0.4, 0.6},
	theta = 0.6 + 0.15*cos(t/15), phi = 0.6, ty = 0.05);
  isosurface ("l2", min (-l2s.stddev, val));
  translate (y = -L0/4)
    cells (n = {0,1,0});
  draw_string (str, pos = 2, lw = 3);
  save ("movl2.mp4");*/
}

event adapt (i++)
  adapt_wavelet ((scalar*){u}, (double[]){ue, ue, ue}, maxlevel);

event stop (t = tend);
