// same as poiseuille.c but with periodic boundary conditions

#include "grid/multigrid.h"
#include "navier-stokes/centered.h"

int main() {
  origin (-0.5, 0);
  periodic (top);
  dimensions (nx = 1);
  
  stokes = true;
  TOLERANCE = 1e-5;
  
  u.t[left] = dirichlet(0);
  u.t[right] = dirichlet(0);
  
  for (N = 8; N <= 64; N *= 2)
    run();
}

scalar un[];

event init (t = 0) {
  // we also check for hydrostatic balance
  const face vector g[] = {1.,1.};
  a = g;
  const face vector muc[] = {1.,1.};
  mu = muc;
  foreach()
    un[] = u.y[];
}

event logfile (t += 0.1; i <= 100) {
  double du = change (u.y, un);
  if (i > 0 && du < 1e-6)
    return 1; /* stop */
}

event profile (t = end) {
  printf ("\n");
  foreach()
    fprintf (stdout, "%g %g %g %g %g\n", x, y, u.x[], u.y[], p[]);
  scalar e[];
  foreach()
    e[] = u.y[] - 0.5*(0.25 - x*x);
  norm n = normf (e);
  fprintf (stderr, "%d %g %g %g\n", N, n.avg, n.rms, n.max);
}
