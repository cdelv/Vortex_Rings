#include "fractions.h"
#include "curvature.h"
#include "utils.h"

int main()
{
  origin (-0.5, -0.5);
  init_grid (64);
  mask (x > 0.35 && y > 0.25 ? top : none);
  scalar c[];
  fraction (c, min(sq(x + 0.1) + sq(y + 0.22) - sq(0.18),
		   sq(x - 0.1) + sq(y - 0.22) - sq(0.18)));
  c.refine = c.prolongation = fraction_refine;
  while (adapt_wavelet ({c}, (double[]){1e-4}, 6).nf)
    ;

  trash ({c});
  fraction (c, min(sq(x + 0.1) + sq(y + 0.22) - sq(0.18),
		   sq(x - 0.1) + sq(y - 0.22) - sq(0.18)));
  output_cells (stdout);
  output_facets (c, stdout);
  vector h[];
  heights (c, h);
  scalar kappa[];
  curvature (c, kappa);
  foreach() {
    if (h.x[] != nodata)
      fprintf (stderr, "%g %g %g %g hx\n", x, y, height(h.x[]), Delta);
    if (h.y[] != nodata)
      fprintf (stderr, "%g %g %g %g hy\n", x, y, height(h.y[]), Delta);
    if (kappa[] != nodata)
      fprintf (stderr, "%g %g %g kappa\n", x, y, kappa[]);
  }
}
