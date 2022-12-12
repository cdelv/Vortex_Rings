#include "utils.h"

int main ()
{
  init_grid (64);

  scalar s[];
  foreach()
    s[] = x + y;

  // statsf() uses reduction operations
  stats stat = statsf (s);
  fprintf (qerr, "%g %g %g\n", stat.min, stat.sum, stat.max);

  // Array reduction
  #define arr_size 10
  int cells[arr_size] = {0};
  foreach (reduction(+:cells[:arr_size])) 
    cells[(int)(10*fabs(x))]++;

  for (int i = 0; i < arr_size; i++) 
    fprintf (qerr, "%d ", cells[i]);
  fputc ('\n', qerr);

  // Coord reduction
  coord p = {0};
  foreach (reduction(+:p))
    foreach_dimension()
      p.x++;
  fprintf (qerr, "%g %g\n", p.x, p.y);
}
