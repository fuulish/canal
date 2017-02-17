#include <math.h>
#include "mol.h"

double get_distance_periodic_sqr ( double x1, double y1, double z1, double x2, double y2, double z2, double cell )
{
  double dst[3];
  double val = 0.;

  dst[0] = x1 - x2;
  dst[1] = y1 - y2;
  dst[2] = z1 - z2;

  dst[0] -= cell * round ( dst[0] / cell );
  dst[1] -= cell * round ( dst[1] / cell );
  dst[2] -= cell * round ( dst[2] / cell );

  val += dst[0]*dst[0];
  val += dst[1]*dst[1];
  val += dst[2]*dst[2];

  return val;
}

double get_distance_periodic ( double x1, double y1, double z1, double x2, double y2, double z2, double cell )
{
  double dsqr;
  dsqr = get_distance_periodic_sqr ( x1, y1, z1, x2, y2, z2, cell );

  return sqrt(dsqr);
}
