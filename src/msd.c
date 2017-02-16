#include <stdlib.h>
#include "msd.h"
#include "tools.h"

void calculate_msd_one ( double *out, double *x, double *xo, int nlns )
{
  int i;

  double *xdlt = (double *) malloc ( nlns * sizeof(double) );
  double *xodlt = (double *) malloc ( nlns * sizeof(double) );

  subtract_array_number ( xdlt, x, x[0], 1, nlns );
  multiply_array_array ( out, xdlt, xodlt, 1, nlns );

  free ( xdlt );
  free ( xodlt );
}

void calculate_msd_xyz ( double *out, double *x, double *y, double *z, int nlns )
{

  double *tmp = (double *) malloc ( nlns * sizeof (double));

  calculate_msd_one ( out, x, x, nlns );

  calculate_msd_one ( tmp, y, y, nlns );
  add_arrays_inplace ( out, tmp, 1, nlns );
  
  calculate_msd_one ( tmp, z, z, nlns );
  add_arrays_inplace ( out, tmp, 1, nlns );

  free (tmp);

}

void get_conductivity ( double *msd, double *x, double *y, double *z, int ncol, int nlns, int nrestart )
{
}
