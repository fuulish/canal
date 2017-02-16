#include <stdio.h>
#include <stdlib.h>
#include "macros.h"
#include "msd.h"
#include "tools.h"

void calculate_msd_one ( double *out, double *x, double *xo, int nlns )
{
  int i;

  double *xdlt = (double *) malloc ( nlns * sizeof(double) );
  double *xodlt = (double *) malloc ( nlns * sizeof(double) );

  subtract_array_number ( xdlt, x, x[0], 1, nlns );
  subtract_array_number ( xodlt, xo, xo[0], 1, nlns );

  multiply_array_array ( out, xdlt, xodlt, 1, nlns );

  free ( xdlt );
  free ( xodlt );
}

void calculate_msd_xyz ( double *out, double *x, double *y, double *z, int nlns, int offset )
{

  double *tmp = (double *) malloc ( nlns * sizeof (double));

  calculate_msd_one ( out, x, x, nlns );

  calculate_msd_one ( tmp, y, y, nlns );
  add_arrays_inplace ( out, tmp, 1, nlns );
  
  calculate_msd_one ( tmp, z, z, nlns );
  add_arrays_inplace ( out, tmp, 1, nlns );

  free (tmp);

}

void calculate_msd_xyz_cross ( double *out, double *xi, double *yi, double *zi, double *xj, double *yj, double *zj, int nlns )
{

  double *tmp = (double *) malloc ( nlns * sizeof (double));

  calculate_msd_one ( out, xi, xj, nlns );

  calculate_msd_one ( tmp, yi, yj, nlns );
  add_arrays_inplace ( out, tmp, 1, nlns );
  
  calculate_msd_one ( tmp, zi, zj, nlns );
  add_arrays_inplace ( out, tmp, 1, nlns );

  free (tmp);

}

void get_qflux ( double *cnd, double *x, double *y, double *z, double *chg, int ncol, int nlns, int nrestart )
{
  int i, j, r;
  int nlns_tmp;
  int offcnt;

  double scale;
  double *xi, *yi, *zi, *xj, *yj, *zj;

  int offset = nlns / nrestart;

  double *tmp = (double *) malloc ( nlns * sizeof(double));
  // double *msd = (double *) calloc ( nlns, sizeof(double));

  // for savety reason assume the worst, i.e., conductivity/output array has not been initialized properly
  for ( i=0; i<nlns; i++ )
    cnd[i] = 0.;

  for ( r=0; r<nrestart; r++ ){

    offcnt = r*offset;
    nlns_tmp = nlns - offcnt;

    printf("OFFCNT: %i NLNS_TMP: %i\n", offcnt, nlns_tmp);

    printf("%i\n", r);
    // different molecules
    for ( i=0; i<ncol; i++ ) {
      for ( j=i; j<ncol; j++ ) {

        if ( i == j )
          scale = 1.;
        else
          scale = 2.;

        // get point to beginning of ncol-th sub-array
        // then use offcnt to start a bit farther down in memory lane
        xi = asub(x, i, nlns, offcnt);
        yi = asub(y, i, nlns, offcnt);
        zi = asub(z, i, nlns, offcnt);

        xj = asub(x, j, nlns, offcnt);
        yj = asub(y, j, nlns, offcnt);
        zj = asub(z, j, nlns, offcnt);

        calculate_msd_xyz_cross ( tmp, xi, yi, zi, xj,  yj, zj, nlns_tmp );
        multiply_array_number_inplace ( tmp, chg[i]*chg[j]*scale, 1, nlns_tmp );

        add_arrays_inplace ( cnd, tmp, 1, nlns_tmp );

      }
    }

  }

  free (tmp);
  // free (msd);
}
