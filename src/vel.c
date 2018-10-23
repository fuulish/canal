/*
canal is C program to calculate electrical conductivities from
molecular simulation data.
Copyright 2017 Frank Uhlig (uhlig.frank@gmail.com)

This file is part of canal.

canal is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

canal is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with canal.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "macros.h"
#include "vel.h"
#include "mol.h"
#include "tools.h"
#ifdef DEBUG
#include "io.h"
#endif

double calculate_dot ( double *xi, double *yi, double *zi, double *xj, double *yj, double *zj, int nlns, int nskip )
{

  int i;
  // double *tmp = (double *) malloc ( nlns * sizeof (double));
  double tmp[nlns];
  double fnl[nlns];

  double out = 0.;

  multiply_array_array ( fnl, xi, xj, 1, nlns, nskip );

  multiply_array_array ( tmp, yi, yj, 1, nlns, nskip );
  add_arrays_inplace ( fnl, tmp, 1, nlns );

  multiply_array_array ( tmp, zi, zj, 1, nlns, nskip );
  add_arrays_inplace ( fnl, tmp, 1, nlns );

  for ( i=0; i<nlns; i++ )
    out += fnl[i];

  return out;

  // free (tmp);

}

double calculate_nrm ( double *xi, double *yi, double *zi, int nlns, int nskip )
{
  double dot = calculate_dot ( xi, yi, zi, xi, yi, zi, nlns, nskip );

  return sqrt(dot);
}

void calculate_velocities ( double *vx, double *vy, double *vz, double *xi, double *yi, double *zi, int nlns, int timestep )
{

  double *ptr_a, *ptr_b, *ptr_v;

  // FUOD| make this more general by using functions for numerical differentiation, to be created in tools.c
  // x component
  ptr_v = &(vx[1]);
  ptr_a = &(xi[2]);
  ptr_b = &(xi[0]);

  subtract_array_array ( ptr_v, ptr_a, ptr_b, 1, nlns-2 );
  divide_array_number_inplace ( ptr_v, 2.*timestep, 1, nlns-1 );

  vx[0] = ( xi[1] - xi[0] ) / timestep;
  vx[nlns-1] = ( xi[nlns-1] - xi[nlns-2] ) / timestep;

  // y component
  ptr_v = &(vy[1]);
  ptr_a = &(yi[2]);
  ptr_b = &(yi[0]);

  subtract_array_array ( ptr_v, ptr_a, ptr_b, 1, nlns-2 );
  divide_array_number_inplace ( ptr_v, 2.*timestep, 1, nlns-1 );

  vy[0] = ( yi[1] - yi[0] ) / timestep;
  vy[nlns-1] = ( yi[nlns-1] - yi[nlns-2] ) / timestep;

  // z component
  ptr_v = &(vz[1]);
  ptr_a = &(zi[2]);
  ptr_b = &(zi[0]);

  subtract_array_array ( ptr_v, ptr_a, ptr_b, 1, nlns-2 );
  divide_array_number_inplace ( ptr_v, 2.*timestep, 1, nlns-1 );

  vz[0] = ( zi[1] - zi[0] ) / timestep;
  vz[nlns-1] = ( zi[nlns-1] - zi[nlns-2] ) / timestep;

}

//FUDO| we shouldn't need neinst, because it'll always be 1
void get_vflux_locl ( double *neinaa, double *neincc, double *cnd_cc, double *cnd_ac, double *cnd_aa, double *x, double *y, double *z, double *chg, int ncol, int nlns, int nrestart, double dr, double rstart, int rnum, double *cell, double timestep, int nskip )
{
  int i, j, n;
  int nlns_tmp;
  int rndx;

  double *vx, *vy, *vz;
  double *ptr_vxi, *ptr_vyi, *ptr_vzi, *ptr_vxj, *ptr_vyj, *ptr_vzj;
  double *xi, *yi, *zi, *xj, *yj, *zj;

  double scale;

  int arrlen;

  //FUDO| each array has nlns time entries, and ncol molecule entries
  arrlen = nlns * ncol;

  vx = (double *) malloc(arrlen * sizeof(double));
  vy = (double *) malloc(arrlen * sizeof(double));
  vz = (double *) malloc(arrlen * sizeof(double));

  //FUDO| need to transform x,y,z data to vx,vy,vz data
  //FUDO| we need to this with respect to each ncol subarrays, otherwise we're using the beginning and ending of the wrong timeseries/molecues

  for ( n=0; n<ncol; n++)
    calculate_velocities ( asub(vx, nlns, n, 0), asub(vy, nlns, n, 0), asub(vz, nlns, n, 0), asub(x, nlns, n, 0), asub(y, nlns, n, 0), asub(z, nlns, n, 0), nlns, timestep );

  int nrmlen;

  //FUDO| for each kind of correlation only a distance dependence, we will average over all molecules/kinds
  nrmlen = rnum;

  double  *nrm_na = calloc(nrmlen, sizeof(double));
  double  *nrm_nc = calloc(nrmlen, sizeof(double));
  double  *nrm_cc = calloc(nrmlen, sizeof(double));
  double  *nrm_ac = calloc(nrmlen, sizeof(double));
  double  *nrm_aa = calloc(nrmlen, sizeof(double));

  for ( i=0; i<rnum; i++ ) {
    cnd_cc[i] = 0.;
    cnd_ac[i] = 0.;
    cnd_aa[i] = 0.;
    neinaa[i] = 0.;
    neincc[i] = 0.;
  }

  for ( n=0; n<nlns; n++ ) {

    for ( i=0; i<ncol; i++ ) {
      for ( j=i; j<ncol; j++ ) {

        if ( i == j )
          scale = 1.;
        else
          scale = 2.;

        ptr_vxi = asub(vx, nlns, i, n);
        ptr_vyi = asub(vy, nlns, i, n);
        ptr_vzi = asub(vz, nlns, i, n);
                               
        ptr_vxj = asub(vx, nlns, j, n);
        ptr_vyj = asub(vy, nlns, j, n);
        ptr_vzj = asub(vz, nlns, j, n);

        xi = asub(x, nlns, i, n);
        yi = asub(y, nlns, i, n);
        zi = asub(z, nlns, i, n);
                         
        xj = asub(x, nlns, j, n);
        yj = asub(y, nlns, j, n);
        zj = asub(z, nlns, j, n);

        nlns_tmp = 1;
        double dot = calculate_dot ( ptr_vxi, ptr_vyi, ptr_vzi, ptr_vxj, ptr_vyj, ptr_vzj, nlns_tmp, nskip );
        double nrm_i = calculate_nrm ( ptr_vxi, ptr_vyi, ptr_vzi, nlns_tmp, nskip );
        double nrm_j = calculate_nrm ( ptr_vxj, ptr_vyj, ptr_vzj, nlns_tmp, nskip );

        dot /= (nrm_i * nrm_j);

        double dst = get_distance_periodic(xi[0], yi[0], zi[0], xj[0], yj[0], zj[0], cell[n]);

        if ( dst >= rstart )
          rndx = (int) floor ( ( dst - rstart ) / dr );
        else
          rndx = 0;

        if ( rndx >= rnum )
          continue;

        if ( i == j ) {
          if ( chg[i] < 0 ) {
            neinaa[rndx] += dot*scale;
            nrm_na[rndx] += scale;
          }
          else if( chg[i] > 0) {
            neincc[rndx] += dot*scale;
            nrm_nc[rndx] += scale;
          }
        }
        else {
          if ( (chg[i] > 0) && (chg[j] > 0) ) {
            cnd_cc[rndx] += dot*scale;
            nrm_cc[rndx] += scale;
          }
          else if ( (chg[i] < 0) && (chg[j] < 0) ) {
            cnd_aa[rndx] += dot*scale;
            nrm_aa[rndx] += scale;
          }
          else {
            cnd_ac[rndx] += dot*scale;
            nrm_ac[rndx] += scale;
          }
        }

      }

    }

    printf(" %i / %i \r", (n+1), nlns);
    fflush(stdout);

  }
  printf("\n");

  //FUDO| division needs to be done for all array elements
  //FUDO| problemativ if nrm == 0 somewhere, which is not unlikely

  int safely = 1;

  divide_array_array_inplace ( neinaa, nrm_na, 1, rnum, safely );
  divide_array_array_inplace ( neincc, nrm_nc, 1, rnum, safely );
  divide_array_array_inplace ( cnd_cc, nrm_cc, 1, rnum, safely );
  divide_array_array_inplace ( cnd_ac, nrm_ac, 1, rnum, safely );
  divide_array_array_inplace ( cnd_aa, nrm_aa, 1, rnum, safely );

#ifdef DEBUG
  write_array_to_file ( "norm_neinaa.out", nrm_na, 1, rnum );
  write_array_to_file ( "norm_neincc.out", nrm_nc, 1, rnum );
  write_array_to_file ( "norm_catcat.out", nrm_cc, 1, rnum );
  write_array_to_file ( "norm_anicat.out", nrm_ac, 1, rnum );
  write_array_to_file ( "norm_aniani.out", nrm_aa, 1, rnum );
#endif

  // free (nrm_neinst);
  free (nrm_na);
  free (nrm_nc);
  free (nrm_cc);
  free (nrm_ac);
  free (nrm_aa);

  free ( vx );
  free ( vy );
  free ( vz );

}
