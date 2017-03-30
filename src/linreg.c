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

#include "linreg.h"
#include "constants.h"
#include "tools.h"
#include "io.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int linear_regression(int n, const double x[], const double y[], double* m, double* b, double* r)
{
  double sum_x = 0.0;
  double sum_x2 = 0.0;
  double sum_xy = 0.0;
  double sum_y = 0.0;
  double sum_y2 = 0.0;

  int i;
  for ( i=0; i<n ; ++i) {
    sum_x  += x[i];
    sum_x2 += x[i] * x[i];
    sum_xy += x[i] * y[i];
    sum_y  += y[i];
    sum_y2 += y[i] * y[i];
  }

  double denom = (n * sum_x2 - sum_x * sum_x);
  if (denom == 0) {
    // singular matrix. can't solve the problem.
    *m = 0;
    *b = 0;
    if ( r != NULL ) *r = 0;
    return 1;
  }

  *m = (n * sum_xy  -  sum_x * sum_y) / denom;
  *b = (sum_y * sum_x2  -  sum_x * sum_xy) / denom;

  if ( r != NULL )
    *r = (sum_xy - sum_x * sum_y / n) / sqrt((sum_x2 - sum_x * sum_x / n) * (sum_y2 - sum_y * sum_y / n));

   return 0;
}

void get_linear_regression ( double *data, int len, double temp, double vol, double timestep )
{
  int i;
  double m, b, r, cond;

  double *time = (double *) malloc(len * sizeof(double));
  for ( i=0; i<len; i++ )
    time[i] = i*timestep;

  linear_regression(len, time, data, &m, &b, &r);

  cond = m / (6. * vol * KBOLTZ * temp );
  cond *= E2C*E2C / A2M / FS2S;

  double m1, b1, r1, cond1;
  int hlen = len/2;
  linear_regression(hlen, time, data, &m1, &b1, &r1);

  cond1 = m1 / (6. * vol * KBOLTZ * temp );
  cond1 *= E2C*E2C / A2M / FS2S;

  double m2, b2, r2, cond2;
  linear_regression(hlen, &(time[hlen]), &(data[hlen]), &m2, &b2, &r2);

  cond2 = m2 / (6. * vol * KBOLTZ * temp );
  cond2 *= E2C*E2C / A2M / FS2S;

  printf("FIT PARAMETER: %14.8f %14.8f\n", m, b);
  printf("CONDUCTIVITY IS: %14.8f +/- %14.8f S/m\n", cond, fabsf(cond2-cond1));

#ifdef DEBUG
  multiply_array_number_inplace ( time, m, 1, len );
  add_array_number_inplace ( time, b, 1, len );
  write_array_to_file ( "fit_cond.out", time, 1, len );
#endif
}
