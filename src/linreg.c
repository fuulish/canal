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
#include <gsl/gsl_fit.h>
#include <gsl/gsl_cdf.h>
// #include <gsl/gsl_randist.h>

double get_linear_regression(size_t n, const double x[], const double y[], double* m, double* b, double* r)
{
  size_t xstride = 1;
  size_t ystride = 1;

  double cov00, cov01, cov11, sumsq;

  int ret = gsl_fit_linear (x, xstride, y, ystride, n, b, m, &cov00, &cov01, &cov11, &sumsq);

  // printf("CHISQ: %e, M: %e +/- %e, B: %e +/- %e\n", sumsq, *m, sqrt(cov11), *b, sqrt(cov00));

  // double stdev0=sqrt(cov00);
  // double t0=*b/stdev0;
  // double pv0=t0<0?2*(1-gsl_cdf_tdist_P(-t0,n-2)):2*(1-gsl_cdf_tdist_P(t0,n-2));//This is the p-value of the constant term
  // cout<<"Intercept\t"<<*b<<"\t"<<stdev0<<"\t"<<t0<<"\t"<<pv0<<endl;

  double stdev1=sqrt(cov11);
  // double t1=*m/stdev1;
  // double pv1=t1<0?2*(1-gsl_cdf_tdist_P(-t1,n-2)):2*(1-gsl_cdf_tdist_P(t1,n-2));//This is the p-value of the linear term
  // cout<<"x\t"<<*m<<"\t"<<stdev1<<"\t"<<t1<<"\t"<<pv1<<endl;

  return stdev1;
  // return pv1;
}

double calculate_conductivity ( double *data, size_t len, double temp, double vol, double timestep, int fitstrt, char *outprefix, char *units )
{
  int i;
  double m, b, r, cond;

  double *time = (double *) malloc(len * sizeof(double));
  for ( i=0; i<len; ++i )
    time[i] = (i + fitstrt) * timestep;

  double cond_err = get_linear_regression(len, time, data, &m, &b, &r);

  cond = m / (6. * vol * KBOLTZ * temp );
  cond *= E2C*E2C / A2M / FS2S;

  double m1, b1, r1, cond1;
  int hlen = len/2;
  get_linear_regression(hlen, time, data, &m1, &b1, &r1);

  cond1 = m1 / (6. * vol * KBOLTZ * temp );
  cond1 *= E2C*E2C / A2M / FS2S;

  double m2, b2, r2, cond2;
  get_linear_regression(hlen, &(time[hlen]), &(data[hlen]), &m2, &b2, &r2);

  cond2 = m2 / (6. * vol * KBOLTZ * temp );
  cond2 *= E2C*E2C / A2M / FS2S;

  cond_err /= (6. * vol * KBOLTZ * temp );
  cond_err *= E2C*E2C / A2M / FS2S;

  printf("FIT PARAMETER: %14.8f %14.8f\n", m, b);
  printf("%s: %e +/- %e %s\n", outprefix, cond, cond_err, units);
  // printf("%s: %e +/- %14.8f %s\n", outprefix, cond, fabs(cond2-cond1), units);
  // printf("CONDUCTIVITY IS: %14.8f +/- %14.8f S/m\n", cond, fabs(cond2-cond1));

  free( time );

#ifdef DEBUG
  multiply_array_number_inplace ( time, m, 1, len );
  add_array_number_inplace ( time, b, 1, len );
  write_array_to_file ( "fit_cond.out", time, 1, len );
#endif

  return cond;
}
