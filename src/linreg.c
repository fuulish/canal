#include "linreg.h"
#include "constants.h"
#include "tools.h"
#include "io.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>                           /* math functions                */

inline static REAL sqr(REAL x) {
  return x*x;
}

int linreg(int n, const REAL x[], const REAL y[], REAL* m, REAL* b, REAL* r)
{
  REAL sumx = 0.0;                        /* sum of x                      */
  REAL sumx2 = 0.0;                       /* sum of x**2                   */
  REAL sumxy = 0.0;                       /* sum of x * y                  */
  REAL sumy = 0.0;                        /* sum of y                      */
  REAL sumy2 = 0.0;                       /* sum of y**2                   */

  int i;
  for ( i=0;i<n;i++)
  {
  sumx  += x[i];
  sumx2 += sqr(x[i]);
  sumxy += x[i] * y[i];
  sumy  += y[i];
  sumy2 += sqr(y[i]);
  }

  REAL denom = (n * sumx2 - sqr(sumx));
  if (denom == 0) {
    // singular matrix. can't solve the problem.
    *m = 0;
    *b = 0;
    if (r) *r = 0;
    return 1;
  }

  *m = (n * sumxy  -  sumx * sumy) / denom;
  *b = (sumy * sumx2  -  sumx * sumxy) / denom;
  if (r!=NULL) {
    *r = (sumxy - sumx * sumy / n) /          /* compute correlation coeff     */
    sqrt((sumx2 - sqr(sumx)/n) *
    (sumy2 - sqr(sumy)/n));
  }

   return 0;
}

void get_linear_regression ( double *data, int len, double temp, double vol, double timestep )
{
  int i;
  double m, b, r, cond;

  double *time = (double *) malloc(len * sizeof(double));
  for ( i=0; i<len; i++ )
    time[i] = i*timestep;

  linreg(len, time, data, &m, &b, &r);

  cond = m / (6. * vol * KBOLTZ * temp );
  cond *= E2C*E2C / A2M / FS2S;

  double m1, b1, r1, cond1;
  int hlen = len/2;
  linreg(hlen, time, data, &m1, &b1, &r1);

  cond1 = m1 / (6. * vol * KBOLTZ * temp );
  cond1 *= E2C*E2C / A2M / FS2S;

  double m2, b2, r2, cond2;
  linreg(hlen, &(time[hlen]), &(data[hlen]), &m2, &b2, &r2);

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
