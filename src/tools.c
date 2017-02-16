#include "tools.h"
#include "macros.h"

void add_array_array ( double *out, double *a, double *b, int ncol, int nlns )
{
  int i, j;

  for ( i=0; i<ncol; i++ )
    for ( j=0; j<nlns; j++ )
      ael(out, nlns, i, j) = ael(a, nlns, i, j) + ael(b, nlns, i, j);

}

void add_arrays_inplace ( double *b, double *a, int ncol, int nlns )
{
  add_array_array ( b, a, b, ncol, nlns );
}

void multiply_array_array ( double *out, double *a, double *b, int ncol, int nlns )
{
  int i, j;

  for ( i=0; i<ncol; i++ )
    for ( j=0; j<nlns; j++ )
      ael(out, nlns, i, j) = ael(a, nlns, i, j) * ael(b, nlns, i, j);

}

void multiply_array_number ( double *out, double *a, double f, int ncol, int nlns )
{
  int i, j;

  for ( i=0; i<ncol; i++ )
    for ( j=0; j<nlns; j++ )
      ael(out, nlns, i, j) = ael(a, nlns, i, j) * f;
}

void subtract_array_number ( double *out, double *a, double f, int ncol, int nlns )
{
  add_array_number ( out, a, -f, ncol, nlns );
}

void add_array_number ( double *out, double *a, double f, int ncol, int nlns )
{
  int i, j;

  for ( i=0; i<ncol; i++ )
    for ( j=0; j<nlns; j++ )
      ael(out, nlns, i, j) = ael(a, nlns, i, j) + f;
}
