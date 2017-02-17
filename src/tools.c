#include "tools.h"
#include "macros.h"

void add_array_array ( double *out, double *a, double *b, int ncol, int nlns )
{
  int i, j;

#ifdef OPENMP
#pragma omp parallel for default(none) \
  private(i,j) shared(out, nlns, ncol, a, b)
#endif
  for ( i=0; i<ncol; i++ )
    for ( j=0; j<nlns; j++ )
      ael(out, nlns, i, j) = ael(a, nlns, i, j) + ael(b, nlns, i, j);

}

void add_arrays_inplace ( double *b, double *a, int ncol, int nlns )
{
  add_array_array ( b, a, b, ncol, nlns );
}

void subtract_array_array ( double *out, double *a, double *b, int ncol, int nlns )
{
  int i, j;

#ifdef OPENMP
#pragma omp parallel for default(none) \
  private(i,j) shared(out, nlns, ncol, a, b)
#endif
  for ( i=0; i<ncol; i++ )
    for ( j=0; j<nlns; j++ )
      ael(out, nlns, i, j) = ael(a, nlns, i, j) - ael(b, nlns, i, j);

}

void subtract_arrays_inplace ( double *a, double *b, int ncol, int nlns )
{
  subtract_array_array ( a, a, b, ncol, nlns );
}

void divide_array_array ( double *out, double *a, double *b, int ncol, int nlns )
{
  int i, j;

#ifdef OPENMP
#pragma omp parallel for default(none) \
  private(i,j) shared(out, nlns, ncol, a, b)
#endif
  for ( i=0; i<ncol; i++ )
    for ( j=0; j<nlns; j++ )
      ael(out, nlns, i, j) = ael(a, nlns, i, j) / ael(b, nlns, i, j);

}

void multiply_array_array_inplace ( double *a, double *b, int ncol, int nlns )
{
  multiply_array_array ( a, a, b, ncol, nlns );
}

void divide_array_array_inplace ( double *a, double *b, int ncol, int nlns )
{
  divide_array_array ( a, a, b, ncol, nlns );
}

void multiply_array_array ( double *out, double *a, double *b, int ncol, int nlns )
{
  int i, j;

#ifdef OPENMP
#pragma omp parallel for default(none) \
  private(i,j) shared(out, nlns, ncol, a, b)
#endif
  for ( i=0; i<ncol; i++ )
    for ( j=0; j<nlns; j++ )
      ael(out, nlns, i, j) = ael(a, nlns, i, j) * ael(b, nlns, i, j);

}

void divide_array_number ( double *out, double *a, double f, int ncol, int nlns )
{
  int i, j;

#ifdef OPENMP
#pragma omp parallel for default(none) \
  private(i,j) shared(out, nlns, ncol, a, f)
#endif
  for ( i=0; i<ncol; i++ )
    for ( j=0; j<nlns; j++ )
      ael(out, nlns, i, j) = ael(a, nlns, i, j) / f;
}

void divide_array_number_inplace (double *out, double f, int ncol, int nlns)
{
  multiply_array_number ( out, out, f, ncol, nlns );
}


void multiply_array_number ( double *out, double *a, double f, int ncol, int nlns )
{
  int i, j;

#ifdef OPENMP
#pragma omp parallel for default(none) \
  private(i,j) shared(out, nlns, ncol, a, f)
#endif
  for ( i=0; i<ncol; i++ )
    for ( j=0; j<nlns; j++ )
      ael(out, nlns, i, j) = ael(a, nlns, i, j) * f;
}

void multiply_array_number_inplace (double *out, double f, int ncol, int nlns)
{
  multiply_array_number ( out, out, f, ncol, nlns );
}

void subtract_array_number ( double *out, double *a, double f, int ncol, int nlns )
{
  add_array_number ( out, a, -f, ncol, nlns );
}
void add_array_number_inplace ( double *out, double f, int ncol, int nlns )
{
  add_array_number ( out, out, f, ncol, nlns );
}

void add_array_number ( double *out, double *a, double f, int ncol, int nlns )
{
  int i, j;

#ifdef OPENMP
#pragma omp parallel for default(none) \
  private(i,j) shared(out, nlns, ncol, a, f)
#endif
  for ( i=0; i<ncol; i++ )
    for ( j=0; j<nlns; j++ )
      ael(out, nlns, i, j) = ael(a, nlns, i, j) + f;
}
