#include "tools.h"
#include "macros.h"

void add_array_array ( double *out, double *a, double *b, int ncol, int nlns )
{
  int i, j;
  double *ptr_a, *ptr_b, *ptr_out;

#ifdef OPENMP
#pragma omp parallel for default(none) \
  private(i,j) shared(out, nlns, ncol, a, b, ptr_a, ptr_b, ptr_out)
#endif
  for ( i=0; i<ncol; i++ ) {
    ptr_out = asub(out, nlns, i, 0);
    ptr_a = asub(a, nlns, i, 0);
    ptr_b = asub(b, nlns, i, 0);

    for ( j=0; j<nlns; j++ )
      ptr_out[j] = ptr_a[j] + ptr_b[j];
  }

}

void add_arrays_inplace ( double *b, double *a, int ncol, int nlns )
{
  add_array_array ( b, a, b, ncol, nlns );
}

void divide_array_array ( double *out, double *a, double *b, int ncol, int nlns )
{
  int i, j;
  double *ptr_a, *ptr_b, *ptr_out;

#ifdef OPENMP
#pragma omp parallel for default(none) \
  private(i,j) shared(out, nlns, ncol, a, b, ptr_a, ptr_b, ptr_out)
#endif
  for ( i=0; i<ncol; i++ ) {
    ptr_out = asub(out, nlns, i, 0);
    ptr_a = asub(a, nlns, i, 0);
    ptr_b = asub(b, nlns, i, 0);

    for ( j=0; j<nlns; j++ )
      ptr_out[j] = ptr_a[j] / ptr_b[j];
  }

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
  double *ptr_a, *ptr_b, *ptr_out;

#ifdef OPENMP
#pragma omp parallel for default(none) \
  private(i,j) shared(out, nlns, ncol, a, b, ptr_a, ptr_b, ptr_out)
#endif
  for ( i=0; i<ncol; i++ ) {
    ptr_out = asub(out, nlns, i, 0);
    ptr_a = asub(a, nlns, i, 0);
    ptr_b = asub(b, nlns, i, 0);

    for ( j=0; j<nlns; j++ )
      ptr_out[j] = ptr_a[j] * ptr_b[j];
  }

}

//FUDO| remove this function, by reusing multiply version
void divide_array_number ( double *out, double *a, double f, int ncol, int nlns )
{
  int i, j;
  double *ptr_a, *ptr_out;

#ifdef OPENMP
#pragma omp parallel for default(none) \
  private(i,j) shared(out, nlns, ncol, a, f, ptr_a, ptr_out)
#endif
  for ( i=0; i<ncol; i++ ) {
    ptr_out = asub(out, nlns, i, 0);
    ptr_a = asub(a, nlns, i, 0);

    for ( j=0; j<nlns; j++ )
      ptr_out[j] = ptr_a[j] / f;
  }
}

void divide_array_number_inplace (double *out, double f, int ncol, int nlns)
{
  multiply_array_number ( out, out, f, ncol, nlns );
}


void multiply_array_number ( double *out, double *a, double f, int ncol, int nlns )
{
  int i, j;
  double *ptr_a, *ptr_out;

#ifdef OPENMP
#pragma omp parallel for default(none) \
  private(i,j) shared(out, nlns, ncol, a, f, ptr_a, ptr_out)
#endif
  for ( i=0; i<ncol; i++ ) {
    ptr_out = asub(out, nlns, i, 0);
    ptr_a = asub(a, nlns, i, 0);

    for ( j=0; j<nlns; j++ )
      ptr_out[j] = ptr_a[j] * f;
  }
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
  double *ptr_a, *ptr_out;

#ifdef OPENMP
#pragma omp parallel for default(none) \
  private(i,j) shared(out, nlns, ncol, a, f)
#endif
  for ( i=0; i<ncol; i++ ) {
    ptr_out = asub(out, nlns, i, 0);
    ptr_a = asub(a, nlns, i, 0);

    for ( j=0; j<nlns; j++ )
      ptr_out[j] = ptr_a[j] + f;
  }
}
