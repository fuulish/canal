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

#include "tools.h"
#include "macros.h"

void add_array_array ( double * restrict out, double * restrict a, double * restrict b, int ncol, int nlns )
{
  int i, j;

  for ( i=0; i<ncol; ++i )
    for ( j=0; j<nlns; ++j )
      ael(out, nlns, i, j) = ael(a, nlns, i, j) + ael(b, nlns, i, j);

}

void add_arrays_inplace ( double * restrict a, double * restrict b, int ncol, int nlns )
{
  int i, j;

  for ( i=0; i<ncol; ++i )
    for ( j=0; j<nlns; ++j )
      ael(a, nlns, i, j) += ael(b, nlns, i, j);
}

void subtract_array_array ( double * restrict out, double * restrict a, double * restrict b, int ncol, int nlns )
{
  int i, j;

  for ( i=0; i<ncol; i++ )
    for ( j=0; j<nlns; j++ )
      ael(out, nlns, i, j) = ael(a, nlns, i, j) - ael(b, nlns, i, j);

}

void subtract_arrays_inplace ( double * restrict a, double * restrict b, int ncol, int nlns )
{
  int i, j;

  for ( i=0; i<ncol; i++ )
    for ( j=0; j<nlns; j++ )
      ael(a, nlns, i, j) -= ael(b, nlns, i, j);

}

void divide_array_array ( double * restrict out, double * restrict a, double * restrict b, int ncol, int nlns )
{
  int i, j;

  for ( i=0; i<ncol; ++i )
    for ( j=0; j<nlns; ++j )
      ael(out, nlns, i, j) = ael(a, nlns, i, j) / ael(b, nlns, i, j);

}

void multiply_array_array_inplace ( double * restrict a, double * restrict b, int ncol, int nlns )
{
  int i, j;

  for ( i=0; i<ncol; ++i )
    for ( j=0; j<nlns; ++j )
      ael(a, nlns, i, j) *= ael(b, nlns, i, j);
}

void divide_array_array_inplace ( double * restrict a, double * restrict b, int ncol, int nlns )
{
  int i, j;

  for ( i=0; i<ncol; ++i )
    for ( j=0; j<nlns; ++j )
      ael(a, nlns, i, j) /= ael(b, nlns, i, j);

}

void multiply_array_array ( double * restrict out, double * restrict a, double * restrict b, int ncol, int nlns )
{
  int i, j;

  for ( i=0; i<ncol; ++i )
    for ( j=0; j<nlns; ++j )
      ael(out, nlns, i, j) = ael(a, nlns, i, j) * ael(b, nlns, i, j);

}

void divide_array_number ( double * restrict out, double * restrict a, double f, int ncol, int nlns )
{
  multiply_array_number( out, a, 1./f, ncol, nlns );
}

void divide_array_number_inplace (double * restrict out, double f, int ncol, int nlns)
{
  multiply_array_number_inplace ( out, 1./f, ncol, nlns );
}

void multiply_array_number ( double * restrict out, double * restrict a, double f, int ncol, int nlns )
{
  int i, j;

  for ( i=0; i<ncol; ++i )
    for ( j=0; j<nlns; ++j )
      ael(out, nlns, i, j) = ael(a, nlns, i, j) * f;
}

void multiply_array_number_inplace (double * restrict out, double f, int ncol, int nlns)
{
  int i, j;

  for ( i=0; i<ncol; ++i )
    for ( j=0; j<nlns; ++j )
      ael(out, nlns, i, j) *= f;
}

void subtract_array_number ( double * restrict out, double * restrict a, double f, int ncol, int nlns )
{
  add_array_number ( out, a, -f, ncol, nlns );
}

void subtract_array_number_inplace ( double * restrict out, double f, int ncol, int nlns )
{
  add_array_number_inplace( out, -f, ncol, nlns );
}

void add_array_number_inplace ( double * restrict out, double f, int ncol, int nlns )
{
  int i, j;

  for ( i=0; i<ncol; ++i )
    for ( j=0; j<nlns; ++j )
      ael(out, nlns, i, j) += f;
}

void add_array_number ( double * restrict out, double * restrict a, double f, int ncol, int nlns )
{
  int i, j;

  for ( i=0; i<ncol; ++i )
    for ( j=0; j<nlns; ++j )
      ael(out, nlns, i, j) = ael(a, nlns, i, j) + f;
}
