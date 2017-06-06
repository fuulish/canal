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


void multiply_array_array ( double *out, double *a, double *b, int ncol, int nlns );
void multiply_array_array_inplace ( double *a, double *b, int ncol, int nlns );
void divide_array_number ( double *out, double *a, double f, int ncol, int nlns );
void divide_array_number_inplace (double *out, double f, int ncol, int nlns);
void multiply_array_number ( double *out, double *a, double f, int ncol, int nlns );
void multiply_array_number_inplace (double *out, double f, int ncol, int nlns);
void add_array_number ( double *out, double *a, double f, int ncol, int nlns );
void add_array_number_inplace ( double *out, double f, int ncol, int nlns );
void subtract_array_number ( double *out, double *a, double f, int ncol, int nlns );
void add_array_array ( double *out, double *a, double *b, int ncol, int nlns );
void add_arrays_inplace ( double *b, double *a, int ncol, int nlns );
void divide_array_array ( double *out, double *a, double *b, int ncol, int nlns );
void divide_array_array_inplace ( double *a, double *b, int ncol, int nlns );
void subtract_array_array ( double *out, double *a, double *b, int ncol, int nlns );
void subtract_arrays_inplace ( double *a, double *b, int ncol, int nlns );
void copy_number_to_array( double* out, double f, int ncol, int nlns );
void zero_out_array( double *out, int ncol, int nlns );
