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

enum tasks_e {
  NTTN=0,
  COND,
  VELP,
};

int analyze_file ( char *fname, int *ncol, int *nlns, char *delim);
double *read_file_double(char *fname, int nlns, int ncol, char *delim);
void write_array_to_file ( char *fname, double *a, int ncol, int nlns );
void read_input( char *fname, int *nrestart, double *avvol, double *temp, double *timestep, int *split, int *spatial, int *rnum, double *rstart, double *dr, char *xcom_fn, char *ycom_fn, char *zcom_fn, char *chgs_fn, char *cell_fn, int *task, double *fitoffset, double *fitlength );
