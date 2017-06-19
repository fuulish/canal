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


enum diff_e {
  FWD=1,
  BWD,
  CDD,
};

double calculate_dot ( double *xi, double *yi, double *zi, double *xj, double *yj, double *zj, int nlns );
double calculate_nrm ( double *xi, double *yi, double *zi, int nlns );
void calculate_velocities ( double *vx, double *vy, double *vz, double *xi, double *yi, double *zi, int nlns, int timestep );
void get_vflux_locl ( double *neinaa, double *neincc, double *cnd_cc, double *cnd_ac, double *cnd_aa, double *x, double *y, double *z, double *chg, int ncol, int nlns, int nrestart, double dr, double rstart, int rnum, double *cell, double timestep );
