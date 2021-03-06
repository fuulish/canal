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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "errors.h"
#include "io.h"
#include "macros.h"
#include "msd.h"
#include "vel.h"
#include "constants.h"
#include "linreg.h"

int main(int argc, char *argv[]) {

  int ncol, nlns;
  int qcol, nchg;
  int ccol, ncll;
  char delim[2] = " ";

  int nrestart = 1;
  double avvol = 1.;
  double temp = 300.;
  double timestep = 0.5;
  int split = 0;
  int spatial = 0;
  int rnum = 1;
  double rstart = 5.;
  double dr = 2.;
  int task = NTTN;
  double fitoffset = 0.1;
  double fitlength = 0.9;

  char xcom_fn[100];
  char ycom_fn[100];
  char zcom_fn[100];

  char chgs_fn[100];
  char cell_fn[100];

  if ( argc > 1 )
    read_input( argv[1], &nrestart, &avvol, &temp, &timestep, &split, &spatial, &rnum, &rstart, &dr, xcom_fn, ycom_fn, zcom_fn, chgs_fn, cell_fn, &task, &fitoffset, &fitlength );
  else
    print_error ( INCOMPLETE_INPUT, "Input file is missing", __FILE__, __LINE__ );

  double *xcom, *ycom, *zcom, *chgs, *cell;

  analyze_file ( xcom_fn, &ncol, &nlns, delim );
  xcom = read_file_double ( xcom_fn, nlns, ncol, delim );

  analyze_file ( ycom_fn, &ncol, &nlns, delim );
  ycom = read_file_double ( ycom_fn, nlns, ncol, delim );

  analyze_file ( zcom_fn, &ncol, &nlns, delim );
  zcom = read_file_double ( zcom_fn, nlns, ncol, delim );

  analyze_file ( chgs_fn, &qcol, &nchg, delim );
  chgs = read_file_double ( chgs_fn, nchg, qcol, delim );

  printf("NCHGS: %i, NCOL: %i\n", nchg, ncol);
  if ( nchg != ncol ) {
    print_error ( FATAL, "Dimensions of charge array and position data don't match", __FILE__, __LINE__ );
    return FATAL;
  }

  //FUDO| this needs fixing, cell allocation below only to avoid issues with maybe-unitialized warning
  if ( rnum > 1 ) {

    analyze_file ( cell_fn, &ccol, &ncll, delim );
    cell = read_file_double ( cell_fn, ncll, ccol, delim );

    if ( ncll != nlns ) {
      print_error ( FATAL, "Dimensions of cell array and position data don't match", __FILE__, __LINE__ );
      return FATAL;
    }
  }
  else {
    cell = (double *) calloc ( nlns, sizeof (double));
  }

  /* data I looks fine
  for ( i=1; i<2; ++i ) {
    for ( j=0; j<nlns; ++j ) {
      printf("%14.8f\n", ael(xcom, nlns, i, j));
    }
  }
  */

  switch ( task ) {
    case COND:
      {

      if ( split ) {

        int arrlen = nlns * rnum;

        double *qflux_neinaa = (double *) calloc ( arrlen, sizeof(double));
        double *qflux_neincc = (double *) calloc ( arrlen, sizeof(double));
        double *qflux_catcat = (double *) calloc ( arrlen, sizeof(double));
        double *qflux_anicat = (double *) calloc ( arrlen, sizeof(double));
        double *qflux_aniani = (double *) calloc ( arrlen, sizeof(double));

        get_qflux_srtd ( qflux_neinaa, qflux_neincc, qflux_catcat, qflux_anicat, qflux_aniani, xcom, ycom, zcom, chgs, ncol, nlns, nrestart, dr, rstart, rnum, cell );
        // get_qflux_srtd ( qflux_neinst, qflux_catcat, qflux_anicat, qflux_aniani, xcom, ycom, zcom, chgs, ncol, nlns, nrestart);

        write_array_to_file ( "cond_neinaa.out", qflux_neinaa, 1, nlns );
        write_array_to_file ( "cond_neincc.out", qflux_neincc, 1, nlns );
        write_array_to_file ( "cond_catcat.out", qflux_catcat, rnum, nlns );
        write_array_to_file ( "cond_anicat.out", qflux_anicat, rnum, nlns );
        write_array_to_file ( "cond_aniani.out", qflux_aniani, rnum, nlns );

        printf("NUMBER OF COLUMNS: %i AND NUMBER OF LINES: %i\n", ncol, nlns);

        int fitstrt, nlns_fit;
        fitstrt = fitoffset * nlns;
        nlns_fit = fitlength * nlns;

        printf("NEINAA: ");
        calculate_conductivity(&(qflux_neinaa[fitstrt]), nlns_fit, temp, avvol, timestep, fitstrt, "CONDUCTIVITY IS:", "S/m");
        printf("NEINCC: ");
        calculate_conductivity(&(qflux_neincc[fitstrt]), nlns_fit, temp, avvol, timestep, fitstrt, "CONDUCTIVITY IS:", "S/m");
        printf("CATCAT: ");
        calculate_conductivity(&(qflux_catcat[fitstrt]), nlns_fit, temp, avvol, timestep, fitstrt, "CONDUCTIVITY IS:", "S/m");
        printf("ANIANI: ");
        calculate_conductivity(&(qflux_aniani[fitstrt]), nlns_fit, temp, avvol, timestep, fitstrt, "CONDUCTIVITY IS:", "S/m");
        printf("ANICAT: ");
        calculate_conductivity(&(qflux_anicat[fitstrt]), nlns_fit, temp, avvol, timestep, fitstrt, "CONDUCTIVITY IS:", "S/m");

        free ( qflux_neinaa );
        free ( qflux_neincc );
        free ( qflux_catcat );
        free ( qflux_anicat );
        free ( qflux_aniani );
      }
      else {
        double *qflux = (double *) calloc ( nlns, sizeof(double));

        get_qflux ( qflux, xcom, ycom, zcom, chgs, ncol, nlns, nrestart);

        int fitstrt, nlns_fit;
        fitstrt = fitoffset * nlns;
        nlns_fit = fitlength * nlns;
        printf("NLNS_FIT: %i %i %i\n", nlns_fit, nlns, fitstrt);

        printf("FULL: ");
        calculate_conductivity(&(qflux[fitstrt]), nlns_fit, temp, avvol, timestep, fitstrt, "CONDUCTIVITY IS:", "S/m");

        write_array_to_file ( "cond_all.out", qflux, 1, nlns );

        printf("NUMBER OF COLUMNS: %i AND NUMBER OF LINES: %i\n", ncol, nlns);

        free ( qflux );
      }
      }
      break;
    case VELP:
      {
      int arrlen = rnum;

      double *qflux_neinaa = (double *) calloc ( arrlen, sizeof(double));
      double *qflux_neincc = (double *) calloc ( arrlen, sizeof(double));
      double *qflux_catcat = (double *) calloc ( arrlen, sizeof(double));
      double *qflux_anicat = (double *) calloc ( arrlen, sizeof(double));
      double *qflux_aniani = (double *) calloc ( arrlen, sizeof(double));

      get_vflux_locl ( qflux_neinaa, qflux_neincc, qflux_catcat, qflux_anicat, qflux_aniani, xcom, ycom, zcom, chgs, ncol, nlns, nrestart, dr, rstart, rnum, cell, timestep );

      write_array_to_file ( "velp_neinaa.out", qflux_neinaa, 1, rnum );
      write_array_to_file ( "velp_neincc.out", qflux_neincc, 1, rnum );
      write_array_to_file ( "velp_catcat.out", qflux_catcat, 1, rnum );
      write_array_to_file ( "velp_anicat.out", qflux_anicat, 1, rnum );
      write_array_to_file ( "velp_aniani.out", qflux_aniani, 1, rnum );

      printf("NUMBER OF COLUMNS: %i AND NUMBER OF LINES: %i\n", ncol, nlns);

      free ( qflux_neinaa );
      free ( qflux_neincc );
      free ( qflux_catcat );
      free ( qflux_anicat );
      free ( qflux_aniani );
      }
      break;
    case ELMO:
      {
        // check a couple of things, and maybe warn
        int arrlen = nlns * rnum;

        double *qflux_neinaa = (double *) calloc ( arrlen, sizeof(double));
        double *qflux_neincc = (double *) calloc ( arrlen, sizeof(double));
        double *qflux_catcat = (double *) calloc ( arrlen, sizeof(double));
        double *qflux_anicat = (double *) calloc ( arrlen, sizeof(double));
        double *qflux_aniani = (double *) calloc ( arrlen, sizeof(double));

        get_mobil_srtd ( qflux_neinaa, qflux_neincc, qflux_catcat, qflux_anicat, qflux_aniani, xcom, ycom, zcom, chgs, ncol, nlns, nrestart, dr, rstart, rnum, cell );

        write_array_to_file ( "mobil_neinaa.out", qflux_neinaa, 1, nlns );
        write_array_to_file ( "mobil_neincc.out", qflux_neincc, 1, nlns );
        write_array_to_file ( "mobil_catcat.out", qflux_catcat, rnum, nlns );
        write_array_to_file ( "mobil_anicat.out", qflux_anicat, rnum, nlns );
        write_array_to_file ( "mobil_aniani.out", qflux_aniani, rnum, nlns );

        printf("NUMBER OF COLUMNS: %i AND NUMBER OF LINES: %i\n", ncol, nlns);

        int fitstrt, nlns_fit;
        fitstrt = fitoffset * nlns;
        nlns_fit = fitlength * nlns;

        //FUX| for now we'll make use of the conductivity calculation, which is essentially the same, except a couple of pre-factors
        //FUX| we use (E2C * NCHG / A2M^3) instead of the volume, need to count anions and cations

        int i;
        int ncat = 0;
        int nani = 0;
        // int ntot = 0;

        for ( i=0; i<ncol; ++i ) {
          if( chgs[i] < 0 )
            nani += 1;
          else if( chgs[i] > 0 )
            ncat += 1;
        }

        // ntot = nani + ncat;
        double uniconv = E2C / (A2M * A2M * A2M);
        double elmo_ani, elmo_cat;

        //FUX| neinst and anicat don't make much sense right now

        // printf("NEINST: ");
        // calculate_conductivity(&(qflux_neinst[fitstrt]), nlns_fit, temp, ntot*uniconv, timestep, fitstrt);
        printf("CATCAT: ");
        elmo_cat = calculate_conductivity(&(qflux_catcat[fitstrt]), nlns_fit, temp, ncat*uniconv, timestep, fitstrt, "MOBILITY IS:", "m^2/Vs");
        printf("ANIANI: ");
        elmo_ani = calculate_conductivity(&(qflux_aniani[fitstrt]), nlns_fit, temp, nani*uniconv, timestep, fitstrt, "MOBILITY IS:", "m^2/Vs");
        // printf("ANICAT: ");
        // calculate_conductivity(&(qflux_anicat[fitstrt]), nlns_fit, temp, ntot*uniconv, timestep, fitstrt, "MOBILITY IS:", "m^2/Vs");

        printf("actual conductivity is: %14.8f\n", -1.*elmo_ani / avvol * nani * uniconv + 1.*elmo_cat / avvol * ncat * uniconv );

        free ( qflux_neinaa );
        free ( qflux_neincc );
        free ( qflux_catcat );
        free ( qflux_anicat );
        free ( qflux_aniani );
      }
      break;
    case DIFF:
      {
      int arrlen = nlns * rnum;

      double *qflux_neinaa = (double *) calloc ( arrlen, sizeof(double));
      double *qflux_neincc = (double *) calloc ( arrlen, sizeof(double));

      get_diff ( qflux_neinaa, qflux_neincc, xcom, ycom, zcom, chgs, ncol, nlns, nrestart, dr, rstart, rnum, cell);

      write_array_to_file ( "diff_neinaa.out", qflux_neinaa, 1, nlns );
      write_array_to_file ( "diff_neincc.out", qflux_neincc, 1, nlns );

      printf("NUMBER OF COLUMNS: %i AND NUMBER OF LINES: %i\n", ncol, nlns);

      int fitstrt, nlns_fit;
      fitstrt = fitoffset * nlns;
      nlns_fit = fitlength * nlns;

      // first terms counter-acting conductivity stuff, last term units conversion angstrom^2/fs -> m^2/s
      //double uniconv = E2C*E2C / A2M / FS2S / KBOLTZ * 1.e5;
      double uniconv = E2C*E2C / A2M / FS2S / KBOLTZ * 1.e1;

      int i;
      int ncat = 0;
      int nani = 0;
      // int ntot = 0;

      for ( i=0; i<ncol; ++i ) {
        if( chgs[i] < 0 )
          nani += 1;
        else if( chgs[i] > 0 )
          ncat += 1;
      }

      printf("NEINAA: ");
      calculate_conductivity(&(qflux_neinaa[fitstrt]), nlns_fit, nani, uniconv, timestep, fitstrt, "DIFFUSION IS:", "cm^2/s");
      printf("NEINCC: ");
      calculate_conductivity(&(qflux_neincc[fitstrt]), nlns_fit, ncat, uniconv, timestep, fitstrt, "DIFFUSION IS:", "cm^2/s");

      free ( qflux_neinaa );
      free ( qflux_neincc );
      }
      break;
    default:
      break;
  }

  free ( cell );

  free ( xcom );
  free ( ycom );
  free ( zcom );

  free ( chgs );

  return 0;

}

