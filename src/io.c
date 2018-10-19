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
#include <limits.h>
#include "errors.h"
#include "io.h"
#include "macros.h"

int analyze_file ( char *fname, int *ncol, int *nlns, char *delim)
{
  FILE *datei;
  char *txt;
  char *tok;

  size_t nbytes = 0;
  int ncolumn;
  *nlns = 0;

  datei = fopen(fname, "r");

  //FUDO| catch error like in read_input
  if ( datei != NULL ) {
    while ( getline ( &txt, &nbytes, datei ) != -1 ) {
      tok = strtok ( txt, delim );

      ncolumn = 0;
      while ( tok !=  NULL ) {
        ++ncolumn;
        tok = strtok ( NULL, delim );
      }

      if ( !(*nlns == 0) && (ncolumn != *ncol) ) {
        print_error(IO_ERROR, "Mismatch in number of columns in file.", __FILE__, __LINE__);
        return IO_ERROR;
      }
      else {
        *ncol = ncolumn;
      }

      ++(*nlns);

    }

    fclose ( datei );
    free ( txt );
  }

  return 0;
}

//FUDO| turn around ncol, nlns
double *read_file_double(char *fname, int nlns, int ncol, char *delim, int stride)
{
  FILE *datei;
  char *txt;
  char *tok;
  int i;
  int j =0;
  int cnt = 0;

  size_t nbytes = 0;

  int localNlns = nlns / stride; // savety/rounding margin
  int ndatpt = ncol * localNlns;

  double *data = (double *) malloc ( ndatpt * sizeof(double) );
  datei = fopen(fname, "r");

  if ( datei != NULL ) {
    while ( (getline ( &txt, &nbytes, datei ) != -1) && (j < localNlns) ) {
      if( cnt % stride == 0 ) {
        tok = strtok ( txt, delim );

        i = 0;

        while ( tok != NULL ) {
          ael(data, localNlns, i, j) = atof ( tok );
          if( aindex(localNlns, i, j) >= ndatpt )
          {
            printf("SOMETHING's WRONG\n");
            printf("%i %i\n", aindex(localNlns, i, j), ndatpt);
          }

          tok = strtok ( NULL, delim );
          ++i;
        }
        ++j;
      }
      ++cnt;
    }
    fclose ( datei );
    free ( txt );
  }

  return data;

}

void write_array_to_file ( char *fname, double *a, int ncol, int nlns )
{
  int i, j;

  FILE *datei;
  datei = fopen(fname, "w");

  //FUDO| catch error like in read_input
  if ( datei != NULL ) {

    for ( j=0; j<nlns; ++j ) {
      for ( i=0; i<ncol; ++i ) {
        fprintf(datei, " %e", ael(a, nlns, i, j));
      }
      fprintf(datei, "\n");
    }

    fclose ( datei );
  }
}

inpArg_t read_input( char *fname )
{
    FILE *datei;
    char *txt;
    char *variable;
    char *value;
    char *buf;
    size_t nbytes = 0;

    inpArg_t inputArgs = {.nrestart = 1, .avvol = 1., .temp = 300., .timestep = 0.5, .split = 0, .spatial = 0, .rnum = 1, .rstart = 5., .dr = 2., .task = NTTN, .fitoffset = 0.1, .fitlength = 0.9, .datastride = 1, .nmaxlns = -1, };
    enum varset_e
    {
      NRESTART = INT_MIN,
      AVVOL,
      TEMP,
      TIMESTEP,
      SPLIT,
      SPATIAL,
      XCOM,
      YCOM,
      ZCOM,
      CHGS,
      CELL,
      TASK,
      FITOFFSET,
      FITLENGTH,
      DATASTRIDE,
      NMAXLINES,
    };

    int def_nrestart = NRESTART;
    int def_avvol = AVVOL;
    int def_temp = TEMP;
    int def_timestep = TIMESTEP;
    int def_split = SPLIT;
    int def_spatial = SPATIAL;
    int def_xcom = XCOM;
    int def_ycom = YCOM;
    int def_zcom = ZCOM;
    int def_chgs = CHGS;
    int def_cell = CELL;
    int def_task = TASK;
    int def_fitoffset = FITOFFSET;
    int def_fitlength = FITLENGTH;
    int def_datastride = DATASTRIDE;
    int def_nmaxlines = NMAXLINES;

    datei = fopen(fname, "r");

    if ( datei == NULL )
      print_error ( FILE_NOT_FOUND, fname, __FILE__, __LINE__ );

    while ( getline ( &txt, &nbytes, datei ) != -1 ) {
      variable = strtok (txt, "=");
      value = strtok (NULL, "=");

      printf("VAR: %s VAL: %s\n", variable, value);

      if ( strstr(variable, "nrestart") != NULL )
      {
        inputArgs.nrestart = atol(value);
        def_nrestart = 0;
      }
      else if ( strstr(variable, "avvol") != NULL )
      {
        inputArgs.avvol = atof(value);
        def_avvol = 0;
      }
      else if ( strstr(variable, "temp") != NULL )
      {
        inputArgs.temp= atof(value);
        def_temp = 0;
      }
      else if ( strstr(variable, "timestep") != NULL )
      {
        inputArgs.timestep = atof(value);
        def_timestep = 0;
      }
      else if ( strstr(variable, "split") != NULL )
      {
        inputArgs.split = atoi(value);
        def_split = 0;
      }
      else if ( strstr(variable, "spatial") != NULL )
      {
        inputArgs.spatial = 1;
        def_spatial = 0;

        buf = strtok (value, " ");
        inputArgs.rnum = atoi(buf);

        if ( inputArgs.rnum < 2 ) {
        // if ( (*rnum == 0) || (*rnum == 1) ) {
          inputArgs.rnum = 1;
          inputArgs.spatial = 0;
          continue;
        }

        buf = strtok (NULL, " ");
        inputArgs.rstart = atof(buf);

        buf = strtok (NULL, " ");
        inputArgs.dr = atof(buf);

      }
      else if ( strstr(variable, "xcom") != NULL )
      {
        buf = strtok (value, " \n");
        strcpy ( inputArgs.xcom_fn, buf );
        def_xcom = 0;
      }
      else if ( strstr(variable, "ycom") != NULL )
      {
        buf = strtok (value, " \n");
        strcpy ( inputArgs.ycom_fn, buf );
        def_ycom = 0;
      }
      else if ( strstr(variable, "zcom") != NULL )
      {
        buf = strtok (value, " \n");
        strcpy ( inputArgs.zcom_fn, buf );
        def_zcom = 0;
      }
      else if ( strstr(variable, "chgs") != NULL )
      {
        buf = strtok (value, " \n");
        strcpy ( inputArgs.chgs_fn, buf );
        def_chgs = 0;
      }
      else if ( strstr(variable, "cell") != NULL )
      {
        buf = strtok (value, " \n");
        strcpy ( inputArgs.cell_fn, buf );
        def_cell = 0;
      }
      else if ( strstr(variable, "task") != NULL )
      {
        if ( strstr ( value, "cond" ) != NULL )
          inputArgs.task = COND;
        else if ( strstr ( value, "velp" ) != NULL )
          inputArgs.task = VELP;
        else if ( strstr ( value, "elmo" ) != NULL )
          inputArgs.task = ELMO;
        else if ( strstr ( value, "diff" ) != NULL )
          inputArgs.task = DIFF;

        def_task = 0;
      }
      else if ( strstr(variable, "fitoffset" ) != NULL )
      {
        buf = strtok (value, " \n");
        inputArgs.fitoffset = atof(buf);
        def_fitoffset = 0;
      }
      else if ( strstr(variable, "fitlength" ) != NULL )
      {
        buf = strtok (value, " \n");
        inputArgs.fitlength = atof(buf);
        def_fitlength = 0;
      }
      else if ( strstr(variable, "datastride" ) != NULL )
      {
        buf = strtok (value, " \n");
        inputArgs.datastride = atol(buf);
        def_datastride = 0;
      }
      else if ( strstr(variable, "nmaxlines" ) != NULL )
      {
        buf = strtok (value, " \n");
        inputArgs.nmaxlns = atol(buf);
        def_nmaxlines = 0;
      }
    }

    if ( def_nrestart == NRESTART )
      print_warning ( YOU_KNOW_WHAT, "Using defaults for NRESTART");

    if ( def_avvol == AVVOL )
      print_warning ( YOU_KNOW_WHAT, "Using defaults for AVVOL");

    if ( def_temp == TEMP )
      print_warning ( YOU_KNOW_WHAT, "Using defaults for TEMP");

    if ( def_timestep == TIMESTEP )
      print_warning ( YOU_KNOW_WHAT, "Using defaults for TIMESTEP");

    if ( def_split == SPLIT )
      print_warning ( YOU_KNOW_WHAT, "Using defaults for SPLIT");

    if ( def_spatial == SPATIAL )
      print_warning ( YOU_KNOW_WHAT, "Using defaults for SPATIAL");

    if ( def_xcom == XCOM )
      print_error ( INCOMPLETE_INPUT, "XCOM file is missing", __FILE__, __LINE__);

    if ( def_ycom == YCOM )
      print_error ( INCOMPLETE_INPUT, "YCOM file is missing", __FILE__, __LINE__);

    if ( def_zcom == ZCOM )
      print_error ( INCOMPLETE_INPUT, "ZCOM file is missing", __FILE__, __LINE__);

    if ( def_chgs == CHGS )
      print_error ( INCOMPLETE_INPUT, "CHGS file is missing", __FILE__, __LINE__);

    if ( def_cell == CELL )
      print_error ( INCOMPLETE_INPUT, "CELL file is missing", __FILE__, __LINE__);

    if ( def_task == TASK )
      print_error ( INCOMPLETE_INPUT, "TASK is missing", __FILE__, __LINE__);

    if ( def_fitoffset == FITOFFSET )
      print_warning ( YOU_KNOW_WHAT, "Using defaults for FITOFFSET");

    if ( def_fitlength == FITLENGTH )
      print_warning ( YOU_KNOW_WHAT, "Using defaults for FITLENGTH");

    if ( def_datastride == DATASTRIDE )
      print_warning ( YOU_KNOW_WHAT, "Using defaults for DATASTRIDE");
    else
      inputArgs.timestep *= inputArgs.datastride;

    if ( def_nmaxlines == NMAXLINES )
      print_warning ( YOU_KNOW_WHAT, "Using defaults for NMAXLINES");


    fclose(datei);
    free(txt);

    return inputArgs;

}

