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
        ncolumn++;
        tok = strtok ( NULL, delim );
      }

      if ( !(*nlns == 0) && (ncolumn != *ncol) ) {
        print_error(IO_ERROR, "Mismatch in number of columns in file.", __FILE__, __LINE__);
        return IO_ERROR;
      }
      else {
        *ncol = ncolumn;
      }

      (*nlns)++;

    }

    fclose ( datei );
  }

  free ( txt );

  return 0;
}

//FUDO| turn around ncol, nlns
double *read_file_double(char *fname, int nlns, int ncol, char *delim)
{
  FILE *datei;
  char *txt;
  char *tok;
  int i;
  int j =0;

  size_t nbytes = 0;
  int ndatpt = ncol * nlns;

  double *data = (double *) malloc ( ndatpt * sizeof(double) );
  datei = fopen(fname, "r");

  if ( datei != NULL ) {
    while ( getline ( &txt, &nbytes, datei ) != -1 ) {
      tok = strtok ( txt, delim );

      i = 0;

      while ( tok != NULL ) {
        ael(data, nlns, i, j) = atof ( tok );

        tok = strtok ( NULL, delim );
        i++;
      }

      j++;

    }
    fclose ( datei );
  }

  free ( txt );
  return data;

}

void write_array_to_file ( char *fname, double *a, int ncol, int nlns )
{
  int i, j;

  FILE *datei;
  datei = fopen(fname, "w");

  //FUDO| catch error like in read_input
  if ( datei != NULL ) {

    for ( j=0; j<nlns; j++ ) {
      for ( i=0; i<ncol; i++ ) {
        fprintf(datei, " %e", ael(a, nlns, i, j));
      }
      fprintf(datei, "\n");
    }

    fclose ( datei );
  }
}

void read_input( char *fname, int *nrestart, double *avvol, double *temp, double *timestep, int *split, int *spatial, int *rnum, double *rstart, double *dr, char *xcom_fn, char *ycom_fn, char *zcom_fn, char *chgs_fn, char *cell_fn, int *task )
{
    FILE *datei;
    char *txt;
    char *variable;
    char *value;
    char *buf;
    size_t nbytes = 0;

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

    datei = fopen(fname, "r");

    if ( datei == NULL )
      print_error ( FILE_NOT_FOUND, "control", __FILE__, __LINE__ );

    while ( getline ( &txt, &nbytes, datei ) != -1 ) {
      variable = strtok (txt, "=");
      value = strtok (NULL, "=");

      printf("VAR: %s VAL: %s\n", variable, value);

      if ( strstr(variable, "nrestart") != NULL )
      {
        *nrestart = atol(value);
        def_nrestart = 0;
      }
      else if ( strstr(variable, "avvol") != NULL )
      {
        *avvol = atof(value);
        def_avvol = 0;
      }
      else if ( strstr(variable, "temp") != NULL )
      {
        *temp= atof(value);
        def_temp = 0;
      }
      else if ( strstr(variable, "timestep") != NULL )
      {
        *timestep = atof(value);
        def_timestep = 0;
      }
      else if ( strstr(variable, "split") != NULL )
      {
        *split = atoi(value);
        def_split = 0;
      }
      else if ( strstr(variable, "spatial") != NULL )
      {
        *spatial = 1;
        def_spatial = 0;

        buf = strtok (value, " ");
        *rnum = atoi(buf);

        buf = strtok (NULL, " ");
        *rstart = atof(buf);

        buf = strtok (NULL, " ");
        *dr = atof(buf);

      }
      else if ( strstr(variable, "xcom") != NULL )
      {
        buf = strtok (value, " \n");
        strcpy ( xcom_fn, buf );
        def_xcom = 0;
      }
      else if ( strstr(variable, "ycom") != NULL )
      {
        buf = strtok (value, " \n");
        strcpy ( ycom_fn, buf );
        def_ycom = 0;
      }
      else if ( strstr(variable, "zcom") != NULL )
      {
        buf = strtok (value, " \n");
        strcpy ( zcom_fn, buf );
        def_zcom = 0;
      }
      else if ( strstr(variable, "chgs") != NULL )
      {
        buf = strtok (value, " \n");
        strcpy ( chgs_fn, buf );
        def_chgs = 0;
      }
      else if ( strstr(variable, "cell") != NULL )
      {
        buf = strtok (value, " \n");
        strcpy ( cell_fn, buf );
        def_cell = 0;
      }
      else if ( strstr(variable, "task") != NULL )
      {
        if ( strstr ( value, "cond" ) != NULL )
          *task = COND;
        else if ( strstr ( value, "velp" ) != NULL )
          *task = VELP;

        def_task = 0;
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

    fclose(datei);
    free(txt);

}

