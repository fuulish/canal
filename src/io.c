#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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

void read_input( char *fname, int *nrestart, double *avvol, double *temp, double *timestep, int *split, int *spatial, int *rnum, double *rstart, double *dr )
{
    FILE *datei;
    char *txt;
    char *variable;
    char *value;
    char *buf;
    size_t nbytes = 0;

    enum varset_e
    {
      NRESTART = 1,
      AVVOL,
      TEMP,
      TIMESTEP,
      SPLIT,
      SPATIAL,
    };

    int set_nrestart = NRESTART;
    int set_avvol = AVVOL;
    int set_temp = TEMP;
    int set_timestep = TIMESTEP;
    int set_split = SPLIT;
    int set_spatial = SPATIAL;

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
        set_nrestart = 0;
      }
      else if ( strstr(variable, "avvol") != NULL )
      {
        *avvol = atof(value);
        set_avvol = 0;
      }
      else if ( strstr(variable, "temp") != NULL )
      {
        *temp= atof(value);
        set_temp = 0;
      }
      else if ( strstr(variable, "timestep") != NULL )
      {
        *timestep = atof(value);
        set_timestep = 0;
      }
      else if ( strstr(variable, "split") != NULL )
      {
        *split = atoi(value);
        set_split = 0;
      }
      else if ( strstr(variable, "spatial") != NULL )
      {
        *spatial = 1;
        set_spatial = 0;

        buf = strtok (value, " ");
        *rnum = atoi(buf);

        buf = strtok (NULL, " ");
        *rstart = atof(buf);

        buf = strtok (NULL, " ");
        *dr = atof(buf);

      }
    }

    if ( set_nrestart == NRESTART )
      print_warning ( YOU_KNOW_WHAT, "Using defaults for NRESTART");

    if ( set_avvol == AVVOL )
      print_warning ( YOU_KNOW_WHAT, "Using defaults for AVVOL");

    if ( set_temp == TEMP )
      print_warning ( YOU_KNOW_WHAT, "Using defaults for TEMP");

    if ( set_timestep == TIMESTEP )
      print_warning ( YOU_KNOW_WHAT, "Using defaults for TIMESTEP");

    if ( set_split == SPLIT )
      print_warning ( YOU_KNOW_WHAT, "Using defaults for SPLIT");

    if ( set_spatial == SPATIAL )
      print_warning ( YOU_KNOW_WHAT, "Using defaults for SPATIAL");

    fclose(datei);
    free(txt);

}

