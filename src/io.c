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

  if ( datei = fopen(fname, "r") ) {
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

double *read_file_double(char *fname, int nlns, int ncol, char *delim)
{
  FILE *datei;
  char *txt;
  char *tok;
  int i;
  int j =0;

  size_t nbytes = 0;
  int nlines = 0;
  int ndatpt = ncol * nlns;

  double *data = (double *) malloc ( ndatpt * sizeof(double) );

  if ( datei = fopen(fname, "r") ) {
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

  if ( datei = fopen(fname, "w") ) {

    for ( j=0; j<nlns; j++ ) {
      for ( i=0; i<ncol; i++ ) {
        fprintf(datei, " %e", ael(a, nlns, i, j));
      }
      fprintf(datei, "\n");
    }

    fclose ( datei );
  }
}
