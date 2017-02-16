#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "errors.h"
#include "io.h"
#include "macros.h"

int main(int argc, char *argv[]) {

  int ncol, nlns;
  char *delim = " ";

  double *xcom, *ycom, *zcom;

  analyze_file ( "../data/xcom.dat", &ncol, &nlns, delim );
  xcom = read_file_double ( "../data/xcom.dat", nlns, ncol, delim );

  analyze_file ( "../data/ycom.dat", &ncol, &nlns, delim );
  ycom = read_file_double ( "../data/ycom.dat", nlns, ncol, delim );

  analyze_file ( "../data/zcom.dat", &ncol, &nlns, delim );
  zcom = read_file_double ( "../data/zcom.dat", nlns, ncol, delim );

  int i, j;

  for ( i=0; i<ncol; i++ ) {
    for ( j=0; j<nlns; j++ ) {
      ael(xcom, nlns, i, j) += 1.;
      // printf("%14.8f\n", ael(xcom, nlns, i, j));
    }
  }

  printf("NUMBER OF COLUMNS: %i AND NUMBER OF LINES: %i\n", ncol, nlns);

  free ( xcom );
  free ( ycom );
  free ( zcom );

  return 0;

}

