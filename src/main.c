#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "errors.h"
#include "io.h"
#include "macros.h"
#include "msd.h"
#include "data.h"

int main(int argc, char *argv[]) {

  int ncol, nlns;
  int qcol, nchg;
  char *delim = " ";

  double *xcom, *ycom, *zcom, *chgs;

  analyze_file ( "../data/xcom.dat", &ncol, &nlns, delim );
  xcom = read_file_double ( "../data/xcom.dat", nlns, ncol, delim );

  analyze_file ( "../data/ycom.dat", &ncol, &nlns, delim );
  ycom = read_file_double ( "../data/ycom.dat", nlns, ncol, delim );

  analyze_file ( "../data/zcom.dat", &ncol, &nlns, delim );
  zcom = read_file_double ( "../data/zcom.dat", nlns, ncol, delim );

  analyze_file ( "../data/charges.dat", &qcol, &nchg, delim );
  chgs = read_file_double ( "../data/charges.dat", nchg, qcol, delim );

  printf("NCHGS: %i, NCOL: %i\n", nchg, ncol);
  if ( nchg != ncol ) {
    print_error ( FATAL, "Dimensions of charge array and position data don't match", __FILE__, __LINE__ );
    return FATAL;
  }

  int i, j;

  /* data I looks fine
  for ( i=1; i<2; i++ ) {
    for ( j=0; j<nlns; j++ ) {
      printf("%14.8f\n", ael(xcom, nlns, i, j));
    }
  }
  */

  double *qflux = (double *) calloc ( nlns, sizeof(double));

  get_qflux ( qflux, xcom, ycom, zcom, chgs, ncol, nlns, NRESTART);

  write_array_to_file ( "cond.out", qflux, 1, nlns );

  printf("NUMBER OF COLUMNS: %i AND NUMBER OF LINES: %i\n", ncol, nlns);

  free ( xcom );
  free ( ycom );
  free ( zcom );

  return 0;

}

