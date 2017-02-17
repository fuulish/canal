#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "errors.h"
#include "io.h"
#include "macros.h"
#include "msd.h"
#include "constants.h"

int main(int argc, char *argv[]) {

  int ncol, nlns;
  int qcol, nchg;
  int ccol, ncll;
  char *delim = " ";

  int nrestart = 1;
  double avvol = 1.;
  double temp = 300.;
  double timestep = 0.5;
  int split = 0;
  int spatial = 0;
  int rnum = 1;
  double rstart = 5.;
  double dr = 2.;

  char xcom_fn[100];
  char ycom_fn[100];
  char zcom_fn[100];

  char chgs_fn[100];
  char cell_fn[100];

  if ( argc > 1 )
    read_input( argv[1], &nrestart, &avvol, &temp, &timestep, &split, &spatial, &rnum, &rstart, &dr, xcom_fn, ycom_fn, zcom_fn, chgs_fn, cell_fn );

  double *xcom, *ycom, *zcom, *chgs, *cell;

  analyze_file ( xcom_fn, &ncol, &nlns, delim );
  xcom = read_file_double ( "../data/xcom.dat", nlns, ncol, delim );

  analyze_file ( ycom_fn, &ncol, &nlns, delim );
  ycom = read_file_double ( "../data/ycom.dat", nlns, ncol, delim );

  analyze_file ( zcom_fn, &ncol, &nlns, delim );
  zcom = read_file_double ( "../data/zcom.dat", nlns, ncol, delim );

  analyze_file ( chgs_fn, &qcol, &nchg, delim );
  chgs = read_file_double ( "../data/charges.dat", nchg, qcol, delim );

  analyze_file ( cell_fn, &ccol, &ncll, delim );
  cell = read_file_double ( "../data/cell.dat", ncll, ccol, delim );

  printf("NCHGS: %i, NCOL: %i\n", nchg, ncol);
  if ( nchg != ncol ) {
    print_error ( FATAL, "Dimensions of charge array and position data don't match", __FILE__, __LINE__ );
    return FATAL;
  }

  if ( ncll != nlns ) {
    print_error ( FATAL, "Dimensions of cell array and position data don't match", __FILE__, __LINE__ );
    return FATAL;
  }

  /* data I looks fine
  for ( i=1; i<2; i++ ) {
    for ( j=0; j<nlns; j++ ) {
      printf("%14.8f\n", ael(xcom, nlns, i, j));
    }
  }
  */

  if ( split ) {

    int arrlen = nlns * rnum;

    double *qflux_neinst = (double *) calloc ( arrlen, sizeof(double));
    double *qflux_catcat = (double *) calloc ( arrlen, sizeof(double));
    double *qflux_anicat = (double *) calloc ( arrlen, sizeof(double));
    double *qflux_aniani = (double *) calloc ( arrlen, sizeof(double));

    get_qflux_srtd ( qflux_neinst, qflux_catcat, qflux_anicat, qflux_aniani, xcom, ycom, zcom, chgs, ncol, nlns, nrestart, dr, rstart, rnum, cell );
    // get_qflux_srtd ( qflux_neinst, qflux_catcat, qflux_anicat, qflux_aniani, xcom, ycom, zcom, chgs, ncol, nlns, nrestart);

    write_array_to_file ( "cond_neinst.out", qflux_neinst, 1, nlns );
    write_array_to_file ( "cond_catcat.out", qflux_catcat, rnum, nlns );
    write_array_to_file ( "cond_anicat.out", qflux_anicat, rnum, nlns );
    write_array_to_file ( "cond_aniani.out", qflux_aniani, rnum, nlns );

    printf("NUMBER OF COLUMNS: %i AND NUMBER OF LINES: %i\n", ncol, nlns);

    free ( qflux_neinst );
    free ( qflux_catcat );
    free ( qflux_anicat );
    free ( qflux_aniani );
  }
  else {
    double *qflux = (double *) calloc ( nlns, sizeof(double));

    get_qflux ( qflux, xcom, ycom, zcom, chgs, ncol, nlns, nrestart);

    write_array_to_file ( "cond_all.out", qflux, 1, nlns );

    printf("NUMBER OF COLUMNS: %i AND NUMBER OF LINES: %i\n", ncol, nlns);

    free ( qflux );
  }

  free ( cell );

  free ( xcom );
  free ( ycom );
  free ( zcom );

  return 0;

}

