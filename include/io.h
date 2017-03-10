int analyze_file ( char *fname, int *ncol, int *nlns, char *delim);
double *read_file_double(char *fname, int nlns, int ncol, char *delim);
void write_array_to_file ( char *fname, double *a, int ncol, int nlns );
void read_input( char *fname, int *nrestart, double *avvol, double *temp, double *timestep, int *split, int *spatial, int *rnum, double *rstart, double *dr, char *xcom_fn, char *ycom_fn, char *zcom_fn, char *chgs_fn, char *cell_fn, double *fitoffset, double *fitlength );
