int analyze_file ( char *fname, int *ncol, int *nlns, char *delim);
double *read_file_double(char *fname, int nlns, int ncol, char *delim);
void write_array_to_file ( char *fname, double *a, int ncol, int nlns );
