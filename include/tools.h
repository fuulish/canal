
void multiply_array_array ( double *out, double *a, double *b, int ncol, int nlns );
void multiply_array_array_inplace ( double *a, double *b, int ncol, int nlns );
void divide_array_number ( double *out, double *a, double f, int ncol, int nlns );
void divide_array_number_inplace (double *out, double f, int ncol, int nlns);
void multiply_array_number ( double *out, double *a, double f, int ncol, int nlns );
void multiply_array_number_inplace (double *out, double f, int ncol, int nlns);
void add_array_number ( double *out, double *a, double f, int ncol, int nlns );
void add_array_number_inplace ( double *out, double f, int ncol, int nlns );
void subtract_array_number ( double *out, double *a, double f, int ncol, int nlns );
void add_array_array ( double *out, double *a, double *b, int ncol, int nlns );
void add_arrays_inplace ( double *b, double *a, int ncol, int nlns );
void divide_array_array ( double *out, double *a, double *b, int ncol, int nlns );
void divide_array_array_inplace ( double *a, double *b, int ncol, int nlns );
void subtract_array_array ( double *out, double *a, double *b, int ncol, int nlns );
void subtract_arrays_inplace ( double *a, double *b, int ncol, int nlns );
