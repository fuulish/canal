
void calculate_msd_one ( double *out, double *x, double *xo, int nlns );
void calculate_msd_xyz ( double *out, double *x, double *y, double *z, int nlns );
void calculate_msd_xyz_cross ( double *out, double *xi, double *yi, double *zi, double *xj, double *yj, double *zj, int nlns);
void get_qflux ( double *cnd, double *x, double *y, double *z, double *chg, int ncol, int nlns, int nrestart );
void get_qflux_srtd ( double *neinst, double *cnd_cc, double *cnd_ca, double *cnd_aa, double *x, double *y, double *z, double *chg, int ncol, int nlns, int nrestart, double dr, double rstart, int rnum, double *cell );
