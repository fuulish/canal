
enum diff_e {
  FWD=1,
  BWD,
  CDD,
};

double calculate_dot ( double *xi, double *yi, double *zi, double *xj, double *yj, double *zj, int nlns );
double calculate_nrm ( double *xi, double *yi, double *zi, int nlns );
void calculate_velocities ( double *vx, double *vy, double *vz, double *xi, double *yi, double *zi, int nlns, int timestep );
void get_vflux_locl ( double *neinst, double *cnd_cc, double *cnd_ac, double *cnd_aa, double *x, double *y, double *z, double *chg, int ncol, int nlns, int nrestart, double dr, double rstart, int rnum, double *cell );
