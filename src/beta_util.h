#ifndef __betautilh__
#define __betautilh__

static const char *expec_apl_name[] = { "CHL1" };
static const double expec_apl[] = { 40.0 };
static const double default_apl = 65.0;

double util_volume( double z1, double z2, double r1, double r2 );
double util_volume( double z1, double z2, double r1, double r2, double rc1, double rc2 );
double smoother( double r, double r_low, double r_high );
void buildData( struct parameter_block *blocks, double *betaz, int N_Z_BINS, double bz_width, double max_r_cutoff, int N_MAX_R_BINS, double *water_rho_average);
void read_residue_bz( const char *fileName, int nres, const char **names, double *max_z, int *nb, double **betaz_out );

void disk_disk_hist( double *beta1z, double *beta2z, double dzbin, int nzbins, double r_hist, double *hist, double bin_width, int nrbins );
void r2_disk_disk_hist( double *beta1z, double *beta2z, double dzbin, int nzbins, double r1_disk1, double r_disk2, double *hist, double bin_width, int nrbins );

#endif
