// by alex sodt
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "util.h"
#include "pdb.h"
#include "dcd.h"
#include "ctype.h"
#include "alignSet.h"
#include "sys/time.h"
#include "input.h"
#include "hs.h"
#include "beta_util.h"

#include "fftw3.h"
#define SUB_PRESS
#define LEN 100
#define SKIPH

#define WATER_SAMPLE (2.0)

#define SELECT_Q

//#define DONT_SUB


double smooth_width = 3.0;
double h_bond_cutoff = 2.3;
double cutoff = 30.0;
double comp_cutoff = 3.0;
double dr_width = 0.5;

//#define DEBUG_POLYHEDRON
extern double *FracINV;

/* HBOND pairs:
 
O31 O32 O21 O22  glycerol
O11 O12 O13 O14  phospho
HN1 HN2 HN3      ethanolamine.

 */

double util_volume( double z1, double z2, double r1, double r2 );
double smoother( double r, double r_low, double r_high );
double max_r_cutoff = 150;
double sum_beta = 0;
int N_Z_BINS;
int N_R_BINS;
int N_MAX_R_BINS;
double dr_integration;
static double av_Lx,av_Ly,av_Lz;
double water_rho_average = 0;
void setupFit( double *fitQ, int fitnq, double *fitExp, double *fitExpErr, double *fitav_rho, double *fitc_beta_real, double **fitav_rho_res, pfft_instruction *fitinstructions, const char **fitresNames, int fitnres_used, double scale_factor );
void fit( void );

int main( int argc, char **argv )
{
          struct timeval tp;
 
          gettimeofday( &tp, NULL );
 
          srand((long)tp.tv_usec);

	char *buffer = (char *)malloc( sizeof(char) * 100000 );

	/*
 *	each parameter file can contain filenames for histograms, etc.
 *
 * 	*/

	if( argc < 2 )
	{
		printf("Syntax: betaPostProcess paramFile1 [paramFile2]\n");	
		return 0;
	}
		

	int mode = 0;
	

	struct parameter_block *blocks = NULL;
	double weight_sum = 0;

	for( int c = 1; c < argc; c++ )
	{
		struct parameter_block *ablock = (struct parameter_block *)malloc( sizeof(struct parameter_block) );
		ablock->post_process = 1;
		getInput(argv[c], ablock );
		weight_sum += ablock->weight;
		ablock->next = blocks;
		blocks = ablock;	
	}

	char *qvals_file = NULL;

	int not_set = 1;
	for( parameter_block *aBlock = blocks; aBlock; aBlock = aBlock->next )
	{
		if( not_set && aBlock->cutoff >= 0 )
		{
			cutoff = aBlock->cutoff;
			not_set = 0;
		}
		else if( aBlock->cutoff >= 0 && fabs(cutoff-aBlock->cutoff) > 1e-3 )
		{
			printf("Cutoff specification error: Different cutoffs in different input files.\n");
			exit(1);
		}

		if( aBlock->qvals && qvals_file == NULL )
			qvals_file = aBlock->qvals;
		else if( aBlock->qvals )
		{
			if( strcmp( aBlock->qvals, qvals_file ) )
			{
				printf("At least two different qvals files are specified, '%s' and '%s'.\n", aBlock->qvals, qvals_file );
				exit(1); 
			}
		}
	}

	for( parameter_block *ablock = blocks; ablock; ablock = ablock->next )
		ablock->weight /= weight_sum;

	double bz_width = 0.5;
//	double bz_width = 0.305864732;
	int N_Z_BINS = 200;
	N_R_BINS = ceil(cutoff/dr_width);
	N_MAX_R_BINS = max_r_cutoff / dr_width;
	dr_integration = max_r_cutoff / N_MAX_R_BINS;

	double *av_rho = (double *)malloc( sizeof(double) * N_Z_BINS );
	double water_rho_average;
	buildData( blocks, av_rho, N_Z_BINS, bz_width, max_r_cutoff, N_MAX_R_BINS, &water_rho_average );
	av_Lz = N_Z_BINS * bz_width;

#ifdef VERBOSE	
	printf("water_rho_average: %le\n", water_rho_average );
#endif
	double cyl_cutoff = -1;
	int cutoff_mode = -1;

	for( parameter_block *aBlock = blocks; aBlock; aBlock = aBlock->next )
	{
		if( cutoff_mode == -1 )
		{
			cyl_cutoff = aBlock->cyl_cutoff;
			cutoff_mode = aBlock->cutoff_mode;
		}

		if( aBlock->cutoff_mode != cutoff_mode || (aBlock->cyl_cutoff-cyl_cutoff) > 1e-10 )
		{
			printf("Cylindrical cutoff inconsistency.\n");
			exit(1);
		}
	}
	
	for( int iz = 0; iz < N_Z_BINS; iz++ )
	{
		sum_beta += av_rho[iz];
	}
	sum_beta /= N_Z_BINS;
	sum_beta *= av_Lz;

	fftw_complex *h_in;
	fftw_complex *h_out;
	fftw_plan p;

	int inflation = 40;
	int inflated_z = N_Z_BINS * inflation;

	h_in = (fftw_complex *)fftw_malloc( sizeof(fftw_complex) * inflated_z );
	h_out = (fftw_complex *)fftw_malloc( sizeof(fftw_complex) * inflated_z );

	p = fftw_plan_dft_1d( inflated_z, h_in, h_out, FFTW_FORWARD, FFTW_ESTIMATE ); 

	for( int j = 0; j < inflated_z; j++ )
	{
		h_in[j][0] = 0;
		h_in[j][1] = 0;
	}

	int xj = 0;

	for( int j = inflated_z/2-N_Z_BINS/2; j < inflated_z/2+N_Z_BINS/2; j++ )
	{
		h_in[j][0] = av_rho[xj];
		h_in[j][1] = 0;
		xj++;
	}
	
	fftw_execute(p);
	
//	FILE *fzFile = fopen("ft.z", "w");
	double max_q = 0.5;
	int n_q = max_q / (2 * M_PI / (inflation*av_Lz) );

	double val_q[n_q];
	double val_iz[n_q];

	int use_q = 0;

	double dz = av_Lz / N_Z_BINS;

	for( int j = 1; j < inflated_z; j++ )
	{
		double q = j * 2 * M_PI / (inflation*av_Lz);

		if( q < max_q && use_q < n_q )
		{
			val_q[use_q] = q;
			val_iz[use_q] = M_PI*dz*dz*(h_out[j][0]*h_out[j][0] + h_out[j][1]*h_out[j][1]) / q / q;
			use_q++;
		}
	}


	double *c_beta_tot = (double *)malloc( sizeof(double) * N_MAX_R_BINS );
	memset( c_beta_tot, 0, sizeof(double) * N_MAX_R_BINS );

	double av_Lx = 1.0;
	double av_Ly = 1.0;

	double av_area = blocks->av_area;
//	av_area = 10479.71439;

#ifdef VERBOSE
	printf("beta_per_unit_area: %lf\n", sum_beta );
#endif
	double total_area = 0;
	double total_vol = 0;
	double nw = 0;
	double total_weight = 0;
		
	for( parameter_block *block = blocks; block; block = block->next )
	{
//		these don't work, we get them from the betaz file.
//		total_area += block->av_bz_area * block->weight;
//		total_vol  += block->av_bz_vol * block->weight;

//		these will be over-ridden by the histogram:
		total_area += block->av_area * block->weight;
		total_vol += block->av_area * block->height * block->weight;

		nw += block->num_water * block->weight;
		total_weight += block->weight;
	}



	nw /= total_weight;

	double slope_sub = (total_area) * 2 * M_PI * sum_beta * sum_beta;

	for( int ir = 0; ir < N_MAX_R_BINS; ir++ )
	{
		double r1 = (ir) * max_r_cutoff / N_MAX_R_BINS;
		double r2 = (ir+1) * max_r_cutoff / N_MAX_R_BINS;
		double dr = max_r_cutoff / N_MAX_R_BINS;

		double continuum_value = 0, real_value = 0, continuum_value_compensated = 0, continuum_value_no_smooth=0, sph_continuum_value = 0;
		
		real_value = 0;

		for( parameter_block *block = blocks; block; block = block->next )
		{
			if( ir < block->max_bin )
			{ 	
				double vol = (4.0/3.0) * M_PI * (r2*r2*r2-r1*r1*r1);
				double l_real = block->beta_hist[ir];// * ( block->av_area / block->av_bz_area); 
				
/*				double layer_vol = block->av_bz_area * av_Lz / N_Z_BINS;

				double bg_bg = 0;
				double corr_cross = 0;

				for( int iz = 0; iz < N_Z_BINS; iz++ )
				{
					double zv = iz * av_Lz / N_Z_BINS;
	
					double bg_vol = util_volume( -zv, av_Lz-zv, r1, r2 ); 
					
					corr_cross += water_rho_average * layer_vol * (vol-bg_vol) * water_rho_average;

					bg_bg += water_rho_average * layer_vol * vol * water_rho_average;
				}

				l_real -= corr_cross;
				l_real += bg_bg;
*/
				real_value += l_real * block->weight;
			}
		}
	
		for( int z = 0; z < N_Z_BINS; z++ )
		{
			double zv = -av_Lz/2 + (z+0.5) * av_Lz / N_Z_BINS;
			double zvol = total_area * (av_Lz/N_Z_BINS);
		
			double bv = av_rho[z];

			for( int z2 = 0; z2 < N_Z_BINS; z2++ )
			{
				double bin_z1 = -av_Lz/2 + z2 * av_Lz / N_Z_BINS - zv;	
				double bin_z2 = -av_Lz/2 + (z2+1) * av_Lz / N_Z_BINS - zv;	
	
				double r1s = r1;;
				double r2s = r2;;

				double approx_vol;
				int use_z2 = z2;

				if( cutoff_mode )
				{
					approx_vol = util_volume( bin_z1, bin_z2, r1s, r2s,  0, cyl_cutoff );
					double approx_vol_c = util_volume( bin_z1, bin_z2, r1s, r2s, cyl_cutoff, r2s );
//					double approx_vol_check = util_volume( bin_z1, bin_z2, r1s, r2s );
//					if( approx_vol_check > 1e-3 )
//					printf(" vol in %lf vol out %lf vol sum %lf vol check %lf\n", approx_vol, approx_vol_c, approx_vol+approx_vol_c, approx_vol_check );
					continuum_value += bv * av_rho[use_z2] * approx_vol * zvol;	
					continuum_value_no_smooth += bv * av_rho[use_z2] * approx_vol_c * zvol;	
					sph_continuum_value += bv * av_rho[use_z2] * approx_vol * zvol;	
				}
				else
				{
					approx_vol = util_volume( bin_z1, bin_z2, r1s, r2s );
					continuum_value += bv * av_rho[use_z2] * approx_vol * zvol;	
				}
			}
		}

		// solvent has already been subtracted from the continuum model, there is no problem.

		double rmean = (r1+r2)/2;

		double alpha = smoother( rmean, cutoff-smooth_width, cutoff );	
	
		if( cutoff_mode )
		{
			alpha = smoother( rmean, cutoff-smooth_width, cutoff );	
				
		}
		else if( rmean >= cutoff )
			alpha = 0.0;
		double del = slope_sub * rmean;
//#define VERBOSE
	//	c_beta_tot[ir] = alpha * real_value/(2*dr_integration) + (1-alpha) * continuum_value/(2*dr_integration) - del; //continuum_value;
		if( cutoff_mode )
			c_beta_tot[ir] = real_value/(dr_integration) + continuum_value_no_smooth / (dr_integration) - del; //continuum_value;
		else
			c_beta_tot[ir] = alpha * real_value/(dr_integration) + (1-alpha) * continuum_value/(dr_integration) + continuum_value_no_smooth / (dr_integration) - del; //continuum_value;
#ifdef VERBOSE
		if( cutoff_mode )
			printf("r: %lf real: %lf continuum: %lf continuum_out: %lf asymptote-continuum_out: %lf val: %lf\n", rmean,
			real_value/(dr_integration), (continuum_value)/(dr_integration), (continuum_value_no_smooth)/(dr_integration), del-  (continuum_value_no_smooth)/(dr_integration), c_beta_tot[ir] );
		else
			printf("r: %lf real: %lf continuum: %lf asymptote: %lf val: %lf\n", rmean,
		real_value/(dr_integration), (continuum_value+continuum_value_no_smooth)/(dr_integration), del, c_beta_tot[ir] );
#endif
	}		
		
	// we have the scattering length per unit area, including contrast
	// changing the fraction d2o changes the "solvent subtraction"
	// it also changes the b/area
	
	double b_per_d2o = 2 * 0.6671 + 0.5804;
	double b_per_h2o = 2 * -0.3742 + 0.5804; 

	double beta_sub_d2o =   0.065;
	double beta_sub_df  = - 0.072; // per fraction.
	double b_area_f0 = sum_beta;


	double f_opt = 1- total_area * sum_beta / ( b_per_d2o * nw - b_per_h2o * nw + beta_sub_df * total_vol);

#ifdef VERBOSE	
	printf("Perfect cancellation at f=%lf D2O.\n", f_opt );
#endif

//	for( double q = 0.001; q < 0.5; q += 0.001)
#ifdef VERBOSE
	printf("q algorithm 1d\n");
#endif
	double eps = 1e-3;

double dq = 0.002;
#ifdef SELECT_Q
	FILE *qFile = NULL;
	if( qvals_file )
	{
		qFile = fopen(qvals_file,"r");
		if( !qFile )
		{
			printf("Couldn't open qvals file '%s'.\n", qvals_file );
			exit(1);
		}
	}
	double *select_q;
	if( qFile )
	{
		int nvals = 0;
		while( !feof(qFile) )
		{
			getLine(qFile, buffer);
			if( feof(qFile) ) 
				break;
			double tval;
			int nr =sscanf( buffer, "%lf", &tval );			
			if( nr == 1 )
				nvals++;
		}		
		rewind(qFile);
		use_q = 0;
		select_q = (double *)malloc( sizeof(double) * nvals );
		while( !feof(qFile) )
		{
			getLine(qFile, buffer);
			if( feof(qFile) ) 
				break;
			double tval;
			int nr =sscanf( buffer, "%lf", &tval );			
			if( nr == 1 )
			{
				select_q[use_q] = tval;
				use_q++;
			}
		}	
		fclose(qFile);	
	}
	else 
	{
		use_q = 0.5 / dq; 
		select_q = (double *)malloc( sizeof(double) * use_q );
		for( int iq = 0; iq < use_q; iq++ )
			select_q[iq] = (1+iq)*dq;
	}

#endif

#ifdef VERBOSE
	printf("q0 %lf %le\n", 0.0, slope_sub/total_area );
#endif
	for( int iq = 0; iq < use_q; iq++ )
	{
		double expVal = 1.0;
#ifdef SELECT_Q
		double q = select_q[iq];
#else
		double q = val_q[iq];
#endif
		double val = 0;

		double dr = max_r_cutoff / N_MAX_R_BINS;
		for( int ir = 0; ir < N_MAX_R_BINS; ir++ )
		{
			double r;

			if( ir == 0 )
				r = (ir+eps) * max_r_cutoff / N_MAX_R_BINS;
			else
				r = (ir+0.5) * max_r_cutoff / N_MAX_R_BINS;

			val += dr * c_beta_tot[ir] * sin( q * r ) / (q * r);	
		}

		printf("%lf %lf\n", q, (slope_sub/q/q+val)/total_area  ); 
//			printf("%lf %lf %lf %lf\n", q, val_iz[iq] + val/av_Lx/av_Ly, val_iz[iq] );
//		else	
//			printf("%lf %lf %lf\n", q, (slope_sub/q/q + val)/av_Lx/av_Ly, val_iz[iq] );
	}
//	printf("FFT %d\n", n_q);
	if( blocks && blocks->print_fft )
	{
		for( int iq = 0; iq < n_q; iq++ )
			printf("%le %le\n", val_q[iq], val_iz[iq] );
	}	
}

