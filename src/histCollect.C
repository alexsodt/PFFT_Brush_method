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
double dr_width = 0.25;

//#define DEBUG_POLYHEDRON
extern double *FracINV;

/* HBOND pairs:
 
O31 O32 O21 O22  glycerol
O11 O12 O13 O14  phospho
HN1 HN2 HN3      ethanolamine.

 */

double util_volume( double z1, double z2, double r1, double r2 );
double smoother( double r, double r_low, double r_high );
double max_r_cutoff = 100;
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

	if( argc < 3 )
	{
		printf("Syntax: histCollect output.name gofrFileName1 [gofrFileName2 ...]\n");	
		return 0;
	}
	

	double av_Lx=0, av_Ly=0, av_Lz=0;
	int prev_l_N_R_BINS = -1;
	double prev_l_cutoff = -1;
	
	double n_g_av = 0;

	int natoms_use = -1;
	int MAX_HIST = 100;

	double *IV = NULL;
	double *bi_hist[MAX_HIST];
	double *bi_hist_cn = NULL;
	double bi_length[MAX_HIST];
	memset( bi_length, 0, sizeof(double) * MAX_HIST );
	double n_hist[MAX_HIST];
	memset( n_hist, 0, sizeof(double) * MAX_HIST );
	int init = 1;
	char atom_names[MAX_HIST];
		
	int cutoff_mode = 0;
	double cyl_cutoff=0;

	int nIV = 0;

	for( int c = 2; c < argc; c++ )
	{
		FILE *hist = fopen(argv[c],"r");

		int l_N_R_BINS;
		double l_cutoff, l_av_Lx, l_av_Ly, l_av_Lz;

		getLine( hist, buffer );
		int nr = sscanf(buffer, "%d %lf %lf %lf %lf cylinder %lf\n", &l_N_R_BINS, &l_cutoff, &l_av_Lx, &l_av_Ly, &l_av_Lz, &cyl_cutoff );	

		if( nr == 6 )
			cutoff_mode = 1;

		prev_l_cutoff = l_cutoff;
	
		av_Lx += l_av_Lx;
		av_Ly += l_av_Ly;
		av_Lz += l_av_Lz;
		
		if( prev_l_N_R_BINS > 0 && l_N_R_BINS != prev_l_N_R_BINS )
		{
			printf("Bin schemes of input histograms don't match.\n");
			return -1;
		}

		prev_l_N_R_BINS = l_N_R_BINS;

		getLine( hist, buffer );

		double l_nhist[MAX_HIST];
		int natoms  = readNDoubles( buffer+strlen("nhist "), l_nhist, MAX_HIST );
		natoms_use = natoms;
		for( int i = 0; i < natoms; i++ )
			n_hist[i] += l_nhist[i];

		getLine( hist, buffer );
		readNDoubles( buffer+strlen("blength "), bi_length, MAX_HIST );

		double temp[l_N_R_BINS];
		if( init )
		{
			bi_hist_cn = (double *)malloc( sizeof(double) * l_N_R_BINS );
			memset( bi_hist_cn, 0, sizeof(double) * l_N_R_BINS );
			for( int x = 0; x < natoms; x++ )
			{
				bi_hist[x] = (double*)malloc( sizeof(double) * l_N_R_BINS );
				memset( bi_hist[x], 0, sizeof(double) * l_N_R_BINS );
			}
			IV = (double *)malloc( sizeof(double) * l_N_R_BINS );
			memset( IV, 0, sizeof(double) * l_N_R_BINS );
			init = 0;
		}

		getLine( hist, buffer );
		readNDoubles( buffer + strlen("cn "), temp, l_N_R_BINS );

		for( int i = 0; i < l_N_R_BINS; i++ )
			bi_hist_cn[i] += temp[i];

		for( int x = 0; x < natoms; x++ )
		{
			getLine( hist, buffer );
			atom_names[x] = buffer[0];
			readNDoubles( buffer + 2, temp, l_N_R_BINS );
			for( int i = 0; i < l_N_R_BINS; i++ )
				bi_hist[x][i] += temp[i];
		}

		getLine( hist, buffer );

		if( !strncasecmp( buffer, "IV", 2) )
		{
			readNDoubles( buffer + 2, temp, l_N_R_BINS );
			for( int i = 0; i < l_N_R_BINS; i++ )
				IV[i] += temp[i];
			nIV++;
		}

		fclose(hist);

		n_g_av += 1;	
	}
	
	FILE *output = fopen(argv[1],"w");

	if( !output )
	{
		printf("Couldn't open file '%s' to write.\n", argv[1] );
		return -1;
	}
	
	if( cutoff_mode )	
		fprintf(output, "%d %lf %lf %lf %lf cylinder %lf\n", prev_l_N_R_BINS, prev_l_cutoff, av_Lx/n_g_av, av_Ly/n_g_av, av_Lz/n_g_av, cyl_cutoff );	
	else
		fprintf(output, "%d %lf %lf %lf %lf\n", prev_l_N_R_BINS, prev_l_cutoff, av_Lx/n_g_av, av_Ly/n_g_av, av_Lz/n_g_av);	

		double l_nhist[MAX_HIST];
		int natoms  = readNDoubles( buffer+strlen("nhist "), l_nhist, MAX_HIST );
	
	fprintf(output, "nhist");
	for( int i = 0; i < natoms_use; i++ )
		fprintf(output, " %lf", n_hist[i] / n_g_av );
	fprintf(output,"\n");

	fprintf(output, "blength");
	for( int i = 0; i < natoms_use; i++ )
		fprintf(output, " %lf", bi_length[i] );
	fprintf(output,"\n");


	fprintf(output, "cn");
	for( int i = 0; i < prev_l_N_R_BINS; i++ )
		fprintf(output, " %lf", bi_hist_cn[i] / n_g_av );
	fprintf(output, "\n");	
		
	for( int x = 0; x < natoms_use; x++ )
	{
		fprintf(output, "%c", atom_names[x] );

		for( int i = 0; i < prev_l_N_R_BINS; i++ )
			fprintf(output, " %lf", bi_hist[x][i] / n_g_av );
		fprintf(output, "\n");	
	}

	if( nIV > 0 )
	{
		fprintf(output, "IV");
		for( int i = 0; i < prev_l_N_R_BINS; i++ )
			fprintf(output, " %lf", IV[i] / nIV );	
		fprintf(output, "\n");
	}

	return 0;
}
