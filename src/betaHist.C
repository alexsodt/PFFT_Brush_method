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
#include "spline.h"
#ifdef FFTW
#include <fftw3.h>
#endif
#include "beta_util.h"
#include "input.h"

#define DISABLE_INTRA
#define N_Z_BINS 150
#define MAX_RES 100
//#define DONT_SUB
double smooth_width = 3.0;
	double dr_width = 0.5;

/* HBOND pairs:
 
O31 O32 O21 O22  glycerol
O11 O12 O13 O14  phospho
HN1 HN2 HN3      ethanolamine.

 */

double smoother( double r, double r_low, double r_high );
double EffVol( char atname );

int main( int argc, char **argv )
{
	if( argc < 4  )
	{
		printf("Syntax: betaHist input_file psf dcd1 [dcd2 ...]\n");	
		return 0;
	}
          
	struct timeval tp;
 
          gettimeofday( &tp, NULL );
 
          srand((long)tp.tv_usec);

	char *buffer = (char *)malloc( sizeof(char) * 100000 ); // buffer for reading in scattering file
		
	double water_rho_average = 0;
	double av_rho[N_Z_BINS];
	memset( av_rho, 0, sizeof(double) * N_Z_BINS );
	
	double av_rho_get[N_Z_BINS];
	memset( av_rho_get, 0, sizeof(double) * N_Z_BINS );
	
	double max_r_cutoff = 150;
	int N_MAX_R_BINS = max_r_cutoff / dr_width;

	int mode = 0;
	double cutoff = 10.0;
		

//Handle input file

	parameter_block *block = (parameter_block *)malloc( sizeof(parameter_block) );
	block->post_process = 0;
	getInput(argv[1], block );
	block->next =NULL;
	
	int cutoff_mode = 1;
	
	printf("cutoff: %le\n", block->cutoff );

	if( block->cutoff > 0 )
	{
		cutoff = fabs(block->cutoff);
		double bz_width = 0.5;
		buildData( block, av_rho, N_Z_BINS, bz_width, max_r_cutoff, N_MAX_R_BINS, &water_rho_average );

		mode = 0;
	}
	else
		mode = 1;

	int disable_histogramming = 0;

	if( cutoff < 0.01 )
	{
		printf("Disabling particle-particle interactions: cutoff very small.\n");
		disable_histogramming = 1;	
	}

	const char *uniqueName = block->name;
	double strain = block->apply_strain;
	double lateral_scale = block->lateral_scale;
//Done with input file


	printf("water_rho_average: %le\n", water_rho_average );
	fflush(stdout);

	double outside_z_box_rho = water_rho_average;


	FILE *psfFile = fopen(argv[2], "r" );

	if( ! psfFile )
	{
		printf("Couldn't open PSF file '%s'.\n", argv[2] );
		return 0;
	}

	if( !strcasecmp( argv[2] + strlen(argv[2])-3, "pdb" ) )
		loadPSFfromPDB( psfFile );	
	else
		loadPSF( psfFile );

	fclose(psfFile);	
	
// set the exchangeable proton's scattering length based on solvent.


	int nat = curNAtoms();

//FIXME: Needs to be freed!
	struct atom_rec *at = (struct atom_rec *)malloc( sizeof(struct atom_rec) * curNAtoms() );

	
	int nx = -1;
	int ny = -1;
	int nz = -1;

struct binnit
{
	int *id;
	int npos;
	int nposSpace;
};
	binnit *bins;

	int init_done = 0;

//FIXME: Needs to be freed!
	double *b_i = (double *)malloc( sizeof(double) * curNAtoms() );
//FIXME: Needs to be freed!
	int    *b_ind = (int *)malloc( sizeof(int) * curNAtoms() );
//FIXME: Needs to be freed!
	int    *b_res = (int *)malloc( sizeof(int ) * curNAtoms() );		
				
	// Y: exchangeable proton takes the average value from water.
	const char *bi_atoms=  "HDCNOPSXY123456789";// MUST go H then D for 'Z' flagged atom hard coded below.
	int natoms = strlen(bi_atoms);

// the indexes of proton and deuteron.
#define PROTON_B 0
#define DEUTERON_B 1

	int nlength = 0;
	double *bi_length = (double *)malloc( sizeof(double) * curNAtoms() );
	memset( bi_length, 0, sizeof(double) * curNAtoms() );
#ifdef DO_SAXS
//	double bi_length[17] = { 1, 1, 6, 7, 8, 15, 16, 1, 1,1,1,1,1,1,1,1,1 };
#else
//	double bi_length[17] = { -0.3742, 0.6671, 0.6651, 0.940, 0.5804, 0.517, 0.2847, 0.531731 /* 87% D 13% H */, 0.0,0,0,0,0,0,0,0,0 /* Z1-9, special flagged hydrogen will be overwritten*/ };

	bi_length[0] = -0.3742;
	bi_length[1] = 0.6671;
	bi_length[2] = 0.6651;
	bi_length[3] = 0.940;
	bi_length[4] = 0.5804;
	bi_length[5] = 0.517;
	bi_length[6] = 0.2847;
	bi_length[7] = 0.531731;
	bi_length[8] = 0.6671;

	for( parameter_block *bl = block; bl; bl = bl->next )
	{
		if( fabs(bl->X_scattering_length) < PFFT_HUGE )
		{
			bi_length[7] = bl->X_scattering_length;
		}
	}

	double av_hd_scat = 0;
	double nav_hd_scat = 0;

	for( int a = 0; a < curNAtoms(); a++ )
	{
		char *atname, *resname;

		getAtName( a, &atname);
		getResName( a, &resname );

		if( !strcasecmp( resname, "TIP3" ) )
		{
			if( atname[0] == 'D' )
			{
				av_hd_scat += bi_length[1];
				nav_hd_scat++; 
			}
			if( atname[0] == 'H' )
			{
				av_hd_scat += bi_length[0];
				nav_hd_scat++; 
			}
		}

		free(atname);
		free(resname);	
	}

	av_hd_scat /= nav_hd_scat;
	bi_length[8] = av_hd_scat;

	printf("Exchangeable proton scattering set to %lf (%lf%% D2O).\n", av_hd_scat, 100 * (av_hd_scat + 0.3742) * 0.960338 );
#endif
//	double bi_length[7] = { 1, 1, 1, 1, 1, 1, 1 };
	double *use_bi_length = bi_length;
	const char *martini_atoms[] =
	{
		"NC3",
		"PO4",
		"GLD",
		"GL",
		"C",
		"D",
		"W",
		"WD",
		"ND3",
		"ROH",
		"R"
	};

	double martini_lengths[] =
	{	
		-0.5158, // choline
		2.672, // phosphate
		4.4913,	// ester deuterated
		1.88805, // ester
		-0.3332, // methylene, saturated
		0.4152, // methylene, flexible
		-0.672, // water
		7.6584, // water, deuterated
	 	9.6597, // choline, deuterated
		-0.3332,  // cholesterol ROH group, arb
		-0.3332,  // cholesterol R group, arb
	};
	
	int natoms_martini = sizeof(martini_lengths)/sizeof(double);
	int max_natoms = natoms;
	if( natoms_martini > max_natoms )
		max_natoms = natoms_martini;
	
	// to handle generic deuteration schemes.
	max_natoms = 1000;

	double  n_hist[MAX_RES][max_natoms];

	double ***b_hist = NULL;
	double ***b_rho = NULL;

//	double *b_hist[MAX_RES][max_natoms];
//	double *b_rho[MAX_RES][max_natoms];
	double nframe_average = 0;

	int N_R_BINS = ceil(cutoff / dr_width);

	
	if( cutoff_mode )
		N_R_BINS = N_MAX_R_BINS;

//FIXME: Needs to be freed!
	double *beta_tot = (double *)malloc( sizeof(double) * N_MAX_R_BINS );
	memset( beta_tot, 0, sizeof(double) * N_MAX_R_BINS );
	
	
//FIXME: Needs to be freed!
	double *termIV = ( double *)malloc( sizeof(double) * N_R_BINS );
	memset( termIV, 0, sizeof(double) * N_R_BINS );
	
	double **b_hist_cn;

/*
	for( int r = 0; r < MAX_RES; r++ )
	{
//FIXME: Needs to be freed!
		b_hist_cn[r] = ( double *)malloc( sizeof(double) * N_R_BINS );
		memset( b_hist_cn[r], 0, sizeof(double) * N_R_BINS );

		for( int a = 0; a < max_natoms; a++ )
		{
//FIXME: Needs to be freed!
			b_rho[r][a] = (double *)malloc( sizeof(double) * N_Z_BINS );
			memset( b_rho[r][a], 0, sizeof(double) * N_Z_BINS );

//FIXME: Needs to be freed!
			b_hist[r][a] = (double *)malloc( sizeof(double) * N_R_BINS );
			memset( b_hist[r][a], 0, sizeof(double) * N_R_BINS );

			n_hist[r][a] = 0;
		}
	}
*/
	double av_Lx = 0;
	double nav_Lx = 0;
	
	double av_Ly = 0;
	double nav_Ly = 0;

	double av_Lz = 0;
	double nav_Lz = 0;

//FIXME: Needs to be freed!
        double *cur_set = (double *)malloc( sizeof(double) * 3 * nat );
//FIXME: Needs to be freed!
        double *last_align = (double *)malloc( sizeof(double) * 3 * nat );
        int aligned = 0;

	double nframes_gofr = 0;
	int gofr_skip = 1;

	double sub_value = 0;
						
	double alt_cn[N_R_BINS];
	memset( alt_cn, 0, sizeof(double) * N_R_BINS );

	int pick_zmin = 0;///2 -50;
	int pick_zmax = N_Z_BINS;
	double pick_beta = 0;

	char *resNames[MAX_RES];
	int nres_used = 0;
	int n_res[MAX_RES];
	memset( n_res, 0, sizeof(int) * MAX_RES );

	int martini_detected = 0;

	int all_f_for_skip = 0;

	// average box dimensions for area compensation, if requested.
	double cLx=0, cLy=0, cLz=0;

	if( block->areaCompensate )
	{
		// number of frames averaged over in compensation.
		double nComp = 0;

		for( int c = 3; c < argc; c++ )
		{
			FILE *dcdFile = fopen(argv[c], "r");
	
			if( ! dcdFile )
			{
				printf("Couldn't open dcd file '%s'.\n", argv[c] );
				return 0;
			}
		
			printf("Processing %s to compute area compensation.\n", argv[c]);
			fflush(stdout);
		
			readDCDHeader(dcdFile);
			setFractional();
			int nframes = curNFrames();
			for( int f = 0; f < nframes; f++ )
			{
				printf("frame %d/%d.\n", f, nframes );
				fflush(stdout);
				loadFrame( dcdFile, at );

				if( !DCDsuccess() )
				{
					nframes = f;
					break;
				}
			
				double La,Lb,Lc,alpha,beta,gamma;
				PBCD( &La, &Lb, &Lc, &alpha, &beta, &gamma );
			
				cLx += La;
				cLy += Lb;
				cLz += Lc;
				nComp += 1;
	
				for( int a = 0; a < curNAtoms(); a++ )
					at[a].zap();
			}

			fclose(dcdFile);
		}

		cLx /= nComp;
		cLy /= nComp;
		cLz /= nComp;
	}

	int warning_trigger = 0;

	int *lipid_start = (int *)malloc( sizeof(int) * curNAtoms() );
	int *lipid_stop = (int *)malloc( sizeof(int) * curNAtoms() );
	int nlipids = 0;
	int *upper_leaflet_lipids = (int *)malloc( sizeof(int) * curNAtoms() );
	int *lower_leaflet_lipids = (int *)malloc( sizeof(int) * curNAtoms() );
	int nupper = 0;
	int nlower = 0;


	for( int c = 3; c < argc; c++ )
	{
		FILE *dcdFile = fopen(argv[c], "r");
	
		if( ! dcdFile )
		{
			printf("Couldn't open dcd file '%s'.\n", argv[c] );
			return 0;
		}
	
		printf("Processing %s.\n", argv[c]);
		fflush(stdout);

		readDCDHeader(dcdFile);
		setFractional();
		int nframes = curNFrames();
		for( int f = 0; f < nframes; f++, all_f_for_skip++ )
		{
			printf("frame %d/%d.\n", f, nframes );
			fflush(stdout);
			loadFrame( dcdFile, at );	

			if( !DCDsuccess() )
			{
				nframes = f;
				break;
			}


			
			double La,Lb,Lc,alpha,beta,gamma;
			PBCD( &La, &Lb, &Lc, &alpha, &beta, &gamma );

			if( block->areaCompensate )
			{
				double alpha = log((cLx*cLy)/(La*Lb));
	
				apply_strain( at, nat, &La, &Lb, &Lc, alpha );
			}
			else
				apply_strain( at, nat, &La, &Lb, &Lc, strain );

			apply_lateral_scale( at, nat, &La, &Lb, &Lc, lateral_scale );

/*			for( int x = 0; x < nat; x++ )
			{
				at[x].x *= lateral_scale;
				at[x].y *= lateral_scale;
			}
*/
			printf("La/b/c: %lf %lf %lf\n", La, Lb, Lc );

			av_Lx += La;
			nav_Lx += 1;
			av_Ly += Lb;
			nav_Ly += 1;
			av_Lz += Lc;
			nav_Lz += 1;

			double vol_scale = 1.0;
			double effvol = 0;
			for( int a = 0; a < curNAtoms(); a++ )
				effvol += EffVol(at[a].atname[0]);
			vol_scale = La*Lb*Lc/effvol;
			printf("vol_scale: %le\n", vol_scale );
			for( int ir = 0; ir < N_R_BINS; ir++ )
			{ 
				double r1 = ir * max_r_cutoff / N_R_BINS;	
				double r2 = (ir+1)* max_r_cutoff / N_R_BINS;

				double layer_vol = La * Lb * (Lc/N_Z_BINS);
	
				for( int iz = 0; iz < N_Z_BINS; iz++ )
				{
					double z = iz * Lc / N_Z_BINS;
				
					double bg_vol;
	
					if( cutoff_mode )
						bg_vol = util_volume( -z, Lc-z, r1, r2, 0, cutoff ); 
					else
						bg_vol = util_volume( -z, Lc-z, r1, r2 ); 

					termIV[ir] += bg_vol * layer_vol;
				}
			}

                        for( int x = 0; x < nat; x++ )
                        {
                                cur_set[x*3+0] = at[x].x;
                                cur_set[x*3+1] = at[x].y;
                                cur_set[x*3+2] = at[x].z;
                        }

                        if( aligned )
                        {
                                for( int p = 0; p < nat; p++ )
                                {
                                        while( cur_set[p*3+0] - last_align[p*3+0] < -0.5 ) cur_set[3*p+0] += 1;
                                        while( cur_set[p*3+0] - last_align[p*3+0] >  0.5 ) cur_set[3*p+0] -= 1;
                                        while( cur_set[p*3+1] - last_align[p*3+1] < -0.5 ) cur_set[3*p+1] += 1;
                                        while( cur_set[p*3+1] - last_align[p*3+1] >  0.5 ) cur_set[3*p+1] -= 1;
                                        while( cur_set[p*3+2] - last_align[p*3+2] < -0.5 ) cur_set[3*p+2] += 1;
                                        while( cur_set[p*3+2] - last_align[p*3+2] >  0.5 ) cur_set[3*p+2] -= 1;
                                }
                        }
                        
			for( int x = 0; x < nat; x++ )
                        {
                               at[x].x=cur_set[x*3+0];
                               at[x].y=cur_set[x*3+1];
                               at[x].z=cur_set[x*3+2];
                        }
			

			double carbon_com[3] = { 0,0,0};
			double ncarbon = 0;		
			double wrap_to[3] = { 0,0,0 };	

			for( int a = 0; a < nat; a++ )
			{
				if( !strcasecmp( at[a].resname, "CLA" ) ) continue;
				
				if( martini_detected )
				{
#ifdef DO_SAXS
					printf("PFFT-MARTINI not currently coded for SAXS.\n");
					exit(1);		
#endif
					if( !strcasecmp( at[a].atname, "D4A" ) || !strcasecmp( at[a].atname, "C4A" ) )
					{
						wrap_to[0] = at[a].x;
						wrap_to[1] = at[a].y;
						wrap_to[2] = at[a].z;
						break;
					}
				}
					
				if( at[a].atname[0] != 'C' ) continue;
				if( !strcasecmp( at[a].atname, "C210" ) || !strcasecmp( at[a].atname, "C12S") )
				{
					wrap_to[0] = at[a].x;
					wrap_to[1] = at[a].y;
					wrap_to[2] = at[a].z;
					break;
				}
			}
			
			for( int a = 0; a < nat; a++ )
			{
				while( at[a].z - wrap_to[2] < -0.5 ) at[a].z += 1.0;
				while( at[a].z - wrap_to[2] > 0.5 ) at[a].z -= 1.0;
			}

			double bilayer_center = 0;
		
#define N_BINS_MOLDIST 1000
			double moldist[N_BINS_MOLDIST];
			memset( moldist, 0, sizeof(double) * N_BINS_MOLDIST );
	
			for( int a = 0; a < curNAtoms(); a++ )
			{ 		
				char *tat = at[a].atname;
				char *trt = at[a].resname;
				
				double az = at[a].z;
				
				if( !strncasecmp( trt, "POT", 3 ) ) continue;
				if( !strncasecmp( trt, "CLA", 3 ) ) continue;

				int doit = 0;
				if( (tat[0] == 'C' && (tat[1] == '3'||tat[1]=='2') && tat[2] != '\0') || !strncasecmp( trt, "MYR", 3) || !strncasecmp(trt, "OLE",3) )
				{
					doit = 1;
				}

				if( strlen(tat) > 2 && (tat[0] == 'C' || tat[0] == 'D') && (tat[2] == 'A' || tat[2] == 'B')  )
				{
//					printf("Using residue %s to average.\n", trt );
					doit = 1;
				}
					
				if( doit )
				{	
					double tz = at[a].z;
					while( tz < 0 ) tz += 1.0;
					while( tz > 1.0 ) tz -= 1.0;

					int zb = N_BINS_MOLDIST * tz; // this is right
					if( zb < 0 ) zb = 0;
					if( zb >= N_BINS_MOLDIST ) zb = N_BINS_MOLDIST-1; 
					moldist[zb] += 1;
				}
			} 
			
			double best_cen = 0;
			double best_chi2 = 1e10;

			int nbins = N_BINS_MOLDIST;
			for( int zb = 0; zb < nbins; zb++ )
			{
				double zv = (zb+0.5) / (double)N_BINS_MOLDIST;
 
	                         int zlow  = zb- nbins/2;
	                         int zhigh = zlow + nbins;
	 
	                         double lchi2 = 0;
	                         for( int iz = zlow; iz < zhigh; iz++ )
	                         {
	                                 double dz = (iz+0.5) / nbins - zv;
	 
	                                 int iiz = iz;
	                                 while( iiz < 0 ) iiz += nbins;
	                                 while( iiz >= nbins ) iiz -= nbins;
	 
	                                 lchi2 += moldist[iiz] * (dz) * (dz);
	                         }
	
	                         if( lchi2 < best_chi2 )
	                         {
	                                 best_chi2 = lchi2;
	                                 best_cen = zv;
	                         }
			}

			if( best_cen > 0.5 ) best_cen -= 1.0;
			if( best_cen < -0.5 ) best_cen += 1.0;

			printf("best_cen: %le\n", best_cen );

			for( int a = 0; a < nat; a++ )
			{
				while( at[a].z - best_cen < -0.5 ) at[a].z += 1.0;
				while( at[a].z - best_cen >  0.5 ) at[a].z -= 1.0;

				at[a].z -= best_cen;
			}

			for( int a = 0; a < nat; a++ )
			{
				if( !strcasecmp( at[a].resname, "CLA" ) ) continue;

				if( martini_detected )
				{
					if( at[a].atname[0] == 'C' || at[a].atname[0] == 'D' )
					{
						carbon_com[0] += at[a].x;
						carbon_com[1] += at[a].y;
						carbon_com[2] += at[a].z;
					
						ncarbon += 1;	
					}		
				}
				else
				{
					if( at[a].atname[0] != 'C' ) continue;
	
					carbon_com[0] += at[a].x;
					carbon_com[1] += at[a].y;
					carbon_com[2] += at[a].z;
				
					ncarbon += 1;			
				}
			}
 
			carbon_com[0] /= ncarbon+1e-10;
			carbon_com[1] /= ncarbon+1e-10;
			carbon_com[2] /= ncarbon+1e-10;

			printf("carbon_com: %lf %lf %lf\n", carbon_com[0], carbon_com[1], carbon_com[2] );
			fflush(stdout);

                        memcpy( last_align, cur_set, sizeof(double) * 3 * nat );
                        aligned = 1;

                        for( int a = 0; a < curNAtoms(); a++ )
                        {
				at[a].x -= carbon_com[0];
				at[a].y -= carbon_com[1];
				at[a].z -= carbon_com[2];

				at[a].z += 0.5; // zero is at the middle.

                                while( at[a].x  < 0 ) at[a].x += 1;
                                while( at[a].y  < 0 ) at[a].y += 1;
                                while( at[a].z  < 0 ) at[a].z += 1;
                                while( at[a].x  > 1 ) at[a].x -= 1;
                                while( at[a].y  > 1 ) at[a].y -= 1;
                                while( at[a].z  > 1 ) at[a].z -= 1;
                        }

#ifndef DISABLE_MARTINI
			for( int a = 0; a < curNAtoms(); a++ )
			{
				if( !strcasecmp( at[a].resname, "W" ) )
					martini_detected = 1;
				if( !strcasecmp( at[a].atname, "NC3") )
					martini_detected = 1;				
				if( !strcasecmp( at[a].atname, "NH3") )
					martini_detected = 1;				
			}
#endif

			if( martini_detected && ! block->scatterFile )
			{
				natoms = natoms_martini;
				use_bi_length = martini_lengths;
			}
			if( !init_done )
			{
				int pres = at[0].res;
				int at_start = 0;
				for( int a = 0; a <= curNAtoms(); a++ )
				{
					if( a == curNAtoms() || at[a].res != pres )
					{
						if( 	
							strcasecmp( at[at_start].resname, "TIP3" )  &&
							strcasecmp( at[at_start].resname, "POT" )   &&
							strcasecmp( at[at_start].resname, "CLA" ) 
						  )
						{ 
							lipid_start[nlipids] = at_start;
							lipid_stop[nlipids] = a;

							double com[3]={0,0,0};
							for( int ax = at_start; ax < a; ax++ )
								com[2] += at[ax].z;
							com[2]/= (a-at_start);

							if( com[2] > 0.5 )
							{
								upper_leaflet_lipids[nupper] = nlipids;
								nupper++;
							}
							else
							{
								lower_leaflet_lipids[nlower] = nlipids;
								nlower++;
							}

							nlipids++;
						}
						at_start = a;
					}
					if( a < curNAtoms() )
						pres = at[a].res;
				} 

				nx = floor(3*La / cutoff);
				ny = floor(3*Lb / cutoff);
				nz = floor(3*Lc / cutoff);
		
				if( cutoff_mode ) // cylindrical
					nz = 1;
			

				if( nx < 1 ) nx = 1;
				if( ny < 1 ) ny = 1;
				if( nz < 1 ) nz = 1;

				if( nx >= 100 ) nx = 100;
				if( ny >= 100 ) ny = 100;
				if( nz >= 100 ) nz = 100;

//FIXME: Needs to be freed!
				bins = (binnit *)malloc( sizeof(binnit) * nx * ny * nz );
	
				for( int x = 0; x < nx*ny*nz; x++ )
				{
					bins[x].nposSpace = 10;
//FIXME: Needs to be freed!
					bins[x].id = (int*)malloc( sizeof(int) * bins[x].nposSpace );
					bins[x].npos = 0; 
				}
			
				pres = -1;
				char pseg[256];
				FILE *scatterFile=NULL;	
				if( block->scatterFile )
				{
					scatterFile = fopen(block->scatterFile, "r");
					if( !scatterFile )
					{
						printf("Couldn't open scattering length file '%s'.\n", block->scatterFile );
						exit(1);
					}
				}
					
				for( int a = 0; a < curNAtoms(); a++ )
				{
					b_res[a] = -1;
					b_ind[a] = -1;
					
					int is_prot = 0;
					if( !strcasecmp( at[a].resname, "ALA") ) is_prot=1;
					if( !strcasecmp( at[a].resname, "LEU") ) is_prot=1;
					if( !strcasecmp( at[a].resname, "VAL") ) is_prot=1;
					if( !strcasecmp( at[a].resname, "GLY") ) is_prot=1;
					if( !strcasecmp( at[a].resname, "TRP") ) is_prot=1;
					if( !strcasecmp( at[a].resname, "ILE") ) is_prot=1;
					if( !strcasecmp( at[a].resname, "PHE") ) is_prot=1;
					if( !strcasecmp( at[a].resname, "TYR") ) is_prot=1;
					if( !strcasecmp( at[a].resname, "THR") ) is_prot=1;
					if( !strcasecmp( at[a].resname, "SER") ) is_prot=1;
					if( !strcasecmp( at[a].resname, "PRO") ) is_prot=1;
					
					if( !strcasecmp( at[a].resname, "MET") ) is_prot=1;
					if( !strcasecmp( at[a].resname, "CYS") ) is_prot=1;
					
					if( !strcasecmp( at[a].resname, "ASN") ) is_prot=1;
					if( !strcasecmp( at[a].resname, "GLN") ) is_prot=1;
					if( !strcasecmp( at[a].resname, "HIS") ) is_prot=1;

					if( !strcasecmp( at[a].resname, "GLU") ) is_prot=1;
					if( !strcasecmp( at[a].resname, "ASP") ) is_prot=1;

					if( !strcasecmp( at[a].resname, "LYS") ) is_prot=1;
					if( !strcasecmp( at[a].resname, "ARG") ) is_prot=1;
					

					if( !strcasecmp( at[a].resname, "HIE") ) is_prot=1;
					if( !strcasecmp( at[a].resname, "HSE") ) is_prot=1;
					if( !strcasecmp( at[a].resname, "CHO") ) is_prot=1;

					if( !strcasecmp( at[a].resname, "DLE") ) is_prot=1;
					if( !strcasecmp( at[a].resname, "DVA") ) is_prot=1;
					if( !strcasecmp( at[a].resname, "EAM") ) is_prot=1;

					if( !strcasecmp( at[a].resname, "CLA" ) || !strcasecmp( at[a].resname, "POT") || !strcasecmp(at[a].resname, "SOD" ) )
					{
						pres = at[a].res;
						strcpy( pseg, at[a].segid );
						continue;
					}		
					for( int u = 0; u < nres_used; u++ )
					{
						if( is_prot && !strcasecmp( resNames[u], "PROT" ) )
							b_res[a] = u;	
						else if( !strcasecmp( at[a].resname, resNames[u] ) )
							b_res[a] = u; 
					}

					if( b_res[a] == -1 )
					{
						if( nres_used == MAX_RES )
						{
							printf("Ran out of space to store residue names.\n");
							exit(1);
						}
						if( is_prot )
						{
//FIXME: Needs to be freed!
							resNames[nres_used] = (char *)malloc( sizeof(char) * (1  +strlen("PROT") ) );
							strcpy( resNames[nres_used], "PROT" );
						}
						else
						{
//FIXME: Needs to be freed!
							resNames[nres_used] = (char *)malloc( sizeof(char) * (1  +strlen(at[a].resname) ) );
							strcpy( resNames[nres_used], at[a].resname );
						}
						b_res[a] = nres_used;
						nres_used++;
					}


					if( scatterFile )
					{
						getLine( scatterFile, buffer );
						int nr = sscanf( buffer, "%lf", b_i+a );
			
						if( nr == 0 )
						{
							printf("Couldn't read scattering length from line %d of '%s'.\n", a+1, block->scatterFile );
							exit(1);
						}

						b_ind[a] = -1;
						int gotit = 0;

						for( int xa = 0; xa < nlength; xa++ )
						{
							if( bi_length[xa] == b_i[a] )
							{
								b_ind[a] = xa;
								gotit = 1;
								break;
							}
						}
						if(!gotit)
						{
							bi_length[nlength] = b_i[a];
							b_ind[a] = nlength;
							nlength++;
						}
					}
					else
					{
						if( martini_detected )
						{
							for( int p = 0; p < natoms; p++ )
							{
								if( !strncasecmp( at[a].atname, martini_atoms[p], strlen(martini_atoms[p]) ) )
								{
									b_ind[a] = p;
									b_i[a] = martini_lengths[p];	
									break;
								}
							}
						}
						else
						{
							for( int p = 0; p < natoms; p++ )
							{
								if( at[a].atname[0] == bi_atoms[p] || (at[a].atname[0] == 'Z' && at[a].atname[1] == bi_atoms[p]) )
								{
									b_ind[a] = p;
#ifdef DO_SAXS
									b_i[a] = bi_length[p] - at[a].charge;
#else
									if( at[a].atname[0] == 'Z' )
									{
										double fr = (1.0 - at[a].deut_fraction) * bi_length[PROTON_B] + at[a].deut_fraction * bi_length[DEUTERON_B]; 
	
										if( fabs(bi_length[p]) > 1e-10 && fabs(bi_length[p]-fr) > 1e-10 )
										{
											printf("ERROR: two fractional deuteration schemes with Z%c?\n", bi_atoms[p] );
											exit(1);
										}
	
										bi_length[p] = fr;
										b_i[a] = fr;
									}
									else
										b_i[a] = bi_length[p];
#endif
								}
							}
						}
					}


					if( pres != at[a].res || strcasecmp( pseg, at[a].segid ) )
						n_res[b_res[a]] += 1;
		
					pres = at[a].res;
					strcpy( pseg, at[a].segid );

					if(  b_ind[a] == -1 )
					{
						printf("Scattering length not set for atom '%s'.\n", at[a].atname );
						exit(1);
					}
				}
			
				if( block->scatterFile )
				{		
					if( nlength > max_natoms )
					{
						printf("Increase max_natoms to at least %d in betaHist.C\n", nlength );
						exit(1);
					}
	
					natoms = nlength;
				}

				for( int a = 0; a < natoms; a++ )
					bi_length[a] *= lateral_scale * lateral_scale;
				for( int a = 0; a < nat; a++ )
					b_i[a] *= lateral_scale * lateral_scale;
				b_rho = (double ***)malloc( sizeof(double **) * nres_used );
				b_hist = (double ***)malloc( sizeof(double **) * nres_used );
				b_hist_cn = (double **)malloc( sizeof(double *) * nres_used );

				for( int r = 0; r < nres_used; r++ )
				{		
					b_rho[r] = (double **)malloc( sizeof(double *) * natoms );
					b_hist[r] = (double **)malloc( sizeof(double *) * natoms );
					b_hist_cn[r] = (double *)malloc( sizeof(double) * N_R_BINS );
					memset( b_hist_cn[r], 0, sizeof(double) * N_R_BINS );
			
					for( int a = 0; a < natoms; a++ )
					{
						b_rho[r][a]  = (double *)malloc( sizeof(double) * N_Z_BINS );
						memset( b_rho[r][a], 0, sizeof(double) * N_Z_BINS );
						b_hist[r][a] = (double *)malloc( sizeof(double) * N_R_BINS );
						memset( b_hist[r][a], 0, sizeof(double) * N_R_BINS );

						n_hist[r][a] = 0;
					}
				}
				init_done = 1;
			}

#if 0
			// this intramolecular correlation correction did not work. I meant it to subtract off the continuum molecule-molecule and add in the explicit molecule-molecule
			// I think it doesn't work because the result does not resemble a squared Fourier transform. Now I am going to try cycling through molecular positions.
	
			if( block->intra )
			{
				// for in-residue atoms outside the cutoff, add in intramolecular correlation contribution.

				int pres = at[0].res;
				int at_start = 0;
				for( int a = 0; a <= curNAtoms(); a++ )
				{
					if( a == curNAtoms() || at[a].res != pres )
					{
						for( int ax1 = at_start; ax1 < a; ax1++ )
						for( int ax2 = at_start; ax2 < a; ax2++ )
						{
							double dr[3] = { at[ax1].x - at[ax2].x,
									 at[ax1].y - at[ax2].y, 	
									 at[ax1].z - at[ax2].z };	

							while(dr[0] < -0.5 ) dr[0] += 1.0;
							while(dr[0] >  0.5 ) dr[0] -= 1.0;
							while(dr[1] < -0.5 ) dr[1] += 1.0;
							while(dr[1] >  0.5 ) dr[1] -= 1.0;
							while(dr[2] < -0.5 ) dr[2] += 1.0;
							while(dr[2] >  0.5 ) dr[2] -= 1.0;
				
							double p1[3] = {at[ax1].x,at[ax1].y,at[ax1].z};
							double p2[3] = {at[ax2].x,at[ax2].y,at[ax2].z};
							TransformFractional(dr);
							TransformFractional(p1);
							TransformFractional(p2);
							double rxy = sqrt(dr[0]*dr[0]+dr[1]*dr[1]);
							double r = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);

							if( disable_histogramming || (cutoff_mode && rxy > cutoff) || (!cutoff_mode && r > cutoff) )
							{	
								double rho_at_1 = evaluateSpline( p2[2] - Lc/2,0 ); 
								double rho_at_2 = evaluateSpline( p2[2] - Lc/2,0 ); 
								
								int rb;

								if( cutoff_mode )	
									rb = (r / max_r_cutoff) * N_R_BINS;
								else
								{
									printf("ERROR: intramolecular correlation doesn't currently work with the radial cutoff.\n");
									exit(1);
									rb = (r / cutoff) * N_R_BINS; // hm this isn't going to work I can see. we need max_r_cutoff here.
								}
	
								b_hist[b_res[ax1]][b_ind[ax1]][rb] += b_i[ax2];

								double length_at_1 = rho_at_1 * vol_scale * EffVol(at[ax1].atname[0]);
								double length_at_2 = rho_at_2 * vol_scale * EffVol(at[ax2].atname[0]);

								b_hist[b_res[ax1]][b_ind[ax1]][rb] += b_i[ax2];
								b_hist[b_res[ax1]][b_ind[ax1]][rb] -= length_at_2 * length_at_1 / b_i[ax1];
				//				printf("r: %lf add %le sub %le del %le z: %lf\n", r, b_i[ax2], rho_at_2 * vol_scale * EffVol(at[ax2].atname[0]),  b_i[ax2] - rho_at_2 * vol_scale * EffVol(at[ax2].atname[0]), p2[2]-Lc/2 ); 
							} 
						}
						at_start = a;
					}
					if( a < curNAtoms() )
						pres = at[a].res;
				} 

			}
#endif
			if( block->intra )
			{
				double lipid_com[3*nlipids];
				memset( lipid_com, 0, sizeof(double) * 3 * nlipids );

				for( int l = 0; l < nlipids; l++ )
				{
					lipid_com[3*l+0] = 0;
					lipid_com[3*l+1] = 0;
					lipid_com[3*l+2] = 0;

					for( int a = lipid_start[l]; a < lipid_stop[l]; a++ )
					{
						double dr[3] = { at[a].x - lipid_com[3*l+0], at[a].y - lipid_com[3*l+1], at[a].z - lipid_com[3*l+2] };

						if( a != lipid_start[l] )
						{
							// wrap, add, maintain com in place.
							while( dr[0] < -0.5 ) dr[0] += 1.0;
							while( dr[0] > 0.5 ) dr[0] -= 1.0;
							while( dr[1] < -0.5 ) dr[1] += 1.0;
							while( dr[1] > 0.5 ) dr[1] -= 1.0;
							while( dr[2] < -0.5 ) dr[2] += 1.0;
							while( dr[2] > 0.5 ) dr[2] -= 1.0;
						}
						dr[0] += lipid_com[3*l+0];
						dr[1] += lipid_com[3*l+1];
						dr[2] += lipid_com[3*l+2];
	
						lipid_com[3*l+0] *= (a-lipid_start[l]);
						lipid_com[3*l+1] *= (a-lipid_start[l]);
						lipid_com[3*l+2] *= (a-lipid_start[l]);

						lipid_com[3*l+0] += dr[0];
						lipid_com[3*l+1] += dr[1];
						lipid_com[3*l+2] += dr[2];

						lipid_com[3*l+0] /= (a-lipid_start[l]+1);
						lipid_com[3*l+1] /= (a-lipid_start[l]+1);
						lipid_com[3*l+2] /= (a-lipid_start[l]+1);
					}
				}

				// interchange COMs.

				int upper_sorter[nupper];
				int lower_sorter[nlower];

				for( int u = 0; u < nupper; u++ )
					upper_sorter[u] = u;
				for( int t = 0; t < nupper; t++ )
				{
					int get = rand() % (nupper-t);

					int old = upper_sorter[t];
					upper_sorter[t] = upper_sorter[t+get]; 
					upper_sorter[t+get] = old;
				}
				
				for( int u = 0; u < nlower; u++ )
					lower_sorter[u] = u;
				for( int t = 0; t < nlower; t++ )
				{
					int get = rand() % (nlower-t);

					int old = lower_sorter[t];
					lower_sorter[t] = lower_sorter[t+get]; 
					lower_sorter[t+get] = old;
				}

				for( int u = 0; u < nupper; u++ )
				{
					int l = upper_leaflet_lipids[u];
					int l_x = upper_leaflet_lipids[upper_sorter[u]];

					for( int a = lipid_start[l]; a < lipid_stop[l]; a++ )
					{
						at[a].x += lipid_com[3*l_x+0] - lipid_com[3*l+0];
						at[a].y += lipid_com[3*l_x+1] - lipid_com[3*l+1];
//						at[a].z += lipid_com[3*l_x+2] - lipid_com[3*l+2];
						while(at[a].x < 0 ) at[a].x += 1.0;
						while(at[a].x > 1 ) at[a].x -= 1.0;
						while(at[a].y < 0 ) at[a].y += 1.0;
						while(at[a].y > 1 ) at[a].y -= 1.0;
						while(at[a].z < 0 ) at[a].z += 1.0;
						while(at[a].z > 1 ) at[a].z -= 1.0;
					}
				}
				
				for( int u = 0; u < nlower; u++ )
				{
					int l = lower_leaflet_lipids[u];
					int l_x = lower_leaflet_lipids[lower_sorter[u]];

					for( int a = lipid_start[l]; a < lipid_stop[l]; a++ )
					{
						at[a].x += lipid_com[3*l_x+0] - lipid_com[3*l+0];
						at[a].y += lipid_com[3*l_x+1] - lipid_com[3*l+1];
//						at[a].z += lipid_com[3*l_x+2] - lipid_com[3*l+2];
	
						while(at[a].x < 0 ) at[a].x += 1.0;
						while(at[a].x > 1 ) at[a].x -= 1.0;
						while(at[a].y < 0 ) at[a].y += 1.0;
						while(at[a].y > 1 ) at[a].y -= 1.0;
						while(at[a].z < 0 ) at[a].z += 1.0;
						while(at[a].z > 1 ) at[a].z -= 1.0;
					}
				}
			}			


			for( int b = 0; b < nx*ny*nz; b++ )
				bins[b].npos = 0;
	
			for( int a = 0; a < curNAtoms(); a++ )
			{

				int bx = nx * (at[a].x );
	 			while( bx < 0 ) bx += nx;
	 			while( bx >= nx ) bx -= nx;
				int by = ny * (at[a].y );
	 			while( by < 0 ) by += ny;
	 			while( by >= ny ) by -= ny;
				int bz = nz * (at[a].z );
	 			while( bz < 0 ) bz += nz;
	 			while( bz >= nz ) bz -= nz;
	
				int tbin = bz+nz*(by+ny*bx);
	
				if( bins[tbin].nposSpace == bins[tbin].npos )
				{
					bins[tbin].nposSpace *= 2;
					bins[tbin].id = (int *)realloc( bins[tbin].id, sizeof(int) * bins[tbin].nposSpace );
				}
	
				bins[tbin].id[bins[tbin].npos] = a;
	
				bins[tbin].npos += 1;
			}
			
	
			double rsearch = 5.0;
	
			int nbox_x = 1;
			int nbox_y = 1;
			int nbox_z = 1;

			double bw_x = La/nx;
			double bw_y = Lb/ny;
			double bw_z = Lc/nz;
	
			int s_mdx = -ceil(cutoff/bw_x);
			int s_pdx =  ceil(cutoff/bw_x);
			int s_mdy = -ceil(cutoff/bw_y);
			int s_pdy =  ceil(cutoff/bw_y);
			int s_mdz = -ceil(cutoff/bw_z);
			int s_pdz =  ceil(cutoff/bw_z);

			if( cutoff_mode ) 
			{ // cylindrical
				s_mdz = -1;
				s_pdz = 1;
			}
/*	
			if( nx == 1 ) { s_mdx = 0; s_pdx = 0; }
			else if( nx == 2 ) { s_mdx = -1; s_pdx = 0; }
			if( ny == 1 ) { s_mdy = 0; s_pdy = 0; }
			else if( ny == 2 ) { s_mdy = -1; s_pdy = 0; }
			if( nz == 1 ) { s_mdz = 0; s_pdz = 0; }
			else if( nz == 2 ) { s_mdz = -1; s_pdz = 0; }
*/
			if( mode == 0 && all_f_for_skip % gofr_skip == 0 ) // create histogram for mode 0 
			{
				for( int a = 0; a < nat; a++ )
				{
					if(b_ind[a] == -1 ) continue;

					int zbin = at[a].z * N_Z_BINS;

					if( zbin < pick_zmin || zbin >= pick_zmax ) continue;

					n_hist[b_res[a]][b_ind[a]] += 1;
				}
		
				if( !disable_histogramming )
				{

				// constructing beta_hist_cn, the continuum background interacting with the all-atom system.
				for( int a = 0; a < nat; a++ )
				{
					if( b_ind[a] == -1 ) continue;
					
					// the particle z coordinate.
					double pz = at[a].z*Lc;
					// min should be zero, max should be Lz. this is the solvent background inside the simulation cell, the simulation cell sum. it is interacting with the infinite exterior all-atom sum.
					double z_bg_max = Lc * pick_zmax / N_Z_BINS;
					double z_bg_min = Lc * pick_zmin / N_Z_BINS;

/*
 *   Attention! calculation of term III is here!
 *
 *
 * */

					if( cutoff_mode )
					{
						for( int ir = 0; ir < N_R_BINS; ir++ )
						{
							double r1 = ir * max_r_cutoff / N_R_BINS;	
							double r2 = (ir+1)* max_r_cutoff / N_R_BINS;	
	
							double vol = util_volume( z_bg_min-pz, z_bg_max-pz, r1, r2, 0, cutoff );

							// will be multiplied by an arbitrary background later.
							b_hist_cn[b_res[a]][ir] += vol * b_i[a];						
						}
					}
					else
					{
						for( int ir = 0; ir < N_R_BINS; ir++ )
						{
							double r1 = ir * cutoff / N_R_BINS;	
							double r2 = (ir+1)* cutoff / N_R_BINS;	
	
							double vol = util_volume( z_bg_min-pz, z_bg_max-pz, r1, r2 );

							// will be multiplied by an arbitrary background later.
							b_hist_cn[b_res[a]][ir] += vol * b_i[a];						
						}
					}				
				} 

				for( int bin_x = 0; bin_x < nx; bin_x++ )
				for( int bin_y = 0; bin_y < ny; bin_y++ )
				for( int bin_z = 0; bin_z < nz; bin_z++ )
				{
					int bin = bin_x*ny*nz + bin_y*nz + bin_z;
		 
					printf("bin: %d/%d %d/%d %d/%d\n", bin_x, nx, bin_y, ny, bin_z, nz );
					fflush(stdout);
					for( int a1x = 0; a1x < bins[bin].npos; a1x++ )
					{
						int a = bins[bin].id[a1x];

						int zbin = at[a].z * N_Z_BINS;
						if( zbin < pick_zmin || zbin >= pick_zmax ) continue;
						if( b_ind[a] == -1) continue;

						pick_beta += b_i[a] / ((Lc/N_Z_BINS)*(pick_zmax-pick_zmin)*La*Lb);
	
						double debug_br[N_R_BINS];
						memset( debug_br, 0, sizeof(double) * N_R_BINS );
						for( int dbx = s_mdx; dbx <= s_pdx; dbx++ )			
						for( int dby = s_mdy; dby <= s_pdy; dby++ )			
						for( int dbz = s_mdz; dbz <= s_pdz; dbz++ )		
						{
							int tbin_x = bin_x + dbx;
							int tbin_y = bin_y + dby;
							int tbin_z = bin_z + dbz;
//#define USE_Z_PBC				
							double del[3] = { 0, 0, 0 };

							while( tbin_x < 0 ) { tbin_x += nx; del[0] -= 1.0; }
							while( tbin_y < 0 ) { tbin_y += ny; del[1] -= 1.0; }
#ifdef USE_Z_PBC
							while( tbin_z < 0 ) { tbin_z += ny; del[2] -= 1.0; }
#else
							if( tbin_z < 0 ) continue; //tbin_z += nz;
#endif					
							while( tbin_x >= nx ) { tbin_x -= nx; del[0] += 1.0; }
							while( tbin_y >= ny ) { tbin_y -= ny; del[1] += 1.0; }
#ifdef USE_Z_PBC
							while( tbin_z >= nz ) { tbin_z -= ny; del[2] += 1.0; }
#else
							if( tbin_z >= nz ) continue; // tbin_z -= nz;
#endif
		
							int near_bin = tbin_x * ny * nz + tbin_y * nz + tbin_z;
				
							for( int a2x = 0; a2x < bins[near_bin].npos; a2x++ )
							{
								int a2 = bins[near_bin].id[a2x];
								if( b_ind[a2] == -1 ) continue;
	
								double dr[3] = { at[a2].x+del[0] - at[a].x, at[a2].y+del[1] - at[a].y, at[a2].z + del[2] - at[a].z };

								
				
//								while( dr[0] < -0.5 ) dr[0] += 1.0;
//								while( dr[0] > 0.5 ) dr[0] -= 1.0;
//								while( dr[1] < -0.5 ) dr[1] += 1.0;
//								while( dr[1] > 0.5 ) dr[1] -= 1.0;

								// DO NOT go outside PBC in the z direction, use the continuum approximation. 

//								if( dr[2] < -0.5 || dr[2] > 0.5 )
//									continue;

//								while( dr[2] < -0.5 ) dr[2] += 1.0;
//								while( dr[2] > 0.5 ) dr[2] -= 1.0;
				
								double op[3] = { at[a].x, at[a].y, at[a].z };
//								double odr[3] = { dr[0], dr[1], dr[2] };
								TransformFractional( dr );
//								TransformFractional( op );
				
								double r = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
								int rb;								
								if( cutoff_mode )
								{
									double rc = sqrt(dr[0]*dr[0]+dr[1]*dr[1]);
									if( rc > cutoff ) 
										continue;

									if( r > 225 )
									{
										printf("okay huh");
									}

									rb = (r / max_r_cutoff) * N_R_BINS;
								}
								else
								{
									rb = (r / cutoff) * N_R_BINS;
								}
								if( rb >= N_R_BINS ) continue;

								if( rb == 0 && a != a2 && !warning_trigger )
								{
									printf("WARNING: ATOMS WERE FOUND IN THE SAME ZERO BIN. THEY ARE TOO CLOSE TOGETHER.\n");
									warning_trigger = 1;
	//								exit(1);
								}

#ifndef DISABLE_INTRA
								if( at[a2].atname[0] == 'Z' && at[a].atname[0] == 'Z' && at[a].res == at[a2].res && !strcasecmp( at[a].segid, at[a2].segid) )
								{
									// Special intra-molecular deuteration fraction hack.

									if( fabs(b_i[a]) < 1e-10 )
									{
										printf("Divide by zero/lack of precision error here.\n");
										exit(1);
									}

									if( !strcasecmp( at[a].atname, at[a2].atname) )
									{
										b_hist[b_res[a]][b_ind[a]][rb] += at[a].deut_fraction * bi_length[DEUTERON_B] * bi_length[DEUTERON_B] / b_i[a]; // here is the divide by zero.
										b_hist[b_res[a]][b_ind[a]][rb] += (1-at[a].deut_fraction) * bi_length[PROTON_B] * bi_length[PROTON_B] / b_i[a];
									}
									else
									{ 
										// when we are a deuteron it's definitely hydrogen
										double f_run = at[a].deut_fraction;
										b_hist[b_res[a]][b_ind[a]][rb] += at[a].deut_fraction * bi_length[PROTON_B] * bi_length[DEUTERON_B] / b_i[a];
										// when it is a deuteron we are definitely hydrogen
										f_run += at[a2].deut_fraction;
										b_hist[b_res[a]][b_ind[a]][rb] += at[a2].deut_fraction * bi_length[PROTON_B] * bi_length[DEUTERON_B] / b_i[a];	
										// otherwise we are both hydrogen
										b_hist[b_res[a]][b_ind[a]][rb] += (1-f_run) * bi_length[PROTON_B] * bi_length[PROTON_B] / b_i[a];	
									}
								}
								else
#endif
								{
					/* Attention! This is where term Ia is calculated: */
									b_hist[b_res[a]][b_ind[a]][rb] += b_i[a2];
								}
							}
						}

						double dz_bot = at[a].z * Lc;
						double dz_top = Lc - at[a].z * Lc;

#ifndef USE_Z_PBC
						if( dz_bot < cutoff || (cutoff_mode && dz_bot < max_r_cutoff ) )
						{
							for( int ir = 0; ir < N_R_BINS; ir++ )
							{
								double r1,r2;

								if( cutoff_mode )
								{
									r1 = ir * max_r_cutoff / (double)N_R_BINS;
									r2 = (ir+1) * max_r_cutoff / (double)N_R_BINS;
								}				
								else
								{
									r1 = ir * cutoff / (double)N_R_BINS;
									r2 = (ir+1) * cutoff / (double)N_R_BINS;
								}				
								// the volume that is within the histogram as well as the cylinder.
								double vol;

								if( cutoff_mode )
									vol = util_volume( -1e10, -dz_bot,  r1, r2, 0, cutoff );
								else
									vol = util_volume( -1e10, -dz_bot,  r1, r2 );
								
								/* Technically this is part of term Ia, it is within the PFFT radius */
								b_hist[b_res[a]][b_ind[a]][ir] += vol * (outside_z_box_rho);
							}				
						}
					
						if( dz_top < cutoff || (cutoff_mode && dz_top < max_r_cutoff ) )
						{
							for( int ir = 0; ir < N_R_BINS; ir++ )
							{
								double r1,r2;

								if( cutoff_mode )
								{
									r1 = ir * max_r_cutoff / (double)N_R_BINS;
									r2 = (ir+1) * max_r_cutoff / (double)N_R_BINS;
								}				
								else
								{
									r1 = ir * cutoff / (double)N_R_BINS;
									r2 = (ir+1) * cutoff / (double)N_R_BINS;
								}				
								// the volume that is within the histogram as well as the cylinder.
								double vol;

								if( cutoff_mode )
									vol = util_volume( dz_top, 1e10, r1, r2, 0, cutoff );
								else
									vol = util_volume( dz_top, 1e10, r1, r2 );
							
								/* Technically this is part of term Ia, it is within the PFFT radius */
								b_hist[b_res[a]][b_ind[a]][ir] += vol * (outside_z_box_rho);
							}				
						}
#endif

					}
				}
				}	

				nframes_gofr += 1;
			}	

			// always create beta(z)
			
			for( int a = 0; a < curNAtoms(); a++ )
			{
				if( b_ind[a] == -1) continue;

				int bz = at[a].z * N_Z_BINS;
				while( bz < 0 ) bz += N_Z_BINS;
				while( bz >= N_Z_BINS ) bz -= N_Z_BINS;
				
				av_rho_get[bz] += b_i[a];
				b_rho[b_res[a]][b_ind[a]][bz] += 1;
			}

			for( int a = 0; a < curNAtoms(); a++ )
				at[a].zap();

			nframe_average += 1;
		} 

		fclose(dcdFile);
	}


	av_Lx /= nav_Lx;
	av_Ly /= nav_Ly;
	av_Lz /= nav_Lz;
		

	if( mode == 1 )
	{
/*		double nwat_average = 0;
		double water_rho_average = 0;
		double len_av = 5;
		for( int z = 0 ;z < N_Z_BINS; z ++ )
		{
			double val = 0;
			double vol = av_Lx * av_Ly * (av_Lz / N_Z_BINS );
			double zbin = -av_Lz/2 + (z+0.5) * (av_Lz/N_Z_BINS);
			for( int xx = 0; xx < natoms; xx++ )
				val += b_rho[xx][z] * bi_length[xx] / nframe_average / vol;
			
			if( fabs(fabs(zbin)-av_Lz/2) < len_av )
			{
				water_rho_average += val;
				nwat_average += 1;
			}
		}
	
		water_rho_average /= nwat_average;
*/
		char fileName[256];
		sprintf(fileName, "%s.bz", uniqueName );
		FILE *bzFile = fopen(fileName, "w");
		fprintf(bzFile, "%d %lf %lf %lf\n", N_Z_BINS, av_Lx, av_Ly, av_Lz );

		for( int ir = 0; ir < nres_used; ir++ )
			fprintf(bzFile, "%d ", n_res[ir] );
		fprintf(bzFile, "\n");

		fprintf(bzFile, "z tot");
		for( int ir = 0; ir < nres_used; ir++ )
			fprintf(bzFile, " %s", resNames[ir] );
		fprintf(bzFile, "\n");
	

		for( int iz = 0; iz < N_Z_BINS; iz++ )
		{
			double vol = av_Lx * av_Ly * (av_Lz / N_Z_BINS );
			double zbin = -av_Lz/2 + (iz+0.5) * (av_Lz/N_Z_BINS);
			double val = 0;

			for( int r = 0; r < nres_used; r++ )
			for( int xx = 0; xx < natoms; xx++ )
				val += b_rho[r][xx][iz] * use_bi_length[xx] / nframe_average / vol;
			fprintf(bzFile, "%lf %lf", zbin, av_rho_get[iz] / nframe_average / vol );

			for( int r = 0; r < nres_used; r++ )
			{
				double valr = 0;
				for( int xx = 0; xx < natoms; xx++ )
					valr += b_rho[r][xx][iz] * use_bi_length[xx] / nframe_average / vol;
				fprintf(bzFile, " %lf", valr );
			}
			fprintf(bzFile, "\n");

			av_rho[iz] = val;
		}

#ifdef FFTW
		fftw_complex *h_in;
		fftw_complex *h_out;
		fftw_plan p;

		int inflation = 10;
		int inflated_z = N_Z_BINS * inflation;

//FIXME: Needs to be freed!
		h_in = (fftw_complex *)fftw_malloc( sizeof(fftw_complex) * inflated_z );
//FIXME: Needs to be freed!
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
		
		FILE *fzFile = fopen("ft.z", "w");
		for( int j = 0; j < N_Z_BINS; j++ )
		{
			double q = j * 2 * M_PI / (inflation*av_Lz);
			fprintf(fzFile, "q %lf %d %le\n", q, j, h_out[j][0]*h_out[j][0] + h_out[j][1]*h_out[j][1] );
		}
#endif

	}
	else 
	{
		// write our histogram.
	
		char fileName[256];
		sprintf(fileName, "%s.br", uniqueName );
		FILE *outBR = fopen( fileName, "w" );
	
		if( cutoff_mode )
			fprintf(outBR, "%d %lf %lf %lf %lf cylinder %lf\n", N_R_BINS, max_r_cutoff, av_Lx, av_Ly, av_Lz, cutoff );	
		else
			fprintf(outBR, "%d %lf %lf %lf %lf\n", N_R_BINS, cutoff, av_Lx, av_Ly, av_Lz );	
		fprintf(outBR, "nhist" );
		for( int x = 0; x < natoms; x++ )
		{
			double num = 0;
			for( int r = 0; r < nres_used; r++ )
				num += n_hist[r][x]/nframes_gofr;

			if( num > 0 )
				fprintf(outBR, " %lf", num );
		}	 
		fprintf(outBR, "\n");
		fprintf(outBR, "blength" );
		for( int x = 0; x < natoms; x++ )
		{
			double num = 0;
			for( int r = 0; r < nres_used; r++ )
				num += n_hist[r][x];

			if( num > 0 )
				fprintf(outBR, " %lf", use_bi_length[x] ); 
		}
		fprintf(outBR, "\n");

			
		fprintf(outBR, "cn");
		for( int ir = 0; ir < N_R_BINS; ir++ )
		{
			double val = 0;
			for( int r = 0; r < nres_used; r++ )
				val += b_hist_cn[r][ir] / nframes_gofr;
			fprintf(outBR, " %lf", val );
		}
		fprintf(outBR, "\n");
		
	
		for( int x = 0; x < natoms; x++ )
		{
			double num = 0;
			for( int r = 0; r < nres_used; r++ )
				num += n_hist[r][x];
			if( num > 0 )
			{
				if( block->scatterFile )
					fprintf(outBR, "%d", x );
				else if( martini_detected )
					fprintf(outBR, "%s", martini_atoms[x] );
				else
					fprintf(outBR, "%c", bi_atoms[x] );
				for( int ir = 0; ir < N_R_BINS; ir++ )
				{
					double val = 0;
					for( int r = 0; r < nres_used; r++ )
						val += b_hist[r][x][ir]/nframes_gofr;
					if( ! ((val <0 ) || (val > -1)) )
					{	
						printf("nan.\n");
					}	
					fprintf(outBR, " %lf", val);
				}
				fprintf(outBR, "\n"); 
			}
		}
		
		fprintf(outBR, "IV");
		for( int ir = 0; ir < N_R_BINS; ir++ )
		{
			fprintf(outBR, " %lf", termIV[ir]/nframes_gofr );
		}
		fprintf(outBR, "\n");
	/*	
		for( int r = 0; r < nres_used; r++ )
		{
			for( int x = 0; x < natoms; x++ )
			{
				double num = 0;
				for( int ir = 0; ir < nres_used; ir++ )
					num += n_hist[ir][x];
				if( num > 0 )
				{
					fprintf(outBR, "%c %s", bi_atoms[x], resNames[r] );
					for( int ir = 0; ir < N_R_BINS; ir++ )
						fprintf(outBR, " %lf", b_hist[r][x][ir]/nframes_gofr );
					fprintf(outBR, "\n"); 
				}
			}
		}

		for( int r = 0; r < nres_used; r++ )
		{
			fprintf(outBR, "cn %s", resNames[r] );
			for( int ir = 0; ir < N_R_BINS; ir++ )
				fprintf(outBR, " %lf", b_hist_cn[r][ir] / nframes_gofr );
			fprintf(outBR, "\n");
		}
*/
		fclose(outBR);

#define SKIP_BH_REPORT
#ifndef SKIP_BH_REPORT
		// read in beta_z
		//
		double real_sub = water_rho_average;
		
		double av_picked_cn = 0;
		for( int z = pick_zmin; z < pick_zmax; z++ )
			av_picked_cn += av_rho[z];
		av_picked_cn /= (pick_zmax-pick_zmin);
#ifndef DONT_SUB
		for( int z = 0; z < N_Z_BINS; z++ )
			av_rho[z] -= water_rho_average;
#else
		real_sub = 0; // DONT_SUB
#endif
    //FIXME: Needs to be freed!
		double *c_beta_tot = (double *)malloc( sizeof(double) * N_MAX_R_BINS );
		memset( c_beta_tot, 0, sizeof(double) * N_MAX_R_BINS );
		
		double sum_beta = 0;
		for( int iz = 0; iz < N_Z_BINS; iz++ )
			sum_beta += av_rho[iz];
		sum_beta /= N_Z_BINS;
		sum_beta *= av_Lz;
		double slope_sub = (av_Lx*av_Ly) * M_PI * sum_beta * sum_beta;
		



			
		for( int ir = 0; ir < N_MAX_R_BINS; ir++ )
		{
			double r1 = (ir) * max_r_cutoff / N_MAX_R_BINS;
			double r2 = (ir+1) * max_r_cutoff / N_MAX_R_BINS;
			double dr = max_r_cutoff / N_MAX_R_BINS;

			double continuum_value = 0, real_value = 0;

			if( ir < N_R_BINS )
			{
				double vol = (4.0/3.0) * M_PI * (r2*r2*r2-r1*r1*r1);


				double pick_vol = av_Lx *av_Ly * av_Lz * (pick_zmax-pick_zmin)/(double)(N_Z_BINS);
//				double bg = av_Lx * av_Ly * av_Lz * real_sub;
				double bg = pick_vol * real_sub;


				double rho_i_total_bg = 0;

				for( int xx = 0; xx < natoms; xx++ )
				{
					if( n_hist[xx] < 1 ) continue;

					// the histogram of the all-atom <-> all-atom interactions
					real_value += use_bi_length[xx] * b_hist[xx][ir] / nframes_gofr;
					rho_i_total_bg += use_bi_length[xx] * (n_hist[xx]/nframes_gofr) * vol * real_sub;
// DO NOT USE:
	
//					real_value += bi_length[xx] * (n_hist[xx]/nframes_gofr) * (b_hist[xx][ir] / n_hist[xx] - vol * real_sub); // beta rho_0 
// OLD					real_value += bi_length[xx] * (n_hist[xx]/nframes_gofr) * (b_hist[xx][ir] / n_hist[xx] - 2 * vol * real_sub ); 
//					alt_cn[ir] += bi_length[xx] * (n_hist[xx]/nframes_gofr) * (-2*vol*real_sub);
//					real_value -= bg * (b_hist[xx][ir] / n_hist[xx]); 
				}
				
				// subtract off the all-atom interacting with the background.
				real_value -= b_hist_cn[ir] * real_sub;
				alt_cn[ir] -= b_hist_cn[ir] * real_sub;
	
				real_value -= rho_i_total_bg;
				alt_cn[ir] -= rho_i_total_bg;

				double layer_vol = av_Lx * av_Ly * av_Lz / N_Z_BINS;

				double bg_bg = 0;
				double corr_cross = 0;
				for( int iz = pick_zmin; iz < pick_zmax; iz++ )
				{
					double zv = iz * av_Lz / N_Z_BINS;
	
					double bg_vol = util_volume( -zv, av_Lz-zv, r1, r2 ); 
					
					corr_cross += real_sub * layer_vol * (vol-bg_vol) * real_sub;

					bg_bg += real_sub * real_sub * layer_vol * vol;
				}

				real_value -= corr_cross;
				alt_cn[ir] -= corr_cross;
		
				real_value += bg_bg;
				alt_cn[ir] += bg_bg;	

//				alt_cn[ir] += bg * vol * real_sub;
//				real_value += bg * vol * real_sub; 
			} 
			
	
			for( int z = 0; z < N_Z_BINS; z++ )
			{
				//if( z != pick_zbin ) continue;
				if( z < pick_zmin || z >= pick_zmax ) continue;

				double zv = -av_Lz/2 + (z+0.5) * av_Lz / N_Z_BINS;
				double zvol = av_Lx * av_Ly * (av_Lz/N_Z_BINS);
	
				double bv = av_rho[z];
				for( int z2 = -2*N_Z_BINS; z2 < N_Z_BINS+2*N_Z_BINS; z2++ )
				{

#ifndef DONT_SUB
					if( z2 < 0 || z2 >= N_Z_BINS ) continue;
#endif
					double bin_z1 = -av_Lz/2 + z2 * av_Lz / N_Z_BINS - zv;	
					double bin_z2 = -av_Lz/2 + (z2+1) * av_Lz / N_Z_BINS - zv;	

					double approx_vol = util_volume( bin_z1, bin_z2, r1, r2 );

					int use_z2 = z2;
#ifdef DONT_SUB
					while( use_z2 < 0 ) use_z2 += N_Z_BINS;
					while( use_z2 >= N_Z_BINS ) use_z2 -= N_Z_BINS;
#endif
					if( z2 < 0 || z2 >= N_Z_BINS )
						continuum_value += bv * outside_z_box_rho * approx_vol * zvol;	
					else			
						continuum_value += bv * av_rho[use_z2] * approx_vol * zvol;	
				}

			}
			double rmean = (r1+r2)/2;	
			double alpha = smoother( rmean, cutoff-smooth_width, cutoff );	
			double del = slope_sub * rmean;

			printf("real: %lf continuum: %lf asymptote: %lf alt_con %lf\n",
				real_value, continuum_value, del, (ir < N_R_BINS ? alt_cn[ir] : 0.0) );
			c_beta_tot[ir] = alpha * real_value + (1-alpha) * continuum_value - del;
		}		
		printf("Pick beta: %lf av_picked_cn: %lf\n", pick_beta, av_picked_cn );
#endif // SKIP_BH_REPORT
	}
  
	free(block);
}

double EffVol( char atname )
{
	double rad = 1;
	switch( atname )
	{
		case 'H':
		case 'D':
		case 'X':
			rad = 1.1;
			break;
		case 'C':
			rad = 1.7;
			break;
		case 'O':
			rad = 1.52;
			break;
		case 'N':
			rad = 1.55;
			break;
		case 'P':
			rad = 1.8;
			break;
//		default:
//			printf("Vol for atom '%c' not defined.\n", atname );
	}

	return (4.0/3.0)*M_PI*rad*rad*rad;
}

