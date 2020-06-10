#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "input.h"
#include "spline.h"
#include "util.h"
#include "beta_util.h"


double Power( double a, double b )
{
	return pow(a,b);
}

// builds beta(z) used in continuum region.
//
void buildData( struct parameter_block *blocks, double *betaz, int N_Z_BINS, double bz_width, double max_r_cutoff, int N_MAX_R_BINS, double *water_rho)
{
	char *buffer = (char *)malloc( sizeof(char) * 100000 );
	memset( betaz, 0, sizeof(double) * N_Z_BINS );

	double water_rho_average = 0;
	double nwater_rho_average = 1e-15;

	double MIN_Z = -(N_Z_BINS*bz_width)/2;
	double MAX_Z =  (N_Z_BINS*bz_width)/2;
	
	for( parameter_block *block = blocks; block; block = block->next )
	{

		FILE *bzFile = fopen(block->bzFile, "r");
	
		if( !bzFile )
		{
			printf("Requested beta.z file \"%s\" not found.\n", block->bzFile );
			exit(1);
		}

		getLine( bzFile, buffer );
		int l_N_Z_BINS;
		double t_av_Lx, t_av_Ly;
		double l_av_Lx, l_av_Ly, l_av_Lz;
		char code;
		
		int nr = sscanf(buffer, "%d %lf %lf %lf\n", &l_N_Z_BINS, &t_av_Lx, &t_av_Ly, &l_av_Lz );
	
		block->av_bz_area   += t_av_Lx * t_av_Ly;
		block->av_bz_vol = t_av_Lx * t_av_Ly * l_av_Lz;
		
		block->av_area = block->av_bz_area;
		block->height = l_av_Lz;


		getLine( bzFile, buffer );
		getLine( bzFile, buffer );


		for( int z = 0; z < l_N_Z_BINS; z ++ )
		{
			double val = 0;
			double zbin;
	
	
			getLine( bzFile, buffer );
			double zv, rv;
			sscanf( buffer, "%lf %lf", &zv, &rv ); 

//			printf("zv: %lf rv: %lf lavz: %lf wat: %lf \n", zv, rv, l_av_Lz, block->water_average_thickness  );	
			if( (zv < -l_av_Lz/2 + block->water_average_thickness ) ||
			    (zv >  l_av_Lz/2 - block->water_average_thickness ) )
			{
				water_rho_average += rv;
				nwater_rho_average += 1;
			}	
		}
		
		fclose(bzFile);
	}

	water_rho_average /= nwater_rho_average;
	int forcing = 0;	
	for( parameter_block *block = blocks; block; block = block->next )
	{
		if( block->forceSolventSub )
		{
			water_rho_average = block->solvent_sub;
			forcing = 1;
		}
	}

	if( forcing )
	{
		printf("Forcing solvent subtraction to %lf.\n", water_rho_average );
	}

//	printf("Water sub: %lf\n", water_rho_average );

	for( parameter_block *block = blocks; block; block = block->next )
	{
		double smooth_width = block->smooth;

		FILE *bzFile = fopen(block->bzFile, "r");
	
		if( !bzFile )
		{
			printf("Requested beta.z file \"%s\" not found.\n", block->bzFile );
			exit(1);
		}

		getLine( bzFile, buffer );
		int l_N_Z_BINS;
		double t_av_Lx, t_av_Ly;
		double l_av_Lx, l_av_Ly, l_av_Lz;
		char code;
		
		int nr = sscanf(buffer, "%d %lf %lf %lf\n", &l_N_Z_BINS, &t_av_Lx, &t_av_Ly, &l_av_Lz );
		l_av_Lx = t_av_Lx;	
		l_av_Ly = t_av_Ly;	
	
		int nres_used = 0;
		int read = 10;
		int *res_counts = (int *)malloc( sizeof(int) * read );
		int done_read = 0;
	
		getLine( bzFile, buffer );
	
		while( !done_read )
		{
			nres_used = readNInts( buffer, res_counts, read );
		
			if( nres_used == read ) 
			{
				read *= 2; 
				res_counts = (int*)realloc( res_counts, sizeof(int) * read );
			} 
			else
				done_read = 1;
		}
	
		char *resNames[nres_used];
		for( int ir = 0; ir < nres_used; ir++ )
			resNames[ir] = (char *)malloc( sizeof(char) * (strlen(buffer )+1) );
		double av_rho[l_N_Z_BINS];
		double **av_rho_res;
	
		getLine( bzFile, buffer );
		char *read_buffer = buffer + strlen("z tot ");

		int names_read = 0;
		int clen = 0;

		while( *read_buffer && names_read < nres_used )
		{	
			while( *read_buffer && *read_buffer != ' ' && *read_buffer != '\t' )
			{
				resNames[names_read][clen] = *read_buffer;
				read_buffer += 1;
				clen += 1;
				resNames[names_read][clen] = '\0';
			}
			int len = (1+strlen(resNames[names_read]) );
			resNames[names_read] = (char *)realloc( resNames[names_read], sizeof(char) * len );

//			printf("Structure has %d residues named '%s'.\n", res_counts[names_read], resNames[names_read] );
		 		clen = 0;

			if( !strcasecmp( resNames[names_read], "TIP3" ) )
				block->num_water = res_counts[names_read];
			if( !strcasecmp( resNames[names_read], "W" ) )
				block->num_water = res_counts[names_read];
	
			names_read++;	
			while( *read_buffer && (*read_buffer == ' ' || *read_buffer == '\t') )
				read_buffer += 1;
		}


		av_rho_res = (double **)malloc( sizeof(double *) * nres_used );
		for( int ir = 0; ir < nres_used; ir++ )
			av_rho_res[ir] = (double *)malloc( sizeof(double) * l_N_Z_BINS );	

		int nwra = 0;

		setupSpline( -l_av_Lz/2, l_av_Lz/2, l_N_Z_BINS, 0, 1 );

		for( int z = 0; z < l_N_Z_BINS; z ++ )
		{
			double val = 0;
			double zbin;
	
			double t_vals[2+nres_used];
	
			getLine( bzFile, buffer );
			readNDoubles( buffer, t_vals, 2+nres_used );
	
			zbin = t_vals[0];
			av_rho[z] = t_vals[1];
//			betaz[z] = t_vals[1];			
			
			AddPointToSpline( t_vals[0], av_rho[z], 0 );	

			for( int ir = 0; ir < nres_used; ir++ )
				av_rho_res[ir][z] = t_vals[2+ir];
		}

		SolveSpline(0);

		for( int iz = 0; iz < N_Z_BINS; iz++ )
		{
			double zv = MIN_Z + (iz+0.5) * (MAX_Z-MIN_Z)/N_Z_BINS;

			if( zv < -l_av_Lz/2 || zv >= l_av_Lz/2 )
				betaz[iz] += water_rho_average * block->weight;	
			else
				betaz[iz] += evaluateSpline( zv, 0 ) * block->weight;
		}
			
		block->beta_hist = NULL;
		block->max_bin = 0;

		if( block->brFile ) 
		{
			FILE *hist = fopen(block->brFile,"r");
			
			int l_N_R_BINS;
			double l_cutoff;
			double l_cyl_cutoff=0;
			int cutoff_mode = 0;
			int MAX_HIST = 100;
			double n_hist[MAX_HIST];
			double h0;
			double bi_length[MAX_HIST];	
			double *bi_hist[MAX_HIST];
			double *bi_hist_cn;
			int natoms;


			if( hist )
			{	
		
				getLine( hist, buffer );
				nr = sscanf(buffer, "%d %lf %lf %lf %lf cylinder %lf\n", &l_N_R_BINS, &l_cutoff, &l_av_Lx, &l_av_Ly, &l_av_Lz, &l_cyl_cutoff );	
		
				if( nr == 6 )
				{
					// cylindrical model
					cutoff_mode = 1;
				}
		
				block->max_bin = l_N_R_BINS;
				block->cutoff_mode = cutoff_mode;
				block->cyl_cutoff = l_cyl_cutoff;
				if( fabs(l_cutoff/l_N_R_BINS-max_r_cutoff/N_MAX_R_BINS) > 1e-3 )
				{
					printf("Graining of beta(r) file %s (%lf) does not match large-scale scheme (%lf).\n", block->brFile, l_cutoff/l_N_R_BINS, max_r_cutoff/N_MAX_R_BINS  );
				}
				getLine( hist, buffer );
				natoms  = readNDoubles( buffer+strlen("nhist "), n_hist, MAX_HIST );
		
				h0 = n_hist[0];
		//		printf("h0: %lf\n", h0 );
		
				getLine( hist, buffer );
				readNDoubles( buffer+strlen("blength "), bi_length, MAX_HIST );
		
				bi_hist_cn = (double *)malloc( sizeof(double) * l_N_R_BINS );
				
				getLine( hist, buffer );
				readNDoubles( buffer + strlen("cn "), bi_hist_cn, l_N_R_BINS );
				
				for( int x = 0; x < natoms; x++ )
				{
					bi_hist[x] = (double*)malloc( sizeof(double) * l_N_R_BINS );
					getLine( hist, buffer );
	
					int init_len = 1;
	
					while( init_len < strlen(buffer) && buffer[init_len] != ' ' )
						init_len++;
	
					readNDoubles( buffer + 1+init_len, bi_hist[x], l_N_R_BINS );
				}

				int fileIV = 0;

				getLine( hist, buffer );

				double *termIV = NULL;

				if( !feof(hist) )
				{
					termIV = (double *)malloc( sizeof(double) * l_N_R_BINS );
					
					int init_len = 1;
	
					while( init_len < strlen(buffer) && buffer[init_len] != ' ' && buffer[init_len] != '\t' )
						init_len++;
					int nr = readNDoubles( buffer+1+init_len, termIV, l_N_R_BINS );

					if( nr == l_N_R_BINS )
						fileIV = 1;
				}

		
				fclose(hist);
	
				double sub_rho  = water_rho_average;	
		
				double *c_beta_real = (double *)malloc( sizeof(double) * l_N_R_BINS );
				memset( c_beta_real, 0, sizeof(double) * l_N_R_BINS );
		
				for( int ir = 0; ir < l_N_R_BINS; ir++ )
				{
					double r1 = (ir) * l_cutoff / l_N_R_BINS;
					double r2 = (ir+1) * l_cutoff / l_N_R_BINS;
					double dr = l_cutoff / l_N_R_BINS;
		
					if( ir < l_N_R_BINS )
					{
						double Vb;
					
						if( cutoff_mode )
							Vb = util_volume( -r2, r2, r1, r2, 0, l_cyl_cutoff );
						else
							Vb = (4.0/3.0) * M_PI * (r2*r2*r2-r1*r1*r1);
		
						double rho_i_total_bg = 0;
						double real_value = 0;
						for( int xx = 0; xx < natoms; xx++ )
						{
							if( n_hist[xx] < 1 ) continue;
		
							// the histogram of the all-atom <-> all-atom interactions, including solvent background above the periodic boundary.
							// This is term Ia
							real_value += bi_length[xx] * bi_hist[xx][ir]; 

							// the interaction of the atom with the background everywhere, subtracted.

							// Attention! This is the calculation of term II:
							rho_i_total_bg += bi_length[xx] * (n_hist[xx]) * Vb * sub_rho;
							// the particle histogram goes up to the cutoff (for spherical) or for the entire PFFT histogram (for cylindrical)
							// both of the two preceding terms are not "missing" anything, their interactions complete the particle histogram
						}
//						printf("ir: %d real_value %lf  total_bg %lf hist_cn %lf ", ir, real_value, rho_i_total_bg, bi_hist_cn[ir]*sub_rho );	

						// subtract off the background interacting with all-atom

						// this is term III, computed in betaHist:
						real_value -= bi_hist_cn[ir] * sub_rho;
						// term II:
						real_value -= rho_i_total_bg;
						
		
		
						double layer_vol = l_av_Lx * l_av_Ly * l_av_Lz / N_Z_BINS;
		
						double corr_cross = 0;
						for( int iz = 0; iz < N_Z_BINS; iz++ )
						{
							double zv = iz * l_av_Lz / N_Z_BINS;
			
							double bg_vol;
		
							if( cutoff_mode )
								bg_vol = util_volume( -zv, l_av_Lz-zv, r1, r2, 0, l_cyl_cutoff ); 
							else
								bg_vol = util_volume( -zv, l_av_Lz-zv, r1, r2 ); 

							// we are computing background-background interactions here.
							
							//corr_cross += sub_rho * layer_vol * (vol-bg_vol) * sub_rho;
							corr_cross += sub_rho * layer_vol * bg_vol * sub_rho;
						}

//						printf("corr_cross %lf\n", - corr_cross );
						// this is term IV:
						if( fileIV )
							real_value += termIV[ir] * sub_rho * sub_rho;
						else
							real_value += corr_cross;
				
						c_beta_real[ir] += real_value;
					}
				}	
		
				block->beta_hist = (double *)malloc( sizeof(double) * N_MAX_R_BINS  );
				memset( block->beta_hist, 0, sizeof(double) * N_MAX_R_BINS );
		
				// make any user requested modifications.
		
				for( int ir = 0; ir < N_MAX_R_BINS; ir++ )
				{
					double r1 = (ir) * max_r_cutoff / N_MAX_R_BINS;
					double r2 = (ir+1) * max_r_cutoff / N_MAX_R_BINS;
					double dr = max_r_cutoff / N_MAX_R_BINS;
			
					double continuum_mod = 0;
					double real_value = 0;
			
					if( ir < l_N_R_BINS )
						real_value = c_beta_real[ir]; 
				
					double rmean = (r1+r2)/2;	
			
					for( pfft_instruction *inst = block->instructions; inst; inst = inst->next )
					{
						if( cutoff_mode )
						{
							printf("Instructions NYI in cylindrical cutoff mode.\n");
							exit(1);
						}
		
						if( rmean < inst->R1-smooth_width || rmean > inst->R2 )
							continue;
						double alpha = 1.0;
			
						if( rmean < inst->R1 )
							alpha = smoother( rmean, inst->R1-smooth_width, inst->R1 );
						else if( rmean > inst->R2-smooth_width )
							alpha = 1-smoother( rmean, inst->R2-smooth_width, inst->R2 );
			//			printf("PFFT!! at r=%lf alpha=%lf.\n", rmean, alpha);
			
						// modify continuum value.
				
						double *select_rho = NULL;
						double *sub_rho = NULL;
						double *water_rho = NULL;
			
						for( int ir = 0; ir < nres_used; ir++ )
						{
							if( !strcasecmp( resNames[ir], inst->forRes ) )
								select_rho = av_rho_res[ir];
							if( !strcasecmp( resNames[ir], inst->depleteRes ) )
								sub_rho = av_rho_res[ir];
							if( !strcasecmp( resNames[ir], "TIP3" ) )
								water_rho = av_rho_res[ir];
						}
					
						for( int z = 0; z < N_Z_BINS; z++ )
						{
							double zv = -l_av_Lz/2 + (z+0.5) * l_av_Lz / N_Z_BINS;
							double zvol = l_av_Lx * l_av_Ly * (l_av_Lz/N_Z_BINS);
				
							double bv = select_rho[z];
				
							for( int z2 = 0; z2 < N_Z_BINS; z2++ )
							{
								double bin_z1 = -l_av_Lz/2 + z2 * l_av_Lz / N_Z_BINS - zv;	
								double bin_z2 = -l_av_Lz/2 + (z2+1) * l_av_Lz / N_Z_BINS - zv;	
				
								double approx_vol = util_volume( bin_z1, bin_z2, r1, r2 );
				
								int use_z2 = z2;
								
								continuum_mod += bv * sub_rho[use_z2] * approx_vol * zvol * alpha * inst->mod;
							}
						}
					}
		
					block->beta_hist[ir] = (real_value + continuum_mod);
				}
		
				if( termIV )
					free(termIV);

				free(bi_hist_cn);
				free(c_beta_real);		
			}
		}

//		printf("lx: %lf ly: %lf\n", l_av_Lx, l_av_Ly );
		block->av_area = l_av_Lx*l_av_Ly;
		block->height = l_av_Lz;
	}
	free(buffer);

	for( int iz = 0; iz < N_Z_BINS; iz++ )
	{
		betaz[iz] -= water_rho_average;
//		printf("BZ %lf %lf\n", MIN_Z + (iz+0.5) * (MAX_Z-MIN_Z)/N_Z_BINS, betaz[iz] );
	}
	
	// construct beta(z).

	*water_rho = water_rho_average;
}

double alt_vol_below_within( double z, double R, double rc1, double rc2 )
{
	double phi2 = 0;
	double phi1 = 0;

	if( rc2 < R )
		phi2 = acos(rc2/R);
	if( rc1 < R )
		phi1 = acos(rc1/R);
	double phi1N = -phi1;
	double phi2N = -phi2;

	double phiZ = -M_PI/2;
	if( z > -R  )
	{
		if( z < R )
			phiZ = asin(z/R);
		else
			phiZ = M_PI/2;
	}

	if( phi1 > phiZ ) phi1 = phiZ;
	if( phi2 > phiZ ) phi2 = phiZ;
	if( phi1N > phiZ ) phi1N = phiZ;
	if( phi2N > phiZ ) phi2N = phiZ;

	// integrate from phi2 to phi1
	double Vol1 = (M_PI*R*R*R/12.0) * (9 * sin(phi1) + sin(3*phi1) - 9*sin(phi2) - sin(3*phi2) );
	// integrate from -phi1 to -phi2
	double Vol2 = (M_PI*R*R*R/12.0) * (9 * sin(phi2N) + sin(3*phi2N) - 9*sin(phi1N) - sin(3*phi1N) );
	// integrate the rc2 cylinder
	double dZCyl2 = R * (sin(phi2) - sin(phi2N));
	double VolC2 = dZCyl2 * M_PI * rc2*rc2;
	// integrate the rc2 cylinder
	double dZCyl1 = R * (sin(phi1) - sin(phi1N));
	double VolC1 = dZCyl1 * M_PI * rc1*rc1;

	return Vol1+Vol2+VolC2-VolC1;	
}

// returns the volume, within a shell of a sphere, below a particular z value, but within rc1 and rc2, cylindrical radii.
double vol_below_within( double z, double R, double rc1, double rc2 )
{
	if( rc1 > R ) return 0;

	double phi_stop;

	if( z <= -R )
		phi_stop = 0;
	else if( z >= R )
		phi_stop = M_PI;
	else
		phi_stop = M_PI - acos( z / R );

	double vol = 0;
	double phi1=0,phi2=0;

	phi1 = asin(rc1/R);
	
	if( rc2 < R )
	{	// here we integrate from phi1 to phi2, using R as the outer limit.
		phi2 = asin(rc2/R);
		
		if( phi_stop > phi1 )
		{
			if( phi_stop < phi2 )
				phi2 = phi_stop;
		
			vol += -(M_PI*R*Power(rc1,2)*(cos(phi1) - cos(phi2))) +    (M_PI*Power(R,3)*(9*cos(phi1) - cos(3*phi1) - 9*cos(phi2) + cos(3*phi2)))/12.;
		}
	
		// now we integrate with rc2 < r, we get the inner cylinder. 
		
		phi1 = asin(rc2/R);
		phi2 = M_PI - asin(rc2/R);
		
		if( phi_stop > phi1 )
		{
			if( phi_stop < phi2 )
				phi2 = phi_stop;
		
			vol += M_PI*R*(-Power(rc1,2) + Power(rc2,2))*(cos(phi1) - cos(phi2));
		}

		// now we finish the integral with the last bit that is nearly the same as the first. 
		
		phi1 = M_PI-asin(rc2/R);
		phi2 = M_PI-asin(rc1/R);
		
		if( phi_stop > phi1 )
		{
			if( phi_stop < phi2 )
				phi2 = phi_stop;
		
			vol += -(M_PI*R*Power(rc1,2)*(cos(phi1) - cos(phi2))) +    (M_PI*Power(R,3)*(9*cos(phi1) - cos(3*phi1) - 9*cos(phi2) + cos(3*phi2)))/12.;
		}
	}
	else
	{
		// here the outer radius is all outside, we always use R as the outer limit over the whole arc.
		
		phi2 = M_PI - asin(rc1/R);
		
		if( phi_stop > phi1 )
		{
			if( phi_stop < phi2 )
				phi2 = phi_stop;
		
			vol += -(M_PI*R*Power(rc1,2)*(cos(phi1) - cos(phi2))) +    (M_PI*Power(R,3)*(9*cos(phi1) - cos(3*phi1) - 9*cos(phi2) + cos(3*phi2)))/12.;
		}
	}

	return vol;
}


double util_volume( double z1, double z2, double r1, double r2 )
{
/*
	the volume within a spherical shell (defined between r1 and r2) that is also between z1 and z2.
	Vol1 = sphere(r2) - sphere(r1)
	Vol2 = surface below (z2) - surface below(z1)
	returns volume of (intersection( Vol1, Vol2 ))
*/
	if( z1 > z2 )
	{
		double t = z1;
		z1 = z2;
		z2 = t;	
	}
	
	double vol2 =(4*M_PI/3.0)*r2*r2*r2;
	double vol1 =(4*M_PI/3.0)*r1*r1*r1;

	double vol_2_below_z2 = vol2;
	double vol_2_below_z1 = vol2;
	
	double vol_1_below_z2 = vol1;
	double vol_1_below_z1 = vol1;

	if( z2 < -r2 )
		vol_2_below_z2 = 0;
	else if( z2 < r2 )
		vol_2_below_z2 = (1.0/3.0)*M_PI*(2*r2-z2)*(r2+z2)*(r2+z2);
	
	if( z2 < -r1 )
		vol_1_below_z2 = 0;
	else if( z2 < r1 )
		vol_1_below_z2 = (1.0/3.0)*M_PI*(2*r1-z2)*(r1+z2)*(r1+z2);
	
	if( z1 < -r2 )
		vol_2_below_z1 = 0;
	else if( z1 < r2 )
		vol_2_below_z1 = (1.0/3.0)*M_PI*(2*r2-z1)*(r2+z1)*(r2+z1);
	
	if( z1 < -r1 )
		vol_1_below_z1 = 0;
	else if( z1 < r1 )
		vol_1_below_z1 = (1.0/3.0)*M_PI*(2*r1-z1)*(r1+z1)*(r1+z1);

	double vol = (vol_2_below_z2-vol_2_below_z1) - (vol_1_below_z2-vol_1_below_z1);	

	if( vol < -1e-8 )
	{
		printf("vol: %le r1: %lf r2: %lf z1: %lf z2: %lf\n", vol, r1, r2, z1, z2 ); 
		printf("CRITICAL ERROR.\n");
		exit(1);
	}
	
	return vol;
} 

double util_volume( double z1, double z2, double r1, double r2, double rc1, double rc2 )
{
	if( z1 > z2 )
	{
		double t = z1;
		z1 = z2;
		z2 = t;	
	}
	
	double vol_2_below_z2 = vol_below_within( z2, r2, rc1, rc2 );
	double vol_2_below_z1 = vol_below_within( z1, r2, rc1, rc2 );
	double vol_1_below_z2 = vol_below_within( z2, r1, rc1, rc2 );
	double vol_1_below_z1 = vol_below_within( z1, r1, rc1, rc2 );
	
	double alt_vol_2_below_z2 = alt_vol_below_within( z2, r2, rc1, rc2 );
	double alt_vol_2_below_z1 = alt_vol_below_within( z1, r2, rc1, rc2 );
	double alt_vol_1_below_z2 = alt_vol_below_within( z2, r1, rc1, rc2 );
	double alt_vol_1_below_z1 = alt_vol_below_within( z1, r1, rc1, rc2 );

	double del = fabs(vol_2_below_z2-alt_vol_2_below_z2) +
		     fabs(vol_2_below_z1-alt_vol_2_below_z1) +
		     fabs(vol_1_below_z2-alt_vol_1_below_z2) +
		     fabs(vol_1_below_z1-alt_vol_1_below_z1);

	if( del > 1e-5 )
	{
		printf("FORMULA MISMATCH\n");
		exit(1);
	}
	 
	double vol = (vol_2_below_z2-vol_2_below_z1) - (vol_1_below_z2-vol_1_below_z1);	

	if( vol < -1e-8 )
	{
		printf("vol: %le r1: %lf r2: %lf z1: %lf z2: %lf rc1: %lf rc2: %lf\n", vol, r1, r2, z1, z2, rc1, rc2); 
		printf("CRITICAL ERROR.\n");
		exit(1);
	}
	
	return vol;

} 

double smoother( double r, double r_low, double r_high )
{
	if( r < r_low )
	        return 1;

	if( r >= r_high )
	        return 0;

	return 1.0 / ( 1 + exp(- ((r_high-r_low)/(r_high-r) - (r_high-r_low)/(r-r_low))));
}

void read_residue_bz( const char *fileName, int nres, const char **names, double *max_z, int *nb, double **betaz_out )
{
	char *buffer = (char *)malloc( sizeof(char) * 100000 );

	double water_rho_average = 0;
	double nwater_rho_average = 0;

	
	FILE *bzFile = fopen(fileName, "r");

	if( !bzFile )
	{
		printf("Requested beta.z file \"%s\" not found.\n", fileName);
		exit(1);
	}

	getLine( bzFile, buffer );
	int N_Z_BINS;
	double t_av_Lx, t_av_Ly;
	double l_av_Lx, l_av_Ly, l_av_Lz;
	char code;
	
	int nr = sscanf(buffer, "%d %lf %lf %lf\n", &N_Z_BINS, &t_av_Lx, &t_av_Ly, &l_av_Lz );

	*max_z = l_av_Lz;
	*nb = N_Z_BINS;

	double *betaz = (double *)malloc( sizeof(double) * N_Z_BINS * nres ); 
	*betaz_out = betaz;

	double nres_types[nres+1];

	getLine( bzFile, buffer );
	nr = readNDoubles( buffer, nres_types, nres+1 );
	getLine( bzFile, buffer );
	char *resNames = buffer + strlen("z tot ");

	int map[nres+1];
	double fac[nres+1];

	for( int x = 0; x < nres+1; x++ )
	{
		map[x] = -2;

		for( int y = 0; y < nres; y++ )
		{
			if( !strncasecmp( resNames, names[y], strlen(names[y]) ) )
				map[x] = y;
		}

		if( !strncasecmp( resNames, "TIP3", 4 ) ) 
			map[x] = -1; 
		if( !strcasecmp( resNames, "W" ) || !strcasecmp( resNames, "WD" ) ) 
			map[x] = -1; 
		int n_apl = sizeof(expec_apl) / sizeof(double);
		fac[x] = default_apl;
		for( int t = 0; t < n_apl; t++ )
		{
			if( !strncasecmp( expec_apl_name[t], resNames, strlen(expec_apl_name[t]) ) )
				fac[x] = expec_apl[t];			
		}

		while( *resNames && !( *resNames == ' ' || *resNames == '\t') ) resNames++ ;
		while( *resNames == ' ' || *resNames == '\t' ) resNames++ ;

		if( *resNames == '\0' )
			break; 
	}	
	
	printf("map: ");
	for( int i = 0; i < nres+1; i++ )
		printf(" %d", map[i] );
	printf("\n");

	double sum = 0;
	for( int x = 0; x < nres+1; x++ )
	{
		if( map[x] >= 0 )
			sum += nres_types[x] * fac[x];
	}
	for( int x = 0; x < nres+1; x++ )
	{
		if( map[x] >= 0 )
			nres_types[x] *= fac[x] / sum;

		printf("nres_types[%d] %lf\n", x, (double) nres_types[x] );
	}
	for( int z = 0; z < N_Z_BINS; z ++ )
	{
		double val = 0;
		double zbin;

		double t_vals[2+nres+1];
		getLine( bzFile, buffer );	
		readNDoubles( buffer, t_vals, 2+nres+1 );

		double zv = t_vals[0];
		double rv = t_vals[1];

		if( (zv < -l_av_Lz/2 + 2.0 ) ||
		    (zv >  l_av_Lz/2 - 2.0 ) )
		{
			water_rho_average += rv;
			nwater_rho_average += 1;
		}	

	
		for( int y = 0; y < nres+1; y++ )
		{
			if( map[y] >= 0 )
				betaz[map[y]*N_Z_BINS+z] += t_vals[2+y] / nres_types[y];	

			if( map[y] == -1 )
			{
				for( int x = 0; x < nres; x++ )
					betaz[x*N_Z_BINS+z] += t_vals[2+y]; // add in the water beta.
			}
		}
	}

	water_rho_average /= nwater_rho_average;
		
	for( int x = 0; x < nres; x++ )
	for( int z = 0; z < N_Z_BINS; z++ )
		betaz[x*N_Z_BINS+z] -= water_rho_average;
	
	for( int z = 0; z < N_Z_BINS; z++ )
	{
		printf("%d ", z );
		for( int x = 0; x < nres; x++ )
			printf(" %lf", betaz[x*N_Z_BINS+z] );
		printf("\n");
	}		
		
	fclose(bzFile);
}

void disk_disk_hist( double *beta1z, double *beta2z, double dzbin, int nzbins, double r_disk, double *hist, double bin_width, int nrbins )
{
	double r_hist = nrbins * bin_width;
	double dr_int = 1.0;
	int n_rbins = r_disk / dr_int + 1;

	double dr_use = r_disk / n_rbins;

	double sum_prof = 0;

	for( int iz = 0; iz < nzbins; iz++ )
		sum_prof += beta1z[iz] * dzbin;
 
	for( int z1 = 0; z1 < nzbins; z1++ )
	for( int z2 = 0; z2 < nzbins; z2++ )
	{
		double dz = (z2-z1) * dzbin;

		for( int i_r1 = 0; i_r1 < n_rbins; i_r1++ )
		{
			double r1 = (i_r1+0.5)*dr_use;

			for( int i_r2 = 0; i_r2 < n_rbins; i_r2++ )
			{
				double r2 = (i_r2+0.5)*dr_use;
						
				double the_factor =  2 * M_PI * r1 * r2 * 2 * dzbin * dzbin * beta1z[z1] * beta2z[z2] * dr_use * dr_use;

				double del2 = r1*r1+r2*r2+dz*dz;
				double rmin = sqrt(1e-14+del2 - 2 * r1 * r2);
				double rmax = sqrt(1e-14+del2 + 2 * r1 * r2);

				int bin_min = (rmin/r_hist) * nrbins;
				int bin_max = (rmax/r_hist) * nrbins;

				if( bin_min == bin_max && bin_min < nrbins)
					hist[bin_min] += the_factor * 2*M_PI; 
				else
				{
					double th_cur = 0;
					
					double sum = 0;

					int broken = 0;
					for( int ib = bin_min; ib <= bin_max; ib++ )
					{
						if( ib >= nrbins ) { broken = 1; break; }
								
						double next_bin = (ib+1) * bin_width;

						double next_theta;

						double TOL = 1e-8;

						double val = (del2 - next_bin*next_bin) / (2*r1*r2) ;

						if( val <= -1.0 && val >= -1.0 - TOL )
							val = -1.0 + TOL;
						
						if( val >= 1.0 && val <=  1.0 + TOL )
							val = 1.0 - TOL;

						if( ib == bin_max )
							next_theta = M_PI;
						else
							next_theta = acos( val );
						
						double dth = next_theta - th_cur;
						
						if( dth > M_PI ) dth -= 2*M_PI;
						if( dth < -M_PI ) dth += 2*M_PI;

						// last factor of two because we are only integrating from 0 to Pi and doubling.
						hist[ib] += the_factor * dth;
				
						th_cur = next_theta;
						sum += dth;
					}
				}
			}
		}
	}

	double tsum = 0;
	for( int ir = 0; ir < nrbins; ir++ )
	{
		double r = (ir+0.5) * bin_width;
		tsum += hist[ir];
	}

	printf("val: %le debug_val: %le\n", tsum, pow( sum_prof * r_disk * r_disk * M_PI, 2.0 ) ); 

}


void r2_disk_disk_hist( double *beta1z, double *beta2z, double dzbin, int nzbins, double r_disk1, double r_disk2, double *hist, double bin_width, int nrbins )
{
	double r_hist = nrbins * bin_width;
	double dr_int = 1.0;
	int n_rbins1 = r_disk1 / dr_int + 1;
	int n_rbins2 = (r_disk2-r_disk1) / dr_int + 1;

	double dr_use1 = r_disk1 / n_rbins1;
	double dr_use2 = (r_disk2-r_disk1) / n_rbins2;

	double sum_prof = 0;

	for( int iz = 0; iz < nzbins; iz++ )
		sum_prof += beta1z[iz] * dzbin;
 
	for( int z1 = 0; z1 < nzbins; z1++ )
	for( int z2 = 0; z2 < nzbins; z2++ )
	{
		double dz = (z2-z1) * dzbin;

		for( int i_r1 = 0; i_r1 < n_rbins1; i_r1++ )
		{
			double r1 = (i_r1+0.5)*dr_use1;

			for( int i_r2 = 0; i_r2 < n_rbins2; i_r2++ )
			{
				double r2 = r_disk1 + (i_r2+0.5)*dr_use2;
						
				double the_factor =  2 * M_PI * r1 * r2 * 2 * dzbin * dzbin * beta1z[z1] * beta2z[z2] * dr_use1 * dr_use2;

				double del2 = r1*r1+r2*r2+dz*dz;
				double rmin = sqrt(1e-14+del2 - 2 * r1 * r2);
				double rmax = sqrt(1e-14+del2 + 2 * r1 * r2);

				int bin_min = (rmin/r_hist) * nrbins;
				int bin_max = (rmax/r_hist) * nrbins;

				if( bin_min == bin_max && bin_min < nrbins)
					hist[bin_min] += the_factor * 2*M_PI; 
				else
				{
					double th_cur = 0;
					
					double sum = 0;

					int broken = 0;
					for( int ib = bin_min; ib <= bin_max; ib++ )
					{
						if( ib >= nrbins ) { broken = 1; break; }
								
						double next_bin = (ib+1) * bin_width;

						double next_theta;

						double TOL = 1e-8;

						double val = (del2 - next_bin*next_bin) / (2*r1*r2) ;

						if( val <= -1.0 && val >= -1.0 - TOL )
							val = -1.0 + TOL;
						
						if( val >= 1.0 && val <=  1.0 + TOL )
							val = 1.0 - TOL;

						if( ib == bin_max )
							next_theta = M_PI;
						else
							next_theta = acos( val );
						
						double dth = next_theta - th_cur;
						
						if( dth > M_PI ) dth -= 2*M_PI;
						if( dth < -M_PI ) dth += 2*M_PI;

						// last factor of two because we are only integrating from 0 to Pi and doubling.
						hist[ib] += the_factor * dth;
				
						th_cur = next_theta;
						sum += dth;
					}
				}
			}
		}
	}

	double tsum = 0;
	for( int ir = 0; ir < nrbins; ir++ )
	{
		double r = (ir+0.5) * bin_width;
		tsum += hist[ir];
	}


}
