#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "util.h"
#include <math.h>
#include "fit_control.h"

#define SCALE_Q

void NNLS( double *A, double *b, double *sol, int n, int m, double eps );
extern "C" void dgesv_( int *, int *, double *, int *, int *, double *, int *, int *);
int main( int argc, char **argv )
{
	if( argc < 5 ) 
	{
		printf("Syntax: bestFitScaling.exe data1 cut Rvesicle srcData1 [srcData2 ...]\n");
		printf("\tIncluding Rvesicle computes I(0) for comparison to a vesicle sample.\n");
		exit(1);
	}

	int do_q0 = 0;

	double pow_coeff = 0.0;

#ifdef SCALE_Q
	pow_coeff = 2.0;
#elif defined(SCALE_Q2)
	pow_coeff = 4.0;
#endif

	double Rves=500;

	Rves = atof(argv[3]);

	if( Rves >= 0 )
		do_q0 = 1;
	else
		Rves = 500;

	int nq = 0;
	int nqs = 10;

	double *qv = (double *)malloc( sizeof(double) * nqs );

	int nfiles = 2 + argc - 5;

	double *fvs = (double *)malloc( sizeof(double) * nqs * nfiles );

	double *fv1 = (double *)malloc( sizeof(double) * nqs );
	double *fv2 = (double *)malloc( sizeof(double) * nqs );
	double low_cut = atof(argv[2]);
	double fit_cut = 0.08;
	char *buffer= (char *)malloc( sizeof(char) * 100000 );

	double sim_q0 = 0;
	int got_sim_q0 = 0;

	int arg_list[nfiles];
	arg_list[0] = 1;
	int narg = 1;
	for( int t = 4; t < argc; t++, narg++ )
		arg_list[narg] = t;

	for( int cx = 0; cx < narg; cx++ )
	{
		int c = arg_list[cx];

		FILE *theFile = fopen(argv[c],"r");
	
		int done = 1;

		double last_val = 0;
		int init =1;

		if( done )
		{
			int nql = 0;
			while( !feof(theFile ) )
			{
				getLine( theFile, buffer );
				if( feof(theFile) ) break;	
				if( buffer[0] == '#' ) continue;
				if( strlen(buffer) < 3 ) continue;
				if( !strncasecmp( buffer, "r: ", 3 ) ) continue;
				if( !strncasecmp( buffer, "water_rho: ", 9 ) ) continue;
				if( !strncasecmp( buffer, "beta_per", 8 ) ) continue;
				if( !strncasecmp( buffer, "Perfect", 7 ) ) continue;
				if( !strncasecmp( buffer, "q algo", 6 ) ) continue;
				if( !strncasecmp( buffer, "slope", 5 ) ) continue;
				
				if( !strncasecmp( buffer, "q0 ", 3 ) )
				{
					got_sim_q0 = 1;
					double qval0=0;
					int nr = sscanf( buffer, "q0 %lf %lf", &qval0, &sim_q0 );
					continue;
				}


				double q, v;
				int nr = sscanf( buffer, "%lf %lf", &q, &v );

				if( q < 1e-10 ) continue;

				if( nql >= nqs )
				{
					nqs *= 2;

					fvs = (double *)realloc( fvs, sizeof(double) * nqs * nfiles );
					memset( fvs + nql * nfiles, 0, sizeof(double) * nfiles * (nqs-nql) );
					fv1 = (double *)realloc( fv1, sizeof(double) * nqs );
					memset( fv1+nql, 0, sizeof(double) * (nqs-nql) );
					fv2 = (double *)realloc( fv2, sizeof(double) * nqs );
					memset( fv2+nql, 0, sizeof(double) * (nqs-nql) );
					qv = (double *)realloc( qv, sizeof(double) * nqs );
				}
				
				if( nql < nq && fabs(q-qv[nql]) > 1e-4 ) 
				{
					printf("Q mismatch for file %s, %lf != %lf.\n", argv[c], q, qv[nql] );
					exit(1);
				}

				if( c == 2 && fabs(qv[nql] -q) > 1e-4 )
				{
					printf("q mismatched at %le %le\n", qv[nql], q );
					exit(1);
				}
				qv[nql] = q;

				fvs[nql*nfiles+cx] = v;	

				if( cx == 0 )
					fv1[nql] = v;
				if( cx == 1 )
					fv2[nql] = v;
				nql++;
			}

			if( nql > nq ) nq = nql;
		}	
	}
	
	// scipy hack
	int start_iq = 0;
	for( int iq = 0; iq < nq; iq++ )
	{
		if( qv[iq] > low_cut )
		{
			start_iq = iq;	
			break;
		}
	}
	double *sol = (double *)malloc( sizeof(double) * nfiles );
	memset( sol, 0, sizeof(double) * nfiles );

#ifdef USE_SCIPY 
	double *Amn = (double *)malloc( sizeof(double) * nfiles * nq );
	double *bm = (double *)malloc( sizeof(double) * nq );

	memcpy( Amn, fvs, sizeof(double) * nfiles * nq );

	for( int i = 0; i < nq; i++ )
	{
		bm[i] = fvs[i*nfiles+0];
		Amn[i*nfiles+0] = 1.0;
	}
	
	
	NNLS( Amn+start_iq*nfiles, bm+start_iq, sol, nfiles, nq-start_iq, 1e-10 );

#else
	// best fit, matrix style.

	double *A = (double *)malloc( sizeof(double) * nfiles * nfiles ); // nfiles adds in incoherent scattering and subtracts out the target.
	memset( A, 0, sizeof(double) * nfiles * nfiles );
	double *b = (double *)malloc( sizeof(double) * nfiles );	
	memset( b, 0, sizeof(double) * nfiles );
	for( int iq = 0; iq < nq; iq++ )
	{
		if( qv[iq] > low_cut )
		{
			// b row.
			for( int f2 = 1; f2 < nfiles; f2++ )
			{
				A[0*nfiles+f2] += fvs[iq*nfiles+f2]*pow(qv[iq],pow_coeff);  // x
				A[f2*nfiles+0] += fvs[iq*nfiles+f2]*pow(qv[iq],pow_coeff);  // x
			}

			A[0] += pow(qv[iq],pow_coeff);
			b[0] += fvs[iq*nfiles+0]*pow(qv[iq],pow_coeff); // y
		}
	}

	for( int f = 1; f < nfiles; f++ )
	for( int iq = 0; iq < nq; iq++ )
	{
		if( qv[iq] > low_cut )
		{
			for( int f2 = 1; f2 < nfiles; f2++ )
				A[f*nfiles+f2] += fvs[iq*nfiles+f2] * fvs[iq*nfiles+f]*pow(qv[iq],pow_coeff);

			b[f] += fvs[iq*nfiles+0] * fvs[iq*nfiles+f]*pow(qv[iq],pow_coeff);
		}	
	}

	NNLS( A, b, sol, nfiles, nq-start_iq, 1e-10 );
#endif

	printf("sol:\n");
	for( int x = 0; x < nfiles; x++ )
		printf("%.14le\n", sol[x] );
/*
	int nrhs = 1;
	int ipiv[nfiles];
	int info = 0;

	dgesv_( &nfiles, &nrhs, A, &nfiles, ipiv, sol, &nfiles, &info );

	if( info != 0 )
	{
		printf("ERROR solving matrix equation.\n");
		exit(1);
	}
*/	
	if( nfiles == 2 )
	{
		printf("b: %lf k: %lf\n", sol[0], sol[1] );
	}

	if( nfiles == 2 )
	{
		double b = sol[0];
		double k = sol[1];
		double Asample = 1.0 / k;
			
		double chi2 = 0;
		int chi2_num = 0;
		
		double chi2_full = 0;
		int chi2_full_num = 0;
		for( int iq = 0; iq < nq; iq++ )
		{
			if( qv[iq] > fit_cut )
			{
				chi2_full += pow(fv2[iq] * k + b - fv1[iq],2.0);
				chi2_full_num++;
			}
			if( qv[iq] > low_cut )
			{
				chi2 += pow(fv2[iq] * k + b - fv1[iq],2.0);
				chi2_num++;
			}
		}
		chi2_full /= chi2_full_num;
		chi2/=chi2_num;
	
		if( do_q0 && got_sim_q0 )
		{
			double bbar2 = sim_q0 / (2*M_PI);
			double A_v = 4 * M_PI * Rves * Rves;
			double N_v = Asample / A_v;
			printf("bbar2: %le\n", bbar2 );		
	
			printf("%le %le XX (Vesicle transform)\n", 0.0, bbar2 * Asample / M_PI );
			
		}
		else if( do_q0 )
		{
			printf("Failed to read q0 value from file. Cannot compute q=0 scattering of simulation.\n");
		}
	
		for( int iq = 0; iq < nq; iq++ )
		{
			printf("%le %le %le", qv[iq], fv2[iq] * k + b, fv1[iq] );
	
			if( iq == nq-1 )
				printf(" chi2 %le %le", chi2, chi2_full );
	
			printf("\n");
		}
	}
	else
	{
		double b = sol[0];	
		double chi2 = 0;
		int chi2_num = 0;
		
		double chi2_full = 0;
		int chi2_full_num = 0;
		for( int f = 1; f < nfiles; f++ )
		for( int iq = 0; iq < nq; iq++ )
		{
			if( qv[iq] > fit_cut )
			{
				chi2_full += pow(fvs[iq*nfiles+f] * sol[f] + b - fvs[iq*nfiles],2.0);
				chi2_full_num++;
			}
			if( qv[iq] > low_cut )
			{
				chi2 += pow(fvs[iq*nfiles+f] * sol[f] + b - fvs[iq*nfiles],2.0);
				chi2_num++;
			}
		}
		chi2_full /= chi2_full_num;
		chi2/=chi2_num;
	
		for( int iq = 0; iq < nq; iq++ )
		{
			double signal = b;

			for( int f = 1; f < nfiles; f++ )
				signal += sol[f] * fvs[iq*nfiles+f];

			printf("%le %le %le %le", qv[iq], signal, fvs[iq*nfiles+0], fvs[iq*nfiles+0]-signal );
	
			if( iq == nq-1 )
				printf(" chi2 %le %le", chi2, chi2_full );
	
			printf("\n");
		}

		printf("Simulation weights:\n");
		int sim = 1;
		for( int f = 4; f< argc; f++, sim++ )
			printf("%s %le\n", argv[f], sol[sim] );
	}
}
