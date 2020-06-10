#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "util.h"
#define TOL (1e-8)

#include "fit_control.h"

#ifdef USE_SCIPY

void NNLS( double *A, double *b,double *sol,  int n, int m, double eps )
{
	// write to disk.
	FILE *interface= fopen("interface.py","w");

	fprintf(interface,
"import numpy as np\n"
"import scipy.optimize\n"
"A = np.loadtxt('scipy_A.inp')\n"
"b = np.loadtxt('scipy_b.inp')\n"
"x, norm = scipy.optimize.nnls(A,b)\n"
"print(x)\n");
	fclose(interface);
	
	FILE *inputFile = fopen("scipy_A.inp","w");
	for( int j = 0; j < m; j++ )
	{
		for( int i = 0; i < n; i++ )
			fprintf(inputFile, "%.14le ", A[j*n+i] );
		fprintf(inputFile,"\n");
	}
	fclose(inputFile);

	inputFile = fopen("scipy_b.inp","w");
	for( int i = 0; i < m; i++ )
		fprintf(inputFile, "%.14le ", b[i] );
	fprintf(inputFile,"\n");
	fclose(inputFile);

	system("python interface.py > interface.out");
	
	FILE *outFile = fopen("interface.out","r");
	char *buffer = (char *)malloc( sizeof(char) * 100000 );
	getLine(outFile,buffer);

	int nr = readNDoubles( buffer+2, sol, n );

		
}

#else

#define dgetrf dgetrf_
#define dgetri dgetri_

extern "C" void dgetrf( int *, int *, double *, int *, int *, int * );
extern "C" void dgetri( int *, double *, int *, int *, double *, int *, int * );

void getAPinv( double *Ap, double *A, int *set_P, int np, int n );
void getSP( double *s, double *Ap, double *A, double *y, int *set_P, int np, int n ); 
void getSimpleSP( double *s, double *A, double * y, int *set_P, int np, int n );

void NNLS( double *A, double *b, double *sol, int n, int m, double eps )
{
	// condition the problem to deal with our tolerances.

	double maxa = 0;
#ifdef SCALE_IT

	for( int i = 0; i < n; i++ )
	{
		if( fabs(A[i*n+i]) > maxa )
			maxa = fabs(A[i*n+i]);
	}

	for( int j = 0; j < n; j++ )
	{
		for( int k = 0; k < n; k++ )
		{
			A[j*n+k] /= maxa;
		}
		b[j] /= maxa;
	}
#endif
	int set_P[n];
	int set_R[n];
	int np = 0;
	int nr = n;

	for( int i = 0; i < n; i++ )
		set_R[i] = i;

	double x[n];
	memset(x,0,sizeof(double)*n);
	double w[n];
	double tw[n];

	// tw = y - A x
	for( int i = 0; i < n; i++ )
	{
		w[i] = b[i];

		for( int j = 0; j < n; j++ )
			w[i] -= A[i*n+j] * x[j];
	}

	printf("Ap={\n");
	for( int k = 0; k < n; k++ )
	{
		if( k == 0 ) 
			printf("{");
		else
			printf(",\n{");
		for( int j = 0; j < n; j++ )
		{
			if( j == 0 )
				printf("%le", A[k*n+j] );
			else
				printf(", %le", A[k*n+j] );
		}
		printf("}");
	}
	printf("};\n");
	
	printf("b={\n");
	for( int k = 0; k < n; k++ )
	{
		if( k == 0 )
			printf("%lf", b[k] );
		else
			printf(", %lf", b[k] );
	}
	printf("};\n");
	double *s = (double *)malloc( sizeof(double) * n );
	double *ts = (double *)malloc( sizeof(double) * n );
	double *Ap = (double *)malloc( sizeof(double) * n * n );
	double *work = (double *)malloc( sizeof(double) * n * n );
	int imaxw=0;
	double maxw=w[0];
	do {
		imaxw=0;
		maxw=w[0];
		for( int i = 1; i < n; i++ )
		{
			if( w[i] > maxw )
			{
				maxw = w[i];
				imaxw = i;
			}
		}	

		if( maxw < TOL ) break;

		int imaxwr = 0;
		double maxwr = w[set_R[0]];
		
		for( int i = 1; i < nr; i++ )
		{
			if( w[set_R[i]] > maxwr )
			{
				maxwr = w[set_R[i]];
				imaxwr = i;
			}
		}	

		set_P[np] = set_R[imaxwr]; np++;
		set_R[imaxwr] = set_R[nr-1]; nr--;

//		getAPinv( Ap, A, set_P, np, n );
		memset(s,0,sizeof(double)*n);
//		getSP( s, Ap, A, b, set_P, np, n ); 
		getSimpleSP( s, A, b, set_P, np, n );	
		double min_s;

		do{ 
			int ainit =0;
			double alpha = 0;
			min_s = s[set_P[0]];

			for( int i = 0; i < np; i++ )
			{
				if( s[set_P[i]] < -TOL )
				{
					if( s[set_P[i]] < min_s )
						min_s = s[set_P[i]];

					double alpha_check = x[set_P[i]] / ( x[set_P[i]] - s[set_P[i]]);

					if( alpha_check < alpha || !ainit )
					{
						alpha = alpha_check;
						ainit = 1;
					}
				}
			}

			printf("min_s: %le\n", min_s );
			if( min_s < -TOL )
			{
				for( int i = 0; i < n; i++ )
					x[i] += alpha * (s[i] - x[i]);
				for( int i = 0; i < np; i++ )
				{
					if( fabs(x[set_P[i]]) < TOL )
					{
						set_R[nr] = set_P[i]; nr++;
						set_P[i] = set_P[np-1]; np--;
					}
				}

				if( np > 0 )
				{
	//				getAPinv( Ap, A, set_P, np, n );
					memset(s,0,sizeof(double)*n);
	//				getSP( s, Ap, A, b, set_P, np, n ); 					
					getSimpleSP( s, A, b, set_P, np, n );	
				}
			}	
		} while( min_s < -TOL && np > 0 );

		memcpy( x, s, sizeof(double) * n );
	
		// tw = y - A x
		for( int i = 0; i < n; i++ )
		{
			w[i] = b[i];
	
			for( int j = 0; j < n; j++ )
				w[i] -= A[i*n+j] * x[j];
		}

		printf("w:\n");
		for( int i = 0; i < n; i++ )
			printf(" %le\n", w[i] );
		printf("\n");
	} while( maxw > eps );

	memcpy( sol, x, sizeof(double) * n );
}

void getSimpleSP( double *s, double *A, double * y, int *set_P, int np, int n )
{
	double *Ap = (double *)malloc( sizeof(double) * np * np );

	for( int i = 0; i < np; i++ )	
	for( int j = 0; j < np; j++ )
		Ap[i*np+j] = A[set_P[i]*n+set_P[j]];
	
	double *work = (double *)malloc( sizeof(double) * np * np );
	int info;
	int ipiv[n];
	dgetrf( &np, &np, Ap, &np, ipiv, &info );
	int lwork = np*np;
	dgetri( &np, Ap, &np, ipiv, work, &lwork, &info );
	free(work);

	for( int i = 0; i < np; i++ )
	for( int j = 0; j < np; j++ )
		s[set_P[i]] += Ap[i*np+j] * y[set_P[j]];
}

void getAPinv( double *Ap, double *A, int *set_P, int np, int n)
{
	double *work = (double *)malloc( sizeof(double) * np * np );

	// Ap inv
	for( int i = 0; i < np; i++ )
	for( int j = 0; j < np; j++ )
	{
		Ap[i*np+j] = 0;

		for( int k = 0; k < np; k++ )
			Ap[i*np+j] += A[set_P[k]*n+set_P[i]] * A[set_P[k]*n+set_P[j]];
	}

	int info;
	int ipiv[n];
	dgetrf( &np, &np, Ap, &np, ipiv, &info );
	int lwork = np*np;
	dgetri( &np, Ap, &np, ipiv, work, &lwork, &info );

	free(work);
} 

void getSP( double *s, double *Ap, double *A, double *y, int *set_P, int np, int n ) 
{
	double *ts = (double *)malloc( sizeof(double) * np );

	// A^P^T . y
	for( int i = 0; i < np; i++ )
	{
		ts[i] = 0;

		for( int j = 0; j < np; j++ )
			ts[i] += A[set_P[j]*n+set_P[i]] * y[set_P[j]];
	}		

	for( int i = 0; i < np; i++ )
	{
		s[set_P[i]] = 0;
		for( int j = 0; j < np; j++ )
			s[set_P[i]] += Ap[i*np+j] * ts[j];
	}

	free(ts);
}
#endif
