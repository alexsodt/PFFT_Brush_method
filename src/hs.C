#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int ndiv = 20;

int hs_get_ndiv( void )
{
	return ndiv;
}

double hs_transform( double z_new, double z_old, double target_spacing, double src_spacing, double x, int isel, int weight  )
{
	x = fabs(x);
	
	double min_target = z_new - target_spacing/2;
	double max_target = z_new + target_spacing/2;
	
	double min_hs = z_old - src_spacing/2;
	double max_hs = z_old + src_spacing/2;

	double exp_target = 5;
	double dlambda = 2*sqrt(exp_target/x) / ndiv; 	

//	printf("x: %le\n", x );
	double max_lambda = dlambda * ndiv / 2;
//	printf("max_lambda: %le\n", max_lambda );
	if( max_lambda > 0.3 )
	{
		max_lambda = 0.3;
		dlambda = 2 * max_lambda / ndiv;
	}
	double nrm = sqrt(M_PI)*erf(max_lambda*sqrt(x))/sqrt(x);
	double val = 0;


	for( int id = 0; id <= ndiv; id++ )
	{
		if( isel >= 0 && isel != id ) 
			continue;

		double lam = -max_lambda + id * dlambda;

		double fac = 1.0;
		if( id == 0 || id == ndiv )
			fac = 0.5;
	
		double ampl;

		if( weight )
			ampl = exp(-lam*lam*x ) * (1-lam)*(1-lam) * dlambda * fac / nrm;
		else
			ampl = (1-lam)*(1-lam);


		// transformed heaviside.
		double min_z_t = min_hs * (1+2*lam);
		double max_z_t = max_hs * (1+2*lam);

		if( max_z_t < min_target ) 
			continue;
		if( max_target < min_z_t )
			continue;
		if( min_z_t >= min_target && max_z_t <= max_target )
			val += ampl;			
		else if( min_z_t <= min_target && max_z_t >= max_target )
			val += ampl * (max_target-min_target) / (max_z_t-min_z_t);
		else if( min_z_t <= min_target )
			val += ampl * (max_z_t-min_target) / ( max_z_t-min_z_t);
		else if( max_z_t >= max_target )
			val += ampl * (max_target-min_z_t) / ( max_z_t-min_z_t);
		else
		{
			printf("Logical error.\n");
			exit(1);
		}
	}

	return val;
	
}
