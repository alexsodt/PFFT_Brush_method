#ifndef __hsh__
#define __hsh__

int hs_get_ndiv( void );
double hs_transform( double z_new, double z_old, double target_spacing, double src_spacing, double x, int isel=-1, int weight=1 );


#endif
