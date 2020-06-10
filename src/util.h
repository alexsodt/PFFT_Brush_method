#include <stdio.h>

#ifndef __utilh__
#define __utilh__


void getLine( FILE *theFile, char *theBuffer );
int readNDoubles( char *buffer, double *vals, int nvalues );
int readNInts( char *buffer, int *vals, int nvalues );
void print5( int val, char *str );

#endif
