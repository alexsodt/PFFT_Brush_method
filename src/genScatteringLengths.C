// by alex sodt
//#define DEBUG_BR 
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

int main( int argc, char **argv )
{
          struct timeval tp;
 
          gettimeofday( &tp, NULL );
 
          srand((long)tp.tv_usec);

	char *buffer = (char *)malloc( sizeof(char) * 100000 );

	if( argc != 2)
	{
		printf("Syntax: genScatteringLengths.exe psf\n");	
		return 0;
	}

	FILE *psfFile = fopen(argv[1], "r" );

	if( ! psfFile )
	{
		printf("Couldn't open PSF file '%s'.\n", argv[1] );
		return 0;
	}

	if( !strcasecmp( argv[1] + strlen(argv[1])-3, "pdb" ) )
		loadPSFfromPDB( psfFile );	
	else
		loadPSF( psfFile );
	
	rewind(psfFile);

	struct deuteration
	{
		char resname[256];
		char atname[256];
		double fraction;
	};

	int ndeutSpace = 10;
	int ndeut = 0;
	deuteration *theDeuteration = (deuteration *)malloc( sizeof(deuteration) * ndeutSpace );

	while(!feof(psfFile))
	{
		getLine( psfFile, buffer );

		if( feof(psfFile) ) break;

		if( !strncasecmp( buffer, " REMARKS deuterate", strlen(" REMARKS deuterate") ) )
		{
			if( ndeut == ndeutSpace )
			{
				ndeutSpace *= 2;
				theDeuteration = (deuteration *)realloc( theDeuteration, sizeof(deuteration) * ndeutSpace );
			}

			sscanf( buffer, " REMARKS deuterate %s %s fraction %lf\n", 
					theDeuteration[ndeut].resname, theDeuteration[ndeut].atname, &(theDeuteration[ndeut].fraction) );

			ndeut++;
		}
	}

	fclose(psfFile);	
	
	int nat = curNAtoms();

	struct atom_rec *at = (struct atom_rec *)malloc( sizeof(struct atom_rec) * curNAtoms() );

	loadInfo(at);

	double *b_i = (double *)malloc( sizeof(double) * curNAtoms() );
	int    *b_ind = (int *)malloc( sizeof(int) * curNAtoms() );
	int    *b_res = (int *)malloc( sizeof(int ) * curNAtoms() );		
				
	const char *bi_atoms=  "HDCNOPSX";
	int natoms = strlen(bi_atoms);
#ifdef DO_SAXS
	double bi_length[8] = { 1, 1, 6, 7, 8, 15, 16, 1 };
#else
	double bi_length[8] = { -0.3742, 0.6671, 0.6651, 0.940, 0.5804, 0.517, 0.2847, 0.64 /* 87% D 13% H */ };
	//double bi_length[8] = { -0.3742, 0.6671, 0.6651, 0.940, 0.5804, 0.517, 0.2847, 0.531731 /* 87% D 13% H */ };
#endif
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
		-0.3332, // methylene, flexible
		0.4152, // methylene, saturated
		-0.672, // water
		7.6584, // water, deuterated
	 	9.6597, // choline, deuterated
		-0.3332,  // cholesterol ROH group, arb
		-0.3332  // cholesterol R group, arb
	};
	
	int natoms_martini = sizeof(martini_lengths)/sizeof(double);
	int max_natoms = natoms;

	if( natoms_martini > max_natoms )
		max_natoms = natoms_martini;
	int martini_detected = 0;

	for( int a = 0; a < curNAtoms() && ! martini_detected; a++ )
	{
		if( !strcasecmp( at[a].resname, "W" ) )
			martini_detected = 1;
		if( !strcasecmp( at[a].atname, "NC3") )
			martini_detected = 1;				
		if( !strcasecmp( at[a].atname, "NH3") )
			martini_detected = 1;				
	}

	for( int a = 0; a < curNAtoms(); a++ )
	{				
		double b = 0;

		int deut = 0;	
		for( int d = 0; d < ndeut; d++ )
		{
			if( !strcasecmp( at[a].resname, theDeuteration[d].resname ) &&
		            !strcasecmp( at[a].atname,  theDeuteration[d].atname ) )
			{
				deut=1;
				b = bi_length[1] * (theDeuteration[d].fraction) + bi_length[0] * (1-theDeuteration[d].fraction);
				break;
			}
		}

		if( deut )
		{
		}
		else if( martini_detected )
		{
			for( int p = 0; p < natoms; p++ )
			{
				if( !strncasecmp( at[a].atname, martini_atoms[p], strlen(martini_atoms[p]) ) )
				{
					b = martini_lengths[p];	
					break;
				}
			}
		}
		else
		{
			for( int p = 0; p < natoms; p++ )
			{
				if( at[a].atname[0] == bi_atoms[p] )
				{
					b= bi_length[p];
				}
			}
		}
		printf("%lf\n", b );
	}
}


