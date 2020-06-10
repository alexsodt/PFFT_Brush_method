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
#include "deuteration.h"
#include "input.h"

/*
 * the program currently breaks deuteration into two blocks assuming the leaflets are symmetric.
 *
 * */


int main( int argc, char **argv )
{
	struct timeval tp;
	
	gettimeofday( &tp, NULL );
	
	srand((long)tp.tv_usec);

	char buffer[4096];

	if( argc != 3 )
	{
		printf("Syntax: deuterate psf input.file\n");
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

	fclose(psfFile);	
	
	int nat = curNAtoms();

	struct atom_rec *at = (struct atom_rec *)malloc( sizeof(struct atom_rec) * curNAtoms() );
	loadInfo(at);

	struct parameter_block *ablock = (struct parameter_block *)malloc( sizeof(struct parameter_block) );
	getInput(argv[2], ablock );
	
	int nspecies = 0;
	int nspecies_space = 10;
	char **species = (char **)malloc( sizeof(char *) * nspecies_space );

	deuterate_request *deut = NULL;

	for( deut = ablock->deuterations; deut; deut = deut->next )
	{
		int gotit = 0;

		for( int x = 0; x < nspecies; x++ )
		{
			if( !strcasecmp( deut->resName, species[x] ) )
				gotit = 1;
		}

		if( !gotit )
		{
			if( nspecies == nspecies_space )
			{
				nspecies_space *= 2;
				species = (char **)realloc( species, sizeof(char *) * nspecies_space );
			}
			
			species[nspecies] = (char *)malloc( sizeof(char) * (1 + strlen(deut->resName) ) );
			strcpy( species[nspecies], deut->resName );
			nspecies++;
		}
	} 
		
	int *species_count = ( int*) malloc( sizeof(int) * nspecies );
	memset( species_count, 0, sizeof(int) * nspecies );
	int pres = -1;
	char presName[256];
	sprintf(presName, "INIT" );

	for( int a = 0; a < curNAtoms(); a++ )
	{
		if( at[a].res == pres && !strcasecmp( at[a].resname, presName ) )
			continue;

		for( int x = 0; x < nspecies; x++ )
		{
			if( !strcasecmp( at[a].resname, species[x] ) )
				species_count[x] += 1;
		}

		pres = at[a].res;
		strcpy( presName, at[a].resname );
	}

	int ncomments = 0;

	for( deut = ablock->deuterations; deut; deut = deut->next )
	{
		if( tolower(deut->scheme[0]) == 'z' )
		{
			ncomments++;
			deut->write_fraction = deut->fraction;
			deut->fraction = 1.0;
		}
	}

	char **comments = NULL;

	if( ncomments > 0 )
	{ 
		comments = (char **)malloc( sizeof(char *) * ncomments );
		ncomments = 0;

		for( int spec = 0; spec < nspecies; spec++ )
		{
			int cntr = 0;
			for( deut = ablock->deuterations; deut; deut = deut->next )
			{
				if( !strcasecmp( deut->resName, species[spec]) )
				{
					cntr++;

					if( tolower(deut->scheme[0]) == 'z' )
					{
						int max_len = strlen(deut->resName) + 16 + 32 + strlen("deuterate Z fraction ") + 256;
						char tstring[max_len];
						sprintf(tstring, "deuterate %s Z%d fraction %lf", deut->resName, cntr, deut->write_fraction );
						comments[ncomments] = (char *)malloc( sizeof(char) * ( 1 + strlen(tstring) ) );
						sprintf(comments[ncomments], tstring ); 
						ncomments++;
						deut->write_fraction = deut->fraction;
						deut->fraction = 1.0;
					}
				}
			}
		}
	}


	for( int spec = 0; spec < nspecies; spec++ )
	{
		int start = 0;

		for( int leaflet = 0; leaflet < 2; leaflet++ )
		{
			// assign sequential residues to each scheme.	
	
			// has this lipid been assigned a special circumstance?
		
			int is_z_deut = 0;
	
			for( deut = ablock->deuterations; deut; deut = deut->next )		
			{
				if( deut->scheme[0] == 'z' && !strcasecmp( deut->resName, species[spec] ) )
					is_z_deut = 1;
			}
	
			if( is_z_deut )
			{
				for( int a = 0; a < curNAtoms(); a++ )
				{
					int cntr = 0;
	
					for( deuterate_request * theScheme = ablock->deuterations; theScheme; theScheme = theScheme->next )
					{
						if( strcasecmp( theScheme->resName, at[a].resname) ) continue;
						if( theScheme->scheme[0] == 'z' )
						{
							if( !strcasecmp( theScheme->resName, species[spec] ) )
							{
								cntr++;
								int query = should_we_deuterate( species[spec], theScheme->scheme, at[a].atname );
				
								if( toupper(at[a].atname[0]) == 'H') 
								{	
									if( query == 1 )
										at[a].atname[0] = 'D';	 
									else if( query == 2 )
										at[a].atname[0] = 'X';	 
									else if( query == 3 )
									{
										free(at[a].atname);
										at[a].atname = (char *)malloc( sizeof(char) * 16 );
										sprintf(at[a].atname, "Z%d", cntr );
									}
								}
							}
						}
					}
					pres = at[a].res;
					strcpy( presName, at[a].resname );
				}
			}
			else
			{
				int deuterated_count = 0;
		
				for( deut = ablock->deuterations; deut; deut = deut->next )
				{
					if( strcasecmp( deut->resName, species[spec] ) ) continue;
		
					if( strcasecmp( deut->scheme, "prot" ) )
					{
						deut->num = (int)lround( deut->fraction * (species_count[spec]/2) ); 
						deuterated_count += deut->num;
					}			
				}
		
				int protiated_count = (species_count[spec]/2) - deuterated_count;
		
				while( protiated_count < 0 )
				{
					double max_off = 0;
					deuterate_request *worst_deut = NULL;
					for( deut = ablock->deuterations; deut; deut = deut->next )
					{
						if( strcasecmp( deut->resName, species[spec] ) ) continue;
		
						double off = deut->num - deut->fraction * (species_count[spec]/2);
		
						if( off > max_off )
						{
							max_off = off;
							worst_deut = deut;
						}
					}
		
					if( worst_deut )
					{
						worst_deut->num -= 1;
						protiated_count += 1;
					}
					else
					{
						printf("Deuteration failure.\n");
						exit(1);
					}
				}
		
				pres = -1;
				sprintf( presName, "INIT" );
				deuterate_request *theScheme = NULL;
		
				int num_remaining = (species_count[spec]/2);
		
				int cur_cntr = 0;

				for( int a = 0; a < curNAtoms(); a++ )
				{
					if( at[a].res != pres || strcasecmp( at[a].resname, presName ) )
					{
						theScheme = NULL;
		
						if( !strcasecmp( at[a].resname, species[spec] ) )
						{
							int skip = 0;
							if( leaflet == 1 && cur_cntr < (species_count[spec]/2) )
							{
								skip = 1;
								pres = at[a].res;
								strcpy( presName, at[a].resname );
							}

							cur_cntr++;
			
							if( skip ) continue;			
								
							// pick a deuteration scheme.
								
							double pcum = 0;
			
							for( deut = ablock->deuterations; deut; deut = deut->next )
							{
								if( !strcasecmp( deut->resName, species[spec] ) )
								{
									deut->temp = pcum + (double)deut->num / (double)num_remaining; 
									pcum = deut->temp;
									//printf("NLeft %d deut scheme %s p %lf cum %lf\n",
								//	num_remaining, deut->scheme, (double)deut->num / (double)num_remaining, pcum );
								}
							}		
			
							double rn = (double)rand() / (double)RAND_MAX;
			
							
							for( deut = ablock->deuterations; deut; deut = deut->next )
							{
								if( !strcasecmp( deut->resName, species[spec] ) )
								{
									if( rn <= deut->temp )	
									{
										theScheme = deut;
										break;
									}
								}
							}		
		
							if( theScheme )	
								theScheme->num -= 1;
		
							num_remaining -= 1;
						}
					}
		
					if( !strcasecmp( at[a].resname, species[spec]) && theScheme )
					{
						int query = should_we_deuterate( species[spec], theScheme->scheme, at[a].atname );
			
						if( query == 1 )
							at[a].atname[0] = 'D';	 
						else if( query == 2 )
							at[a].atname[0] = 'X';	 
						else if( query == 3 )
							at[a].atname[0] = 'Z';	 
					}
		
					pres = at[a].res;
					strcpy( presName, at[a].resname );
				}
			}
		}
	}
	printPSF( stdout, at, comments, ncomments ); 
}








