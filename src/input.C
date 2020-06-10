#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "input.h"
#include "util.h"

void loadBlockDefaults( struct parameter_block *block )
{
	// defaults.
	block->scatterFile = NULL;
	block->cutoff_mode = 0;
	block->weight = 1.0;
//	block->cutoff = 30.0;
	block->cutoff = -1;
	block->smooth = 3.0;
	block->qvals = NULL;
	block->print_fft = 0;
	block->apply_strain = 0;
	block->lateral_scale = 1.0;
	block->water_average_thickness = 2.0;
	block->intra = 0;

	block->name = (char *)malloc( sizeof(char) * 256 );
	sprintf(block->name, "default");

	block->X_scattering_length = PFFT_HUGE*1.2; // 0.531731; // 87% Deut 
	block->brFile = NULL; 
	block->expFit = (char *)malloc( sizeof(char) * 256 );
	sprintf( block->expFit, "exp.fit");
	block->bzFile = (char *)malloc( sizeof(char) * 256 );
	sprintf( block->bzFile, "beta.z");
	block->activateFit = 0;
	block->instructions = NULL;
	block->deuterations = NULL;
	block->forceSolventSub = 0;
	block->solvent_sub=0;
}

void getInput( const char *fileName, struct parameter_block *block) 
{
	FILE *theFile = fopen(fileName,"r");

	if( !theFile )
	{
		printf("Couldn't open input file '%s'\n", fileName );
		exit(1);
	}

	char *buffer = (char *)malloc( sizeof(char) * 100000 ); // getLine stops at character 99999 so this does avoid buffer overruns.
	char *word1 = (char *)malloc( sizeof(char) * 4096 );
	char *word2 = (char *)malloc( sizeof(char) * 4096 );

	loadBlockDefaults( block );

	int name_not_default = 0;
	int betaz_not_default = 0;

	while( !feof(theFile) )
	{
		getLine( theFile, buffer );

		char *ptr = buffer;
		while( *ptr == ' ' || *ptr == '\t' ) ptr += 1;
		if( *ptr == '#' ) continue;

		if( strlen(buffer) > 4095 )
		{
			printf("Line %s too long.\n", buffer );
			exit(1);
		}

		if( feof(theFile) ) break;

		int nr = sscanf( buffer, "%s %s", word1, word2 );

		if( nr != 2 )
		{
			printf("Could not interpret input line '%s'.\n", buffer );
			exit(1);
		}

		if( !strcasecmp( word1, "HIST" ) )
		{
			block->brFile = (char *)malloc( sizeof(char) * (1+strlen(word2) ) );
			strcpy( block->brFile, word2 );
			
		}
		else if( !strcasecmp( word1, "SCATTER" ) )
		{
			block->scatterFile = (char *)malloc( sizeof(char) * (1+strlen(word2) ) );
			strcpy( block->scatterFile, word2 );
			
		}
		else if( !strcasecmp( word1, "intra" ) )
		{
			if( !strcasecmp( word2, "yes") || !strcasecmp( word2, "true") || !strcasecmp( word2, "1") ) 
				block->intra = 1;
			else if( !strcasecmp( word2, "no") || !strcasecmp( word2, "false") || !strcasecmp( word2, "0") ) 
				block->intra = 0;
			else
			{
				printf("Could not interpret input line '%s'.\n", buffer );
				exit(1);
			}
		}
		else if( !strcasecmp( word1, "EXPFIT" ) )
		{
			free(block->expFit);
			block->expFit = (char *)malloc( sizeof(char) * (1+strlen(word2) ) );
			strcpy( block->expFit, word2 );
			block->activateFit = 1;
		}
		else if( !strcasecmp( word1, "QVALS" ) )
		{
			free(block->qvals);
			block->qvals = (char *)malloc( sizeof(char) * (1+strlen(word2) ) );
			strcpy( block->qvals, word2 );
		}
		else if( !strcasecmp( word1, "STRAIN" ) )
		{
			block->apply_strain = atof(word2);
		}
		else if( !strcasecmp( word1, "LATERAL_SCALE" ) )
		{
			block->lateral_scale = atof(word2);
		}
		else if( !strcasecmp( word1, "SOLVENT_SUB" ) )
		{
			block->solvent_sub = atof(word2);
			block->forceSolventSub = 1;
		}
		else if( !strcasecmp( word1, "NAME" ) )
		{
			name_not_default = 1;
			free(block->name);
			block->name = (char *)malloc( sizeof(char) * (1+strlen(word2) ) );
			strcpy( block->name, word2 );
		}
		else if( !strcasecmp( word1, "BETAZ" ) )
		{
			betaz_not_default = 1;
			free(block->bzFile);
			block->bzFile = (char *)malloc( sizeof(char) * (1+strlen(word2) ) );
			strcpy( block->bzFile, word2 );
		}
		else if( !strcasecmp( word1, "DEUTERATE" ) )
		{
			double fr;
			char garbage[256];
			int nr = sscanf( buffer, "%s %s %s %lf", garbage, word1, word2, &fr );

			if( nr != 4 )
			{
				printf("Could not interpret 'DEUTERATE' command:\n%s\n", buffer );
				printf("%d %s %s %lf\n", nr, word1, word2, fr );
				exit(1);
			} 

			deuterate_request *deut = (deuterate_request *)malloc( sizeof( deuterate_request) );

			deut->next = block->deuterations;
			block->deuterations = deut;

			deut->resName = (char *)malloc( sizeof(char) * (1 + strlen(word1) ) );
			strcpy( deut->resName, word1 ); 	

			deut->scheme = (char *)malloc( sizeof(char) * (1 + strlen(word2) ) );
			strcpy( deut->scheme, word2 );

			deut->fraction = fr;
		}
		else if( !strcasecmp( word1, "MODIFY" ) ) 
		{
			double mod, R1, R2;
			int nr = sscanf( buffer, "MODIFY %s %s", word1, word2 ); //%lf %lf %lf", word1, word2, &mod, &R1, &R2  );

			if( nr != 2 )
			{
				printf("Could not interpret 'MODIFY' command:\n%s\n", buffer );
				exit(1);
			}

			char *read = buffer;
			// to the end of MODIFY
			while( *read && *read != ' ' && *read != '\t' ) read++;
			while( *read && (*read == ' ' || *read == '\t') ) read++;
			// to the end of word1
			while( *read && *read != ' ' && *read != '\t' ) read++;
			while( *read && (*read == ' ' || *read == '\t') ) read++;
			// to the end of word2
			while( *read && *read != ' ' && *read != '\t' ) read++;
			while( *read && (*read == ' ' || *read == '\t') ) read++;
			
			if( !*read) { printf("Could not interpret line '%s'.\n", buffer); exit(1); }
			pfft_instruction *inst = (pfft_instruction *)malloc( sizeof(pfft_instruction) );
	
			if( !strncasecmp( read, "X", 1  ) )
			{
				block->activateFit = 1;
				inst->fit_parm[0] = 1;
				inst->mod = atof(read+1);
			}
			else
			{
				inst->fit_parm[0] = 0;
				inst->mod = atof(read);
			}
			while( *read && *read != ' ' && *read != '\t' ) read++;
			while( *read && (*read == ' ' || *read == '\t') ) read++;
			if( !*read) { printf("Could not interpret line '%s'.\n", buffer); exit(1); }
			if( !strncasecmp( read, "X", 1  ) )
			{
				block->activateFit = 1;
				inst->fit_parm[1] = 1;
				inst->R1 = atof(read+1);
			}
			else
			{
				inst->fit_parm[1] = 0;
				inst->R1 = atof(read);
			}
			
			while( *read && *read != ' ' && *read != '\t' ) read++;
			while( *read && (*read == ' ' || *read == '\t') ) read++;

			if( !*read )
			{
				inst->R2 = inst->R1;
			inst->R1 = -1;
				inst->fit_parm[2] = inst->fit_parm[1];
				inst->fit_parm[1] = 0;
			}
			else
			{
				if( !strncasecmp( read, "X", 1  ) )
				{
					block->activateFit = 1;
					inst->fit_parm[2] = 1;
					inst->R2 = atof(read+1);
				}
				else
				{
					inst->fit_parm[2] = 0;
					inst->R2 = atof(read);
				}

			}				
			
			inst->next = block->instructions;
			block->instructions = inst;
		
			inst->forRes = (char *)malloc( sizeof(char) * (1 + strlen(word1) )  );			
			inst->depleteRes = (char *)malloc( sizeof(char) * (1 + strlen(word2) )  );		
	
			strcpy( inst->forRes, word1 );	
			strcpy( inst->depleteRes, word2 );	
			
			if( R1 >= 0 )
			{
				printf("Modifying residue %s for residue %s by %lf between %lf and %lf.\n", inst->depleteRes, inst->forRes,
					inst->mod, inst->R1, inst->R2 );
			}
			else
			{
				printf("Depleting residue %s for residue %s by %lf between the particle cutoff and %lf.\n", inst->depleteRes, inst->forRes,
					inst->mod, inst->R2 );
			}
		}
		else if( !strcasecmp( word1, "WEIGHT" ) )
			block->weight = atof(word2);
		else if( !strcasecmp( word1, "WATER_AVERAGE_THICKNESS" ) )
			block->water_average_thickness = atof(word2);
		else if( !strcasecmp( word1, "X_scattering_length" ) )
			block->X_scattering_length = atof(word2);
		else if( !strcasecmp( word1, "CUTOFF" ) )
			block->cutoff = atof(word2);
		else if( !strcasecmp( word1, "fft") ) 
		{
			if( !strcasecmp( word2, "yes") || !strcasecmp( word2, "true") || !strcasecmp( word2, "1") ) 
				block->print_fft = 1;
			else if( !strcasecmp( word2, "no") || !strcasecmp( word2, "false") || !strcasecmp( word2, "0") ) 
				block->print_fft = 0;
			else
			{
				printf("Could not interpret input line '%s'.\n", buffer );
				exit(1);
			}
		}
		else if( !strcasecmp( word1, "areaCompensate") ) 
		{
			if( !strcasecmp( word2, "yes") || !strcasecmp( word2, "true") || !strcasecmp( word2, "1") ) 
				block->areaCompensate = 1;
			else if( !strcasecmp( word2, "no") || !strcasecmp( word2, "false") || !strcasecmp( word2, "0") ) 
				block->areaCompensate = 0;
			else
			{
				printf("Could not interpret input line '%s'.\n", buffer );
				exit(1);
			}
		}
		else if( !strcasecmp( word1, "SMOOTH" ) )
			block->smooth = atof(word2);
		else
		{
			printf("Could not interpret input line '%s'.\n", buffer );
			exit(1);
		}
	}

	if( name_not_default && !betaz_not_default )
	{
		free(block->bzFile);
		int len = strlen(block->name) + 4 ;
		block->bzFile = (char *)malloc( sizeof(char) * len );
		sprintf(block->bzFile, "%s.bz", block->name );
	}


	if( !block->brFile && 0)
	{
		printf("Parameter file '%s' does not contain a histogram.\n", fileName );
		exit(1);
	}

	if( !block->brFile && block->post_process )
	{
		block->cutoff = 0.0001;
		block->smooth = 0.0001;
	}
}



