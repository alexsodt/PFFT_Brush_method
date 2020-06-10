// by alex sodt
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

int readNInts( char *buffer, int *vals, int nvalues )
{
        int s = 0;

        int slen = strlen(buffer);
	int nread = 0;

        for( int x = 0; x < nvalues; x++ )
        {
                while( s < slen && (buffer[s] == ' ' || buffer[s] == '\t') )// || isalpha(buffer[s])) )
                        s++;

                int nr = sscanf( buffer + s, "%d", vals+x );

		if( nr != 1 )
			return nread;

		nread++;

                while( s < slen && !(buffer[s] == ' ' || buffer[s] == '\t'))// || isalpha(buffer[s])) )
                        s++;
        }

	return nread;
}


int readNDoubles( char *buffer, double *vals, int nvalues )
{
        int s = 0;

        int slen = strlen(buffer);
	int nread = 0;

        for( int x = 0; x < nvalues; x++ )
        {
                while( s < slen && (buffer[s] == ' ' || buffer[s] == '\t') )// || isalpha(buffer[s])) )
                        s++;

                int nr = sscanf( buffer + s, "%lf", vals+x );

		if( nr != 1 )
			return nread;

		nread++;

                while( s < slen && !(buffer[s] == ' ' || buffer[s] == '\t'))// || isalpha(buffer[s])) )
                        s++;
        }

	return nread;
}



void getLine( FILE *theFile, char *theBuffer )
{
        int i = 0;

        while( !feof(theFile) )
        {
                char tc = fgetc(theFile);

                if( tc != '\n' && i < 99999 )
                {
                        theBuffer[i++] = tc;
                }
                else if( tc != '\n' && i >= 99999 )
		{
		}
		else
                        break;
        }

        theBuffer[i] = '\0';
}

void print5( int val, char *str )
{
        if( val < 10 )
                sprintf(str, "0000%d", val );
        else if( val < 100 )
                sprintf(str, "000%d", val );
        else if( val < 1000 )
                sprintf(str, "00%d", val );
        else if( val < 10000 )
                sprintf(str, "0%d", val );
        else
                sprintf(str, "%d", val );
}

