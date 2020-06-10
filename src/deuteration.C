#include "deuteration.h"
#include <stdio.h>
#include <string.h>
int deuterationAvailable( const char *resName, const char *scheme )
{
	if( !strcasecmp( resName, "TIP3" ) )
	{
		if( !strcasecmp( scheme, "d2o" ) )
			return 1;
		if( !strcasecmp( scheme, "z2o" ) )
			return 1;
	}
	if( !strcasecmp( resName, "CHL1" ) )
	{
		if( !strcasecmp( scheme, "prot" ) )
			return 1;
		if( !strcasecmp( scheme, "d41" ) )
			return 1;
		if( !strcasecmp( scheme, "yeast" ) )
			return 1;
	}
	else if( !strcasecmp( resName, "DPPC" ) )
	{
		if( !strcasecmp( scheme, "d13" ) )
			return 1;
		if( !strcasecmp( scheme, "d62" ) )
			return 1;
		if( !strcasecmp( scheme, "z62" ) )
			return 1;
		if( !strcasecmp( scheme, "z13" ) )
			return 1;
		if( !strcasecmp( scheme, "prot") )
			return 1; 
	}
	else if( !strcasecmp( resName, "DPPC" ) )
	{
		if( !strcasecmp( scheme, "d66" ) )
			return 1;
		if( !strcasecmp( scheme, "prot") )
			return 1; 
	}
	else if( !strcasecmp( resName, "PSM" ) )
	{
		if( !strcasecmp( scheme, "d31" ) )
			return 1;
	}
	else if( !strcasecmp( resName, "POPE" ) )
	{
		if( !strcasecmp( scheme, "d31" ) )
			return 1;
	}
	else if( !strcasecmp( resName, "POPS" ) )
	{
		if( !strcasecmp( scheme, "d31" ) )
			return 1;
	}
	else if( !strcasecmp( resName, "POPC" ) )
	{
		if( !strcasecmp( scheme, "d31" ) )
			return 1;
		if( !strcasecmp( scheme, "d82" ) )
			return 1;
	}

	return 0;
}

int should_we_deuterate( const char *resName, const char *scheme, const char *atName )
{
	if( atName[0] != 'H' ) return 0;

	/* Is this atom part of the deuteration scheme? */

	if( !strcasecmp( resName, "TIP3" ) )
	{
		if( !strcasecmp( scheme, "prot" ) )
			return 0;
		else if( !strcasecmp( scheme, "d2o") )
			return 1;
		else if( !strcasecmp( scheme, "z2o") )
			return 3;
	}
	else if( !strcasecmp( resName, "CHL1" ) )
	{	
		if( !strcasecmp( scheme, "prot" ) )
			return 0;
		else if( !strcasecmp( scheme, "d41") ) 
		{
			if( !strcasecmp( atName, "H3'") )
				return 1;
			else
				return 2; // special case
			
		}
		else if( !strcasecmp( scheme, "yeast") ) 
		{
			if( !strcasecmp( atName, "H3'") )
				return 0;
			if( !strcasecmp( atName, "H19A") )
				return 0;
			if( !strcasecmp( atName, "H19B") )
				return 0;
			if( !strcasecmp( atName, "H18A") )
				return 0;
			if( !strcasecmp( atName, "H18B") )
				return 0;
			if( !strcasecmp( atName, "H21A") )
				return 0;
			if( !strcasecmp( atName, "H21B") )
				return 0;
//			if( !strcasecmp( atName, "H26A") )
//				return 0;
			return 1; // special case
		}
	}
	else if( !strcasecmp( resName, "DPPC" ) )
	{
		char last_char = atName[strlen(atName)-1];
		if( !strcasecmp( scheme, "prot" ) )
			return 0;
		else if( !strcasecmp( scheme, "d13") ) 
		{
			if( (last_char == 'A' || last_char == 'B' || last_char == 'C' ) && 
			    (!strncasecmp( atName, "H13",3) ||
			    !strncasecmp( atName, "H14",3) ||
			    !strncasecmp( atName, "H15",3)) )
				return 1;
		}
		else if( !strcasecmp( scheme, "z13") ) 
		{
			if( (last_char == 'A' || last_char == 'B' || last_char == 'C' ) && 
			    (!strncasecmp( atName, "H13",3) ||
			    !strncasecmp( atName, "H14",3) ||
			    !strncasecmp( atName, "H15",3)) )
				return 3;
		}
		else if( !strcasecmp( scheme, "d62") ) 
		{
			if( last_char == 'X' || last_char == 'Y' || last_char == 'Z' ||
			    last_char == 'R' || last_char == 'S' || last_char == 'T' )
			{
				if( strlen(atName) == 2  ) // HS, HX, HY: glycerol 
					return 0;
				return 1;
			}
		}
		else if( !strcasecmp( scheme, "z62") ) 
		{
			if( last_char == 'X' || last_char == 'Y' || last_char == 'Z' ||
			    last_char == 'R' || last_char == 'S' || last_char == 'T' )
			{
				if( strlen(atName) == 2  ) // HS, HX, HY: glycerol 
					return 0;
				return 3;
			}
		}
	}
	else if( !strcasecmp( resName, "DOPC" ) )
	{
		char last_char = atName[strlen(atName)-1];
		if( !strcasecmp( scheme, "prot" ) )
			return 0;
		else if( !strcasecmp( scheme, "d13") ) 
		{
			if( (last_char == 'A' || last_char == 'B' || last_char == 'C' ) && 
			    (!strncasecmp( atName, "H13",3) ||
			    !strncasecmp( atName, "H14",3) ||
			    !strncasecmp( atName, "H15",3)) )
				return 1;
		}
		else if( !strcasecmp( scheme, "d66") ) 
		{
			if( last_char == 'X' || last_char == 'Y' || last_char == 'Z' ||
			    last_char == 'R' || last_char == 'S' || last_char == 'T' )
			{
				if( strlen(atName) == 2  ) // HS, HX, HY: glycerol 
					return 0;
				return 1;
			}
		}
	}
	else if( !strcasecmp( resName, "DOPS" ) )
	{
		char last_char = atName[strlen(atName)-1];
		if( !strcasecmp( scheme, "prot" ) )
			return 0;
		else if( !strcasecmp( scheme, "d62") ) 
		{
			if( last_char == 'X' || last_char == 'Y' || last_char == 'Z' ||
			    last_char == 'R' || last_char == 'S' || last_char == 'T' )
			{
				if( strlen(atName) == 2  ) // HS, HX, HY: glycerol 
					return 0;
				return 1;
			}
		}
	}
	else if( !strcasecmp( resName, "POPS" ) )
	{
		char last_char = atName[strlen(atName)-1];
		if( !strcasecmp( scheme, "prot" ) )
			return 0;
		else if( !strcasecmp( scheme, "d31") ) 
		{
			if( last_char == 'X' || last_char == 'Y' || last_char == 'Z' ) 
			{
				if( strlen(atName) == 2  ) // HS, HX, HY: glycerol 
					return 0;
				return 1;
			}
		}
	}
	else if( !strcasecmp( resName, "POPE" ) )
	{
		char last_char = atName[strlen(atName)-1];
		if( !strcasecmp( scheme, "prot" ) )
			return 0;
		else if( !strcasecmp( scheme, "d31") ) 
		{
			if( last_char == 'X' || last_char == 'Y' || last_char == 'Z' ) 
			{
				if( strlen(atName) == 2  ) // HS, HX, HY: glycerol 
					return 0;
				return 1;
			}
		}
	}
	else if( !strcasecmp( resName, "POPC" ) )
	{
		char last_char = atName[strlen(atName)-1];
		if( !strcasecmp( scheme, "prot" ) )
			return 0;
		else if( !strcasecmp( scheme, "d82") ) 
		{
			return 1;
		}
		else if( !strcasecmp( scheme, "d31") ) 
		{
			if( last_char == 'X' || last_char == 'Y' || last_char == 'Z' ) 
			{
				if( strlen(atName) == 2  ) // HS, HX, HY: glycerol 
					return 0;
				return 1;
			}
		}
	}
	else if( !strcasecmp( resName, "PSM" ) )
	{
		char last_char = atName[strlen(atName)-1];
		if( !strcasecmp( scheme, "prot" ) )
			return 0;
		if( !strcasecmp( scheme, "d31" ) )
		{
			if( last_char == 'F' || last_char == 'G' || last_char == 'H' )
			{
				if( strcasecmp(atName,"HNF") ) 	
					return 1;
			}
		}
	}
	return 0;
}
