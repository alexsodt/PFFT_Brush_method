#include <stdio.h>
#include <assert.h>
#include <stdlib.h>

#include "bfile.h"

void getStructure(char *path,int natom,float *bs) {
  FILE *bfile;
  bfile=fopen(path,"r");
  if( !bfile )
  {
	printf("Couldn't open file '%s'.\n", path );
	exit(-1);
  }
  for(int ai=0; ai<natom; ai++)
    assert(1==fscanf (bfile, "%f", &bs[ai]));
  fclose(bfile);
}
