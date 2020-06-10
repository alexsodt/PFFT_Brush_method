#include <stdlib.h>
#include <assert.h>

#include "dg_dcd.h"

DCD *openDCD(char *path) {
    DCD *d = malloc(sizeof(DCD));

    d->hdl=fopen(path,"r");

    if( ! d->hdl ) {
      printf("Couldn't open file '%s'.\n", path );
      exit(1);
    }

    fseek(d->hdl, 8, SEEK_SET);
    assert(1==fread(&(d->nframes), 4, 1, d->hdl));
    fseek(d->hdl, 96, SEEK_SET);
    uint32_t ts;
    assert(1==fread(&ts, 4, 1, d->hdl));
    fseek(d->hdl, 80*((long int) ts) + 8, SEEK_CUR);
    assert(1==fread(&(d->natoms), 4, 1, d->hdl));
    fseek(d->hdl, 4, SEEK_CUR);
    d->offset = ftell(d->hdl);

    return d;
}

void closeDCD(DCD *d) {
    fclose(d->hdl);
    free(d);
}

void goToFrame(DCD *d,uint32_t f) {
    fseek(d->hdl,((long int) (d->offset)) + ((long int) f)*(12*((long int) (d->natoms)) + 80),SEEK_SET);
}

void nextFrame(DCD *d) {
    fseek(d->hdl,(12*((long int) (d->natoms)) + 80),SEEK_CUR);
}

uint32_t getFrame(DCD *d) {
    return (ftell(d->hdl) - ((long int) (d->offset)))/(12*((long int) (d->natoms)) + 80);
}

void getUnitCell(DCD *d, double *uc) {
    fseek(d->hdl, 4, SEEK_CUR);
    assert(1==fread(&(uc[0]), 8, 1, d->hdl));
    fseek(d->hdl, 8, SEEK_CUR);
    assert(1==fread(&(uc[1]), 8, 1, d->hdl));
    fseek(d->hdl, 16, SEEK_CUR);
    assert(1==fread(&(uc[2]), 8, 1, d->hdl));
    fseek(d->hdl, 4, SEEK_CUR);
    fseek(d->hdl, -56, SEEK_CUR);
}

void getCoords(DCD *d, float *xs, float *ys, float *zs) {
    fseek(d->hdl, 56, SEEK_CUR);
    fseek(d->hdl, 4, SEEK_CUR);
    assert((d->natoms)==fread(xs, 4, (d->natoms), d->hdl));
    fseek(d->hdl, 4, SEEK_CUR);
    fseek(d->hdl, 4, SEEK_CUR);
    assert((d->natoms)==fread(ys, 4, (d->natoms), d->hdl));
    fseek(d->hdl, 4, SEEK_CUR);
    fseek(d->hdl, 4, SEEK_CUR);
    assert((d->natoms)==fread(zs, 4, (d->natoms), d->hdl));
    fseek(d->hdl, 4, SEEK_CUR);
    fseek(d->hdl, (-1)*(12*((long int) (d->natoms)) + 80), SEEK_CUR);
}
