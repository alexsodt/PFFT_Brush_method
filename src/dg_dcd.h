#ifndef DCD_H_
#define DCD_H_

#include <stdio.h>
#include <stdint.h>

typedef struct {
    int is_pdb;
    FILE *hdl;
    uint32_t  nframes;
    uint32_t  natoms;
    long int  offset;
} DCD;

DCD *openDCD(char *);

void closeDCD(DCD *);

void goToFrame(DCD *,uint32_t);

void nextFrame(DCD *);

uint32_t getFrame(DCD *);

void getUnitCell(DCD *, double *);

void getCoords(DCD *, float *, float *, float *);

#endif
