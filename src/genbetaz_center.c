#include "dg_dcd.h"
#include "psf_b.h"

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

//SETTINGS//

//histogram size
double histdz = 0.5;
double histmaxz = 60;

//Skip frames
int firstframe = 0;

//END//

int main(int argc, char *argv[]) {

  //open DCD and read basic values
  DCD *dcd;
  dcd = openDCD(argv[1]);
  int frames = dcd->nframes;
  int natom = dcd->natoms;

  //get scattering lengths from file
  double bs[natom];
  psfToBs(argv[2],bs);

  //get residue names from file
  char resnames[natom][9];
  for(int i=0; i<natom; i++)
    for(int j=0; j<9; j++)
      resnames[i][j]='\0';
  psfToResNames(argv[2],resnames);

  //prepare output array
  int numbins = ceil(histmaxz/histdz)*2+1;
  double hist[numbins];
  for (int i=0; i<numbins; i++) {
    hist[i]=0;
  }

  //prepare frame data arrays
  double uc[6];
  float xs[natom];
  float ys[natom];
  float zs[natom];

  //DEBUG
//  frames=1;
  printf("natom: %d frames: %d\n", natom, frames );
  //DEBUG

  double minboxz=0;

//  goToFrame(dcd,firstframe);

  for (int fi=firstframe; fi<frames; fi++) {
    if(fi%1000==0)
      fprintf(stderr,"frame=%d\n",fi);

    //get frame data
    getUnitCell(dcd,uc);
    getCoords(dcd,xs,ys,zs);

    minboxz=(minboxz==0)?uc[2]:fmin(uc[2],minboxz);

    //get membrane center of scattering density
    double testcenter=0;
    double zcenter=0;
    do {
      double slxmemb=0;
      zcenter = testcenter;
      testcenter=0;
      for (int ai=0; ai<natom; ai++) {
        if(strncmp(resnames[ai],"TIP3",4)){
          double z = zs[ai] - zcenter;
          while(z>=0.5*uc[2])
            z-=uc[2];
          while(z<-0.5*uc[2])
            z+=uc[2];

          testcenter+=bs[ai]*bs[ai]*z;
          slxmemb+=bs[ai]*bs[ai];
        }
      }
      testcenter/=slxmemb;
      testcenter+=zcenter;
      while(testcenter>=0.5*uc[2])
        testcenter-=uc[2];
      while(testcenter<-0.5*uc[2])
        testcenter+=uc[2];
      //fprintf(stderr,"testcenter: %lf, zcenter: %lf, %%diff=%lf\n",testcenter,zcenter,fabs(2*(testcenter-zcenter)/(testcenter+zcenter)));
      fprintf(stderr,".");
    } while(fabs(2*(testcenter-zcenter)/(testcenter+zcenter))>0.0001);
    zcenter = testcenter;
    fprintf(stderr,"/");

    //add atoms to bins
    for (int ai=0; ai<natom; ai++) {
      float z = zs[ai]-zcenter;
      while(z>=0.5*uc[2])
        z-=uc[2];
      while(z<-0.5*uc[2])
        z+=uc[2];
      int index = (floor((z/histdz)+0.5)+ceil(histmaxz/histdz));

      hist[index]+=bs[ai]/(uc[0]*uc[1]*histdz);
    }

    nextFrame(dcd);
  }

  //-36.6 to +36.6
  int minindex = (floor(((-0.5*minboxz)/histdz)+0.5)+ceil(histmaxz/histdz)) + 1;
  //(floor(((0.5*minboxz)/histdz)+0.5)+ceil(histmaxz/histdz)) - 1;
  //765;
  int maxindex = (floor(((0.5*minboxz)/histdz)+0.5)+ceil(histmaxz/histdz)) - 1;

  //normalize and print NSLD profile
  for (int i=minindex; i<=maxindex; i++)
    printf("%.1f\t%e\n",i*histdz-histmaxz,hist[i]/(frames-firstframe));

  closeDCD(dcd);

  return (0);
}
