#include "dg_dcd.h"
#include "bfile.h"

#include <math.h>
#include <complex.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#define SOLVENT_SUBTRACTION

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

void inferCellAlignment(int,float*,float*,float*,int*);

//SETTINGS//

////q stuff//
int numqs = 250;
bool logspacing = false;
float firstq = 0.002;
float lastq = 0.5;

//Solvent Scattering//
double betaSolv=6.545E-6;

//Time average//
int firstFrame = 0;
int frameModulus = 1;

//END//

int main(int argc, char *argv[]) {

  if( argc < 3 )
	{
		printf("Syntax: fast_dg dcd scattering_lengths.bs [betaSolv] [qvals] [lateral_scale]\n");
		return 0;
	}

  //Prepare list of q-values
  //float qs[268] = {0.001587,0.001875,0.002164,0.002452,0.002741,0.003029,0.003318,0.003606,0.003895,0.004184,0.004241,0.004472,0.004645,0.004761,0.005049,0.005338,0.005453,0.005626,0.005857,0.005915,0.006203,0.006261,0.006665,0.007069,0.007473,0.007876,0.00828,0.008684,0.009088,0.009492,0.009896,0.0103,0.0107,0.01111,0.0113,0.01151,0.01192,0.01232,0.01263,0.01272,0.01313,0.01353,0.01393,0.01396,0.01434,0.01474,0.01515,0.01529,0.01555,0.01595,0.01636,0.01662,0.01676,0.01717,0.01757,0.01795,0.01797,0.01838,0.01878,0.01918,0.01928,0.01959,0.01999,0.0204,0.02061,0.0208,0.0212,0.02161,0.02194,0.02201,0.02241,0.02282,0.02322,0.02327,0.02363,0.02403,0.02443,0.0246,0.02484,0.02524,0.02564,0.02593,0.02605,0.02645,0.02685,0.02726,0.02766,0.02807,0.02847,0.02859,0.02887,0.02928,0.02968,0.02991,0.03008,0.03049,0.03124,0.03257,0.0339,0.03399,0.03523,0.03656,0.03788,0.03798,0.03921,0.04054,0.04187,0.04197,0.0432,0.04452,0.04585,0.04596,0.04718,0.0485,0.04983,0.04996,0.05116,0.05248,0.05381,0.05394,0.05514,0.05646,0.05779,0.05793,0.05911,0.06044,0.06176,0.06192,0.06309,0.06441,0.06573,0.0659,0.06706,0.06838,0.06971,0.06988,0.07103,0.07235,0.07367,0.07386,0.075,0.07632,0.07764,0.07783,0.07896,0.08028,0.08161,0.08181,0.08293,0.08425,0.08557,0.08578,0.08689,0.08821,0.08952,0.08975,0.09084,0.09216,0.09348,0.09371,0.0948,0.09612,0.09743,0.09768,0.09875,0.1001,0.1014,0.1016,0.1027,0.104,0.1053,0.1056,0.1066,0.108,0.1095,0.1135,0.1174,0.1214,0.1253,0.1293,0.1332,0.1371,0.141,0.1449,0.1489,0.1528,0.1567,0.1606,0.1645,0.1683,0.1722,0.1761,0.18,0.1838,0.1877,0.1915,0.1954,0.1992,0.2031,0.2069,0.2107,0.2145,0.2184,0.2222,0.226,0.2297,0.2335,0.2373,0.2411,0.2448,0.2486,0.2523,0.2561,0.2598,0.2635,0.2673,0.271,0.2747,0.2784,0.282,0.2857,0.2894,0.293,0.2967,0.3003,0.304,0.3076,0.3112,0.3148,0.3184,0.322,0.3256,0.3292,0.3327,0.3363,0.3398,0.3434,0.3469,0.3504,0.3539,0.3574,0.3609,0.3644,0.3679,0.3713,0.3748,0.3782,0.3817,0.3851,0.3885,0.3919,0.3953,0.3987,0.4021,0.4054,0.4088,0.4121,0.4155,0.4188,0.4221,0.4254,0.4287,0.432,0.4352,0.4385,0.4417,0.445,0.4482};
//  float qs[190] = {0.014395,0.0179898,0.0215845,0.0251792,0.0287739,0.0323686,0.0359634,0.0395581,0.0431528,0.0467475,0.0503422,0.053937,0.0575317,0.0611264,0.0647211,0.0683158,0.0719106,0.0755053,0.0791,0.0826947,0.0862894,0.0898842,0.0934789,0.0970736,0.100668,0.104263,0.107858,0.111452,0.115047,0.118642,0.122237,0.125831,0.129426,0.133021,0.136616,0.14021,0.143805,0.1474,0.150994,0.154589,0.158184,0.161779,0.165373,0.168968,0.172563,0.176157,0.179752,0.183347,0.186942,0.190536,0.194131,0.197726,0.20132,0.204915,0.20851,0.212105,0.215699,0.219294,0.222889,0.226484,0.230078,0.233673,0.237268,0.240862,0.244457,0.248052,0.251647,0.255241,0.258836,0.262431,0.266025,0.26962,0.273215,0.27681,0.280404,0.283999,0.287594,0.291188,0.294783,0.298378,0.301973,0.305567,0.309162,0.312757,0.316352,0.319946,0.323541,0.327136,0.33073,0.334325,0.33792,0.341515,0.345109,0.348704,0.352299,0.355893,0.359488,0.363083,0.366678,0.370272,0.373867,0.377462,0.381056,0.384651,0.388246,0.391841,0.395435,0.39903,0.402625,0.40622,0.409814,0.413409,0.417004,0.420598,0.424193,0.427788,0.431383,0.434977,0.438572,0.442167,0.445761,0.449356,0.452951,0.456546,0.46014,0.463735,0.46733,0.470924,0.474519,0.478114,0.481709,0.485303,0.488898,0.492493,0.496088,0.499682,0.503277,0.506872,0.510466,0.514061,0.517656,0.521251,0.524845,0.52844,0.532035,0.535629,0.539224,0.542819,0.546414,0.550008,0.553603,0.557198,0.560792,0.564387,0.567982,0.571577,0.575171,0.578766,0.582361,0.585956,0.58955,0.593145,0.59674,0.600334,0.603929,0.607524,0.611119,0.614713,0.618308,0.621903,0.625497,0.629092,0.632687,0.636282,0.639876,0.643471,0.647066,0.65066,0.654255,0.65785,0.661445,0.665039,0.668634,0.672229,0.675824,0.679418,0.683013,0.686608,0.690202,0.693797};
//  float qs[257] = {0.00325,0.00375,0.00425,0.00475,0.00525,0.00575,0.00625,0.00675,0.00725,0.00775,0.00825,0.00875,0.00925,0.00975,0.01025,0.01075,0.01125,0.01175,0.01225,0.01275,0.01325,0.01375,0.01425,0.01475,0.01525,0.01575,0.01625,0.01675,0.01725,0.01775,0.01825,0.01875,0.01925,0.01975,0.02025,0.02075,0.02125,0.02175,0.02225,0.02275,0.02325,0.02375,0.02425,0.02475,0.02525,0.02575,0.02625,0.02675,0.02725,0.02775,0.02825,0.02875,0.02925,0.02975,0.03025,0.03075,0.03125,0.03175,0.03225,0.03275,0.03325,0.03375,0.03425,0.03475,0.03525,0.03575,0.03625,0.03675,0.03725,0.0375,0.03775,0.03825,0.03875,0.03925,0.03975,0.04025,0.04075,0.04125,0.04175,0.04225,0.0425,0.04275,0.04325,0.04375,0.04425,0.04475,0.04525,0.04575,0.04625,0.04675,0.04725,0.0475,0.04775,0.04825,0.04875,0.04925,0.04975,0.0525,0.0575,0.0625,0.0675,0.0725,0.0775,0.0825,0.0875,0.0925,0.0975,0.1025,0.1075,0.1125,0.1175,0.1225,0.1275,0.1325,0.1375,0.1425,0.1475,0.1525,0.1575,0.1625,0.1675,0.1725,0.1775,0.1825,0.1875,0.1925,0.1975,0.2025,0.2075,0.2125,0.2175,0.2225,0.2275,0.2325,0.2375,0.2425,0.2475,0.2525,0.2575,0.2625,0.2675,0.2725,0.2775,0.2825,0.2875,0.2925,0.2975,0.3025,0.3075,0.3125,0.3175,0.3225,0.3275,0.3325,0.3375,0.3425,0.3475,0.3525,0.3575,0.3625,0.3675,0.3725,0.3775,0.3825,0.3875,0.3925,0.3975,0.4025,0.4075,0.4125,0.4175,0.4225,0.4275,0.4325,0.4375,0.4425,0.4475,0.4525,0.4575,0.4625,0.4675,0.4725,0.4775,0.4825,0.4875,0.4925,0.4975,0.5025,0.5075,0.5125,0.5175,0.5225,0.5275,0.5325,0.5375,0.5425,0.5475,0.5525,0.5575,0.5625,0.5675,0.5725,0.5775,0.5825,0.5875,0.5925,0.5975,0.6025,0.6075,0.6125,0.6175,0.6225,0.6275,0.6325,0.6375,0.6425,0.6475,0.6525,0.6575,0.6625,0.6675,0.6725,0.6775,0.6825,0.6875,0.6925,0.6975,0.7025,0.7075,0.7125,0.7175,0.7225,0.7275,0.7325,0.7375,0.7425,0.7475,0.7525,0.7575,0.7625,0.7675,0.7725,0.7775,0.7825,0.7875,0.7925,0.7975,0.8025,0.8075,0.8125,0.8175,0.8225,0.8275,0.8325,0.8375,0.8425,0.8475};
  
  float *qs = NULL;

  if( argc > 4 )
	{
		FILE *qfile = fopen(argv[4],"r");
		if( !qfile) 
		{
			printf("Couldn't open q-file '%s'.\n", argv[4] );
			exit(1);
		}

		char *buffer = (char *)malloc( sizeof(char) * 40000 );

		int nq = 0;
		for( int pass = 0; pass < 2; pass++ )
		{
			rewind(qfile);		
			while( !feof(qfile) ) 
			{
				getLine(qfile, buffer );
				if( feof(qfile) ) break;
				double q = 0;
				int nr = sscanf( buffer, "%lf", &q );

				if( nr == 1 )	
				{
					if( pass == 1 ) 
						qs[nq] = q;
					nq++;
				}
			}
  	
			numqs = nq;
			if( pass == 0 )
				qs = (float *)malloc( sizeof(float) * numqs );
			nq = 0;
		}
		free(buffer);
	} 
	else
	{
  qs = (float *)malloc( sizeof(float) * numqs );

  if (logspacing) { //logarithmic spacing
    float qf = pow(lastq/firstq,1.0/(((float) numqs)-1.0));
    for (int qi=0; qi<numqs; qi++)
      qs[qi]=firstq*pow(qf,qi);
  } else { //linear spacing
    float qd = (lastq-firstq)/(((float) numqs)-1.0);
    for (int qi=0; qi<numqs; qi++)
      qs[qi]=firstq + qd*qi;
  }
	}

	printf("numqs: %d q[0]: %e\n", numqs, qs[0] );

  double lateral_scale = 1.0;

  if( argc > 5 )
	lateral_scale = atof(argv[5]);

  // compute average area so we can divide by it in the end -- alex
  double average_area = 0;
  double naverage_area = 0;

  //open DCD and read basic values
  DCD *dcd;
  dcd = openDCD(argv[1]);

  if( !dcd )
	{
		printf("Couldn't open dcd file '%s'.\n", argv[1] );
		return -1;
	}

  int frames = dcd->nframes;
  int natom = dcd->natoms;

  //get scattering lengths from file
  float bs[natom];

  getStructure(argv[2],natom,bs);

  if( argc > 3 )
		betaSolv = atof(argv[3]);	

  //prepare output array for each frame
  int numFrames = 0;
  for (int fi=firstFrame; fi<frames; fi+=frameModulus)
    numFrames++;
  double * fis = malloc(numqs * numFrames * sizeof(double));
  double * fis_lutz = malloc(numqs * numFrames * sizeof(double));
  for (int qi=0; qi<numqs; qi++)
    for (int fi=0; fi<numFrames; fi++)
	{
      fis[fi*numqs + qi]=0;
      fis_lutz[fi*numqs + qi]=0;
	}
  //prepare output means array
  double is[numqs];
  for(int qi=0; qi<numqs; qi++)
    is[qi]=0;
  
  double is_lutz[numqs];
  for(int qi=0; qi<numqs; qi++)
    is_lutz[qi]=0;

  //prepare output sigmas array
  double sis[numqs];
  for(int qi=0; qi<numqs; qi++)
    sis[qi]=0;
  for( int a = 0; a < natom; a++ )
	bs[a] *= lateral_scale*lateral_scale;

  //prepare frame data arrays
  double uc[6];
  float xs[natom];
  float ys[natom];
  float zs[natom];

  //Get cell offset
  goToFrame(dcd,0);
  getCoords(dcd,xs,ys,zs);
  int cv[3]; //vector pointing toward center of cell
  inferCellAlignment(natom,xs,ys,zs,cv);

  int framecount=0;
  for (int fi=firstFrame; fi<frames; fi+=frameModulus) {
    goToFrame(dcd,fi);
    fprintf(stderr,"Frame %d:\n",fi);
      

    //get frame data
    getUnitCell(dcd,uc);

    uc[0] *= lateral_scale;
    uc[1] *= lateral_scale;

    average_area += uc[0] * uc[1];
    naverage_area += 1;
    getCoords(dcd,xs,ys,zs);

    for( int a = 0; a < natom; a++ )
	{	
		xs[a] *= lateral_scale;
		ys[a] *= lateral_scale;
	}
   
 
    //for (int ai=0; ai<natom; ai++) 
    //	zs[ai] -= uc[2]/2;
    #pragma omp parallel for schedule(dynamic)
    for (int qi=numqs-1; qi>=0; qi--) {
	printf("q: %lf\n", qs[qi] );
      double intensity=0;
      double gridX = 2*M_PI/uc[0];
      double gridY = 2*M_PI/uc[1];
      int minX = floor(qs[qi]/gridX);
      int minY = floor(qs[qi]/gridY);
      int bufSize = 2*(2*minX + 1)*(2*minY + 1);
      double buffer[3*bufSize];
      int gi=0;
      for (int gx=-minX; gx<=minX; gx++) {
        for (int gy=-minY; gy<=minY; gy++) {
          double qx = gx*gridX;
          double qy = gy*gridY;
          if(qx*qx + qy*qy <= qs[qi]*qs[qi]) {
            assert(gi+2<=bufSize);
            buffer[3*gi]=qx;
            buffer[3*gi+3]=qx;
            buffer[3*gi+1]=qy;
            buffer[3*gi+4]=qy;
            buffer[3*gi+2]=sqrt(qs[qi]*qs[qi] - qx*qx - qy*qy);
            buffer[3*gi+5]=-sqrt(qs[qi]*qs[qi] - qx*qx - qy*qy);
            gi+=2;
          }
        }
      }
      int numOs = gi;
      complex float * ffo = malloc(numOs * sizeof(complex float));
      for (int oi=0; oi<numOs; oi++) {
        ffo[oi]=0;
      }
      for (int ai=0; ai<natom; ai++) {
        for (int oi=0; oi<numOs; oi++) {
          ffo[oi]+=bs[ai]*cexp((-1) * I * (xs[ai]*buffer[3*oi] + ys[ai]*buffer[3*oi+1] + zs[ai]*buffer[3*oi+2]));
        }
      }
      for (int oi=0; oi<numOs; oi++) {
#ifdef SOLVENT_SUBTRACTION
        ffo[oi]-=(buffer[3*oi+0]==0 ? uc[0] : (2*sin(0.5*buffer[3*oi+0]*uc[0])/buffer[3*oi+0]))*cexp(-0.5*I*uc[0]*buffer[3*oi+0]*cv[0])*
                 (buffer[3*oi+1]==0 ? uc[1] : (2*sin(0.5*buffer[3*oi+1]*uc[1])/buffer[3*oi+1]))*cexp(-0.5*I*uc[1]*buffer[3*oi+1]*cv[1])*
                 (buffer[3*oi+2]==0 ? uc[2] : (2*sin(0.5*buffer[3*oi+2]*uc[2])/buffer[3*oi+2]))*cexp(-0.5*I*uc[2]*buffer[3*oi+2]*cv[2])*
                 betaSolv;
#endif
        fis[framecount*numqs + qi]+=creal(conj(ffo[oi]) * ffo[oi]) * M_PI/(uc[0] * uc[1] * qs[qi] * fabs(buffer[3*oi+2]) );
        fis_lutz[framecount*numqs + qi]+=creal(conj(ffo[oi]) * ffo[oi]) / numOs;
      }
      free(ffo);
      //fprintf(stderr,"Done %d / %d qs   \r",qi+1,numqs);
    }
    //fprintf(stderr,"\n");
    //fprintf(stderr,"%e\n",fis[framecount*numqs]);
    framecount++;
  }
  for (int qi=0; qi<numqs; qi++) {

    for (int fi=0; fi<numFrames; fi++)
	    is_lutz[qi]+=fis_lutz[fi*numqs+qi] / numFrames;

    for (int fi=0; fi<numFrames; fi++) {
      double delta = fis[fi*numqs + qi] - is[qi];
      is[qi]+=delta/((double) fi+1);
      sis[qi]+=delta * (fis[fi*numqs + qi] - is[qi]);
    }
    sis[qi]/=((double) numFrames - 1);
  }
  free(fis);


  average_area /= naverage_area;

  //print intensity
  for (int qi=0; qi<numqs; qi++) {
    printf("%f\t%e\t%e\t\n",qs[qi],is[qi]/average_area,sqrt(sis[qi])/average_area, is_lutz[qi]/average_area );
  }

  closeDCD(dcd);

  return (0);
}

void inferCellAlignment(int natom, float *xs, float *ys, float *zs, int *loc) {
  int xcount=0;
  int ycount=0;
  int zcount=0;
  for(int i=0; i<natom; i++) {
    if(xs[i]>0)
      xcount++;
    else
      xcount--;

    if(ys[i]>0)
      ycount++;
    else
      ycount--;

    if(zs[i]>0)
      zcount++;
    else
      zcount--;
  }

  if(xcount > natom/2)
    loc[0]=1;
  else if(xcount < -1*(natom/2))
    loc[0]=-1;
  else
    loc[0]=0;

  if(ycount > natom/2)
    loc[1]=1;
  else if(ycount < -1*(natom/2))
    loc[1]=-1;
  else
    loc[1]=0;

  if(zcount > natom/2)
    loc[2]=1;
  else if(zcount < -1*(natom/2))
    loc[2]=-1;
  else
    loc[2]=0;
}
