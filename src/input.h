#ifndef __INPUTH__
#define __INPUTH__

#define PFFT_HUGE			(1e10)

#define SPHERICAL_CUTOFF	0
#define CYLINDRICAL_CUTOFF	1

struct pfft_instruction
{
	char *forRes;	
	char *depleteRes;
	double mod;
	double R1,R2;
	int fit_parm[3];
	struct pfft_instruction *next;
};

typedef struct deuterate_request
{
	char *resName;
	double fraction;
	double write_fraction;
	double temp;
	int num;
	char *scheme;
	struct deuterate_request *next;
} deuterate_request;

struct parameter_block
{
	double *beta_hist; // filled in by buildData.
	char *scatterFile; // file of scattering lengths

	double lateral_scale;
	double num_water;
	int areaCompensate;
	int post_process;
	int print_fft;
	int intra;
	int max_bin;
	int cutoff_mode; // set in build data
	double cyl_cutoff;
	char *bzFile;
	char *brFile;
	char *name;
	char *qvals;
	double apply_strain;
	double height;
	double av_area;
	double av_bz_area;
	double av_bz_vol;
	double weight;
	char *expFit;
	double cutoff;
	double smooth;
	double water_average_thickness; // distance from the box edge to average to get the water NSLD
	int activateFit;
	double solvent_sub;
	int forceSolventSub;

	double X_scattering_length;


	deuterate_request *deuterations;
	pfft_instruction *instructions;
	parameter_block *next;
};

void getInput( const char *fileName, struct parameter_block *block  ); 
void loadBlockDefaults( struct parameter_block *block );

#endif
