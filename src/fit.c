
/*********************************************************************************** 
 FIT
************************************************************************************ 
This file is part of HALOGEN (http://arxiv.org/abs/1412.5228).
HALOGEN has been developped by Santiago Avila and Steven Murray 
(https://github.com/savila/HALOGEN).
************************************************************************************ */





/*=============================================================================
 *                              LIBRARIES
 *=============================================================================*/
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <omp.h>

#include "read_snapshot.h"
#include "populate_mass_function.h"
#include "place_halos.h"
#include "correlate.h"

//USE GSL FOR SMOOTHED SPLINE+MINIMIZATION
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_statistics.h>

/*=============================================================================
 *                              COMMON DEFINES
 *=============================================================================*/

#define rho_crit (27.755e10)
#define LINELENGTH 256
#define NParam 26

char ParameterList[NParam][32] = {"Snapshot","GadgetFormat","MassFunctionFile",
		"OutputDir","NCellsLin","rho_ref","Overdensity","Mmin",
		"GadL_Unit","GadM_Unit","GadSwap","GadDouble","GadLong","NbodyFile","recalc_frac",
		"Dmax","Dmin","alpha_ratio","alpha_ratio_1","best_alpha","num_alpha","ntrials",
		"nthreads","nr","minr","maxr"};


int ParameterSet[NParam];
int NParametersSet = 0;


// -------------  INPUT PARAMETERS ---------------------------------------------
long seed;


char Snapshot[LINELENGTH];
char OutputDir[LINELENGTH],MassFunctionFile[LINELENGTH];//, alphaFile[LINELENGTH];
char NbodyFile[LINELENGTH];

int format;
double recalc_frac;
float Dmax, Dmin;

int Nlin;
float Mmin;
float OVD;
char rho_ref[8];

//int Nalpha=0;
//double *alpha;
//double *Malpha;

float alpha_ratio,alpha_ratio_1,best_alpha;
int num_alpha, ntrials, nthreads, Nalpha,nr;

double minr, maxr;
float LUNIT, MUNIT;
int SWP, LGADGET, DGADGET;


// ------------- GLOBAL PARAMETERS ---------------------------------------------
float *x, *y, *z, *vx,*vy,*vz,Lbox, mpart, *HaloMass, om_m, rho;
long Npart, Nhalos, **ListOfParticles, *NPartPerCell;
double *mcuts;
int *ncuts;
double *alphavec;
double *betavec;
double **nbody_2pcf;
int total_nr;
int Nmin,Nmax;

char OutputNbody2PCF[LINELENGTH];
char OutputHalogen2PCF[LINELENGTH];
char OutputHalogenErr[LINELENGTH];
char OutputAlphaM[LINELENGTH];

//--------------- FITTING PARAMETERS -------------------------------------------
gsl_bspline_workspace *bw;
gsl_vector *chi2_alpha,*alpha_grid, *B, *c,*weights;

gsl_matrix *X, *cov;

gsl_multifit_linear_workspace *mw;
double chisq, Rsq, dof, tss;
int ncoeffs;
/*=============================================================================
 *                              PROTOTYPES
 *=============================================================================*/

int read_input_file(char *);
int get_mass_bins(float *, long, int, int);
int find_best_alpha();
long read_nbody(float **, float **, float **, float **, float **, float **, float **);
int get_nbody_2pcf(float *, float *, float *, float *, long);
double minimize(double,double,double);
float * compute_fvel(float *vx,float *vy,float *vz,float *nbvx,float *nbvy,float *nbvz,float *nbm);
float vel_std(float *x, float* y, float*z, long start, long end);
/*=============================================================================
 *                              MAIN()
 *=============================================================================*/

int main(int argc, char **argv){

	fprintf(stderr,"\n*******************************************************************\n");
	fprintf(stderr,"**                                                               **\n");
	fprintf(stderr,"**            =         HALOGEN FIT V0.5        =                **\n");
	fprintf(stderr,"**                                                               **\n");
	fprintf(stderr,"**                                                               **\n");
	fprintf(stderr,"**                           let there be knowledge of the dark  **\n");
	fprintf(stderr,"*******************************************************************\n\n");

	if (argc!=2){
		fprintf(stderr,"usage: %s inputfile\n",argv[0]);	
		return -1;
	}

	char inname[256];
	float *hx, *hy, *hz, *hR;
	strcpy(inname,argv[1]);
	double dr;
	double **halogen_2pcf, **halogen_err;
	long nend;
	unsigned long long *dd;
	long i,j,ii;

#ifdef VERB
	fprintf(stderr,"#def VERB\n");
#endif 
#ifdef DEBUG
	fprintf(stderr,"#def DEBUG \n");
#endif
#ifdef ULTRADEBUG
	fprintf(stderr,"#def ULTRADEBUG \n");
#endif
#ifdef ONLYBIG
	fprintf(stderr,"#def ONLYBIG\n");
#endif
#ifdef NO_EXCLUSION
	fprintf(stderr,"#def NO_EXCLUSION\n");
#endif
#ifdef NO_MASS_CONSERVATION
	fprintf(stderr,"#def NO_MASS_CONSERVATION\n");
#endif
#ifdef MASS_OF_PARTS
	fprintf(stderr,"#def MASS_OF_PARTS\n");
#endif
#ifdef RANKED
	fprintf(stderr,"#def RANKED\n");
#endif
#ifdef NDENS
	fprintf(stderr,"#def NDENS\n");
#endif
#ifdef NO_PROB_PARALEL
	fprintf(stderr,"#def NO_PROB_PARALEL\n");
#endif
#ifdef NO_PAR_FIT
	fprintf(stderr,"#def NO_PAR_FIT\n");
#endif
#ifdef MASS_CUTS_FIT
	fprintf(stderr,"#def MASS_CUTS_FIT\n");
#endif
#ifdef BETA0
	fprintf(stderr,"#def BETA0\n");
#endif

	fprintf(stderr,"\nReading input file...\n");
	if (read_input_file(inname)<0)
		return -1;
	fprintf(stderr,"... file read correctly!\n");

	//SETUP OUTPUT FILES
	strcpy(OutputNbody2PCF,OutputDir);
	strcpy(OutputHalogen2PCF,OutputDir);
	strcpy(OutputHalogenErr,OutputDir);
	strcpy(OutputAlphaM,OutputDir);

	strcat(OutputNbody2PCF,"/nbody.2pcf");
	strcat(OutputHalogen2PCF,"/halogen.2pcf");
	strcat(OutputHalogenErr,"/halogen.err");
	strcat(OutputAlphaM,"/M-alpha.txt");

	fprintf(stderr,"Reading Gadget file(s)...\n");
	if (read_snapshot(Snapshot, format, LUNIT, MUNIT, SWP, LGADGET, DGADGET,Nlin,&x, &y, &z, &vx, &vy, &vz, &Npart, &mpart, &Lbox, &om_m,&ListOfParticles,&NPartPerCell)==0)
		fprintf(stderr,"...Gadget file(s) correctly read!\n");
	else {
		fprintf(stderr,"error: Something went wrong reading the gadget file %s\n",inname);
		return -1;
	}
	


	#ifdef VERB
	fprintf(stderr,"\n\tCheck: Npart=%ld, mpart=%e, Lbox=%f\n",Npart,mpart,Lbox);
	fprintf(stderr,"\tx[0]= %f, y[0]= %f, z[0]= %f\n",(x)[0],(y)[0],(z)[0]);
	fprintf(stderr,"\t\tvx[0]= %f, vy[0]= %f, vz[0]= %f\n",(vx)[0],(vy)[0],(vz)[0]);
	fprintf(stderr,"\t      ...\n");
	fprintf(stderr,"\tx[%ld]= %f, y[%ld]= %f, z[%ld]= %f\n",Npart-1,(x)[Npart-1],Npart-1,(y)[Npart-1],Npart-1,(z)[Npart-1]);
	fprintf(stderr,"\t\tvx[%ld]= %f, vy[%ld]= %f, vz[%ld]= %f\n\n",Npart-1,(vx)[Npart-1],Npart-1,(vy)[Npart-1],Npart-1,(vz)[Npart-1]);
	#endif
	
	seed = time(NULL);

	//Generate the halo masses from the mass function
	fprintf(stderr,"Generating Halo Masses...\n");
	Nhalos = populate_mass_function(MassFunctionFile,Mmin,Lbox,&HaloMass,seed,160);
	if (Nhalos<0){
		fprintf(stderr,"error: Couldnt create HaloMass array\n");	
		return -1;
	}
	fprintf(stderr,"...Halo Masses Generated\n");

	// GET THE ACTUAL NUMBER OF R VALUES IN CORR (LINEARLY SPACED)
	dr = (maxr - minr) / nr;
	total_nr = (int) ceil(maxr / dr) + 2;
	maxr += 2* dr;

	//density at the boundary of a halo
	if (strcmp(rho_ref,"crit")==0)
		rho = OVD*rho_crit;
	else 
		rho = OVD*rho_crit*om_m;

	Nmin = Lbox*Lbox*Lbox*Dmin;
	Nmax = Lbox*Lbox*Lbox*Dmax;
#ifdef VERB
	fprintf(stderr,"Generating Mass Bins...\n");
#endif
	Nalpha = get_mass_bins(HaloMass,Nhalos,Nmax,Nmin);
#ifdef VERB
	fprintf(stderr,"... generated mass bins\n");
	fprintf(stderr,"Nmin=%d, Nmax=%d, Nmassbins=%d (=Nalphas).\n",Nmin,Nmax,Nalpha);
	fprintf(stderr,"Getting NBODY 2PCF...\n");
#endif
	


	float *nbx, *nby, *nbz, *nbvx, *nbvy, *nbvz, *nbm,*fvel;
	long nb_n;

	// Import the Nbody halos
	nb_n = read_nbody(&nbx,&nby,&nbz,&nbvx,&nbvy,&nbvz,&nbm);
	fprintf(stderr,"READ IT IN\n");
	
	fprintf(stderr,"\tComputing velocity bias fvel...\n");
	fvel = compute_fvel(vx,vy,vz,nbvx,nbvy,nbvz,nbm);	
	free(vx);
	free(vy);
	free(vz);
	free(nbvx);
	free(nbvy);
	free(nbvz);

	#ifdef VERB
	for (ii=0;ii<Nalpha;ii++)
		fprintf(stderr,"\t\tfvel[%ld]=%f\n",ii,fvel[ii]);
	#endif
	fprintf(stderr,"\t...done!\n");
	

	// Get the Nbody 2PCF
	get_nbody_2pcf(nbx,nby,nbz,nbm,nb_n);

	free(nbx);
	free(nby);
	free(nbz);
	free(nbm);

	// DO THE FIT
	if(find_best_alpha()==0)
		fprintf(stderr, "... done fitting.\n");
	else {
		fprintf(stderr,"Problem fitting\n");
		return -1;
	}

	hx = (float *) calloc(Nhalos,sizeof(float));
	hy = (float *) calloc(Nhalos,sizeof(float));
	hz = (float *) calloc(Nhalos,sizeof(float));
	hR = (float *) calloc(Nhalos,sizeof(float));


	seed++;

	// Do a final placement with the correct alphavec
	place_halos(Nhalos,HaloMass, Nlin, Npart, x, y, z, NULL,NULL,NULL,Lbox,
				rho,seed,
				mpart,nthreads, alphavec, betavec, mcuts, Nalpha, recalc_frac,hx, hy, hz,
				NULL,NULL,NULL, hR,ListOfParticles,NPartPerCell);
		
	// ALLOCATE the halogen_2pcf array
	halogen_2pcf = (double **) calloc(Nalpha,sizeof(double *));
	halogen_err = (double **) calloc(Nalpha,sizeof(double *));

	nend = 0;
	dd = (unsigned long long *) calloc(total_nr,sizeof(unsigned long long));
	for(ii=0;ii<Nalpha;ii++){
		(halogen_2pcf)[ii] = (double *) calloc(total_nr,sizeof(double));
		(halogen_err)[ii] = (double *) calloc(total_nr,sizeof(double));
		while(HaloMass[nend]>mcuts[ii] && nend<Nhalos){
			nend++;
		}
		if(nend==Nhalos){
			fprintf(stderr,"ERROR: HALOMASSES DON'T REACH MALPHA_MIN\n");
			return -1;
		}
		// Do the correlation
		correlate(nend, Lbox,hx,hy,hz,halogen_2pcf[ii],halogen_err[ii], dd,total_nr,maxr,160);
	}

	// OUTPUT
	FILE* tpcfile;
	FILE* errfile;
	tpcfile = fopen(OutputHalogen2PCF,"w");
	errfile = fopen(OutputHalogenErr,"w");
	if(tpcfile==NULL){
		fprintf(stderr,"UNABLE TO OPEN NBODY 2PCF FILE\n");
		exit(1);
	}
	if(errfile==NULL){
		fprintf(stderr,"UNABLE TO OPEN NBODY 2PCF FILE\n");
		exit(1);
	}

	fprintf(stderr,"OPENED TPC AND ERR\n");
	for(i=0;i<total_nr;i++){
		fprintf(tpcfile,"%e\t",(i+0.5)*dr);
		fprintf(errfile,"%e\t",(i+0.5)*dr);
		for(j=0;j<Nalpha;j++){
			fprintf(tpcfile,"%e\t",halogen_2pcf[j][i]);
			fprintf(errfile,"%e\t",halogen_err[j][i]);
		}
		fprintf(stderr,"\n");
		fprintf(tpcfile,"\n");
		fprintf(errfile,"\n");
	}
	fclose(tpcfile);
	fclose(errfile);

	fprintf(stderr,"outputted tpc and errfile\n");
	FILE* alphafile;
	alphafile = fopen(OutputAlphaM,"w");
	for(i=0;i<Nalpha;i++){
		fprintf(alphafile,"%e\t%e\t%e\n",mcuts[i],alphavec[i],fvel[i]);
	}
	fclose(alphafile);
	for(ii=0;ii<Nalpha;ii++){
		free((halogen_2pcf)[ii]);
		free((halogen_err)[ii]);
	}
	free(halogen_2pcf);
	free(halogen_err);
	free(dd);	
	free(hx);	
	free(hy);	
	free(hz);	
	free(hR);	


	fprintf(stderr,"\n*******************************************************************\n");
	fprintf(stderr,"**                                     ... the fit is complete.  **\n");
	fprintf(stderr,"*******************************************************************\n\n");

	return 0;
}


/*=============================================================================
 *                              GET MASS BINS
 *=============================================================================*/


int get_mass_bins(float *HaloMass, long Nhalos, int Nmax, int Nmin){
	//Get appropriate mass bins.
	
	float dm_max=0.5;
	long halos_used = 0;
	int i,nbins=0;
	float dm;
	float thiscutoff = 0.0;

	if (Nmax==0 && Nmin==0){
		fprintf(stderr,"Using a single alpha\n");
		mcuts = (double *) malloc(1*sizeof(double));
		ncuts = (int *) malloc(1*sizeof(int));
		mcuts[0] = HaloMass[Nhalos-1];
		ncuts[0] = Nhalos;
		return 1;
	}
		
	// First count how many bins there needs to be
#ifdef VERB
	fprintf(stderr,"\tCounting number of bins: ");
#endif
	while (Nhalos - halos_used > Nmax){
		nbins++;
		i = 0;
		dm = 0.0;
		while (dm < dm_max && i * Nmin < Nmax){
			 i++;
			 thiscutoff = HaloMass[halos_used + i * Nmin - 1];
			 dm = log10(HaloMass[halos_used] / thiscutoff);
		}
		
		halos_used += Nmin*i;
	}


	// Allocate memory
	mcuts = (double *) malloc(nbins*sizeof(double));
	ncuts = (int *) malloc(nbins*sizeof(int));
	
#ifdef VERB
	fprintf(stderr,"\tAssigning masses\n");
#endif
	// Assign values to arrays
	nbins = 0;
	halos_used = 0;
	while (Nhalos - halos_used > Nmax){
		i = 0;
		dm = 0.0;
		while (dm < dm_max && i * Nmin < Nmax){
			 i++;
			 thiscutoff = HaloMass[halos_used + i * Nmin - 1];
			 dm = log10(HaloMass[halos_used] / thiscutoff);
		}
		mcuts[nbins] = (double) thiscutoff ;
		halos_used += Nmin*i;
		ncuts[nbins] = Nmin*i;
		nbins++;
	}
	// Make sure the last bin gets the rest.
	mcuts[nbins-1] = HaloMass[Nhalos-1];
	ncuts[nbins-1] += Nhalos - halos_used;
			
#ifdef VERB
			fprintf(stderr,"\tCUTOFFS: ");
			for (i=0;i<nbins;i++){
				fprintf(stderr,"%e, ",mcuts[i]);
			}
			fprintf(stderr,"\n\tNbins: ");
			for (i=0;i<nbins;i++){
				fprintf(stderr,"%d, ",ncuts[i]);
			}
			fprintf(stderr,"\n");
#endif
	return nbins;
}


/*=============================================================================
 *                              READ ASCII NBODY
 *=============================================================================*/

long read_nbody(float **nbx, float **nby, float **nbz, float **nbvx, float **nbvy, float **nbvz, float **nbm){

	long ii;
	long n_nb=0;
	int ch=0;

	FILE* nbodyf;

	nbodyf = fopen(NbodyFile,"r");
	//First get number of lines
	do
	{
	    ch = fgetc(nbodyf);
	    if(ch == '\n')
	    	n_nb++;
	} while (ch != EOF);

	// last line doesn't end with a new line!
	// but there has to be a line at least before the last line
	if(ch != '\n' && n_nb != 0)
	    n_nb++;

	n_nb--;

	fprintf(stderr,"SIZE OF NBODY FILE: %ld\n", n_nb);
	//Allocate
	(*nbx) = (float *) calloc(n_nb,sizeof(float));
	(*nby) = (float *) calloc(n_nb,sizeof(float));
	(*nbz) = (float *) calloc(n_nb,sizeof(float));
	(*nbvx) = (float *) calloc(n_nb,sizeof(float));
	(*nbvy) = (float *) calloc(n_nb,sizeof(float));
	(*nbvz) = (float *) calloc(n_nb,sizeof(float));
	(*nbm) = (float *) calloc(n_nb,sizeof(float));

	rewind(nbodyf);

	// Read in data

	float a,b,c,d,e,f,g;
	char line[1000];

	for(ii=0;ii<n_nb;ii++){
		//fscanf(nbodyf,"%f%f%f%f",&((*nbx)[ii]),&((*nby)[ii]),&((*nbz)[ii]),&((*nbm)[ii]));
		fgets(line,LINELENGTH,nbodyf);
		sscanf(line,"%f %f %f %f %f %f %f",&a,&b,&c,&d,&e,&f,&g);
		if (a==Lbox)
			a=0.;
		(*nbx)[ii]=a;

		if (b==Lbox)
			b=0.;
		(*nby)[ii]=b;

		if (c==Lbox)
			c=0.;
		(*nbz)[ii]=c;
		
		(*nbvx)[ii]=d;
		(*nbvy)[ii]=e;
		(*nbvz)[ii]=f;
		(*nbm)[ii]=g;
		
	}
	fclose(nbodyf);
	return n_nb;
}

/*=============================================================================
 *                              COMPUTE Fvel
 *=============================================================================*/
float * compute_fvel(float *vx,float *vy,float *vz,float *nbvx,float *nbvy,float *nbvz,float *nbm){
	
	long ii,Nstart=0,Nend=0;
	float * vbias;

	vbias = (float *) calloc(Nalpha,sizeof(float));	
        for(ii=0;ii<Nalpha;ii++){
	Nstart=Nend;
	#ifndef MASS_CUTS_FIT
		Nend+=ncuts[ii];
	#else
                while(nbm[nend]>mcuts[ii]){
                        nend++;
                        if(nend==n_nb){
                                fprintf(stderr,"ERROR: NBODY MASSES DON'T REACH MALPHA_MIN\n");
                                return NULL;
                        }
                }

	#endif
		vbias[ii]=vel_std(nbvx,nbvy,nbvz,Nstart,Nend)/vel_std(vx,vy,vz,Nstart,Nend);	
	}
	return vbias;
}


float vel_std(float *x, float* y, float*z, long start, long end){
	int i;
	float mean=0.0, std=0.0;
	for (i=start;i<end;i++){
		mean+=x[i];
		mean+=y[i];
		mean+=z[i];
		std+=(x[i]*x[i]);
		std+=(y[i]*y[i]);
		std+=(z[i]*z[i]);
	}
	mean = mean/(3*(end-start));
	std = std/(3*(end-start));
	std = sqrt(std-mean*mean);
	return std;
}


/*=============================================================================
 *                              GET NBODY_2PCF
 *=============================================================================*/
int get_nbody_2pcf(float *nbx, float *nby, float *nbz, float *nbm, long n_nb){
	int ii,i,j;
	long nend=0;

	double *ercor;
	unsigned long long *dd;

	ercor = (double *) calloc(total_nr,sizeof(double));
	dd = (unsigned long long *) calloc(total_nr,sizeof(unsigned long long));
	// Find the boxsize of the nbody file.
	float nbody_lbox = 0.0;
	for (ii=0;ii<n_nb;ii++){
		if(nbx[ii]>nbody_lbox) nbody_lbox = nbx[ii];
		if(nby[ii]>nbody_lbox) nbody_lbox = nby[ii];
		if(nbz[ii]>nbody_lbox) nbody_lbox = nbz[ii];
	}
	nbody_lbox = (float) ceil(nbody_lbox);
	
	fprintf(stderr,"BOXSIZE OF NBODY: %e\n", nbody_lbox);
	fprintf(stderr,"HALOS in NBODY: %ld\n",n_nb);

	fprintf(stderr,"M[first]=%e, M[last]=%e \n",nbm[0],nbm[n_nb-1]);

	// ALLOCATE the nbody_2pcf array
	nbody_2pcf = (double **) calloc(Nalpha,sizeof(double *));
	for(ii=0;ii<Nalpha;ii++){
//		fprintf(stderr,"\ti_alpha=%d, mcuts= %e",ii,mcuts[ii]);
		nbody_2pcf[ii] = (double *) calloc(total_nr,sizeof(double));
		fprintf(stderr,"    allocated %d ",total_nr);
//		fprintf(stderr,"\tNend=%ld\n",nend);
		fprintf(stderr,"M[first]=%e \n",nbm[0]);
//		fprintf(stderr,"M[last]=%e \n",nbm[n_nb-1]);
//		fprintf(stderr,"\tnbm=%e\n",nbm[nend]);
		#ifdef MASS_CUTS_FIT
		while(nbm[nend]>mcuts[ii]){
			nend++;
			if(nend==n_nb){
				fprintf(stderr,"ERROR: NBODY MASSES DON'T REACH MALPHA_MIN\n");
				return -1;
			}
		//	fprintf(stderr,"\n-%ld\n",nend);
		}
		#else
		nend+=ncuts[ii];
		#endif
		fprintf(stderr,"\tNend=%ld  nbm=%e\n",nend,nbm[nend]);
		// Do the correlation
		correlate(nend, nbody_lbox,nbx,nby,nbz,nbody_2pcf[ii],ercor,dd,total_nr,maxr,160);
		fprintf(stderr,"\tCorrelated\n\n");
	}

	// OUTPUT
	FILE* tpcfile;
	tpcfile = fopen(OutputNbody2PCF,"w");
	if(tpcfile==NULL){
		fprintf(stderr,"UNABLE TO OPEN NBODY 2PCF FILE\n");
		exit(1);
	}
	for(i=0;i<total_nr;i++){
		for(j=0;j<Nalpha;j++){
			fprintf(tpcfile,"%e\t",nbody_2pcf[j][i]);
		}
		fprintf(tpcfile,"\n");
	}
	fclose(tpcfile);
	fprintf(stderr,"NBODY FILE WRITTEN\n");
	return 0;

}


/*=============================================================================
 *                              FITTING CODE
 *=============================================================================*/
int find_best_alpha(){
	
	//double alpha_2pcf[num_alpha][nr], trials_2pcf[num_alpha][ntrials][total_nr], alpha_err[num_alpha][nr];
	double **alpha_2pcf, ***trials_2pcf, **alpha_err;
	double low, high, da;
	long i,nend,k,j,ir,jk;

	double chi2;

	double yi,yerr,xi,minerr;
	int nbreak, iii,ind,min_ind,max_ind;

	long thisseed[Nalpha*num_alpha*ntrials];
#ifdef ALLOUT
	char Output2PCF[LINELENGTH];
	char ThisFile[15];
#endif

	ncoeffs = num_alpha-1;
	nbreak =ncoeffs-2;
	// Allocate alphavec, which is our result
	alphavec = (double *) calloc(Nalpha,sizeof(double));
	betavec = (double *) calloc(Nalpha,sizeof(double));


	alpha_2pcf = (double **) calloc(num_alpha,sizeof(double *));
	alpha_err = (double **) calloc(num_alpha,sizeof(double *));
	trials_2pcf = (double ***) calloc(num_alpha,sizeof(double **));
	for (i=0;i<num_alpha;i++){
		alpha_2pcf[i] = (double *) calloc(nr,sizeof(double));
		alpha_err[i] = (double *) calloc(nr,sizeof(double));
		trials_2pcf[i] = (double **) calloc(ntrials,sizeof(double*));
		for (j=0;j<ntrials;j++)
			trials_2pcf[i][j] = (double *) calloc(total_nr,sizeof(double));

	}

	// allocate a cubic bspline workspace (k = 4) 
	bw = gsl_bspline_alloc(4, nbreak);
	B = gsl_vector_alloc(ncoeffs);
	alpha_grid = gsl_vector_alloc(num_alpha);
	chi2_alpha = gsl_vector_alloc(num_alpha);
	X = gsl_matrix_alloc(num_alpha,ncoeffs);

	c = gsl_vector_alloc(ncoeffs);
	weights = gsl_vector_alloc(num_alpha);
	cov = gsl_matrix_alloc(ncoeffs, ncoeffs);
	mw = gsl_multifit_linear_alloc(num_alpha, ncoeffs);


	// Loop through each mass bin
	nend = 0;
	iii  = 0;

	seed = time(NULL);
	for (i=0;i<Nalpha;i++){
		for(j=0;j<num_alpha;j++){
			for(k=0;k<ntrials;k++){
				thisseed[iii] = seed + iii;
				iii++;
			}
		}
	}
	fprintf(stderr,"STARTING LOOPS...\n");
	iii = 0;
	for (i=0;i<Nalpha;i++){
		// Get where to place halos to
		nend += ncuts[i];
		
		// Calculate a more optimum alpha_grid for this bin
		if (i < 2)
			low = alpha_ratio_1*best_alpha;
		else{
			low = alpha_ratio*best_alpha;
		}
		high = 1.05*best_alpha;
		da = (high-low)/(num_alpha-1);

		fprintf(stderr,"STARTING THREADS\n");
	#ifndef NO_PAR_FIT
		#pragma omp parallel for  num_threads(nthreads) private(jk,j,k) \
		shared(num_alpha,ntrials,i,alpha_grid,low,da,stderr,Nalpha,alphavec,iii,mcuts,ListOfParticles,\
		NPartPerCell,x,y,z,Lbox,Npart,Nlin,HaloMass,nend,mpart,recalc_frac,betavec,thisseed,trials_2pcf,rho,maxr,\
		total_nr,Nhalos) default(none)
	#endif
		for(jk=0;jk<num_alpha*ntrials;jk++){
  			fprintf(stderr,"Thread %d/%d\n",omp_get_thread_num(),omp_get_num_threads());
			
			float *hx,*hy,*hz,*hR;			
			double *thisalphavec;
			double *ercorr;
			unsigned long long *DD;
			hx = (float *) calloc(Nhalos,sizeof(float));
			hy = (float *) calloc(Nhalos,sizeof(float));
			hz = (float *) calloc(Nhalos,sizeof(float));
			hR = (float *) calloc(Nhalos,sizeof(float));
			thisalphavec = (double *) calloc(Nalpha,sizeof(double));
			DD = (unsigned long long *) calloc(total_nr,sizeof(unsigned long long));
			ercorr = (double *) calloc(total_nr,sizeof(double));

			k=jk%ntrials;
			j=jk/ntrials;

			fprintf(stderr,"MOVED MEMORY\n");
			fprintf(stderr,"GOT k,j: %ld, %ld\n",k,j);
			
			fprintf(stderr,"sizeof M : %f MB\n",(Nlin*Nlin*Nlin*sizeof(double)/1024./1024.));
			//define the alpha grid for this mass bin
			gsl_vector_set(alpha_grid,j,low+j*da);
			

			// create a local alphavec, to which we add the current gridded alpha
			int itemp;
			for (itemp=0;itemp<Nalpha;itemp++)
				thisalphavec[itemp] = alphavec[itemp];


			//fprintf(stderr,"mem copied\n");
			thisalphavec[i] = gsl_vector_get(alpha_grid,j);

			// Place the halos

			#ifdef NO_PAR_FIT
			place_halos(nend,HaloMass, Nlin, Npart, x, y, z, NULL,NULL,NULL,Lbox,
					rho,thisseed[jk+i*num_alpha*ntrials],
					mpart,nthreads,thisalphavec, betavec, mcuts, Nalpha, recalc_frac,hx, hy, hz,
					NULL,NULL,NULL,hR,ListOfParticles,NPartPerCell);
			#else 
			place_halos(nend,HaloMass, Nlin, Npart, x, y, z, NULL,NULL,NULL,Lbox,
					rho,thisseed[jk+i*num_alpha*ntrials],
					mpart,1,thisalphavec, betavec, mcuts, Nalpha, recalc_frac,hx, hy, hz,
					NULL,NULL,NULL,hR,ListOfParticles,NPartPerCell);

			#endif			
			fprintf(stderr,"correlating...\n");
			//Get the 2PCF
			correlate(nend, Lbox,hx,hy,hz,trials_2pcf[j][k], ercorr, DD,total_nr,maxr,1);
			fprintf(stderr,"...correlated\n");
			
			free(hx);
			free(hy);
			free(hz);
			free(hR);
			free(thisalphavec);
			free(ercorr);
			free(DD);

		}

		for (jk=0;jk<num_alpha*ntrials;jk++){
			k=jk%ntrials;
			j=jk/ntrials;
			fprintf(stderr,"RAWCORR %ld, %ld: %e\n",j,k,trials_2pcf[j][k][0]);
		}

	//	#pragma omp parallel for private(ir,j,chi2,k) shared(num_alpha,stderr,i,nr,alpha_2pcf) default(none)
		fprintf(stderr,"GOT CORRELATIONS , GETTING CHI2...\n");
		for(j=0;j<num_alpha;j++){
			//Get mean and stdev of trials_2pcf
			chi2 = 0.0;
			for (ir=0;ir<nr;ir++){
				alpha_2pcf[j][ir] = 0.0;
				for(k=0;k<ntrials;k++){
					alpha_2pcf[j][ir] += trials_2pcf[j][k][ir+total_nr-nr-2];
				}
				alpha_2pcf[j][ir] = alpha_2pcf[j][ir]/ntrials; 
				alpha_err[j][ir] = 0.0;
				for(k=0;k<ntrials;k++){
					alpha_err[j][ir] += pow((trials_2pcf[j][k][ir+total_nr-nr-2]-alpha_2pcf[j][ir]),2);
				}
				alpha_err[j][ir] = pow((alpha_err[j][ir]/(ntrials-1)),0.5);

				// Now get chi^2 values
				#ifdef REL_CHI2
				chi2 += pow(((alpha_2pcf[j][ir]-nbody_2pcf[i][ir+total_nr-nr-2])/nbody_2pcf[i][ir+total_nr-nr-2]),2);
				#else
				chi2 += pow(((alpha_2pcf[j][ir]-nbody_2pcf[i][ir+total_nr-nr-2])/alpha_err[j][ir]),2);
				#endif

				fprintf(stderr,"%ld, %ld, %e, %e, %e, %e\n",j,ir,alpha_2pcf[j][ir],nbody_2pcf[i][ir+total_nr-nr-2],alpha_err[j][ir],chi2);
			}
			gsl_vector_set(chi2_alpha,j,chi2/nr);
			gsl_vector_set(weights,j,nr/chi2);
		}
	//	#endif //NO_PARALLEL_FIT
//*/

#ifdef ALLOUT
		//OUTPUT SOME THINGS FOR CHECKING
		sprintf(ThisFile,"/raw.2pcf.%ld",i);
		strcpy(Output2PCF,OutputDir);
		strcat(Output2PCF,ThisFile);
		FILE* output_2pcf;
		output_2pcf = fopen(Output2PCF,"w");

		//header
		fprintf(output_2pcf,"# r, ");
		for (k=0;k<num_alpha;k++){
			fprintf(output_2pcf,"%e\t",gsl_vector_get(alpha_grid,k));
		}
		fprintf(output_2pcf,"\n");

		//table
		for (ir=0;ir<nr;ir++){
			fprintf(output_2pcf,"%e\t",(ir+total_nr-nr+0.5)*maxr/total_nr);
			for(k=0;k<num_alpha-1;k++){
				fprintf(output_2pcf,"%e\t",alpha_2pcf[k][ir]);
			}
			fprintf(output_2pcf,"%e\n",alpha_2pcf[num_alpha-1][ir]);
		}
		fclose(output_2pcf);

		sprintf(ThisFile,"/raw.err.%ld",i);
		strcpy(Output2PCF,OutputDir);
		strcat(Output2PCF,ThisFile);
		FILE* output_err;
		output_err = fopen(Output2PCF,"w");

		//header
		fprintf(output_err,"# r, ");
		for (k=0;k<num_alpha;k++){
			fprintf(output_err,"%e\t",gsl_vector_get(alpha_grid,k));
		}
		fprintf(output_err,"\n");

		//table
		for (ir=0;ir<nr;ir++){
			fprintf(output_err,"%e\t",(ir+total_nr-nr+0.5)*maxr/total_nr);
			for(k=0;k<num_alpha-1;k++){
				fprintf(output_err,"%e\t",alpha_err[k][ir]);
			}
			fprintf(output_err,"%e\n",alpha_err[num_alpha-1][ir]);
		}
		fclose(output_err);

		sprintf(ThisFile,"/chi2.%ld",i);
		strcpy(Output2PCF,OutputDir);
		strcat(Output2PCF,ThisFile);
		FILE* chifile;
		chifile = fopen(Output2PCF,"w");
		fprintf(chifile,"# alpha, chi2, weight\n");
		for (k=0;k<num_alpha;k++){
			fprintf(chifile,"%e\t%e\t%e\n",gsl_vector_get(alpha_grid,k),gsl_vector_get(chi2_alpha,k),
					gsl_vector_get(weights,k));
		}
		fclose(chifile);

#endif


		// Check if final value or initial value is the smallest
		minerr = gsl_vector_get(chi2_alpha,0);
		ind = 0;
		for(k=1;k<num_alpha;k++){
			if(minerr>gsl_vector_get(chi2_alpha,k)){
				minerr = gsl_vector_get(chi2_alpha,k);
				ind = k;
			}
		}
		if(ind==0){
			fprintf(stderr,"ERROR: alpha_grid doesn't extend low enough, set alpha_ratio lower");
		}
		if(ind==num_alpha-1){
			fprintf(stderr,"ERROR: alpha_grid doesn't extend high enough, set best_alpha higher");
		}
		if (ind>=2){
			min_ind = ind-2;
		}else{
			min_ind = 0;
		}
		if (ind<=num_alpha-3){
			max_ind = ind+2;
		}else{
			max_ind = num_alpha-1;
		}

		/* use uniform breakpoints on interval */
		gsl_bspline_knots_uniform(gsl_vector_get(alpha_grid,0), gsl_vector_get(alpha_grid,num_alpha-1), bw);

		/* construct the fit matrix X */
		for (k = 0; k < num_alpha; ++k){
			double xi = gsl_vector_get(alpha_grid, k);

			/* compute B_j(xi) for all j */
			gsl_bspline_eval(xi, B, bw);

			/* fill in row i of X */
			for (j = 0; j < ncoeffs; ++j){
				double Bj = gsl_vector_get(B, j);
				gsl_matrix_set(X, k, j, Bj);
		    }
		}

		fprintf(stderr,"Got to the fit part\n");
		/* do the fit */
		gsl_multifit_wlinear(X, weights, chi2_alpha, c, cov, &chisq, mw);

		// PRINT OUT EVERYTHING WE HAVE SO FAR
		fprintf(stderr,"alpha\tchi2_alpha\t\n");
		for (k=0;k<num_alpha;k++){
			fprintf(stderr,"%e\t%e\t\n",gsl_vector_get(alpha_grid,k),
					gsl_vector_get(chi2_alpha,k));
		}

#ifdef VERB
		dof = num_alpha - ncoeffs;
		tss = gsl_stats_wtss(weights->data, 1, chi2_alpha->data, 1, chi2_alpha->size);
		Rsq = 1.0 - chisq / tss;

		fprintf(stderr, "chisq/dof = %e, Rsq = %f\n",chisq / dof, Rsq);
#endif

		for (xi=gsl_vector_get(alpha_grid,0);xi<gsl_vector_get(alpha_grid,num_alpha-1);xi+=0.01){
			gsl_bspline_eval(xi, B, bw);
			gsl_multifit_linear_est(B, c, cov, &yi, &yerr);
			fprintf(stderr,"%e\t%e\t%e\n",xi,yi,gsl_vector_get(B,0));
		}
		//DO THE MINIMIZATION
		alphavec[i] = minimize(gsl_vector_get(alpha_grid,min_ind), gsl_vector_get(alpha_grid,max_ind),
									   gsl_vector_get(alpha_grid,ind));
		best_alpha = alphavec[i];

		fprintf(stderr,"\n Best alpha: %f\n\n",best_alpha);
	}

	return 0;

}
/*=============================================================================
 *                              MINIMIZING FUNCTION
 *=============================================================================*/
double min_func(double x, void *params){

	double y, yerr;
	gsl_vector *B;

	B = gsl_vector_alloc(ncoeffs);


	gsl_bspline_eval(x, B, bw);
	gsl_multifit_linear_est(B, c, cov, &y, &yerr);

	return y;
}
/*=============================================================================
 *                              MINIMIZATION
 *=============================================================================*/
double minimize(double a,double b, double m){
	  int status;
	  int iter = 0, max_iter = 100;
	  const gsl_min_fminimizer_type *T;
	  gsl_min_fminimizer *s;
	  double epsabs= 0.005; //eps of solution
	  double epsrel= 0.0; //eps of solution

	  gsl_function F;

	  F.function = &min_func;
	  F.params = 0;

	  T = gsl_min_fminimizer_brent;
	  s = gsl_min_fminimizer_alloc (T);
	  gsl_min_fminimizer_set (s, &F, m, a, b);

	  printf ("using %s method\n",
	          gsl_min_fminimizer_name (s));

	  printf ("%5s [%9s, %9s] %9s %9s\n",
	          "iter", "lower", "upper", "min",
	          "err(est)");

	  printf ("%5d [%.7f, %.7f] %.7f %.7f\n",
	          iter, a, b,
	          m, b - a);

	  do
	    {
	      iter++;
	      status = gsl_min_fminimizer_iterate (s);

	      m = gsl_min_fminimizer_x_minimum (s);
	      a = gsl_min_fminimizer_x_lower (s);
	      b = gsl_min_fminimizer_x_upper (s);

	      status
	        = gsl_min_test_interval (a, b, epsabs, epsrel);

	      if (status == GSL_SUCCESS)
	        printf ("Converged:\n");

	      printf ("%5d [%.7f, %.7f] "
	              "%.7f %.7f\n",
	              iter, a, b,
	              m, b - a);
	    }
	  while (status == GSL_CONTINUE && iter < max_iter);

	  gsl_min_fminimizer_free (s);

	  return m;
	}

/*=============================================================================
 *                              set seed
 *=============================================================================*/
/*
int set_seed(int index){
	//index is the index of the seed array to be populated

	long this_seed;
	int ii;
	int its_in_there=1;

	this_seed = time(NULL);

	do{
		its_in_there = 0;
		for(ii=0;ii<index;ii++){
			if (seed[ii]==this_seed){
				its_in_there = 1;
				this_seed++;
			}
		}
	} while(its_in_there==1);

	seed[index] = this_seed;
	return 0;
}
*/
/*=============================================================================
 *                              I/O
 *=============================================================================*/
int read_input_file(char *name){
	FILE *f;
	char line[LINELENGTH],word[LINELENGTH];
	int i;
	//float x,y;
	for (i=0;i<NParam;i++) 
		ParameterSet[i]=0;	

	if ((f = fopen(name,"r"))==NULL){
		fprintf(stderr,"Could not open input file %s\n",name);
		return -1;
	}
        while(fgets(line,LINELENGTH,f)!=NULL) {	
                if (line[0] != '#') {				//Omit comment lines
			if (sscanf(line,"%s",word)!=EOF){ 	//Check it is not an empty line
			  for (i=0;i<NParam;i++)		//Find the parameter specified
				if (strcmp(ParameterList[i],word)==0)
					break;
			  switch (i){	
				case 0:
					sscanf(line,"%s %s",word,Snapshot);
					if (ParameterSet[i]>0)
						fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
					else
						NParametersSet++;	
					#ifdef VERB 
					fprintf(stderr,"\t%s: %s\n",word,Snapshot);
					#endif
					ParameterSet[i]++;
					break;
				case 1:
					sscanf(line,"%s %d",word,&format);
					if (format!=1 && format !=2){
						fprintf(stderr,"ERROR: invalid parameter %s: %d. Please select 1 or 2\n",word,format);
						return -1;
					}	
					if (ParameterSet[i]>0)
						fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
					else
						NParametersSet++;	
					#ifdef VERB 
					fprintf(stderr,"\t%s: %d\n",word,format);
					#endif
					ParameterSet[i]++;
					break;
					
				case 2:
					sscanf(line,"%s %s",word,MassFunctionFile);
					if (ParameterSet[i]>0)
						fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
					else
						NParametersSet++;
					#ifdef VERB 
					fprintf(stderr,"\t%s: %s\n",word,MassFunctionFile);
					#endif
					ParameterSet[i]++;
					break;

				case 3:
					sscanf(line,"%s %s",word,OutputDir);
					if (ParameterSet[i]>0)
						fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
					else
						NParametersSet++;
					#ifdef VERB 
					fprintf(stderr,"\t%s: %s\n",word,OutputDir);
					#endif
					ParameterSet[i]++;
					break;

				case 4:
					sscanf(line,"%s %d",word,&Nlin);
					if (ParameterSet[i]>0)
						fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
					else
						NParametersSet++;
					#ifdef VERB 
					fprintf(stderr,"\t%s: %d\n",word,Nlin);
					#endif
					ParameterSet[i]++;
					break;

				case 5:
					sscanf(line,"%s %s",word,rho_ref);
					if (ParameterSet[i]>0)
						fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
					else
						NParametersSet++;
					#ifdef VERB 
					fprintf(stderr,"\t%s: %s\n",word,rho_ref);
					#endif
					if ((strcmp(rho_ref,"crit")!=0) && (strcmp(rho_ref,"matter")!=0)){
							fprintf(stderr,"ERROR: Not valid option for %s: %s.\nPlease select \"crit\" or \"matter\". \n",word,rho_ref);
							return -1;
					}
					ParameterSet[i]++;
					break;

				case 6:
					sscanf(line,"%s %f",word,&OVD);
					if (ParameterSet[i]>0)
						fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
					else
						NParametersSet++;
					#ifdef VERB 
					fprintf(stderr,"\t%s: %f\n",word,OVD);
					#endif
					ParameterSet[i]++;
					break;

				case 7:
					sscanf(line,"%s %f",word,&Mmin);
					if (ParameterSet[i]>0)
						fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
					else
						NParametersSet++;
					#ifdef VERB 
					fprintf(stderr,"\t%s: %f.\n",word,Mmin);
					#endif
					ParameterSet[i]++;
					break;

				case 8:
					sscanf(line,"%s %f",word,&LUNIT);
					if (ParameterSet[i]>0)
						fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
					else
						NParametersSet++;
					ParameterSet[i]++;
					break;

				case 9:
					sscanf(line,"%s %f",word,&MUNIT);
					if (ParameterSet[i]>0)
						fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
					else
						NParametersSet++;
					#ifdef VERB 
					fprintf(stderr,"\t%s: %f\n",word,MUNIT);
					#endif
					ParameterSet[i]++;
					break;

				case 10:
					sscanf(line,"%s %d",word, &SWP);
					if (ParameterSet[i]>0)
						fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
					else
						NParametersSet++;
					#ifdef VERB 
					fprintf(stderr,"\t%s: %d\n",word,SWP);
					#endif
					ParameterSet[i]++;
					break;

				case 11:
					sscanf(line,"%s %d",word,&DGADGET);
					if (ParameterSet[i]>0)
						fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
					else
						NParametersSet++;
					#ifdef VERB 
					fprintf(stderr,"\t%s: %d\n",word,DGADGET);
					#endif
					ParameterSet[i]++;
					break;

				case 12:
					sscanf(line,"%s %d",word,&LGADGET);
					if (ParameterSet[i]>0)
						fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
					else
						NParametersSet++;
					#ifdef VERB
					fprintf(stderr,"\t%s: %d\n",word,LGADGET);
					#endif
					ParameterSet[i]++;
					break;
				case 13:
					sscanf(line,"%s %s",word,NbodyFile);
					if (ParameterSet[i]>0)
						fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
					else
						NParametersSet++;
					#ifdef VERB
					fprintf(stderr,"\t%s: %s\n",word,NbodyFile);
					#endif
					ParameterSet[i]++;
					break;
				case 14:
					sscanf(line,"%s %lf",word,&recalc_frac);
					if (ParameterSet[i]>0)
						fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
					else
						NParametersSet++;
					#ifdef VERB
					fprintf(stderr,"\t%s: %lf\n",word,recalc_frac);
					#endif
					ParameterSet[i]++;
					break;
				case 15:
					sscanf(line,"%s %f",word,&Dmax);
					if (ParameterSet[i]>0)
						fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
					else
						NParametersSet++;
					#ifdef VERB
					fprintf(stderr,"\t%s: %f\n",word,Dmax);
					#endif
					ParameterSet[i]++;
					break;
				case 16:
					sscanf(line,"%s %f",word,&Dmin);
					if (ParameterSet[i]>0)
						fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
					else
						NParametersSet++;
					#ifdef VERB
					fprintf(stderr,"\t%s: %f\n",word,Dmin);
					#endif
					ParameterSet[i]++;
					break;

				case 17:
					sscanf(line,"%s %f",word,&alpha_ratio);
					if (ParameterSet[i]>0)
						fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
					else
						NParametersSet++;
					#ifdef VERB
					fprintf(stderr,"\t%s: %f\n",word,alpha_ratio);
					#endif
					ParameterSet[i]++;
					break;
				case 18:
					sscanf(line,"%s %f",word,&alpha_ratio_1);
					if (ParameterSet[i]>0)
						fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
					else
						NParametersSet++;
					#ifdef VERB
					fprintf(stderr,"\t%s: %f\n",word,alpha_ratio_1);
					#endif
					ParameterSet[i]++;
					break;
				case 19:
					sscanf(line,"%s %f",word,&best_alpha);
					if (ParameterSet[i]>0)
						fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
					else
						NParametersSet++;
					#ifdef VERB
					fprintf(stderr,"\t%s: %f\n",word,best_alpha);
					#endif
					ParameterSet[i]++;
					break;
				case 20:
					sscanf(line,"%s %d",word,&num_alpha);
					if (ParameterSet[i]>0)
						fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
					else
						NParametersSet++;
					#ifdef VERB
					fprintf(stderr,"\t%s: %d\n",word,num_alpha);
					#endif
					ParameterSet[i]++;
					break;
				case 21:
					sscanf(line,"%s %d",word,&ntrials);
					if (ParameterSet[i]>0)
						fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
					else
						NParametersSet++;
					#ifdef VERB
					fprintf(stderr,"\t%s: %d\n",word,ntrials);
					#endif
					ParameterSet[i]++;
					break;
				case 22:
					sscanf(line,"%s %d",word,&nthreads);
					if (ParameterSet[i]>0)
						fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
					else
						NParametersSet++;
					#ifdef VERB
					fprintf(stderr,"\t%s: %d\n",word,nthreads);
					#endif
					ParameterSet[i]++;
					break;
				case 23:
					sscanf(line,"%s %d",word,&nr);
					if (ParameterSet[i]>0)
						fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
					else
						NParametersSet++;
					#ifdef VERB
					fprintf(stderr,"\t%s: %d\n",word,nr);
					#endif
					ParameterSet[i]++;
					break;
				case 24:
					sscanf(line,"%s %lf",word,&minr);
					if (ParameterSet[i]>0)
						fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
					else
						NParametersSet++;
					#ifdef VERB
					fprintf(stderr,"\t%s: %lf\n",word,minr);
					#endif
					ParameterSet[i]++;
					break;
				case 25:
					sscanf(line,"%s %lf",word,&maxr);
					if (ParameterSet[i]>0)
						fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
					else
						NParametersSet++;
					#ifdef VERB
					fprintf(stderr,"\t%s: %lf\n",word,maxr);
					#endif
					ParameterSet[i]++;
					break;
				default:
					fprintf(stderr,"WARNING: Unknown parameter %s\n",word);
					break;
			}//switch
		  }// if (non-empty line)
		}//if (no comment line)
	}//while	
	fclose(f);
	if(NParametersSet!=NParam){
		for (i=0;i<NParam;i++){
			if (ParameterSet[i]==0)
				fprintf(stderr,"ERROR: Parameter %s not set in input file\n",ParameterList[i]);
		}
		return -1;
	}
/*
	if ((f = fopen(alphaFile,"r"))==NULL){
		fprintf(stderr,"Could not open input file %s\n",alphaFile);
		return -1;
	}
        while(fgets(line,LINELENGTH,f)!=NULL)
                if (line[0] != '#')
			Nalpha++;	
	fclose(f);
	
	alpha = (double *) calloc(Nalpha,sizeof(double));
	Malpha = (double *) calloc(Nalpha,sizeof(double));

	if ((f = fopen(alphaFile,"r"))==NULL){
		fprintf(stderr,"Could not open alpha file %s\n",alphaFile);
		return -1;
	}

	i=0;
        while(fgets(line,LINELENGTH,f)!=NULL){
                if (line[0] != '#'){
			fgets(line,LINELENGTH,f);
			sscanf(line,"%f %f",&x,&y);
			alpha[i]=x;
			Malpha[i]=y;
			i++;
		}
	}
	fclose(f);
*/
	return 0;
}



