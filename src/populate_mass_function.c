#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>

#include "populate_mass_function.h"

#define LINELENGTH 1024

/*=============================================================================
 *                             PROTOTYPES
 *=============================================================================*/


void get_cubic_splines(double *, double *, int , double *);
double spline_inter(double *, double *, int , double *, double );
int read_mass_function(char *, double **, double **, long *,double,long *);
//static int compare_double(const void * , const void * );
static int compare_float(const void * , const void * );




/*=============================================================================
 *                              populate_mass_funtcion()
 *=============================================================================*/
//Reads the cumulative mass function file, and creates and array of halo_masses which populate that function.
//The return value (if successful) is the number of haloes

long populate_mass_function(char *filename, double Mmin, double Lbox, float **halo_masses,long seed, int nthreads){
	long Npoints,Nshort,Nhalos,i,j;
	double *mass_arr,*mass_inv_arr, *dens_arr, *dens_inv_arr, *y2, d_rand;
	double dens_max;

	if (nthreads<1){
		nthreads = omp_get_max_threads();
	}

	//Reading cumulative mass function
	fprintf(stderr,"\tReading Mass Function... \n");		
	if (read_mass_function(filename,&mass_arr,&dens_arr,&Npoints,Mmin,&Nshort)!=1){
		fprintf(stderr,"ERROR: could not read %s\n",filename);
		return -1;	
	}
	Nshort+=2; //Add 2 points to avoid border effects
	fprintf(stderr,"\t... read\n");		
	srand(seed);
	#ifdef VERB
	fprintf(stderr,"\tSeed used: %ld\n",seed);
	#endif

	//prepare splines for Cumulative density as a funtion of M: n(M)
	y2 = (double *) calloc(Npoints,sizeof(double));

	get_cubic_splines(mass_arr,dens_arr,Npoints,y2);
	
	#ifndef NDENS
	 dens_max = spline_inter(mass_arr,dens_arr,Npoints,y2,Mmin);
	#else
	 dens_max = Mmin;
	#endif

	Nhalos = (long)(dens_max*Lbox*Lbox*Lbox+0.5);	

	fprintf(stderr,"\tNumber of halos: %ld\n",Nhalos);	



	fprintf(stderr,"\tInverting Mass Function... \n");	
	//prepare splines for inverse function M(n)
	mass_inv_arr = (double *) calloc(Nshort,sizeof(double));
	dens_inv_arr = (double *) calloc(Nshort,sizeof(double)); 
	for (i=0;i<Nshort;i++){
		j = Npoints-1-i;
		mass_inv_arr[i]=mass_arr[j];
		dens_inv_arr[i]=dens_arr[j];
	}
	fprintf(stderr,"\t...done!\n");
	*halo_masses = (float *) calloc(Nhalos,sizeof(float));
	get_cubic_splines(dens_inv_arr,mass_inv_arr,Nshort,y2);



	//Generate masses
	fprintf(stderr,"\tGenerating halo masses...\n");	
	#pragma omp parallel for num_threads(nthreads) private(i,d_rand) shared(dens_max,halo_masses,Nhalos,dens_inv_arr,mass_inv_arr,Nshort,y2) default(none)
	for (i=0;i<Nhalos;i++){
		d_rand = (double) rand()*dens_max/RAND_MAX;
		(*halo_masses)[i] = (float) spline_inter(dens_inv_arr,mass_inv_arr,Nshort,y2,d_rand);
	}
	fprintf(stderr,"\t...sorting them...\n");	
	
	qsort(*halo_masses, Nhalos, sizeof(*halo_masses[0]), compare_float);
	
	fprintf(stderr,"\t...done\n");

	return Nhalos;
}


/*=============================================================================
 *                                read_mass_function()
 *=============================================================================*/

int read_mass_function(char *fname, double **x, double **y, long *N, double xmin, long *Nred){
	long Nlines=0, N_gtmin=0;
	FILE *f;
	char line[LINELENGTH];
	double X,Y;

 	if ((f = fopen(fname, "r"))==NULL) {
                perror("fopen:");
                fprintf(stderr,"Error while reading the file %s\n",fname);
                return -1;
        }
	while(fgets(line,LINELENGTH,f)!=NULL) {
                if (line[0] != '#') {
			Nlines++;	
		}	
        }
        fclose(f);

	(*x) = (double *) calloc (Nlines,sizeof(double));
	(*y) = (double *) calloc (Nlines,sizeof(double));

	Nlines=0;
	f = fopen (fname,"r");
	while(fgets(line,LINELENGTH,f)!=NULL) {
                if (line[0] != '#') {
			sscanf(line,"%lf %lf",&X,&Y);	
			(*x)[Nlines]=X;
			(*y)[Nlines]=Y;
			Nlines++;
#ifndef NDENS
			if (X>xmin)
				N_gtmin++;
#else
			if (Y<xmin)
				N_gtmin++;
#endif
        }
        }
        fclose(f);

    if (N_gtmin==Nlines){
    	fprintf(stderr,"Mass Function Doesn't Reach Lowest Mass Needed, %e, %e\n", xmin, (*x)[0]);
    	exit(1);
    }
	*N = Nlines;
	*Nred = N_gtmin;
	return 1;
}

/*=============================================================================
 *                              SPLINES functions
 *=============================================================================*/
// get_cubic_splines(): computes the  coefficients (y2) of the cubic splines for a function y(x) sampled with N points


void get_cubic_splines(double *x, double *y, int N, double *y2){
        int i;
        double *u,s,p;

        u = (double *) calloc(N,sizeof(double));

        y2[0]=0.0; //Natural spline
        u[0]=0.0;
        //fprintf(stderr,"%f\n",y2[0]);


        for (i=1;i<N-1;i++){
                s = (x[i]-x[i-1]) / (x[i+1]-x[i-1]);
                p = s * y2[i-1] + 2.0;
                y2[i] = (s - 1.0) / p;
                u[i] = (y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
                u[i] = (6.0*u[i]/(x[i+1]-x[i-1])-s*u[i-1])/p;
                //fprintf(stderr,"%f %f %f\n",s,p,y2[i]);
        }

        y2[N-1]=0.0; //Natural spline

        for (i=N-2;i>=0;i--){
                y2[i] = y2[i]*y2[i+1]+u[i];
                //fprintf(stderr,"%f\n",y2[i]);

        }
        free(u);
}


// spline_inter(): given a function y(x) sampled with N points and the splines coefficients (y2),
// it returns the yi value at an xi point 
double spline_inter(double *x, double *y, int N, double *y2, double xi){

        double h,b,a;
        int i_low=0,i_high=N-1,i_mid;

        while(i_high-i_low>1){
                i_mid=(i_high+i_low)/2;
                if (x[i_mid]>xi)
                        i_high = i_mid;
                else
                        i_low = i_mid;
        }

        h = x[i_high] - x[i_low];
        a = (x[i_high]-xi)/h;
        b = (xi-x[i_low])/h;
        return (a*y[i_low]+b*y[i_high]+((a*a*a-a)*y2[i_low]+(b*b*b-b)*y2[i_high])*h*h/6.0);
}
/*
static int compare_double(const void * a, const void * b)
{
  if (*(double*)a > *(double*)b) return -1;
  else if (*(double*)a < *(double*)b) return 1;
  else return 0;  
}*/

/*=============================================================================
 *                              COMPARE
 *=============================================================================*/
//compares 2 float variables in the way required by qsort()

static int compare_float(const void * a, const void * b)
{
  if (*(float*)a > *(float*)b) return -1;
  else if (*(float*)a < *(float*)b) return 1;
  else return 0;  
}
