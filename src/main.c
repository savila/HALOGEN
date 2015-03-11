/*********************************************************************************** 
 HALOGEN-2LPT
************************************************************************************ 
This code is a modification of 2LPTic (http://arxiv.org/abs/astro-ph/9711187)
in which we incorporated HALOGEN (http://arxiv.org/abs/1412.5228)
The original code 2LPTic was developped by 
Roman Scoccimarro and Sebastian Pueblas (http://cosmo.nyu.edu/roman/2LPT/)
HALOGEN has been developped by Santiago Avila and Steven Murray 
(https://github.com/savila/HALOGEN).
************************************************************************************ */


#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <drfftw_mpi.h>
#include <mpi.h>
#include <gsl/gsl_rng.h>
// HALOGEN libs //
#include <time.h>
#include <string.h>
#include "read_snapshot.h"
#include "populate_mass_function.h"
#include "place_halos.h"


#include "allvars.h"
#include "proto.h"
/* *********************************************** */

#define rho_crit (27.755e10)
#define LINELENGTH 256
#define NParam 10

int ParameterSet[NParam];
//char ParameterList[NParam][32];
int NParametersSet;
long seed;
int nthreads;

char ParameterList[NParam][32] = {"MassFunctionFile",
                "OutputFile","NCellsLin","alphaFile","rho_ref","Overdensity","Mmin",
                "Seed","recalc_frac",
                "nthreads"};



char Snapshot[LINELENGTH];
char OutputFile[LINELENGTH],MassFunctionFile[LINELENGTH], alphaFile[LINELENGTH];
int format;
float recalc_frac;

int Nlin;
float Mmin;
float OVD;
char rho_ref[8];

int Nalpha=0;
double *alpha_vec,*fvel;
double *Malpha;


float LUNIT, MUNIT;
int SWP, LGADGET, DGADGET;

int read_input_file(char *);

int write_halogen_cat(char *, float *, float *, float *, float *, float *, float *, float *, float *, long);

/*   ********************************************************************************************** */







#define ASSERT_ALLOC(cond) {                                                                                  \
   if(cond)                                                                                                   \
    {                                                                                                         \
      if(ThisTask == 0)                                                                                       \
	fprintf(stderr,"\tallocated %g Mbyte on Task %d\n", bytes / (1024.0 * 1024.0), ThisTask);                     \
    }                                                                                                         \
  else                                                                                                        \
    {                                                                                                         \
      fprintf(stderr,"failed to allocate %g Mbyte on Task %d\n", bytes / (1024.0 * 1024.0), ThisTask);                \
      fprintf(stderr,"bailing out.\n");                                                                               \
      FatalError(1);                                                                                          \
    }                                                                                                         \
}





int main(int argc, char **argv) {

  	MPI_Init(&argc, &argv);
  	MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
  	MPI_Comm_size(MPI_COMM_WORLD, &NTask);

  	float Lbox, mpart, *x, *y, *z, *vx,*vy,*vz,*hx, *hy, *hz, *hvx,*hvy,*hvz,*hR, om_m;
        char tlpt_filename[256],halogen_filename[256];
        long Npart, Nhalos, **ListOfParticles, *NPartPerCell;
        float *HaloMass, rho;

        int number_exclusion=0;


  	if(ThisTask == 0) 
  	{
       		fprintf(stderr,"\n*******************************************************************\n");
        	fprintf(stderr,"**                                                               **\n");
        	fprintf(stderr,"**            =        HALOGEN-2LPT v1.0         =               **\n");
        	fprintf(stderr,"**                                                               **\n");
        	fprintf(stderr,"**                                                               **\n");
        	fprintf(stderr,"**                                            let there be dark  **\n");
        	fprintf(stderr,"*******************************************************************\n\n");
		#ifdef VERB
        	fprintf(stderr,"#def VERB\n");
		#endif
		#ifdef DEBUG
        	fprintf(stderr,"#def DEBUG \n");
		#endif
		#ifdef ULTRADEBUG
        	fprintf(stderr,"#def ULTRADEBUG \n");
		#endif
		#ifdef FULL_EXCLUSION
        	fprintf(stderr,"#def FULL_EXCLUSION\n");
		#endif
		#ifdef HALF_EXCLUSION
        	fprintf(stderr,"#def HALF_EXCLUSION\n");
		#endif
		#ifdef NO_EXCLUSION
	        fprintf(stderr,"#def NO_EXCLUSION\n");
		#endif
		#ifdef RANKED
        	fprintf(stderr,"#def RANKED\n");
		#endif
		#ifdef NDENS
	        fprintf(stderr,"#def NDENS\n");
		#endif
	//FLAGS to be unified
  	}

	#ifdef FULL_EXCLUSION
        number_exclusion++;
	#endif
	#ifdef HALF_EXCLUSION
	number_exclusion++;
	#endif
	#ifdef NO_EXCLUSION
        number_exclusion++;
	#endif

	if (number_exclusion!=1){
      		if(ThisTask == 0)
        		fprintf(stderr,"ERROR: You must select one and only one exclusion criterion in Makefile.def\n");
        		fprintf(stderr,"Task %d, number of exclusion criterion:%d \n",ThisTask,number_exclusion);
  			MPI_Finalize();		/* clean up & finalize MPI */
        		exit(0);
	}

  	if(argc < 3)
    	{
      		if(ThisTask == 0)
		{
	  		fprintf(stdout, "\nParameters are missing.\n");
	  		fprintf(stdout, "Usage: %s <2LPTParameterFile> <Halogen_input_file>\n\n",argv[0]);
		}
      		MPI_Finalize();
      		exit(0);
    	}
        strcpy(tlpt_filename,argv[1]);
        strcpy(halogen_filename,argv[2]);

	if(ThisTask == 0) fprintf(stderr,"\nReading 2LPT input file...\n");
	read_parameterfile(tlpt_filename);
	if(ThisTask == 0) fprintf(stderr,"...2LPT input file read\n\n");
        

	if(ThisTask == 0) fprintf(stderr,"\nReading HALOGEN input file...\n");
        if (read_input_file(halogen_filename)<0){
                MPI_Finalize();
                exit(0);
	}
        if(ThisTask == 0) fprintf(stderr,"... HALOGEN file read correctly!\n\n");

	if(ThisTask == 0) fprintf(stderr,"Initializing 2LPT...\n");
  	set_units();
	initialize_powerspectrum();
  	initialize_ffts();
  	read_glass(GlassFile);
	if(ThisTask == 0) fprintf(stderr,"... initialization done!\n\n");

	if(ThisTask == 0) fprintf(stderr,"Computing displacement fields...\n");
  	displacement_fields();
	if(ThisTask == 0) fprintf(stderr,"...done with displacement fields\n\n");


	 MPI_Barrier(MPI_COMM_WORLD);
	//done with 2LPT, halogen starting

 
        if(ThisTask == 0) fprintf(stderr,"Distributing particles in grid... \n");
        if (distribute_part(Nlin,nthreads,&x, &y, &z, &vx, &vy, &vz, &Npart, &Lbox, &om_m,&ListOfParticles,&NPartPerCell)==0) {
                if(ThisTask == 0) 
			fprintf(stderr,"...particles correctly distributed!\n\n");
	}
        else {
                fprintf(stderr,"ERROR: Something went wrong distributing the particles\n");
                MPI_Finalize();
                exit(0);
        }


  if(ThisTask == 0) 
  {

	
  	#ifdef DEBUG
	int i,j,k,lin_ijk;	
	for (i=0; i<Nlin;i++)
	for (j=0; j<Nlin;j++)
	for (k=0; k<Nlin;k++)	
	{
		lin_ijk = k+j*Nlin+i*Nlin*Nlin;
		if ((i<3 || i>Nlin-4) && (j<3 || j>Nlin-4) && (k<3 || k>Nlin-4))
			fprintf(stderr,"Npart[%d][%d][%d]=%d\n",i,j,k,NPartPerCell[lin_ijk]);
		if ((float)rand()/RAND_MAX<5.0e-8 && ListOfParticles[lin_ijk]>0){
			fprintf(stderr,"(%f,%f,%f) -> [%d,%d,%d]\n",x[ListOfParticles[lin_ijk][0]],y[ListOfParticles[lin_ijk][0]],z[ListOfParticles[lin_ijk][0]],i,j,k);
		}
	}
	#endif
   	mpart = rho_crit * om_m * Lbox*Lbox*Lbox /Npart;

        #ifdef VERB
        fprintf(stderr,"\n\tCheck: Npart=%ld, mpart=%e, Lbox=%f\n",Npart,mpart,Lbox);
        fprintf(stderr,"\tx[0]= %f, y[0]= %f, z[0]= %f\n",(x)[0],(y)[0],(z)[0]);
        fprintf(stderr,"\t      ...\n");
        fprintf(stderr,"\tx[%ld]= %f, y[%ld]= %f, z[%ld]= %f\n\n",Npart-1,(x)[Npart-1],Npart-1,(y)[Npart-1],Npart-1,(z)[Npart-1]);
        #endif

        if (seed<0){
                seed = time(NULL);
                fprintf(stderr,"Seed used: %ld\n",seed);
        }

        //Generate the halo masses from the mass function
        fprintf(stderr,"Generating Halo Masses...\n");
        Nhalos = populate_mass_function(MassFunctionFile,Mmin,Lbox,&HaloMass,seed);
        if (Nhalos<0)
                fprintf(stderr,"error: Couldnt create HaloMass array\n");
        fprintf(stderr,"...Halo Masses Generated\n\n");

        //Allocalte memory for the halo XYZR vector
        hx = (float *) calloc(Nhalos,sizeof(float));
        hy = (float *) calloc(Nhalos,sizeof(float));
        hz = (float *) calloc(Nhalos,sizeof(float));
        hvx = (float *) calloc(Nhalos,sizeof(float));
        hvy = (float *) calloc(Nhalos,sizeof(float));
        hvz = (float *) calloc(Nhalos,sizeof(float));
        hR = (float *) calloc(Nhalos,sizeof(float));

        //density at the boundary of a halo
        if (strcmp(rho_ref,"crit")==0)
                rho = OVD*rho_crit;
        else
                rho = OVD*rho_crit*om_m;

        //place the halos
        fprintf(stderr,"Placing halos down...\n");

        // Check the M-alpha vector against produced halos.
        if( Malpha[Nalpha-1] > HaloMass[Nhalos-1]){
                if (Malpha[Nalpha-1] < 1.2*HaloMass[Nhalos-1]) {
                        Malpha[Nalpha-1] = HaloMass[Nhalos-1]/1.01;
                }
                else{
                        fprintf(stderr,"ERROR: lowest M(alpha) larger than the lowest halo Mass: %e < %e",Malpha[Nalpha-1],HaloMass[Nhalos-1]);
                        exit(0);
                }
        }

        if (place_halos(Nhalos,HaloMass, Nlin, Npart, x, y, z, vx,vy,vz,Lbox, rho,seed,mpart, nthreads,alpha_vec, fvel, Malpha, Nalpha,recalc_frac,hx, hy, hz, hvx,hvy,hvz, hR,ListOfParticles,NPartPerCell)==0){
                fprintf(stderr,"...halos placed correctly\n");
        }
        else {
                fprintf(stderr,"ERROR: Problem placing halos\n");
  		MPI_Finalize();
		exit(0);
        }
        fprintf(stderr,"\n");

        //writting output       
        fprintf(stderr,"Writing Halo catalogue...\n");
        write_halogen_cat(OutputFile,hx,hy,hz,hvx,hvy,hvz,HaloMass,hR,Nhalos);
        fprintf(stderr,"...halo catalogue written in %s\n",OutputFile);

        free(hx);free(hy);free(hz);free(hR);
        free(alpha_vec); free(Malpha);


        fprintf(stderr,"\n*******************************************************************\n");
        fprintf(stderr,"**                        ... and there were dark matter haloes  **\n");
        fprintf(stderr,"*******************************************************************\n\n");

 

  }

  MPI_Barrier(MPI_COMM_WORLD);

  if(NumPart){
    free(partX);
    free(partY);
    free(partZ);
    free(partVX);
    free(partVY);
    free(partVZ);
  }
 

  free_ffts();
  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Finalize();		/* clean up & finalize MPI */
  exit(0);
}





void displacement_fields(void)
{
  MPI_Request request;
  MPI_Status status;
  gsl_rng *random_generator;
  int i, j, k, ii, jj, kk, axes;
  int n;
  int sendTask, recvTask;
  double fac, vel_prefac, vel_prefac2;
  double kvec[3], kmag, kmag2, p_of_k;
  double delta, phase, ampl, hubble_a;
  double u, v, w;
  double f1, f2, f3, f4, f5, f6, f7, f8;
  double dis, dis2, maxdisp, max_disp_glob;
  unsigned int *seedtable;

  unsigned int bytes, nmesh3;
  int coord;
  fftw_complex *(cdisp[3]), *(cdisp2[3]) ; /* ZA and 2nd order displacements */
  fftw_real *(disp[3]), *(disp2[3]) ;

  fftw_complex *(cdigrad[6]);
  fftw_real *(digrad[6]);


#ifdef CORRECT_CIC
  double fx, fy, fz, ff, smth;
#endif

#ifdef VERB
  if(ThisTask == 0)
    {
      printf("\tstart computing displacement fields...\n");
      fflush(stdout);
    }
#endif

  hubble_a =
    Hubble * sqrt(Omega / pow(InitTime, 3) + (1 - Omega - OmegaLambda) / pow(InitTime, 2) + OmegaLambda);

  vel_prefac = InitTime * hubble_a * F_Omega(InitTime);
  vel_prefac2 = InitTime * hubble_a * F2_Omega(InitTime);

  vel_prefac /= sqrt(InitTime);	/* converts to Gadget velocity */
  vel_prefac2 /= sqrt(InitTime);

#ifdef VERB
  if(ThisTask == 0)
    printf("\tvel_prefac= %g, vel_prefac2= %g,  hubble_a=%g fom=%g \n", vel_prefac, vel_prefac2, 
                                                                      hubble_a, F_Omega(InitTime));
#endif

  fac = pow(2 * PI / Box, 1.5);

  maxdisp = 0;

  random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);

  gsl_rng_set(random_generator, Seed);

  if(!(seedtable = malloc(Nmesh * Nmesh * sizeof(unsigned int))))
    FatalError(4);

  for(i = 0; i < Nmesh / 2; i++)
    {
      for(j = 0; j < i; j++)
	seedtable[i * Nmesh + j] = 0x7fffffff * gsl_rng_uniform(random_generator);

      for(j = 0; j < i + 1; j++)
	seedtable[j * Nmesh + i] = 0x7fffffff * gsl_rng_uniform(random_generator);

      for(j = 0; j < i; j++)
	seedtable[(Nmesh - 1 - i) * Nmesh + j] = 0x7fffffff * gsl_rng_uniform(random_generator);

      for(j = 0; j < i + 1; j++)
	seedtable[(Nmesh - 1 - j) * Nmesh + i] = 0x7fffffff * gsl_rng_uniform(random_generator);

      for(j = 0; j < i; j++)
	seedtable[i * Nmesh + (Nmesh - 1 - j)] = 0x7fffffff * gsl_rng_uniform(random_generator);

      for(j = 0; j < i + 1; j++)
	seedtable[j * Nmesh + (Nmesh - 1 - i)] = 0x7fffffff * gsl_rng_uniform(random_generator);

      for(j = 0; j < i; j++)
	seedtable[(Nmesh - 1 - i) * Nmesh + (Nmesh - 1 - j)] = 0x7fffffff * gsl_rng_uniform(random_generator);

      for(j = 0; j < i + 1; j++)
	seedtable[(Nmesh - 1 - j) * Nmesh + (Nmesh - 1 - i)] = 0x7fffffff * gsl_rng_uniform(random_generator);
    }



  for(axes=0,bytes=0; axes < 3; axes++)
    {
      cdisp[axes] = (fftw_complex *) malloc(bytes += sizeof(fftw_real) * TotalSizePlusAdditional);
      disp[axes] = (fftw_real *) cdisp[axes];
    }

#ifdef VERB
  ASSERT_ALLOC(cdisp[0] && cdisp[1] && cdisp[2]);
#endif

#if defined(MULTICOMPONENTGLASSFILE) && defined(DIFFERENT_TRANSFER_FUNC)
  for(Type = MinType; Type <= MaxType; Type++)
#endif
    {
      if(ThisTask == 0)
	{
	#ifdef VERB
	  printf("\tstarting axes=%d...\n", axes);
	  fflush(stdout);
	#endif
	}

      /* first, clean the array */
      for(i = 0; i < Local_nx; i++)
	for(j = 0; j < Nmesh; j++)
	  for(k = 0; k <= Nmesh / 2; k++)
	    for(axes = 0; axes < 3; axes++)
	      {
		cdisp[axes][(i * Nmesh + j) * (Nmesh / 2 + 1) + k].re = 0;
		cdisp[axes][(i * Nmesh + j) * (Nmesh / 2 + 1) + k].im = 0;
	      }

      for(i = 0; i < Nmesh; i++)
	{
	  ii = Nmesh - i;
	  if(ii == Nmesh)
	    ii = 0;
	  if((i >= Local_x_start && i < (Local_x_start + Local_nx)) ||
	     (ii >= Local_x_start && ii < (Local_x_start + Local_nx)))
	    {
	      for(j = 0; j < Nmesh; j++)
		{
		  gsl_rng_set(random_generator, seedtable[i * Nmesh + j]);
		  
		  for(k = 0; k < Nmesh / 2; k++)
		    {
		      phase = gsl_rng_uniform(random_generator) * 2 * PI;
		      do
			ampl = gsl_rng_uniform(random_generator);
		      while(ampl == 0);
		      
		      if(i == Nmesh / 2 || j == Nmesh / 2 || k == Nmesh / 2)
			continue;
		      if(i == 0 && j == 0 && k == 0)
			continue;
		      
		      if(i < Nmesh / 2)
			kvec[0] = i * 2 * PI / Box;
		      else
			kvec[0] = -(Nmesh - i) * 2 * PI / Box;
		      
		      if(j < Nmesh / 2)
			kvec[1] = j * 2 * PI / Box;
		      else
			kvec[1] = -(Nmesh - j) * 2 * PI / Box;
		      
		      if(k < Nmesh / 2)
			kvec[2] = k * 2 * PI / Box;
		      else
			kvec[2] = -(Nmesh - k) * 2 * PI / Box;
		      
		      kmag2 = kvec[0] * kvec[0] + kvec[1] * kvec[1] + kvec[2] * kvec[2];
		      kmag = sqrt(kmag2);
		      
		      if(SphereMode == 1)
			{
			  if(kmag * Box / (2 * PI) > Nsample / 2)	/* select a sphere in k-space */
			    continue;
			}
		      else
			{
			  if(fabs(kvec[0]) * Box / (2 * PI) > Nsample / 2)
			    continue;
			  if(fabs(kvec[1]) * Box / (2 * PI) > Nsample / 2)
			    continue;
			  if(fabs(kvec[2]) * Box / (2 * PI) > Nsample / 2)
			    continue;
			}
		      
		      p_of_k = PowerSpec(kmag);
		      
		      p_of_k *= -log(ampl);
		      
		      delta = fac * sqrt(p_of_k) / Dplus;	/* scale back to starting redshift */
		      
		      if(k > 0)
			{
			  if(i >= Local_x_start && i < (Local_x_start + Local_nx))
			    for(axes = 0; axes < 3; axes++)
			      {
				cdisp[axes][((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k].re =
				  -kvec[axes] / kmag2 * delta * sin(phase);
				cdisp[axes][((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k].im =
				  kvec[axes] / kmag2 * delta * cos(phase);
			      }
			}
		      else	/* k=0 plane needs special treatment */
			{
			  if(i == 0)
			    {
			      if(j >= Nmesh / 2)
				continue;
			      else
				{
				  if(i >= Local_x_start && i < (Local_x_start + Local_nx))
				    {
				      jj = Nmesh - j;	/* note: j!=0 surely holds at this point */
				      
				      for(axes = 0; axes < 3; axes++)
					{
					  cdisp[axes][((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k].re =
					    -kvec[axes] / kmag2 * delta * sin(phase);
					  cdisp[axes][((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k].im =
					    kvec[axes] / kmag2 * delta * cos(phase);
					  
					  cdisp[axes][((i - Local_x_start) * Nmesh + jj) * (Nmesh / 2 + 1) + k].re =
					    -kvec[axes] / kmag2 * delta * sin(phase);
					  cdisp[axes][((i - Local_x_start) * Nmesh + jj) * (Nmesh / 2 + 1) + k].im =
					    -kvec[axes] / kmag2 * delta * cos(phase);
					}
				    }
				}
			    }
			  else	/* here comes i!=0 : conjugate can be on other processor! */
			    {
			      if(i >= Nmesh / 2)
				continue;
			      else
				{
				  ii = Nmesh - i;
				  if(ii == Nmesh)
				    ii = 0;
				  jj = Nmesh - j;
				  if(jj == Nmesh)
				    jj = 0;
				  
				  if(i >= Local_x_start && i < (Local_x_start + Local_nx))
				    for(axes = 0; axes < 3; axes++)
				      {
					cdisp[axes][((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k].re =
					  -kvec[axes] / kmag2 * delta * sin(phase);
					cdisp[axes][((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k].im =
					  kvec[axes] / kmag2 * delta * cos(phase);
				      }
				  
				  if(ii >= Local_x_start && ii < (Local_x_start + Local_nx))
				    for(axes = 0; axes < 3; axes++)
				      {
					cdisp[axes][((ii - Local_x_start) * Nmesh + jj) * (Nmesh / 2 + 1) +
					      k].re = -kvec[axes] / kmag2 * delta * sin(phase);
					cdisp[axes][((ii - Local_x_start) * Nmesh + jj) * (Nmesh / 2 + 1) +
					      k].im = -kvec[axes] / kmag2 * delta * cos(phase);
				      }
				}
			    }
			}
		    }
		}
	    }
	}
      


      /* At this point, cdisp[axes] contains the complex Zeldovich displacement */
     
       if(ThisTask == 0) printf("\tDone Zeldovich.\n");
      
      /* Compute displacement gradient */

      for(i = 0; i < 6; i++)
	{
	  cdigrad[i] = (fftw_complex *) malloc(bytes = sizeof(fftw_real) * TotalSizePlusAdditional);
	  digrad[i] = (fftw_real *) cdigrad[i];

	  #ifdef DEBUG
	  ASSERT_ALLOC(cdigrad[i]);
	  #endif
	}
      
      for(i = 0; i < Local_nx; i++)
	for(j = 0; j < Nmesh; j++)
	  for(k = 0; k <= Nmesh / 2; k++)
	    {
	      coord = (i * Nmesh + j) * (Nmesh / 2 + 1) + k;
	      if((i + Local_x_start) < Nmesh / 2)
		kvec[0] = (i + Local_x_start) * 2 * PI / Box;
	      else
		kvec[0] = -(Nmesh - (i + Local_x_start)) * 2 * PI / Box;
	      
	      if(j < Nmesh / 2)
		kvec[1] = j * 2 * PI / Box;
	      else
		kvec[1] = -(Nmesh - j) * 2 * PI / Box;
	      
	      if(k < Nmesh / 2)
		kvec[2] = k * 2 * PI / Box;
	      else
		kvec[2] = -(Nmesh - k) * 2 * PI / Box;
	      
	      /* Derivatives of ZA displacement  */
	      /* d(dis_i)/d(q_j)  -> sqrt(-1) k_j dis_i */
	      cdigrad[0][coord].re = -cdisp[0][coord].im * kvec[0]; /* disp0,0 */
	      cdigrad[0][coord].im = cdisp[0][coord].re * kvec[0];

	      cdigrad[1][coord].re = -cdisp[0][coord].im * kvec[1]; /* disp0,1 */
	      cdigrad[1][coord].im = cdisp[0][coord].re * kvec[1];

	      cdigrad[2][coord].re = -cdisp[0][coord].im * kvec[2]; /* disp0,2 */
	      cdigrad[2][coord].im = cdisp[0][coord].re * kvec[2];
	      
	      cdigrad[3][coord].re = -cdisp[1][coord].im * kvec[1]; /* disp1,1 */
	      cdigrad[3][coord].im = cdisp[1][coord].re * kvec[1];

	      cdigrad[4][coord].re = -cdisp[1][coord].im * kvec[2]; /* disp1,2 */
	      cdigrad[4][coord].im = cdisp[1][coord].re * kvec[2];

	      cdigrad[5][coord].re = -cdisp[2][coord].im * kvec[2]; /* disp2,2 */
	      cdigrad[5][coord].im = cdisp[2][coord].re * kvec[2];
	    }


      if(ThisTask == 0) fprintf(stderr,"\tFourier transforming displacement gradient...\n");
      for(i = 0; i < 6; i++) rfftwnd_mpi(Inverse_plan, 1, digrad[i], Workspace, FFTW_NORMAL_ORDER);
      if(ThisTask == 0) fprintf(stderr,"\tDone.\n");

      /* Compute second order source and store it in digrad[3]*/

      for(i = 0; i < Local_nx; i++)
	for(j = 0; j < Nmesh; j++)
	  for(k = 0; k < Nmesh; k++)
	    {
	      coord = (i * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + k;

	      digrad[3][coord] =

		digrad[0][coord]*(digrad[3][coord]+digrad[5][coord])+digrad[3][coord]*digrad[5][coord]
                -digrad[1][coord]*digrad[1][coord]-digrad[2][coord]*digrad[2][coord]-digrad[4][coord]*digrad[4][coord];
	    }

      if(ThisTask == 0) fprintf(stderr,"\tFourier transforming second order source...");
      rfftwnd_mpi(Forward_plan, 1, digrad[3], Workspace, FFTW_NORMAL_ORDER);
      if(ThisTask == 0) fprintf(stderr,"\tDone.\n");
      
      /* The memory allocated for cdigrad[0], [1], and [2] will be used for 2nd order displacements */
      /* Freeing the rest. cdigrad[3] still has 2nd order displacement source, free later */

      for(axes = 0; axes < 3; axes++) 
	{
	  cdisp2[axes] = cdigrad[axes]; 
	  disp2[axes] = (fftw_real *) cdisp2[axes];
	}

      free(cdigrad[4]); free(cdigrad[5]); 

      /* Solve Poisson eq. and calculate 2nd order displacements */

      for(i = 0; i < Local_nx; i++)
	for(j = 0; j < Nmesh; j++)
	  for(k = 0; k <= Nmesh / 2; k++)
	    {
	      coord = (i * Nmesh + j) * (Nmesh / 2 + 1) + k;
	      if((i + Local_x_start) < Nmesh / 2)
		kvec[0] = (i + Local_x_start) * 2 * PI / Box;
	      else
		kvec[0] = -(Nmesh - (i + Local_x_start)) * 2 * PI / Box;
	      
	      if(j < Nmesh / 2)
		kvec[1] = j * 2 * PI / Box;
	      else
		kvec[1] = -(Nmesh - j) * 2 * PI / Box;
	      
	      if(k < Nmesh / 2)
		kvec[2] = k * 2 * PI / Box;
	      else
		kvec[2] = -(Nmesh - k) * 2 * PI / Box;

	      kmag2 = kvec[0] * kvec[0] + kvec[1] * kvec[1] + kvec[2] * kvec[2];
#ifdef CORRECT_CIC
	      /* calculate smooth factor for deconvolution of CIC interpolation */
	      fx = fy = fz = 1;
	      if(kvec[0] != 0)
		{
		  fx = (kvec[0] * Box / 2) / Nmesh;
		  fx = sin(fx) / fx;
		}
	      if(kvec[1] != 0)
		{
		  fy = (kvec[1] * Box / 2) / Nmesh;
		  fy = sin(fy) / fy;
		}
	      if(kvec[2] != 0)
		{
		  fz = (kvec[2] * Box / 2) / Nmesh;
		  fz = sin(fz) / fz;
		}
	      ff = 1 / (fx * fy * fz);
	      smth = ff * ff;
	      /*  */
#endif

	      /* cdisp2 = source * k / (sqrt(-1) k^2) */
	      for(axes = 0; axes < 3; axes++)
		{
		  if(kmag2 > 0.0) 
		    {
		      cdisp2[axes][coord].re = cdigrad[3][coord].im * kvec[axes] / kmag2;
		      cdisp2[axes][coord].im = -cdigrad[3][coord].re * kvec[axes] / kmag2;
		    }
		  else cdisp2[axes][coord].re = cdisp2[axes][coord].im = 0.0;
#ifdef CORRECT_CIC
		  cdisp[axes][coord].re *= smth;   cdisp[axes][coord].im *= smth;
		  cdisp2[axes][coord].re *= smth;  cdisp2[axes][coord].im *= smth;
#endif
		}
	    }
      
      /* Free cdigrad[3] */
      free(cdigrad[3]);

      
      /* Now, both cdisp, and cdisp2 have the ZA and 2nd order displacements */

      for(axes = 0; axes < 3; axes++)
	{
          if(ThisTask == 0) fprintf(stderr,"\tFourier transforming displacements, axis %d.\n",axes);

	  rfftwnd_mpi(Inverse_plan, 1, disp[axes], Workspace, FFTW_NORMAL_ORDER);
	  rfftwnd_mpi(Inverse_plan, 1, disp2[axes], Workspace, FFTW_NORMAL_ORDER);

	  /* now get the plane on the right side from neighbour on the right, 
	     and send the left plane */
      
	  recvTask = ThisTask;
	  do
	    {
	      recvTask--;
	      if(recvTask < 0)
		recvTask = NTask - 1;
	    }
	  while(Local_nx_table[recvTask] == 0);
      
	  sendTask = ThisTask;
	  do
	    {
	      sendTask++;
	      if(sendTask >= NTask)
		sendTask = 0;
	    }
	  while(Local_nx_table[sendTask] == 0);
      
	  /* use non-blocking send */
      
	  if(Local_nx > 0)
	    {
	      /* send ZA disp */
	      MPI_Isend(&(disp[axes][0]),
			sizeof(fftw_real) * Nmesh * (2 * (Nmesh / 2 + 1)),
			MPI_BYTE, recvTask, 10, MPI_COMM_WORLD, &request);
	      
	      MPI_Recv(&(disp[axes][(Local_nx * Nmesh) * (2 * (Nmesh / 2 + 1))]),
		       sizeof(fftw_real) * Nmesh * (2 * (Nmesh / 2 + 1)),
		       MPI_BYTE, sendTask, 10, MPI_COMM_WORLD, &status);
	      
	      MPI_Wait(&request, &status);

	      
	      /* send 2nd order disp */
	      MPI_Isend(&(disp2[axes][0]),
			sizeof(fftw_real) * Nmesh * (2 * (Nmesh / 2 + 1)),
			MPI_BYTE, recvTask, 10, MPI_COMM_WORLD, &request);
	      
	      MPI_Recv(&(disp2[axes][(Local_nx * Nmesh) * (2 * (Nmesh / 2 + 1))]),
		       sizeof(fftw_real) * Nmesh * (2 * (Nmesh / 2 + 1)),
		       MPI_BYTE, sendTask, 10, MPI_COMM_WORLD, &status);
	      
	      MPI_Wait(&request, &status);
	    }
	}
      
      /* read-out displacements */

	  nmesh3 = ((unsigned int ) Nmesh ) * ((unsigned int) Nmesh) *  ((unsigned int) Nmesh);
      
      for(n = 0; n < NumPart; n++)
	{
	      /*
	      u = P[n].Pos[0] / Box * Nmesh;
	      v = P[n].Pos[1] / Box * Nmesh;
	      w = P[n].Pos[2] / Box * Nmesh;
	      */
	      u = partX[n] / Box * Nmesh;
	      v = partY[n] / Box * Nmesh;
	      w = partZ[n] / Box * Nmesh;
	      
	      i = (int) u;
	      j = (int) v;
	      k = (int) w;
	      
	      if(i == (Local_x_start + Local_nx))
		i = (Local_x_start + Local_nx) - 1;
	      if(i < Local_x_start)
		i = Local_x_start;
	      if(j == Nmesh)
		j = Nmesh - 1;
	      if(k == Nmesh)
		k = Nmesh - 1;
	      
	      u -= i;
	      v -= j;
	      w -= k;
	      
	      i -= Local_x_start;
	      ii = i + 1;
	      jj = j + 1;
	      kk = k + 1;
	      
	      if(jj >= Nmesh)
		jj -= Nmesh;
	      if(kk >= Nmesh)
		kk -= Nmesh;
	      
	      f1 = (1 - u) * (1 - v) * (1 - w);
	      f2 = (1 - u) * (1 - v) * (w);
	      f3 = (1 - u) * (v) * (1 - w);
	      f4 = (1 - u) * (v) * (w);
	      f5 = (u) * (1 - v) * (1 - w);
	      f6 = (u) * (1 - v) * (w); 
	      f7 = (u) * (v) * (1 - w);
	      f8 = (u) * (v) * (w);
	     
	      for(axes = 0; axes < 3; axes++)
		{
		  dis = disp[axes][(i * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + k] * f1 +
		    disp[axes][(i * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + kk] * f2 +
		    disp[axes][(i * Nmesh + jj) * (2 * (Nmesh / 2 + 1)) + k] * f3 +
		    disp[axes][(i * Nmesh + jj) * (2 * (Nmesh / 2 + 1)) + kk] * f4 +
		    disp[axes][(ii * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + k] * f5 +
		    disp[axes][(ii * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + kk] * f6 +
		    disp[axes][(ii * Nmesh + jj) * (2 * (Nmesh / 2 + 1)) + k] * f7 +
		    disp[axes][(ii * Nmesh + jj) * (2 * (Nmesh / 2 + 1)) + kk] * f8;

		  dis2 = disp2[axes][(i * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + k] * f1 +
		    disp2[axes][(i * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + kk] * f2 +
		    disp2[axes][(i * Nmesh + jj) * (2 * (Nmesh / 2 + 1)) + k] * f3 +
		    disp2[axes][(i * Nmesh + jj) * (2 * (Nmesh / 2 + 1)) + kk] * f4 +
		    disp2[axes][(ii * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + k] * f5 +
		    disp2[axes][(ii * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + kk] * f6 +
		    disp2[axes][(ii * Nmesh + jj) * (2 * (Nmesh / 2 + 1)) + k] * f7 +
		    disp2[axes][(ii * Nmesh + jj) * (2 * (Nmesh / 2 + 1)) + kk] * f8;
		  dis2 /= (float) nmesh3;
	      
/*		  
#ifdef ONLY_ZA 	  
		  P[n].Pos[axes] += dis;
		  P[n].Vel[axes] = dis * vel_prefac;
#else
		  P[n].Pos[axes] += dis - 3./7. * dis2;
		  P[n].Vel[axes] = dis * vel_prefac - 3./7. * dis2 * vel_prefac2;
#endif

		  P[n].Pos[axes] = periodic_wrap(P[n].Pos[axes]);
*/
		  
		  if (axes==0){
			#ifdef ONLY_ZA
		  	partX[n] += dis;
		  	partVX[n] = dis * vel_prefac;
			#else
		  	partX[n] += dis - 3./7. * dis2;
		  	partVX[n] = dis * vel_prefac - 3./7. * dis2 * vel_prefac2;
			#endif
		  	partX[n] = periodic_wrap(partX[n]);
		  }
		  else if (axes==1){
			#ifdef ONLY_ZA
		  	partY[n] += dis;
		  	partVY[n] = dis * vel_prefac;
			#else
		  	partY[n] += dis - 3./7. * dis2;
		  	partVY[n] = dis * vel_prefac - 3./7. * dis2 * vel_prefac2;
			#endif
		  	partY[n] = periodic_wrap(partY[n]);
		  }
		  if (axes==2){
			#ifdef ONLY_ZA
		  	partZ[n] += dis;
		  	partVZ[n] = dis * vel_prefac;
			#else
		  	partZ[n] += dis - 3./7. * dis2;
		  	partVZ[n] = dis * vel_prefac - 3./7. * dis2 * vel_prefac2;
			#endif
		  	partZ[n] = periodic_wrap(partZ[n]);
		  }

		  
		  
		  

		  if(dis - 3./7. * dis2 > maxdisp)
		    maxdisp = dis;
		}
	    
	}
    }
  

  for(axes = 0; axes < 3; axes++) free(cdisp[axes]);
  for(axes = 0; axes < 3; axes++) free(cdisp2[axes]);

  gsl_rng_free(random_generator);

  MPI_Reduce(&maxdisp, &max_disp_glob, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      fprintf(stderr,"\tMaximum displacement: %g kpc/h, in units of the part-spacing= %g\n",
	     max_disp_glob, max_disp_glob / (Box / Nmesh));
    }
}

double periodic_wrap(double x)
{
  while(x >= Box)
    x -= Box;

  while(x < 0)
    x += Box;

  return x;
}


void set_units(void)		/* ... set some units */
{
  UnitTime_in_s = UnitLength_in_cm / UnitVelocity_in_cm_per_s;

  G = GRAVITY / pow(UnitLength_in_cm, 3) * UnitMass_in_g * pow(UnitTime_in_s, 2);
  Hubble = HUBBLE * UnitTime_in_s;
}



void initialize_ffts(void)
{
  int total_size, i, additional;
  int local_ny_after_transpose, local_y_start_after_transpose;
  int *slab_to_task_local;
  size_t bytes;


  Inverse_plan = rfftw3d_mpi_create_plan(MPI_COMM_WORLD,
					 Nmesh, Nmesh, Nmesh, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE);

  Forward_plan = rfftw3d_mpi_create_plan(MPI_COMM_WORLD,
					 Nmesh, Nmesh, Nmesh, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE);

  rfftwnd_mpi_local_sizes(Forward_plan, &Local_nx, &Local_x_start,
			  &local_ny_after_transpose, &local_y_start_after_transpose, &total_size);

  Local_nx_table = malloc(sizeof(int) * NTask);
  MPI_Allgather(&Local_nx, 1, MPI_INT, Local_nx_table, 1, MPI_INT, MPI_COMM_WORLD);

#ifdef DEBUG
  if(ThisTask == 0)
    {
      for(i = 0; i < NTask; i++)
	fprintf(stderr,"Task=%d Local_nx=%d\n", i, Local_nx_table[i]);
      fflush(stdout);
    }
#endif

  Slab_to_task = malloc(sizeof(int) * Nmesh);
  slab_to_task_local = malloc(sizeof(int) * Nmesh);

  for(i = 0; i < Nmesh; i++)
    slab_to_task_local[i] = 0;

  for(i = 0; i < Local_nx; i++)
    slab_to_task_local[Local_x_start + i] = ThisTask;

  MPI_Allreduce(slab_to_task_local, Slab_to_task, Nmesh, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  free(slab_to_task_local);



  additional = (Nmesh) * (2 * (Nmesh / 2 + 1));	/* additional plane on the right side */

  TotalSizePlusAdditional = total_size + additional;

  //Disp = (fftw_real *) malloc(bytes = sizeof(fftw_real) * (total_size + additional));

  Workspace = (fftw_real *) malloc(bytes = sizeof(fftw_real) * total_size);

  #ifdef DEBUG
  ASSERT_ALLOC(Workspace)
  #endif

  //Cdata = (fftw_complex *) Disp;	/* transformed array */
}



void free_ffts(void)
{
  free(Workspace);
  //free(Disp);
  free(Slab_to_task);
  rfftwnd_mpi_destroy_plan(Inverse_plan);
  rfftwnd_mpi_destroy_plan(Forward_plan);
}





static double A, B, alpha, beta, V, gf;

double fnl(double x)		/* Peacock & Dodds formula */
{
  return x * pow((1 + B * beta * x + pow(A * x, alpha * beta)) /
		 (1 + pow(pow(A * x, alpha) * gf * gf * gf / (V * sqrt(x)), beta)), 1 / beta);
}

void print_spec(void)
{
  double k, knl, po, dl, dnl, neff, kf, kstart, kend, po2, po1, DDD;
  char buf[1000];
  FILE *fd;

  if(ThisTask == 0)
    {
      sprintf(buf, "%s/inputspec_%s.txt", OutputDir, FileBase);

      fd = fopen(buf, "w");

      gf = GrowthFactor(0.001, 1.0) / (1.0 / 0.001);

      DDD = GrowthFactor(1.0 / (Redshift + 1), 1.0);

      fprintf(fd, "%12g %12g\n", Redshift, DDD);	/* print actual starting redshift and 
							   linear growth factor for this cosmology */

      kstart = 2 * PI / (1000.0 * (3.085678e24 / UnitLength_in_cm));	/* 1000 Mpc/h */
      kend = 2 * PI / (0.001 * (3.085678e24 / UnitLength_in_cm));	/* 0.001 Mpc/h */

      for(k = kstart; k < kend; k *= 1.025)
	{
	  po = PowerSpec(k);
	  dl = 4.0 * PI * k * k * k * po;

	  kf = 0.5;

	  po2 = PowerSpec(1.001 * k * kf);
	  po1 = PowerSpec(k * kf);

	  if(po != 0 && po1 != 0 && po2 != 0)
	    {
	      neff = (log(po2) - log(po1)) / (log(1.001 * k * kf) - log(k * kf));

	      if(1 + neff / 3 > 0)
		{
		  A = 0.482 * pow(1 + neff / 3, -0.947);
		  B = 0.226 * pow(1 + neff / 3, -1.778);
		  alpha = 3.310 * pow(1 + neff / 3, -0.244);
		  beta = 0.862 * pow(1 + neff / 3, -0.287);
		  V = 11.55 * pow(1 + neff / 3, -0.423) * 1.2;

		  dnl = fnl(dl);
		  knl = k * pow(1 + dnl, 1.0 / 3);
		}
	      else
		{
		  dnl = 0;
		  knl = 0;
		}
	    }
	  else
	    {
	      dnl = 0;
	      knl = 0;
	    }

	  fprintf(fd, "%12g %12g    %12g %12g\n", k, dl, knl, dnl);
	}
      fclose(fd);
    }
}



int write_halogen_cat(char *filename, float *x, float *y, float *z, float *vx, float *vy, float *vz, float *M, float *R,long N){
        FILE *f;
        long i;

        if ((f=fopen(filename,"w") )== NULL){
                fprintf(stderr,"Couldnt open output file %s\n",filename);
                return -1;
        }
        for(i=0;i<N;i++){
                fprintf(f,"%f %f %f %f %f %f %e %f\n",x[i],y[i],z[i],vx[i],vy[i],vz[i],M[i],R[i]);
        }
        fclose(f);
        return 0;
}

int read_input_file(char *name){
        FILE *f;
        char line[LINELENGTH],word[LINELENGTH];
        int i;
        float x,y,z;
	NParametersSet = 0;

        for (i=0;i<NParam;i++)
                ParameterSet[i]=0;

        if ((f = fopen(name,"r"))==NULL){
                fprintf(stderr,"Could not open input file %s\n",name);
                return -1;
        }
        while(fgets(line,LINELENGTH,f)!=NULL) {
                if (line[0] != '#') {                           //Omit comment lines
                        if (sscanf(line,"%s",word)!=EOF){       //Check it is not an empty line
                          for (i=0;i<NParam;i++)                //Find the parameter specified
                                if (strcmp(ParameterList[i],word)==0)
                                        break;
                          switch (i){

                                case 0:
                                        sscanf(line,"%s %s",word,MassFunctionFile);
                                        if (ParameterSet[i]>0)
                                                fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
                                        else
                                                NParametersSet++;
                                        #ifdef VERB 
                                        if(ThisTask == 0) fprintf(stderr,"\t%s: %s\n",word,MassFunctionFile);
                                        #endif
                                        ParameterSet[i]++;
                                        break;

                                case 1:
                                        sscanf(line,"%s %s",word,OutputFile);
                                        if (ParameterSet[i]>0)
                                                fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
                                        else
                                                NParametersSet++;
                                        #ifdef VERB 
                                        if(ThisTask == 0) fprintf(stderr,"\t%s: %s\n",word,OutputFile);
                                        #endif
                                        ParameterSet[i]++;
                                        break;

                                case 2:
                                        sscanf(line,"%s %d",word,&Nlin);
                                        if (ParameterSet[i]>0)
                                                fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
                                        else
                                                NParametersSet++;
                                        #ifdef VERB 
                                        if(ThisTask == 0) fprintf(stderr,"\t%s: %d\n",word,Nlin);
                                        #endif
                                        ParameterSet[i]++;
                                        break;
                                case 3:
                                        sscanf(line,"%s %s",word,alphaFile);
                                        if (ParameterSet[i]>0)
                                                fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
                                        else
                                                NParametersSet++;
                                        #ifdef VERB 
                                        if(ThisTask == 0) fprintf(stderr,"\t%s: %s\n",word,alphaFile);
                                        #endif
                                        ParameterSet[i]++;
                                        break;

                                case 4:
                                        sscanf(line,"%s %s",word,rho_ref);
                                        if (ParameterSet[i]>0)
                                                fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
                                        else
                                                NParametersSet++;
                                        #ifdef VERB 
                                        if(ThisTask == 0) fprintf(stderr,"\t%s: %s\n",word,rho_ref);
                                        #endif
                                        if ((strcmp(rho_ref,"crit")!=0) && (strcmp(rho_ref,"matter")!=0)){
                                                        fprintf(stderr,"ERROR: Not valid option for %s: %s.\nPlease select \"crit\" or \"matter\". \n",word,rho_ref);
                                                        return -1;
                                        }
                                        ParameterSet[i]++;
                                        break;

                                case 5:
                                        sscanf(line,"%s %f",word,&OVD);
                                        if (ParameterSet[i]>0)
                                                fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
                                        else
                                                NParametersSet++;
                                        #ifdef VERB 
                                        if(ThisTask == 0) fprintf(stderr,"\t%s: %f\n",word,OVD);
                                        #endif
                                        ParameterSet[i]++;
                                        break;

                                case 6:
                                        sscanf(line,"%s %f",word,&Mmin);
                                        if (ParameterSet[i]>0)
                                                fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
                                        else
                                                NParametersSet++;
                                        #ifdef VERB
                                                #ifdef NDENS
                                                if (Mmin>1)
                                                        fprintf(stderr,"WARNING: This does not seem a number density in [Mpc/h]^{-3}\n");
                                                if(ThisTask == 0) fprintf(stderr,"\tNdens: %e.\n",Mmin);
                                                #else
                                                if (Mmin<1e5)
                                                        fprintf(stderr,"WARNING: This does not seem a Mass in [Msun/h]\n");
                                                if(ThisTask == 0) fprintf(stderr,"\t%s: %e.\n",word,Mmin);
                                                #endif
                                        #endif
                                        ParameterSet[i]++;
                                        break;

                                case 7:
                                        sscanf(line,"%s %ld",word,&seed);
                                        if (ParameterSet[i]>0)
                                                fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
                                        else
                                                NParametersSet++;
                                        #ifdef VERB
                                        if(ThisTask == 0) fprintf(stderr,"\t%s: %ld\n",word,seed);
                                        #endif
                                        ParameterSet[i]++;
					break;

                                case 8:
                                        sscanf(line,"%s %f",word,&recalc_frac);
                                        if (ParameterSet[i]>0)
                                                fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
                                        else
                                                NParametersSet++;
                                        #ifdef VERB
                                        if(ThisTask == 0) fprintf(stderr,"\t%s: %f\n",word,recalc_frac);
                                        #endif
                                        ParameterSet[i]++;
                                        break;
                                case 9:
                                        sscanf(line,"%s %d",word,&nthreads);
                                        if (ParameterSet[i]>0)
                                                fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
                                        else
                                                NParametersSet++;
                                        #ifdef VERB
                                        if(ThisTask == 0) fprintf(stderr,"\t%s: %d\n",word,nthreads);
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

	#ifdef DEBUG
	fprintf(stderr,"\tNumber of parameters set: %d\n",NParametersSet);
	#endif

        if ((f = fopen(alphaFile,"r"))==NULL){
                fprintf(stderr,"Could not open input file %s\n",alphaFile);
                return -1;
        }
        while(fgets(line,LINELENGTH,f)!=NULL)
                if (line[0] != '#')
                        Nalpha++;
        fclose(f);

        alpha_vec = (double *) calloc(Nalpha,sizeof(double));
        fvel = (double *) calloc(Nalpha,sizeof(double));
        Malpha = (double *) calloc(Nalpha,sizeof(double));

        if ((f = fopen(alphaFile,"r"))==NULL){
                fprintf(stderr,"Could not open alpha file %s\n",alphaFile);
                return -1;
        }

        i=0;
        while(fgets(line,LINELENGTH,f)!=NULL){
                if (line[0] != '#'){
                        //fgets(line,LINELENGTH,f);
                        sscanf(line,"%f %f %f",&x,&y,&z);
                        Malpha[i]= (double) x;
                        alpha_vec[i]=(double) y;
                        fvel[i]= (double) z;
                        i++;
                }
        }
        fclose(f);
        return 0;
}

