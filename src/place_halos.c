
/*********************************************************************************** 
 PLACE_HALOS
************************************************************************************ 
This file is part of HALOGEN (http://arxiv.org/abs/1412.5228).
HALOGEN has been developped by Santiago Avila and Steven Murray 
(https://github.com/savila/HALOGEN).
************************************************************************************ */



#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h> 
#include <omp.h> 
#include <string.h>

#include "place_halos.h"
 

#define MAXTRIALS (50)


/*=============================================================================
 *                             PROTOTYPES
 *=============================================================================*/


//long select_cell_rnd(long *, long *, long *); 
long select_cell(double , double *); 
long select_part_beta_0(long, long **,long *);
int check_HaloR_in_mesh(long,float *, float *, float * , float *,long,long,long,long*,long*,long*,int);
int check_HaloR_in_cell(long ,float *, float *, float * , float *,long ,long,long,long *,long *,long*);
double ComputeCumulative(double,double, double *, double *);

//Global Variables
long NCells,NTotCells;
float Lbox,lcell;
int NTHREADS;


/*=============================================================================
 *                             simple functions
 *=============================================================================*/

//square()
float square(float a){
	return a*a;
}

// R_from_mass():
// takes a halo mass and transforms it to a halo radius, for a given definition 
// in terms of density at the edge (rho, in physical units).
float R_from_mass(float Mass,float rho) {
	return  (float) pow((3./(4.*rho*M_PI)*Mass),(1./3.));
}

//chec_limit(): check that i is between 0 and N-1
long check_limit(long i, long N){
	if (i==N)
		return 0; //Apply boundary conditions
	if (i<0 || i>N){
		fprintf(stderr,"particle assigned to unexisting cell %ld\nExiting...",i);
		exit(0);
	}
	return i;
}


/*=============================================================================
 *                             place_halos()
 *=============================================================================*/

//place_halos():
//
//Takes a list of halo masses (Nhalos, HaloMass), a list of particles 
// (NTotPart,PartX,PartY,PartZ), some simulation parameters (L, mp), and 
// user-defined parameters (Nlin,rho_ref,alpha,Malpha,Nalpha,seed)
//and returns a list of halo positions and radii (HaloX,HaloY,HaloZ,HaloR)
int place_halos(long Nend, float *HaloMass, long Nlin, long NTotPart, float *PartX, 
		float *PartY, float *PartZ, float *PartVX, float *PartVY, float *PartVZ,
		float L, float rho_ref, long seed, float mp, int nthreads, double *alpha, double *fvel, double *Malpha,
		long Nalpha,float recalc_frac, float *HaloX, float *HaloY, float *HaloZ, float *HaloVX,
		float *HaloVY, float *HaloVZ,float *HaloR,long **ListOfPart, 
		long *NPartPerCell){


fprintf(stderr,"\tThis is place_halos.c\n");

//Initiallising -------------------------------------------------
	long i,j,k,lin_ijk, Nmin;
	long *count,trials;
	long ihalo, ipart,i_alpha;
	double invL = 1./L;
	float Mcell,Mhalo,Mchange; 
	float R;
	time_t t0,tI,tII;
	int check;

	double mpart,fvel_i;
	double exponent;
	double TotProb;
	double prob_repicked = 0.0;
	double *MassLeft;
	double *CumulativeProb; 
	long *ListOfHalos,  *NHalosPerCellStart, *NHalosPerCellEnd;
	long Nhalos;
	int recalc;



	float diff;
	time_t t5;
	#ifdef VERB
	time_t t1,t3,t4,t4_5;
	#endif

	long n_recalc =0;
	int use_vel=1;
	
	if (HaloVX==NULL)
		use_vel=0;


	NCells = Nlin;
	Lbox = L;
	
	t0=time(NULL);
	NTotCells = NCells*NCells*NCells;

	
	//Allocate memory for the arrays 
	MassLeft = (double *) calloc(NTotCells,sizeof(double));
  	if(MassLeft == NULL) {
    		fprintf(stderr,"\tplace_halos(): could not allocate %ld array for MassLeft[]\nABORTING",NTotCells);
    		exit(-1);
	}	
	NHalosPerCellStart = (long *) calloc(NTotCells,sizeof(long));
  	if(NHalosPerCellStart == NULL) {
    		fprintf(stderr,"\tplace_halos(): could not allocate %ld array for NHalosPerCell[]\nABORTING",NTotCells);
    		exit(-1);
	}

  	NHalosPerCellEnd = (long *) calloc(NTotCells,sizeof(long));
  	if(NHalosPerCellEnd == NULL) {
    		fprintf(stderr,"\tplace_halos(): could not allocate %ld array for NHalosPerCell[]\nABORTING",NTotCells);
    		exit(-1);
	}

  	count = (long *) calloc(NTotCells,sizeof(long));
  	if(count == NULL) {
    		fprintf(stderr,"\tplace_halos(): could not allocate %ld array for NTotCells[]\nABORTING",NTotCells);
    		exit(-1);
	}

	CumulativeProb = (double *) calloc(NTotCells, sizeof(double));
  	if(CumulativeProb == NULL) {
    		fprintf(stderr,"\tplace_halos(): could not allocate %ld array for CumulativeProb[]\nABORTING",NTotCells);
    		exit(-1);
	}

  	if (nthreads<1){
  		NTHREADS = omp_get_max_threads();
  	}else{
  		NTHREADS = nthreads;
  	}
	#ifdef NO_EXCLUSION
	int *already_chosen;
	already_chosen = (int*) calloc(NTotPart,sizeof(int));
 	if(already_chosen == NULL) {
                fprintf(stderr,"\tplace_halos(): could not allocate %ld array for already_chosen[]\nABORTING",NTotPart);
                exit(-1);
        }
	#endif

        //Initiallise random numbers
	#ifdef VERB
        fprintf(stderr,"\tinput seed: %ld.    time0: %f.",seed, (float) t0);
	#endif
	if (seed>=0){
		srand(seed);
		#ifdef VERB
			fprintf(stderr,"\tUsed: %ld \n",seed);
		#endif
	}
	else {
		srand(t0);
#ifdef VERB
		fprintf(stderr,"\tSeed Used: %ld \n",t0);
#endif
	}


	mpart = (double) mp;
	Nmin = (long)ceil(HaloMass[Nend-1]*0.9/mpart);
	lcell = (float) L/NCells;

	#ifdef VERB
	fprintf(stderr,"\n\tParticles and Halos placed in %ld^3 cells\n",NCells);
	fprintf(stderr,"\tBOX = %f  lcell =%f   rho_ref = %e  invL %f\n",L,L/NCells,rho_ref,invL);
	fprintf(stderr,"\tNhalostart = %d,Nhalosend = %ld,  NPart = %ld\n",0, Nend, NTotPart);
	fprintf(stderr,"\n\tMinimmum mass= %e. Minimum part per halo = %ld. mpart %e\n",HaloMass[Nend-1],Nmin,mpart);
	#endif
	

	#ifdef DEBUG
	fprintf(stderr,"\n\tRAND_MAX=%d\n",RAND_MAX);
	fprintf(stderr,"\tX[0] = %f Y[0] = %f Z[0] = %f\n",PartX[0],PartY[0],PartZ[0]);
	fprintf(stderr,"\tX[1] = %f Y[1] = %f Z[1] = %f\n",PartX[1],PartY[1],PartZ[1]);
	fprintf(stderr,"\tM[0] = %e \n",HaloMass[0]);
	fprintf(stderr,"\tM[1] = %e \n",HaloMass[1]);
	fprintf(stderr,"\t    ... \n");
	fprintf(stderr,"\tM[%ld] = %e \n",Nend-1,HaloMass[Nend-1]);
	fprintf(stderr,"\tX[%ld] = %f Y[%ld] = %f Z[%ld] = %f\n",Nend-1,PartX[Nend-1],Nend-1,PartY[Nend-1],Nend-1,PartZ[Nend-1]);
	#endif	
	
	int r = (int) (R_from_mass(HaloMass[0],rho_ref)/(L/NCells));
	if (L/NCells<R_from_mass(HaloMass[0],rho_ref)){
		fprintf(stderr,"WARNING: cell size is smaller than the radius of the biggest halo. Using r=%i. This may be problematic\n",r);
	}

#ifdef VERB
	fprintf(stderr,"\tR_max=%f, lcell=%f, r=%d\n",R_from_mass(HaloMass[0],rho_ref),(L/NCells),r);
	t1=time(NULL);
 	diff = difftime(t1,t0);
	fprintf(stderr,"\ttime of initialisation %f\n",diff);
#endif
// ------------------------------------------------- Initiallised

	//Alloc Enough Memory
	Nhalos=0;
	for (i=0;i<NCells;i++){
	for (j=0;j<NCells;j++){
	for (k=0;k<NCells;k++){
		lin_ijk = k+j*NCells+i*NCells*NCells;
		NHalosPerCellStart[lin_ijk] = Nhalos;
		NHalosPerCellEnd[lin_ijk] = Nhalos;
		Nhalos += (long) floor(NPartPerCell[lin_ijk]/Nmin+1);
		MassLeft[lin_ijk] = (double) NPartPerCell[lin_ijk]*mpart;
#ifdef ULTRADEBUG
		if (lin_ijk<10 || lin_ijk > (NCells*NCells*NCells) - 10){
			fprintf(stderr,"\tAllocated %ld (longs) in ListOfPart(%ld=[%ld,%ld,%ld])\n",NPartPerCell[lin_ijk],lin_ijk,i,j,k);
		}
#endif		
	}	
	}
	}


	ListOfHalos = (long *) calloc(Nhalos,sizeof(long ));
	if(ListOfHalos == NULL) {
		fprintf(stderr,"\tplace_halos(): could not allocate %ld array for ListOfHalos[]\nABORTING",Nhalos);
		exit(-1);
	}

#ifdef VERB
//	fprintf(stderr,"\tAllocated %ld (longs) in ListOfHalos\n",Nhalos);
	t3=time(NULL);
 	diff = difftime(t3,t1);
	fprintf(stderr,"\t... memory allocated in %f\n",diff);
	fprintf(stderr,"\tComputing probabilities...\n");
#endif 

#ifdef DEBUG
        fprintf(stderr,"\tMass_cell[0]=%e",MassLeft[0]);
	fprintf(stderr,"\t Mass Function\n");
	for (ihalo=0;ihalo<15;ihalo++){
		fprintf(stderr,"\thalo %ld: ",ihalo);
		fprintf(stderr,"M=%e\n",HaloMass[ihalo]);
	}
#endif

//----------------------------------- Particles and haloes assigned to grid



//Computing Cumulative Probability -----------------------------
	
	//find the right alpha
	Mhalo = HaloMass[0];
	i_alpha = 0;
	while(Mhalo<Malpha[i_alpha]) {
		i_alpha++;
		if (i_alpha==Nalpha){
			fprintf(stderr,"\tERROR: No M_alpha low enough found\n");
			fprintf(stderr,"\tERROR: N_alpha = %ld, Mh=%e, Ma= %e\n",Nalpha,Mhalo,Malpha[i_alpha-1]);
			exit(0);
		}
	}	
	Mchange = Malpha[i_alpha];
	exponent = alpha[i_alpha];
	fvel_i = fvel[i_alpha];

	//compute the probability

#ifdef VERB
	fprintf(stderr,"\tUsing OMP with %d threads\n",NTHREADS);
	t4=time(NULL);
#endif
	TotProb = ComputeCumulative(exponent, mpart, MassLeft, CumulativeProb);	
#ifdef VERB
	fprintf(stderr,"\n\tcase 0, TotProb=%e\n",TotProb);
#endif


#ifdef VERB
        fprintf(stderr,"\tNumber of alphas: %ld\n",Nalpha);
        fprintf(stderr,"\tUsing alpha_%ld=%f for M>%e\n",i_alpha,exponent,Mchange);
	t4_5=time(NULL);
 	diff = difftime(t4_5,t4);
	fprintf(stderr,"\tprobabilty computed in %f secods\n",diff);
#endif
// ----------------------------------------- Computed Probability



//Actually placing the haloes----------------------------------- 
#ifdef VERB
	fprintf(stderr,"\n\tPlacing Halos...\n\n");
#endif

	//Place one by one all the haloes (assumed to be ordered from the most massive to the least massive)
	for (ihalo=0;ihalo<Nend;ihalo++){

		#ifdef DEBUG
		fprintf(stderr,"\n\t- Halo %ld ",ihalo);
		#endif
		#ifdef VERB
		if (ihalo%(Nend/10)==0 && ihalo>0){
			//TEMPORARY
			fprintf(stderr,"\t\tFRAC, TOTPROB: %e, %e",(pow(Mcell/mpart,exponent)/TotProb),TotProb);
			fprintf(stderr,"\t%ld%% done\n",(ihalo/(Nend/100)));
		}
		#endif
		//Check whether or not, a change of alpha is needed for this halo mass 		
		Mhalo= HaloMass[ihalo];
		recalc = 0;
		while (Mhalo < Mchange){//if so search the right alpha, and recompute probabilities
			i_alpha++;		
			if (i_alpha==Nalpha){
				fprintf(stderr,"\tERROR: No M_alpha low enough found: %e <%e\n",Mhalo,Malpha[Nalpha-1]);
				exit(0);
			}
			Mchange = Malpha[i_alpha];
			exponent = alpha[i_alpha];
			fvel_i = fvel[i_alpha];


		#ifdef VERB
        		fprintf(stderr,"\n\tUsing alpha_%ld=%f for M>%e\n",i_alpha,exponent,Mchange);
		#endif
        	recalc = 1;
		}
		
		// recalc if different alpha, OR there's a significant chance of choosing the same cell again.
		if(ihalo>0){
		  if(prob_repicked>=recalc_frac){
			recalc = 1;
			n_recalc += 1;
			fprintf(stderr,"RECALCULATING: %ld, %e,    ihalo=%ld\n",n_recalc,prob_repicked,ihalo);

		  }
		}

		if (recalc==1){
			tI=time(NULL);
			fprintf(stderr,"\tcase 1, TotProb_bef=%e",TotProb);
			TotProb=ComputeCumulative(exponent, mpart, MassLeft, CumulativeProb);
			fprintf(stderr,"    TotProb_aft=%e       ihalo=%ld\n\n",TotProb,ihalo);

			prob_repicked=0.0;
#ifdef VERB
			tII=time(NULL);
			diff = difftime(tII,tI);
			fprintf(stderr,"\tProbabilty recomputed in %f secods\n",diff);
#endif

			recalc = 0;
		}


		do {	
		  //First, choose a cell	
		  #ifndef RANKED	
		  trials=0;
		  do{			
			if (trials==MAXTRIALS){
				fprintf(stderr,"MAXTRIALS=%d times picked an empty cell, recomputing Probs...\n",MAXTRIALS);
				fprintf(stderr,"\n\tcase 2, TotProb_bef=%e",TotProb);
				TotProb=ComputeCumulative(exponent, mpart, MassLeft, CumulativeProb);
				fprintf(stderr,"    TotProb_aft=%e       ihalo=%ld\n",TotProb,ihalo);
				prob_repicked = 0.0;
				trials=0;
				
			}
		  	lin_ijk = select_cell(TotProb, CumulativeProb);
			trials++;
		
		  }while (MassLeft[lin_ijk]==0.);



		  k=lin_ijk%(NCells);
		  j=((lin_ijk-k)/NCells)%NCells;
	  	  i=(lin_ijk-k-j*NCells)/(NCells*NCells);

		  #else //RANKED option: deprecated and not optimised
		  lin_ijk=select_heaviest_cell(&i,&j,&k,MassLeft);		  
		  #endif

		  trials=0;


		  //Second, choose a particle in that cell
		  do {

			ipart = select_part_beta_0(lin_ijk,ListOfPart, NPartPerCell);		
			if (ipart<0){
				fprintf(stderr,"WARNING: Picked up an completely empty cell (ihalo %ld) lin_ijk=%ld \n",ihalo,lin_ijk);
				MassLeft[lin_ijk]=0.;
				check=1;  //Choose another cell
				break;
			}

               		HaloX[ihalo] = PartX[ipart];
               		HaloY[ihalo] = PartY[ipart];
               		HaloZ[ihalo] = PartZ[ipart];
			#ifdef DEBUG
			fprintf(stderr,"HaloX=%f PartX=%f\n",HaloX[ihalo],PartX[ipart]);
			#endif

			if (use_vel==1){
               			HaloVX[ihalo] = fvel_i * PartVX[ipart];
               			HaloVY[ihalo] = fvel_i * PartVY[ipart];
               			HaloVZ[ihalo] = fvel_i * PartVZ[ipart];
			}
			R=R_from_mass(HaloMass[ihalo],rho_ref);
			HaloR[ihalo]= R;

			#ifdef NO_EXCLUSION
			  	check = already_chosen[part];
				already_chosen[ipart]=1;
			#else
			//Third, check that is not overlapping a previous halo
			check = check_HaloR_in_mesh(ihalo,HaloX,HaloY,HaloZ,HaloR,i,j,k,ListOfHalos,NHalosPerCellStart,NHalosPerCellEnd,r);
			#endif
			

			if (check==1){
				#ifdef DEBUG
				fprintf(stderr,"Refused part : %ld\n",ipart);
				#endif
				trials++;
			}
			if (trials == MAXTRIALS){
				//in order to avoid infinite loop, we will exit this loop, after MAXTRIALS trials
				#ifdef VERB
				fprintf(stderr,"MAXTRIALS=%d reached, removing cell [%ld,%ld,%ld]\n",MAXTRIALS,i,j,k);
				#endif
				MassLeft[lin_ijk]=0.;
				fprintf(stderr,"\n\tcase 3, TotProb_bef=%e",TotProb);
				TotProb=ComputeCumulative(exponent, mpart, MassLeft, CumulativeProb);
				fprintf(stderr,"    TotProb_aft=%e       ihalo=%ld, R=%f\n",TotProb,ihalo,R);
				prob_repicked=0.0;
				trials=0;
				break;
			}
		  } while (check==1);//If the particle was excluded, try another one in the same cell

	        } while(check==1); //if reached MAXTRIALS, select another cell
		//Particle chosen!
		
		//mass in cell before assignment
                Mcell=MassLeft[lin_ijk];

		
		  #ifndef MASS_OF_PARTS 
                  if (Mcell>HaloMass[ihalo])
			MassLeft[lin_ijk] -= Mhalo; 
                  else
			MassLeft[lin_ijk] = 0.;
		  #else
			exclude(ipart,R,PartX,PartY,PartZ,i,j,k);
		  #endif

		prob_repicked += pow(Mcell/mpart,exponent)/TotProb;


		#ifdef DEBUG
		fprintf(stderr,"\tAfter: Mcell=%e, CProbCell=%e, TotProb=%e.   , Mhalo=%e. CProb[last]=%e\n",MassLeft[lin_ijk],CumulativeProb[lin_ijk],TotProb,Mhalo,CumulativeProb[NTotCells-1]);
		#endif
		#ifdef DEBUG
		fprintf(stderr,"\thalo %ld assigned to particle %ld at [%f,%f,%f]. R= %f, M= %e\n",ihalo,ipart,HaloX[ihalo],HaloY[ihalo],HaloZ[ihalo],R,Mhalo);
		#endif
		#ifdef DEBUG
		fprintf(stderr,"HaloX=%f PartX=%f\n",HaloX[ihalo],PartX[ipart]);
		#endif

		ListOfHalos[NHalosPerCellEnd[lin_ijk]]=ihalo;
		NHalosPerCellEnd[lin_ijk]++;

	}//for(ihalo=Nstart:Nend)
//----------------------------------- Haloes Placed

	fprintf(stderr,"\t... placement Done!\n");
	fprintf(stderr,"\t\tTOTAL NUMBER OF RE-CALCULATIONS: %ld\n",n_recalc);

#ifdef VERB
	t5=time(NULL);
 	diff = difftime(t5,t4_5);
	fprintf(stderr,"\ttime placing %f\n",diff);
	fprintf(stderr,"\tfreeing...\n");
#endif

	free(NHalosPerCellStart);
	free(NHalosPerCellEnd);
        free(count); 
        free(CumulativeProb);
	free(MassLeft);
        free(ListOfHalos);
#ifdef VERB
 	diff = difftime(t5,t0);
	fprintf(stderr,"\ttotal time in place_halos.c %f\n",diff);
	fprintf(stderr,"\tPlacement done!!!\n");
#endif

#ifdef MASS_OF_PARTS
//	free(excluded); free(Nexcluded);
#endif
	return 0;
}
//end of place_halos()



//ComputeCumulative():
//it takes the user-defined alpha (exponent) and (re-)computes the cumulative probability for the cells in the grid
//from the masses in the cells (MassLeft
#ifndef NO_PROB_PARALLEL
double ComputeCumulative(double alpha, double mpart, double *MassLeft, double *CumulativeProb){
        long i,j,k,lin_ijk;
        double *PartProb;
        PartProb = (double *) calloc(NCells,sizeof(double));
        #pragma omp parallel for num_threads(NTHREADS) private(i,j,k,lin_ijk) shared(CumulativeProb,PartProb,NCells,alpha,mpart,MassLeft) default(none)
        for(i=0;i<NCells;i++){
                for(j=0;j<NCells;j++){
                for(k=0;k<NCells;k++){
                        lin_ijk = k+j*NCells+i*NCells*NCells;
                        PartProb[i] += pow(MassLeft[lin_ijk]/mpart,alpha);
                        CumulativeProb[lin_ijk] = PartProb[i];
                }
                }
        }
        for (i=1;i<NCells;i++){
                PartProb[i]+=PartProb[i-1];
        }

        #pragma omp parallel for num_threads(NTHREADS) private(i,j,k,lin_ijk) shared(CumulativeProb,PartProb,NCells) default(none)
        for(i=1;i<NCells;i++){
                for(j=0;j<NCells;j++){
                for(k=0;k<NCells;k++){
                        lin_ijk = k+j*NCells+i*NCells*NCells;
                        CumulativeProb[lin_ijk] += PartProb[i-1];
                }
                }
        }

        return PartProb[NCells-1];
}
#else
double ComputeCumulative(double alpha, double mpart, double *MassLeft, double *CumulativeProb){
        long i;
        double TotProb = 0.;
        for(i=0;i<NTotCells;i++){
                CumulativeProb[i] = TotProb + pow(MassLeft[i]/mpart,alpha);
                TotProb = CumulativeProb[i];
        }
        return TotProb;
}
#endif



//select_part():


//randomly selects a particle from a cell
long select_part_beta_0(long ijk,long **ListOfPart,long *NPartPerCell){
	long i_rnd,ipart;
	if(NPartPerCell[ijk]==0)
		return -1;
	i_rnd = (long) (NPartPerCell[ijk] * ((double)rand()/(RAND_MAX+1.0)));
	#ifdef DEBUG
	fprintf(stderr,"irnd: %ld/%ld ",i_rnd,NPartPerCell[ijk]);
	#endif
        ipart = ListOfPart[ijk][i_rnd];
	return ipart;
}




//select_heaviestcell()
long select_heaviest_cell(long *x, long *y, long *z, double* MassLeft) {
	long i,j,k,lin_ijk, out_ijk=-1;
	float max=0.0;	
	for (i=0;i<NCells;i++){
        for (j=0;j<NCells;j++){
        for (k=0;k<NCells;k++){
		lin_ijk = k+j*NCells+i*NCells*NCells;
		if (max<MassLeft[lin_ijk]){
			max=MassLeft[lin_ijk];
			*x=i;
			*y=j;
			*z=k;
			out_ijk = lin_ijk;	
		}
	}
	}
	}
	return out_ijk;
}



//select_cell():
//select one cell following the weight of the probabilities
long select_cell(double TotProb, double *CumulativeProb) {
        long   i_low, i_up, i_mid;
        double d_rand;

        d_rand = TotProb * ((double)rand()/(RAND_MAX+1.0));

        i_low = 0;
        i_up  = NTotCells-1;
        while(i_low < i_up) {
                i_mid = (i_low+i_up)/2;
                if(d_rand < CumulativeProb[i_mid])
                        i_up  = i_mid;
                else
                        i_low = i_mid + 1;
         }
        return(i_low);
}



//check_HaloR_in_cell():
//checks if there is any collision between the halo just been placed and any previous one, in the cell specified (i,j,k)
//returns 1 if there is any collision, returns 0 for no collision
int check_HaloR_in_cell(long ipart,float *PartX, float *PartY, float *PartZ, float *PartR, long i,long j, long k,long *ListOfHalos,long *NHalosPerCellStart,long *NHalosPerCellEnd){
        long jpart,jj;
        double X=PartX[ipart],Y=PartY[ipart],Z=PartZ[ipart];
	#ifndef HALF_EXCLUSION
	double R=PartR[ipart];
	#endif 

//The cell passed, might not exists, but be a virtual one (for periodic conditions): check, and correct:
//#ifdef PERIODIC

                if (i<0){
                        i += NCells;
                        X += Lbox;
                }
                else if (i>=NCells){
                        i -= NCells;
                        X -= Lbox;
                }
                if (j<0){
                        j += NCells;
                        Y += Lbox;
                }
                else if (j>=NCells){
                        j -= NCells;
                        Y -= Lbox;
                }
                if (k<0){
                        k += NCells;
                        Z += Lbox;
                }
                else if (k>=NCells){
                        k -= NCells;
                        Z -= Lbox;
		}
//#endif

  if (i>=0 && i<NCells && j>=0 && j<NCells && k>=0 && k<NCells){
#ifdef DEBUG
       fprintf(stderr,"Checking cell [%ld,%ld,%ld] for halo %ld",i,j,k,ipart);
#endif
        long lin_ijk = k+j*NCells+i*NCells*NCells;
#ifdef ULTRADEBUG
	fprintf(stderr,"= %ld.",lin_ijk);
#endif
	//loop over all the halos in that cell
	for (jj=NHalosPerCellStart[lin_ijk]; jj<NHalosPerCellEnd[lin_ijk]; jj++){
#ifdef ULTRADEBUG
		fprintf(stderr,"jj=%ld/%ld ",jj-jj=NHalosPerCellStart[lin_ijk],NHalosPerCell[lin_ijk]-jj=NHalosPerCellStart[lin_ijk]);
#endif
		jpart=ListOfHalos[jj];
#ifdef ULTRADEBUG
		fprintf(stderr,"jpart=%ld ",jpart);
#endif	
		//Check for overlapping
		#ifdef HALF_EXCLUSION
		if ((square(X-PartX[jpart])+square(Y-PartY[jpart])+square(Z-PartZ[jpart]))<square(PartR[jpart])) 
		#else
		if ((square(X-PartX[jpart])+square(Y-PartY[jpart])+square(Z-PartZ[jpart]))<square(R+PartR[jpart])) 
		#endif
		{
#ifdef DEBUG
                        fprintf(stderr,"\nChecking cell [%ld,%ld,%ld] for halo %ld",i,j,k,ipart);
                        fprintf(stderr," lin_ijk= %ld.",lin_ijk);
                        fprintf(stderr,"jj=%ld/%ld ",jj,NHalosPerCellEnd[lin_ijk]-NHalosPerCellStart[lin_ijk]);
                        fprintf(stderr,"jhalo=%ld ",jpart);
                        fprintf(stderr,"refused!\n");
#endif
			return 1;
		}
	}
	return 0;
  }
  else {
	fprintf(stderr,"WARNING: Computing distances outside the box. Cell: [%ld,%ld,%ld].\n",i,j,k);
	return 0;
  }
}



//cheack_haloR_in_mesh():
//checks if there is any collision between the halo just been placed and any previous one (going through all the neighbour cells)
//returns 1 if there is any collision, returns 0 for no collision
int check_HaloR_in_mesh(long ihalo,float *X, float *Y, float *Z , float *R,long i,long j,long k,long *ListOfHalos,long *NHalosPerCellStart, long *NHalosPerCellEnd, int r){
	int l,m,n;
	for (l=i-(1+r);l<=i+(1+r);l++){
	for (m=j-(1+r);m<=j+(1+r);m++){
	for (n=k-(1+r);n<=k+(1+r);n++){
		if (check_HaloR_in_cell(ihalo,X,Y,Z,R,l,m,n,ListOfHalos,NHalosPerCellStart,NHalosPerCellEnd)==1)
			return 1;

	}
	}
	}

	return 0;
}

