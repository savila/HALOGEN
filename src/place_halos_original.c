#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h> 
#include <omp.h> 
#include <string.h>

#include "place_halos.h"
 

#define MAXTRIALS (100)



/*=============================================================================
 *                             PROTOTYPES
 *=============================================================================*/


//long select_cell_rnd(long *, long *, long *); 
long select_cell(double , double *); 
//long select_heaviest_cell(long *, long *, long *); 
long select_part(long, long **,long *);
int check_HaloR_in_mesh(long,float *, float *, float * , float *,long,long,long,long**,long*);
int check_HaloR_in_cell(long ,float *, float *, float * , float *,long ,long,long,long **,long *);
double ComputeCumulative(double,double, double *, double *);

//Global Variables
long NCells,NTotCells;
float Lbox,lcell;

#ifdef MASS_OF_PARTS
  int exclude_cell(long,float , float *, float *, float *, long ,long, long);
  void exclude(long,float,float *,float *,float *,long ,long,long);
  int *excluded;
  long *Nexcluded;
#endif

/*=============================================================================
 *                             simple functions
 *=============================================================================*/

//     square()
float square(float a){
	return a*a;
}

// R_from_mass():
// takes a halo mass and transforms it to a halo radius, for a given definition in terms of density at the edge (rho, in physical units).
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
//Takes a list of halo masses (Nhalos, HaloMass), a list of particles (NTotPart,PartX,PartY,PartZ), some simulation parameters (L, mp), and user-defined parameters (Nlin,rho_ref,alpha,Malpha,Nalpha,seed)
//and returns a list of halo positions and radii (HaloX,HaloY,HaloZ,HaloR)
int place_halos(long Nend, float *HaloMass, long Nlin, long NTotPart, float *PartX, float *PartY, float *PartZ, float *PartVX, float *PartVY, float *PartVZ,float L, float rho_ref, long seed, float mp, double *alpha, double *Malpha,long Nalpha,float *HaloX, float *HaloY, float *HaloZ, float *HaloVX, float *HaloVY, float *HaloVZ,float *HaloR,long **ListOfPart, long *NPartPerCell){


fprintf(stderr,"\tThis is place_halos.c v11\n");


//Initiallising -------------------------------------------------
	long i,j,k,lin_ijk, Nmin;
	long *count,trials;
	long ihalo,ilong, ipart,i_alpha;
	double invL = 1./L;
	float Mcell,Mhalo,Mchange; 
	float R;
	time_t t0;
	int check;

	double mpart;
	double exponent;
	double TotProb;
	double *MassLeft;
	double *CumulativeProb; 
	long **ListOfHalos,  *NHalosPerCell;
	long Nstart=0,Nhalos;

	#ifdef VERB
	time_t t1,t2,t3,t4,t4_5,t5;
	float diff;
	#endif

	NCells = Nlin;
	Lbox = L;
	
	t0=time(NULL);
	NTotCells = NCells*NCells*NCells;
	
	MassLeft = malloc(NTotCells*sizeof(double));


	//Allocate memory for the arrays 
	//NPartPerCell = (long *) calloc(NTotCells,sizeof(long));
  /*	if( NPartPerCell == NULL) {
    		fprintf(stderr,"\tplace_halos(): could not allocate %ld array for NPartPerCell[]\nABORTING",NTotCells);
    		exit(-1);
	}*/
	NHalosPerCell = (long *) calloc(NTotCells,sizeof(long));
  	if(NHalosPerCell == NULL) {
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
	fprintf(stderr,"\tUsing OMP with %d threads\n",omp_get_max_threads());
	
	#ifdef MASS_OF_PARTS
	Nexcluded = (long *) calloc(NTotCells,sizeof(long));
  	if(Nexcluded == NULL) {
    		fprintf(stderr,"\tplace_halos(): could not allocate %ld array for Nexcluded[]\nABORTING",NTotCells);
    		exit(-1);
	}
 	excluded  = (int *) calloc(NTotPart, sizeof(long));
  	if(excluded == NULL) {
    		fprintf(stderr,"\tplace_halos(): could not allocate %ld array for excluded[]\nABORTING",NTotPart);
    		exit(-1);
	}
	#endif

        //Initiallise random numbers
	#ifdef VERB
        fprintf(stderr,"\tinput seed: %ld.    time0: %ld.",seed,t0);
	#endif

        if (seed>=0){
                srand(seed);
		#ifdef VERB
        	fprintf(stderr,"Used: %ld \n",seed);
		#endif
	}
        else {
                srand(t0);
        	fprintf(stderr,"Seed Used: %ld \n",t0);
	}

	mpart = (double) mp;
	//Nmin = (long)ceil(HaloMass[Nhalos-1]/mpart);
	Nmin = (long)ceil(HaloMass[Nend-1]*0.8/mpart);

	lcell = (float) L/NCells;
	#ifdef VERB
	fprintf(stderr,"\n\tParticles and Halos placed in %ld^3 cells\n",NCells);
	fprintf(stderr,"\tBOX = %f  lcell =%f   rho_ref = %e  invL %f\n",L,L/NCells,rho_ref,invL);
	fprintf(stderr,"\tNhalostart = %ld,Nhalosend = %ld,  NPart = %ld\n",Nstart, Nend, NTotPart);
	#endif
	

	#ifdef DEBUG
	fprintf(stderr,"\n\tRAND_MAX=%d\n",RAND_MAX);
	fprintf(stderr,"\tX[0] = %f Y[0] = %f Z[0] = %f\n",PartX[0],PartY[0],PartZ[0]);
	fprintf(stderr,"\tX[1] = %f Y[1] = %f Z[1] = %f\n",PartX[1],PartY[1],PartZ[1]);
	fprintf(stderr,"\tM[0] = %e \n",HaloMass[0]);
	fprintf(stderr,"\tM[1] = %e \n",HaloMass[1]);
	fprintf(stderr,"\tM[%ld] = %e \n",Nend-1,HaloMass[Nend-1]);
	fprintf(stderr,"\n\tMinimmum mass= %e. Minimum part per halo = %ld. mpart %e\n",HaloMass[Nend-1],Nmin,mpart);
	#endif	
	
	if (L/NCells<R_from_mass(HaloMass[0],rho_ref)){
		fprintf(stderr,"ERROR: cell size is smaller than the radius of the biggest halo. Please, change the number of cells\n");
		exit(0);
	}
	
	#ifdef VERB
	t1=time(NULL);
 	diff = difftime(t1,t0);
	fprintf(stderr,"\ttime of initialisation %f\n",diff);
	#endif
// ------------------------------------------------- Initiallised



#ifdef VERB
	fprintf(stderr,"\tAssigning particles to grid ...\n");
#endif


//Assign particles to grid ------------------------------------
	//count particles per cell
	/*for (ilong=0;ilong<NTotPart;ilong++) {

                if (PartX[ilong]==Lbox)
                        PartX[ilong]=0.;
                if (PartY[ilong]==Lbox)
                        PartY[ilong]=0.;
                if (PartZ[ilong]==Lbox)
                        PartZ[ilong]=0.;
                i = (long) (invL * PartX[ilong]*NCells);
                j = (long) (invL * PartY[ilong]*NCells);
                k = (long) (invL * PartZ[ilong]*NCells);
                if (i<0 || i>=NCells || j<0 || j>=NCells || k<0 || k>=NCells){
                        fprintf(stderr,"\tERROR: Particle %ld at [%f,%f,%f] seems to be out of the right box interval [0.,%f)",ilong,PartX[ilong],PartY[ilong],PartZ[ilong],L);
		}

		lin_ijk = k+j*NCells+i*NCells*NCells;
//		NPartPerCell[lin_ijk]++;
#ifdef DEBUG
		if(ilong<10 || ilong > NTotPart -10 || ilong==243666)
			fprintf(stderr,"\tipart=%ld  cell: %ld=[%ld,%ld,%ld] Parts in cell=%ld, Pos= [%f,%f,%f]\n",ilong,lin_ijk,i,j,k,NPartPerCell[lin_ijk],PartX[ilong],PartY[ilong],PartZ[ilong]);
#endif
	}*/
#ifdef VERB
	fprintf(stderr,"\t... particles counted ...\n");
	t2=time(NULL);
 	diff = difftime(t2,t1);
	fprintf(stderr,"\ttime counting %f\n",diff);
#endif
	//Alloc Enough Memory
	//ListOfPart = (long **) calloc(NCells*NCells*NCells,sizeof(long *));
	ListOfHalos = (long **) calloc(NCells*NCells*NCells,sizeof(long *));
	for (i=0;i<NCells;i++){
	for (j=0;j<NCells;j++){
	for (k=0;k<NCells;k++){
		lin_ijk = k+j*NCells+i*NCells*NCells;
		//ListOfPart[lin_ijk] = (long *) calloc(NPartPerCell[lin_ijk],sizeof(long));
		Nhalos = (long) (NPartPerCell[lin_ijk]/Nmin);
		ListOfHalos[lin_ijk] = (long *) calloc(Nhalos,sizeof(long));
		if (Nstart==0)
			MassLeft[lin_ijk] = (double) NPartPerCell[lin_ijk]*mpart; 
#ifdef ULTRADEBUG
		if (lin_ijk<10 || lin_ijk > (NCells*NCells*NCells) - 10){
			fprintf(stderr,"\tAllocated %ld (longs) in ListOfPart(%ld=[%ld,%ld,%ld])\n",NPartPerCell[lin_ijk],lin_ijk,i,j,k);
			fprintf(stderr,"\tAllocated %ld (longs) in ListOfHalos(%ld=[%ld,%ld,%ld])\n",Nhalos,lin_ijk,i,j,k);
		}
#endif		
	}	
	}
	}

#ifdef VERB
	fprintf(stderr,"\t... memory allocated ...\n");
	t3=time(NULL);
 	diff = difftime(t3,t2);
	fprintf(stderr,"\ttime allocating %f\n",diff);
#endif 
/*
	#pragma omp parallel for private(ilong,i,k,j,lin_ijk) shared(NTotPart,invL,NCells,ListOfPart,count,PartZ,PartX,PartY,L,stderr) default(none)
	for (ilong=0;ilong<NTotPart;ilong++) {
		i = (long) (invL * PartX[ilong]*NCells);
		j = (long) (invL * PartY[ilong]*NCells);
		k = (long) (invL * PartZ[ilong]*NCells);
                if (i<0 || i>=NCells || j<0 || j>=NCells || k<0 || k>=NCells){
                        fprintf(stderr,"\tERROR: Particle %ld at [%f,%f,%f] seems to be out of the right box interval [0.,%f)",ilong,PartX[ilong],PartY[ilong],PartZ[ilong],L);
		}
		lin_ijk = k+j*NCells+i*NCells*NCells;
		ListOfPart[lin_ijk][count[lin_ijk]] = ilong;
		count[lin_ijk]++;
	}

#ifdef VERB
	fprintf(stderr,"\t...particles assigned, now haloes...\n");
#endif */
	for (ihalo=0;ihalo<Nstart;ihalo++){
		i = (long) (invL * HaloX[ihalo]*NCells);
		j = (long) (invL * HaloY[ihalo]*NCells);
		k = (long) (invL * HaloZ[ihalo]*NCells);
		i=check_limit(i,NCells);
		j=check_limit(j,NCells);
		k=check_limit(k,NCells);
		lin_ijk = k+j*NCells+i*NCells*NCells;
		ListOfHalos[lin_ijk][NHalosPerCell[lin_ijk]] = ihalo;
		NHalosPerCell[lin_ijk]++;
	}

#ifdef DEBUG
        fprintf(stderr,"\tMass_cell[0]=%e",MassLeft[0]);
	fprintf(stderr,"\t Mass Function\n");
	for (ihalo=0;ihalo<15;ihalo++){
		fprintf(stderr,"\thalo %ld: ",ihalo);
		fprintf(stderr,"M=%e\n",HaloMass[ihalo]);
	}
#endif


#ifdef VERB
	fprintf(stderr,"\t ...done\n\n");
	t4=time(NULL);
 	diff = difftime(t4,t3);
	fprintf(stderr,"\ttime of the actual assignment %f\n",diff);
	fprintf(stderr,"\tComputing probabilities...\n");
#endif

//----------------------------------- Particles and haloes assigned to grid



//Computing Cumulative Probability -----------------------------
	
	//find the right alpha
	Mhalo = HaloMass[Nstart];
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
	//compute the probability
	TotProb = ComputeCumulative(exponent, mpart, MassLeft, CumulativeProb);
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
	for (ihalo=Nstart;ihalo<Nend;ihalo++){

		#ifdef DEBUG
		fprintf(stderr,"\n\t- Halo %ld ",ihalo);
		#endif
		#ifdef VERB
		if ((ihalo%1000000)==0)
			fprintf(stderr,"\t%ld million haloes done\n",(ihalo/1000000));
		#endif
		//Check whether or not, a change of alpha is needed for this halo mass 		
		Mhalo= HaloMass[ihalo];
//		if (ihalo<=-1)
//			fprintf(stderr,"\tMhalo=%e\n",Mhalo);

		while (Mhalo < Mchange){//if so search the right alpha, and recompute probabilities
			i_alpha++;		
			if (i_alpha==Nalpha){
				fprintf(stderr,"\tERROR: No M_alpha low enough found\n");
				exit(0);
			}
			Mchange = Malpha[i_alpha];
			exponent = alpha[i_alpha];
			TotProb=ComputeCumulative(exponent, mpart, MassLeft, CumulativeProb);
		#ifdef VERB
        		fprintf(stderr,"\n\tUsing alpha_%ld=%f for M>%e\n",i_alpha,exponent,Mchange);
		#endif
		}
//		if (ihalo<=-1)
//			fprintf(stderr,"\talpha=%f\n",exponent);


		do {	
		  //First, choose a cell	
		  #ifndef RANKED				
		  lin_ijk = select_cell(TotProb, CumulativeProb);
//		  if (ihalo<=-1)
//			fprintf(stderr,"\tlin_ijk=%ld  ",lin_ijk);
		 
		  k=lin_ijk%(NCells);
		  j=((lin_ijk-k)/NCells)%NCells;
	  	  i=(lin_ijk-k-j*NCells)/(NCells*NCells);
		  if (ihalo<=-1)
			fprintf(stderr," = [%ld,%ld,%ld]\n",i,j,k);
		  #else
		  lin_ijk=select_heaviest_cell(&i,&j,&k);		  
		  #endif

		  trials=0;


		  //Second, choose a particle in that cell
		  do {
			ipart = select_part(lin_ijk,ListOfPart, NPartPerCell);		
//			if (ihalo<=-1)
//				fprintf(stderr,"\tipart=%ld\n",ipart);
               		HaloX[ihalo] = PartX[ipart];
               		HaloY[ihalo] = PartY[ipart];
               		HaloZ[ihalo] = PartZ[ipart];
               		HaloVX[ihalo] = PartVX[ipart];
               		HaloVY[ihalo] = PartVY[ipart];
               		HaloVZ[ihalo] = PartVZ[ipart];
			R=R_from_mass(HaloMass[ihalo],rho_ref);
			HaloR[ihalo]= R;
//			if (ihalo<=-1)
//				fprintf(stderr,"\thalo: [[%f,%f,%f],[%f,%f,%f],%f]\n",HaloX[ihalo],HaloY[ihalo],HaloZ[ihalo],HaloVX[ihalo],HaloVY[ihalo],HaloVZ[ihalo],HaloR[ihalo]);
			#ifdef NO_EXCLUSION
			check = 0;
			#else
			//Third, check that is not overlapping a previous halo
			check = check_HaloR_in_mesh(ihalo,HaloX,HaloY,HaloZ,HaloR,i,j,k,ListOfHalos,NHalosPerCell);
			#endif
//			if (ihalo<=-1)
//				fprintf(stderr,"\tCHECK=%d\n",check);
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
				TotProb=ComputeCumulative(exponent, mpart, MassLeft, CumulativeProb);
				break;
			}
		  } while (check==1);//If the particle was excluded, try another one in the same cell

	        } while(check==1); //if reached MAXTRIALS, select another cell
		//Particle chosen!
		
		//mass in cell before assignment
                Mcell=MassLeft[lin_ijk];
//		if (ihalo<=-1)
//			fprintf(stderr,"\tMbefore=%e\n",Mcell);

		
		  #ifndef MASS_OF_PARTS 
                  if (Mcell>HaloMass[ihalo])
			MassLeft[lin_ijk] -= Mhalo; 
                  else
                        MassLeft[lin_ijk] = 0.;
		  #else
			exclude(ipart,R,PartX,PartY,PartZ,i,j,k);
		  #endif
////		if (ihalo<=-1)
//			fprintf(stderr,"\tMafter=%e\n",MassLeft[lin_ijk]);



	#ifndef NO_MASS_CONSERVATION
		if (ihalo<=-1)
			fprintf(stderr,"\tOld version!\n");
		double ProbDiff = pow(MassLeft[lin_ijk]/mpart,exponent)-pow(Mcell/mpart,exponent);

		#ifdef DEBUG
		fprintf(stderr,"\n \tassigned to cell %ld=[%ld,%ld,%ld]\n\t Before: Mcell=%e, CProbCell=%e,  TotProb=%e. ",lin_ijk,i,j,k,Mcell,CumulativeProb[lin_ijk],TotProb);
		#endif
		
		#ifndef MASS_OF_PARTS
		long icell;
		  //Substract the Probability difference from the array (only affected those cells after the selected one)
                  #pragma omp parallel for private(icell) shared(CumulativeProb,robDiff,NTotCells,lin_ijk) default(none)
                  for(icell=lin_ijk;icell<NTotCells;icell++){
                        CumulativeProb[icell]+=ProbDiff;
                  }
                  //TotProb+=ProbDiff;
                  TotProb=CumulativeProb[NCells*NCells*NCells-1];
		#endif
	#endif



		#ifdef DEBUG
		fprintf(stderr," After: Mcell=%e, CProbCell=%e, TotProb=%e.   , Mhalo=%e. CProb[last]=%e\n",MassLeft[lin_ijk],CumulativeProb[lin_ijk],TotProb,Mhalo,CumulativeProb[NTotCells-1]);
		#endif
	
		#ifdef DEBUG
		fprintf(stderr,"\thalo %ld assigned to particle %ld at [%f,%f,%f]. R= %f, M= %e\n",ihalo,ipart,HaloX[ihalo],HaloY[ihalo],HaloZ[ihalo],R,Mhalo);
		#endif

//		if (ihalo<=-1)
//			fprintf(stderr,"\tplace ihalo=%ld",ihalo);
		ListOfHalos[lin_ijk][NHalosPerCell[lin_ijk]]=ihalo;
//		if (ihalo<=-1)
//			fprintf(stderr,"\t check: %ld\n",ListOfHalos[lin_ijk][NHalosPerCell[lin_ijk]]);
		NHalosPerCell[lin_ijk]++;
//		if (ihalo<=-1)
//			fprintf(stderr,"\t Nhalospercell[%ld]= %ld\n",lin_ijk,NHalosPerCell[lin_ijk]);
	}//for(ihalo=Nstart:Nend)
//----------------------------------- Haloes Placed

#ifdef VERB
	t5=time(NULL);
 	diff = difftime(t5,t4_5);
	fprintf(stderr,"\ttime placing %f\n",diff);
 	diff = difftime(t5,t0);
	fprintf(stderr,"\ttotal time in place_halos.c %f\n",diff);
	fprintf(stderr,"\n\tPlacement done!!!\n");
#endif
	free(NHalosPerCell);
        free(count); 
	//free(NPartPerCell);
        free(CumulativeProb);
	free(MassLeft);
        for (i=0;i<NCells;i++){
                for (j=0;j<NCells;j++){
                        for (k=0;k<NCells;k++){
                                lin_ijk = k+j*NCells+i*NCells*NCells;
                                //free(ListOfPart[lin_ijk]);
                                //free(ListOfHalos[lin_ijk]);
                        }
                }
        }

        //free(ListOfPart);
        free(ListOfHalos);
#ifdef MASS_OF_PARTS
	free(excluded); free(Nexcluded);
#endif
		fprintf(stderr," e ");
	return 0;
}
//end of place_halos()



//ComputeCumulative():
//it takes the user-defined alpha (exponent) and (re-)computes the cumulative probability for the cells in the grid
//from the masses in the cells (MassLeft)
double ComputeCumulative(double alpha, double mpart, double *MassLeft, double *CumulativeProb){
	long i;
	double TotProb = 0.;
	for(i=0;i<NTotCells;i++){
		CumulativeProb[i] = TotProb + pow(MassLeft[i]/mpart,alpha);
		TotProb = CumulativeProb[i];
	}
	return TotProb;
}

//select_part():
//randomly selects a particle from a cell
long select_part(long ijk,long **ListOfPart,long *NPartPerCell){
	long i_rnd,ipart;

	i_rnd = (long) (NPartPerCell[ijk] * ((double)rand()/(RAND_MAX+1.0)));
	#ifdef DEBUG
	fprintf(stderr,"irnd: %ld/%ld ",i_rnd,NPartPerCell[ijk]);
	#endif
        ipart = ListOfPart[ijk][i_rnd];

	return ipart;
}

/*
//select_heaviestcell()
long select_heaviest_cell(long *x, long *y, long *z) {
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
*/
/*
long select_cell_rnd(long *x, long *y, long *z) {
	(*x) = rand()%NCells;
	(*y) = rand()%NCells;
	(*z) = rand()%NCells;
	return (*z)+(*y)*NCells+(*x)*NCells*NCells;
}
*/

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
int check_HaloR_in_cell(long ipart,float *PartX, float *PartY, float *PartZ, float *PartR, long i,long j, long k,long **ListOfHalos,long *NHalosPerCell){
        long jpart,jj;
        double X=PartX[ipart],Y=PartY[ipart],Z=PartZ[ipart];
	#ifndef ONLYBIG
	double R=PartR[ipart];
	#endif 

//The cell passed, might not exists, but be a virtual one (for periodic conditions): check, and correct:
//#ifdef PERIODIC
                if (i==-1){
                        i = NCells -1;
                        X = X + Lbox;
                }
                else if (i==NCells){
                        i = 0;
                        X = X - Lbox;
                }

                if (j==-1){
                        j = NCells -1;
                        Y = Y + Lbox;
                }
                else if (j==NCells){
                        j = 0;
                        Y = Y - Lbox;
                }
                if (k==-1){
                        k = NCells -1;
                        Z = Z + Lbox;
                }
                else if (k==NCells){
                        k = 0;
                        Z = Z - Lbox;
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
	for (jj=0; jj<NHalosPerCell[lin_ijk]; jj++){
#ifdef ULTRADEBUG
		fprintf(stderr,"jj=%ld/%ld ",jj,NHalosPerCell[lin_ijk]);
#endif
		jpart=ListOfHalos[lin_ijk][jj];
#ifdef ULTRADEBUG
		fprintf(stderr,"jpart=%ld ",jpart);
#endif	
		//Check for overlapping
		#ifdef ONLYBIG
		if ((square(X-PartX[jpart])+square(Y-PartY[jpart])+square(Z-PartZ[jpart]))<square(PartR[jpart])) 
		#else
		if ((square(X-PartX[jpart])+square(Y-PartY[jpart])+square(Z-PartZ[jpart]))<square(R+PartR[jpart])) 
		#endif
		{
#ifdef DEBUG
                        fprintf(stderr,"\nChecking cell [%ld,%ld,%ld] for halo %ld",i,j,k,ipart);
                        fprintf(stderr," lin_ijk= %ld.",lin_ijk);
                        fprintf(stderr,"jj=%ld/%ld ",jj,NHalosPerCell[lin_ijk]);
                        fprintf(stderr,"jhalo=%ld ",jpart);
                        fprintf(stderr,"refused!\n");
#endif
			return 1;
		}
	}
	return 0;
  }
  else {
	fprintf(stderr,"WARNING: Computing distances outside the box\n");
	return 0;
  }
}



//cheack_haloR_in_mesh():
//checks if there is any collision between the halo just been placed and any previous one (going through all the neighbour cells)
//returns 1 if there is any collision, returns 0 for no collision
int check_HaloR_in_mesh(long ihalo,float *X, float *Y, float *Z , float *R,long i,long j,long k,long **ListOfHalos,long *NHalosPerCell){
	int l,m,n;
	for (l=i-1;l<=i+1;l++){
	for (m=j-1;m<=j+1;m++){
	for (n=k-1;n<=k+1;n++){
		if (check_HaloR_in_cell(ihalo,X,Y,Z,R,l,m,n,ListOfHalos,NHalosPerCell)==1)
			return 1;

	}
	}
	}

	return 0;
}

