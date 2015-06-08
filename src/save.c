/*********************************************************************************** 
This code is a part of 2LPTic (http://arxiv.org/abs/astro-ph/9711187)
adapted to HALOGEN (http://arxiv.org/abs/1412.5228)
The original code 2LPTic was developped by 
Roman Scoccimarro and Sebastian Pueblas (http://cosmo.nyu.edu/roman/2LPT/)
************************************************************************************ */

#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>

#include "allvars.h"
#include "proto.h"

#define ESTIMATED_FRACTION 20.0


//#ifdef ONLY_2LPT
int distribute_part(int Nlin, int Nx, long ***ListOfPart,long **NPartPerCell){
	int NextSlice, PrevSlice, NfromPrevSlice,NfromNextSlice;
	int NPartPrevSlice=0,NPartNextSlice=0; 
	float *NextSliceX, *NextSliceY, *NextSliceZ, *NextSliceVX, *NextSliceVY, *NextSliceVZ, *FromPrevSliceX, *FromPrevSliceY, *FromPrevSliceZ, *FromPrevSliceVX,*FromPrevSliceVY,*FromPrevSliceVZ;
	float *PrevSliceX, *PrevSliceY, *PrevSliceZ, *PrevSliceVX, *PrevSliceVY, *PrevSliceVZ, *FromNextSliceX, *FromNextSliceY, *FromNextSliceZ, *FromNextSliceVX, *FromNextSliceVY, *FromNextSliceVZ;
	long *count;
	double invL = 1.0/Box;
	long ipart=0;
	int i,j,k,ilocal,lin_ijk,ii;
  	MPI_Status status;

	int NTransferPart = (int) (GlassTileFac*GlassTileFac * ESTIMATED_FRACTION);
	
	#ifdef DEBUG
	fprintf(stderr,"\t Note: %f%% of extra memory used\n", (float) NTransferPart*2.0/(GlassTileFac*GlassTileFac*GlassTileFac));
	#endif

	NextSliceX = (float *)malloc(NTransferPart*sizeof(float));
	NextSliceY = (float *)malloc(NTransferPart*sizeof(float));
	NextSliceZ = (float *)malloc(NTransferPart*sizeof(float));
	NextSliceVX = (float *)malloc(NTransferPart*sizeof(float));
	NextSliceVY = (float *)malloc(NTransferPart*sizeof(float));
	NextSliceVZ = (float *)malloc(NTransferPart*sizeof(float));
	
	if (NextSliceX==NULL || NextSliceY==NULL || NextSliceZ==NULL || NextSliceVX==NULL || NextSliceVY==NULL || NextSliceVZ==NULL)  {
		fprintf(stderr,"ERROR: Couldnt allocate NextSliceCoords\n");
		return -1;
	}	
	PrevSliceX = (float *)malloc(NTransferPart*sizeof(float));
	PrevSliceY = (float *)malloc(NTransferPart*sizeof(float));
	PrevSliceZ = (float *)malloc(NTransferPart*sizeof(float));
	PrevSliceVX = (float *)malloc(NTransferPart*sizeof(float));
	PrevSliceVY = (float *)malloc(NTransferPart*sizeof(float));
	PrevSliceVZ = (float *)malloc(NTransferPart*sizeof(float));
	
	if (PrevSliceX==NULL || PrevSliceY==NULL || PrevSliceZ==NULL || PrevSliceVX==NULL || PrevSliceVY==NULL || PrevSliceVZ==NULL)  {
		fprintf(stderr,"ERROR: Couldnt allocate PrevSliceCoords\n");
		return -1;
	}	


	#ifdef DEBUG
	fprintf(stderr,"\tL=%f, invL=%f, Nlin=%d\n",Box,invL,Nlin);
	#endif

	(*NPartPerCell) = (long *) calloc(Nlin*Nlin*Nx,sizeof(long));
	if( *NPartPerCell == NULL) {
		fprintf(stderr,"\tplace_halos(): could not allocate %d array for NPartPerCell[]\nABORTING",Nlin*Nlin*Nx);
		exit(-1);
	}
	#ifdef DEBUG
	fprintf(stderr,"\t...NPartPerCell allocaled, continue...\n");
	#endif
	(*ListOfPart) = (long **) calloc(Nlin*Nlin*Nx,sizeof(long *));
	if( *ListOfPart == NULL) {
		fprintf(stderr,"\tplace_halos(): could not allocate %d array for LisOfPart[]\nABORTING",Nlin*Nlin*Nx);
		exit(-1);
	}
	#ifdef DEBUG
	fprintf(stderr,"\t...ListOFPART allocaled, continue...\n");
	#endif

	count = (long *) calloc(Nlin*Nlin*Nx,sizeof(long));
	if( count == NULL) {
		fprintf(stderr,"\tplace_halos(): could not allocate %d array for count[]\nABORTING",Nlin*Nlin*Nx);
		exit(-1);
	}
	#ifdef DEBUG
	fprintf(stderr,"\t Counters allocated \n\n");
	#endif


	NPartPrevSlice=0;		
	NPartNextSlice=0;		
	//Counting Particles in cells
	//t3=time(NULL);
	//#omp
	for (ipart=0;ipart<NumPart;ipart++) {
		#ifdef ULTRADEBUG
		if(ThisTask==0) fprintf(stderr,"Task %d. part %d/%ld\n",ThisTask, ipart,NumPart);
		#endif
		if (P[ipart].Pos[0]==Box)
			P[ipart].Pos[0]=0.;
		if (P[ipart].Pos[1]==Box)
			P[ipart].Pos[1]=0.;
		if (P[ipart].Pos[2]==Box)
			P[ipart].Pos[2]=0.;
		i = (long) (invL * P[ipart].Pos[0]*Nlin);
		ilocal = (long) i - Nx*ThisTask;
		j = (long) (invL * P[ipart].Pos[1]*Nlin);
		k = (long) (invL * P[ipart].Pos[2]*Nlin);

		if ((j<0) || (j>=Nlin) || (k<0) || (k>=Nlin)){
			fprintf(stderr,"\tERROR: Particle %ld at [%f,%f,%f]->[%d,%d,%d] seems to be out of the right box interval [0.,%f)\n",ipart,P[ipart].Pos[0],P[ipart].Pos[1],P[ipart].Pos[2],i,j,k,Box);
			exit(0);
		}

		if (ilocal>=0 && ilocal<Nx) {
			lin_ijk = k+j*Nlin+ilocal*Nlin*Nlin;
			(*NPartPerCell)[lin_ijk]++;

			//WhichSlice[ipart]=0;		
		}
		else if ((ilocal<0 && ThisTask!=0 && ilocal>=-Nx) || (ThisTask==0 && i<Nlin && i>=(Nlin-Nx))){
			//WhichSlice[ipart]=-1;		
			//PrevSlice[NPartPrevSlice]=ipart;		
			PrevSliceX[NPartPrevSlice]=P[ipart].Pos[0];		
			PrevSliceY[NPartPrevSlice]=P[ipart].Pos[1];		
			PrevSliceZ[NPartPrevSlice]=P[ipart].Pos[2];		
			PrevSliceVX[NPartPrevSlice]=P[ipart].Vel[0];		
			PrevSliceVY[NPartPrevSlice]=P[ipart].Vel[1];		
			PrevSliceVZ[NPartPrevSlice]=P[ipart].Vel[2];		
			NPartPrevSlice++;		
		}
		else if ((ilocal>=Nx && ilocal<2*Nx && ThisTask!=NTask-1) || (ThisTask==NTask-1 && i>=0 && i<Nx)){
			//WhichSlice[ipart]=1;		
			//NextSlice[NPartNextSlice]=ipart;
			NextSliceX[NPartNextSlice]= P[ipart].Pos[0];		
			NextSliceY[NPartNextSlice]= P[ipart].Pos[1];		
			NextSliceZ[NPartNextSlice]= P[ipart].Pos[2];		
			NextSliceVX[NPartNextSlice]=P[ipart].Vel[0];		
			NextSliceVY[NPartNextSlice]=P[ipart].Vel[1];		
			NextSliceVZ[NPartNextSlice]=P[ipart].Vel[2];		
			NPartNextSlice++;
		}
		else{
			fprintf(stderr,"ERROR: Something went wrong, maybe a particle moved accross 2 slices:\n");
			//fprintf(stderr,"\t- Part: %ld/%d. [i,j,k]=[%d,%d,%d], ilocal=%d, [X,Y,Z]=[%f,%f,%f], Task: %d/%d\n",ipart,NumPart,i,j,k,ilocal,partX[ipart],partY[ipart],partZ[ipart], ThisTask, NTask);
			return -1;
		}
	}


	//Send to next slice
	
	if (ThisTask!=NTask-1)
		NextSlice = ThisTask+1;
	else 
		NextSlice = 0;

	if (ThisTask!=0)
		PrevSlice = ThisTask-1;
	else 
		PrevSlice = NTask-1;


	#ifdef VERB
		#ifndef DEBUG
		if (ThisTask==0)
		#endif
			fprintf(stderr,"\tTask %d sending %d parts to task %d, and %d parts to task %d (%f%% stay)\n",ThisTask,NPartPrevSlice,PrevSlice,NPartNextSlice,NextSlice,(float)(NumPart-NPartPrevSlice-NPartNextSlice)/(NumPart)*100.0);
	#endif



	if (ThisTask%2==0){
  		#ifdef DEBUG
  		fprintf(stderr,"\tTask %d, sending %d particles to Task %d\n",ThisTask,NPartNextSlice,NextSlice);
  		#endif
  		MPI_Send( &NPartNextSlice, 1, MPI_INT, NextSlice, 1234, MPI_COMM_WORLD); 
  		#ifdef DEBUG
  		fprintf(stderr,"\tTask %d, sending x[%d]=%f \n",ThisTask,NPartNextSlice-1,NextSliceX[NPartNextSlice-1]);
  		#endif
  		MPI_Send( &NextSliceX[0], NPartNextSlice, MPI_FLOAT, NextSlice, 1111, MPI_COMM_WORLD);
  		#ifdef DEBUG
  		fprintf(stderr,"\tTask %d, sending y[%d]=%f \n",ThisTask,NPartNextSlice-1,NextSliceY[NPartNextSlice-1]);
  		#endif
  		MPI_Send( &NextSliceY[0], NPartNextSlice, MPI_FLOAT, NextSlice, 1112, MPI_COMM_WORLD);
  		#ifdef DEBUG
  		fprintf(stderr,"\tTask %d, sending z[%d]=%f \n",ThisTask,NPartNextSlice-1,NextSliceZ[NPartNextSlice-1]);
  		#endif
  		MPI_Send( &NextSliceZ[0], NPartNextSlice, MPI_FLOAT, NextSlice, 1113, MPI_COMM_WORLD);
  		#ifdef DEBUG
  		fprintf(stderr,"\tTask %d, sending vx[%d]=%f \n",ThisTask,NPartNextSlice-1,NextSliceVX[NPartNextSlice-1]);
  		#endif
  		MPI_Send( &NextSliceVX[0], NPartNextSlice, MPI_FLOAT, NextSlice, 1114, MPI_COMM_WORLD);
  		#ifdef DEBUG
  		fprintf(stderr,"\tTask %d, sending vy[%d]=%f \n",ThisTask,NPartNextSlice-1,NextSliceVY[NPartNextSlice-1]);
  		#endif
  		MPI_Send( &NextSliceVY[0], NPartNextSlice, MPI_FLOAT, NextSlice, 1115, MPI_COMM_WORLD);
  		#ifdef DEBUG
  		fprintf(stderr,"\tTask %d, sending vz[%d]=%f \n",ThisTask,NPartNextSlice-1,NextSliceVZ[NPartNextSlice-1]);
  		#endif
  		MPI_Send( &NextSliceVZ[0], NPartNextSlice, MPI_FLOAT, NextSlice, 1116, MPI_COMM_WORLD);


		MPI_Recv(&NfromPrevSlice, 1, MPI_INT, PrevSlice, 1234, MPI_COMM_WORLD, &status);

  		#ifdef DEBUG
		fprintf(stderr,"\tTask %d receiving %d particles from Task %d\n",ThisTask,NfromPrevSlice,PrevSlice);
  		#endif

		FromPrevSliceX = (float *) malloc(NfromPrevSlice*sizeof(float));
		FromPrevSliceY = (float *) malloc(NfromPrevSlice*sizeof(float));
		FromPrevSliceZ = (float *) malloc(NfromPrevSlice*sizeof(float));
		FromPrevSliceVX = (float *) malloc(NfromPrevSlice*sizeof(float));
		FromPrevSliceVY = (float *) malloc(NfromPrevSlice*sizeof(float));
		FromPrevSliceVZ = (float *) malloc(NfromPrevSlice*sizeof(float));
		MPI_Recv( &FromPrevSliceX[0], NfromPrevSlice, MPI_FLOAT, PrevSlice, 1111, MPI_COMM_WORLD, &status);
		MPI_Recv( &(FromPrevSliceY[0]), NfromPrevSlice, MPI_FLOAT, PrevSlice, 1112, MPI_COMM_WORLD, &status);
		MPI_Recv( &(FromPrevSliceZ[0]), NfromPrevSlice, MPI_FLOAT, PrevSlice, 1113, MPI_COMM_WORLD, &status);
		MPI_Recv( &(FromPrevSliceVX[0]), NfromPrevSlice, MPI_FLOAT, PrevSlice, 1114, MPI_COMM_WORLD, &status);
		MPI_Recv( &(FromPrevSliceVY[0]), NfromPrevSlice, MPI_FLOAT, PrevSlice, 1115, MPI_COMM_WORLD, &status);
		MPI_Recv( &(FromPrevSliceVZ[0]), NfromPrevSlice, MPI_FLOAT, PrevSlice, 1116, MPI_COMM_WORLD, &status);
		#ifdef DEBUG
  		fprintf(stderr,"\tTask %d, received Vy[%d]=%f  from task %d\n",ThisTask,NfromPrevSlice-1,FromPrevSliceVY[NfromPrevSlice-1],PrevSlice);
		#endif
		#ifdef DEBUG
 		fprintf(stderr,"\tTask %d, sending %d particles to task %d\n",ThisTask,NPartPrevSlice,PrevSlice);
  		#endif
  		MPI_Send( &NPartPrevSlice, 1, MPI_INT, PrevSlice, 2234, MPI_COMM_WORLD); 
  		#ifdef DEBUG
  		fprintf(stderr,"\tTask %d, sending x[%d]=%f \n",ThisTask,NPartPrevSlice-1,PrevSliceX[NPartPrevSlice-1]);
  		#endif
  		MPI_Send( &PrevSliceX[0], NPartPrevSlice, MPI_FLOAT, PrevSlice, 2111, MPI_COMM_WORLD);
  		#ifdef DEBUG
  		fprintf(stderr,"\tTask %d, sending y[%d]=%f \n",ThisTask,NPartPrevSlice-1,PrevSliceY[NPartPrevSlice-1]);
  		#endif
  		MPI_Send( &PrevSliceY[0], NPartPrevSlice, MPI_FLOAT, PrevSlice, 2112, MPI_COMM_WORLD);
  		#ifdef DEBUG
  		fprintf(stderr,"\tTask %d, sending z[%d]=%f \n",ThisTask,NPartPrevSlice-1,PrevSliceZ[NPartPrevSlice-1]);
  		#endif
  		MPI_Send( &PrevSliceZ[0], NPartPrevSlice, MPI_FLOAT, PrevSlice, 2113, MPI_COMM_WORLD);
  		#ifdef DEBUG
  		fprintf(stderr,"\tTask %d, sending vx[%d]=%f \n",ThisTask,NPartPrevSlice-1,PrevSliceVX[NPartPrevSlice-1]);
  		#endif
  		MPI_Send( &PrevSliceVX[0], NPartPrevSlice, MPI_FLOAT, PrevSlice, 2114, MPI_COMM_WORLD);
  		#ifdef DEBUG
  		fprintf(stderr,"\tTask %d, sending vy[%d]=%f \n",ThisTask,NPartPrevSlice-1,PrevSliceVY[NPartPrevSlice-1]);
  		#endif
  		MPI_Send( &PrevSliceVY[0], NPartPrevSlice, MPI_FLOAT, PrevSlice, 2115, MPI_COMM_WORLD);
  		#ifdef DEBUG
  		fprintf(stderr,"\tTask %d, sending vz[%d]=%f \n",ThisTask,NPartPrevSlice-1,PrevSliceVZ[NPartPrevSlice-1]);
  		#endif
  		MPI_Send( &PrevSliceVZ[0], NPartPrevSlice, MPI_FLOAT, PrevSlice, 2116, MPI_COMM_WORLD);

		


		MPI_Recv(&NfromNextSlice, 1, MPI_INT, NextSlice, 2234, MPI_COMM_WORLD, &status);
  		#ifdef DEBUG
		fprintf(stderr,"\tTask %d receiving %d particles from Task %d\n",ThisTask,NfromNextSlice,NextSlice);
		#endif
		FromNextSliceX = (float *) malloc(NfromNextSlice*sizeof(float));
		FromNextSliceY = (float *) malloc(NfromNextSlice*sizeof(float));
		FromNextSliceZ = (float *) malloc(NfromNextSlice*sizeof(float));
		FromNextSliceVX = (float *) malloc(NfromNextSlice*sizeof(float));
		FromNextSliceVY = (float *) malloc(NfromNextSlice*sizeof(float));
		FromNextSliceVZ = (float *) malloc(NfromNextSlice*sizeof(float));
		MPI_Recv( &FromNextSliceX[0], NfromNextSlice, MPI_FLOAT, NextSlice, 2111, MPI_COMM_WORLD, &status);
		MPI_Recv( &(FromNextSliceY[0]), NfromNextSlice, MPI_FLOAT, NextSlice, 2112, MPI_COMM_WORLD, &status);
		MPI_Recv( &(FromNextSliceZ[0]), NfromNextSlice, MPI_FLOAT, NextSlice, 2113, MPI_COMM_WORLD, &status);
		MPI_Recv( &(FromNextSliceVX[0]), NfromNextSlice, MPI_FLOAT, NextSlice, 2114, MPI_COMM_WORLD, &status);
		MPI_Recv( &(FromNextSliceVY[0]), NfromNextSlice, MPI_FLOAT, NextSlice, 2115, MPI_COMM_WORLD, &status);
		MPI_Recv( &(FromNextSliceVZ[0]), NfromNextSlice, MPI_FLOAT, NextSlice, 2116, MPI_COMM_WORLD, &status);
		#ifdef DEBUG
  		//sleep(2);
  		fprintf(stderr,"\tTask %d, received Vy[%d]=%f  from task %d\n",ThisTask,NfromNextSlice-1,FromNextSliceVY[NfromNextSlice-1],NextSlice);
  		//sleep(2);
		#endif
	}
	else {
		MPI_Recv(&NfromPrevSlice, 1, MPI_INT, PrevSlice, 1234, MPI_COMM_WORLD, &status);

  		#ifdef DEBUG
		fprintf(stderr,"\tTask %d receiving %d particles from Task %d\n",ThisTask,NfromPrevSlice,PrevSlice);
		#endif
		FromPrevSliceX = (float *) malloc(NfromPrevSlice*sizeof(float));
		FromPrevSliceY = (float *) malloc(NfromPrevSlice*sizeof(float));
		FromPrevSliceZ = (float *) malloc(NfromPrevSlice*sizeof(float));
		FromPrevSliceVX = (float *) malloc(NfromPrevSlice*sizeof(float));
		FromPrevSliceVY = (float *) malloc(NfromPrevSlice*sizeof(float));
		FromPrevSliceVZ = (float *) malloc(NfromPrevSlice*sizeof(float));
		MPI_Recv( &FromPrevSliceX[0], NfromPrevSlice, MPI_FLOAT, PrevSlice, 1111, MPI_COMM_WORLD, &status);
		MPI_Recv( &(FromPrevSliceY[0]), NfromPrevSlice, MPI_FLOAT, PrevSlice, 1112, MPI_COMM_WORLD, &status);
		MPI_Recv( &(FromPrevSliceZ[0]), NfromPrevSlice, MPI_FLOAT, PrevSlice, 1113, MPI_COMM_WORLD, &status);
		MPI_Recv( &(FromPrevSliceVX[0]), NfromPrevSlice, MPI_FLOAT, PrevSlice, 1114, MPI_COMM_WORLD, &status);
		MPI_Recv( &(FromPrevSliceVY[0]), NfromPrevSlice, MPI_FLOAT, PrevSlice, 1115, MPI_COMM_WORLD, &status);
		MPI_Recv( &(FromPrevSliceVZ[0]), NfromPrevSlice, MPI_FLOAT, PrevSlice, 1116, MPI_COMM_WORLD, &status);
		#ifdef DEBUG
  		fprintf(stderr,"\tTask %d, received Vy[%d]=%f  from task %d\n",ThisTask,NfromPrevSlice-1,FromPrevSliceVY[NfromPrevSlice-1],PrevSlice);
		#endif

  		#ifdef DEBUG
  		fprintf(stderr,"\tTask %d, sending %d particles to Task %d\n",ThisTask,NPartNextSlice,NextSlice);
  		#endif
  		MPI_Send( &NPartNextSlice, 1, MPI_INT, NextSlice, 1234, MPI_COMM_WORLD); 
  		#ifdef DEBUG
  		fprintf(stderr,"\tTask %d, sending x[%d]=%f \n",ThisTask,NPartNextSlice-1,NextSliceX[NPartNextSlice-1]);
  		#endif
  		MPI_Send( &NextSliceX[0], NPartNextSlice, MPI_FLOAT, NextSlice, 1111, MPI_COMM_WORLD);
  		#ifdef DEBUG
  		fprintf(stderr,"\tTask %d, sending y[%d]=%f \n",ThisTask,NPartNextSlice-1,NextSliceY[NPartNextSlice-1]);
  		#endif
  		MPI_Send( &NextSliceY[0], NPartNextSlice, MPI_FLOAT, NextSlice, 1112, MPI_COMM_WORLD);
  		#ifdef DEBUG
  		fprintf(stderr,"\tTask %d, sending z[%d]=%f \n",ThisTask,NPartNextSlice-1,NextSliceZ[NPartNextSlice-1]);
  		#endif
  		MPI_Send( &NextSliceZ[0], NPartNextSlice, MPI_FLOAT, NextSlice, 1113, MPI_COMM_WORLD);
  		#ifdef DEBUG
  		fprintf(stderr,"\tTask %d, sending vx[%d]=%f \n",ThisTask,NPartNextSlice-1,NextSliceVX[NPartNextSlice-1]);
  		#endif
  		MPI_Send( &NextSliceVX[0], NPartNextSlice, MPI_FLOAT, NextSlice, 1114, MPI_COMM_WORLD);
  		#ifdef DEBUG
  		fprintf(stderr,"\tTask %d, sending vy[%d]=%f \n",ThisTask,NPartNextSlice-1,NextSliceVY[NPartNextSlice-1]);
  		#endif
  		MPI_Send( &NextSliceVY[0], NPartNextSlice, MPI_FLOAT, NextSlice, 1115, MPI_COMM_WORLD);
  		#ifdef DEBUG
  		fprintf(stderr,"\tTask %d, sending vz[%d]=%f \n",ThisTask,NPartNextSlice-1,NextSliceVZ[NPartNextSlice-1]);
  		#endif
  		MPI_Send( &NextSliceVZ[0], NPartNextSlice, MPI_FLOAT, NextSlice, 1116, MPI_COMM_WORLD);


		MPI_Recv(&NfromNextSlice, 1, MPI_INT, NextSlice, 2234, MPI_COMM_WORLD, &status);
  		#ifdef DEBUG
		fprintf(stderr,"\tTask %d receiving %d particles from Task %d\n",ThisTask,NfromNextSlice,NextSlice);
		#endif
		FromNextSliceX = (float *) malloc(NfromNextSlice*sizeof(float));
		FromNextSliceY = (float *) malloc(NfromNextSlice*sizeof(float));
		FromNextSliceZ = (float *) malloc(NfromNextSlice*sizeof(float));
		FromNextSliceVX = (float *) malloc(NfromNextSlice*sizeof(float));
		FromNextSliceVY = (float *) malloc(NfromNextSlice*sizeof(float));
		FromNextSliceVZ = (float *) malloc(NfromNextSlice*sizeof(float));
		MPI_Recv( &(FromNextSliceX[0]), NfromNextSlice, MPI_FLOAT, NextSlice, 2111, MPI_COMM_WORLD, &status);
		MPI_Recv( &(FromNextSliceY[0]), NfromNextSlice, MPI_FLOAT, NextSlice, 2112, MPI_COMM_WORLD, &status);
		MPI_Recv( &(FromNextSliceZ[0]), NfromNextSlice, MPI_FLOAT, NextSlice, 2113, MPI_COMM_WORLD, &status);
		MPI_Recv( &(FromNextSliceVX[0]), NfromNextSlice, MPI_FLOAT, NextSlice, 2114, MPI_COMM_WORLD, &status);
		MPI_Recv( &(FromNextSliceVY[0]), NfromNextSlice, MPI_FLOAT, NextSlice, 2115, MPI_COMM_WORLD, &status);
		MPI_Recv( &(FromNextSliceVZ[0]), NfromNextSlice, MPI_FLOAT, NextSlice, 2116, MPI_COMM_WORLD, &status);
		#ifdef DEBUG
  		fprintf(stderr,"Task %d, received Vy[%d]=%f  from task %d\n",ThisTask,NfromNextSlice-1,FromNextSliceVY[NfromNextSlice-1],NextSlice);
		#endif


		#ifdef DEBUG
 		fprintf(stderr,"\tTask %d, sending %d particles to task %d\n",ThisTask,NPartPrevSlice,PrevSlice);
  		#endif
  		MPI_Send( &NPartPrevSlice, 1, MPI_INT, PrevSlice, 2234, MPI_COMM_WORLD); 
  		#ifdef DEBUG
  		fprintf(stderr,"\tTask %d, sending x[%d]=%f \n",ThisTask,NPartPrevSlice-1,PrevSliceX[NPartPrevSlice-1]);
  		#endif
  		MPI_Send( &PrevSliceX[0], NPartPrevSlice, MPI_FLOAT, PrevSlice, 2111, MPI_COMM_WORLD);
  		#ifdef DEBUG
  		fprintf(stderr,"\tTask %d, sending y[%d]=%f \n",ThisTask,NPartPrevSlice-1,PrevSliceY[NPartPrevSlice-1]);
  		#endif
  		MPI_Send( &PrevSliceY[0], NPartPrevSlice, MPI_FLOAT, PrevSlice, 2112, MPI_COMM_WORLD);
  		#ifdef DEBUG
  		fprintf(stderr,"\tTask %d, sending z[%d]=%f \n",ThisTask,NPartPrevSlice-1,PrevSliceZ[NPartPrevSlice-1]);
  		#endif
  		MPI_Send( &PrevSliceZ[0], NPartPrevSlice, MPI_FLOAT, PrevSlice, 2113, MPI_COMM_WORLD);
  		#ifdef DEBUG
  		fprintf(stderr,"\tTask %d, sending vx[%d]=%f \n",ThisTask,NPartPrevSlice-1,PrevSliceVX[NPartPrevSlice-1]);
  		#endif
  		MPI_Send( &PrevSliceVX[0], NPartPrevSlice, MPI_FLOAT, PrevSlice, 2114, MPI_COMM_WORLD);
  		#ifdef DEBUG
  		fprintf(stderr,"\tTask %d, sending vy[%d]=%f \n",ThisTask,NPartPrevSlice-1,PrevSliceVY[NPartPrevSlice-1]);
  		#endif
  		MPI_Send( &PrevSliceVY[0], NPartPrevSlice, MPI_FLOAT, PrevSlice, 2115, MPI_COMM_WORLD);
  		#ifdef DEBUG
  		fprintf(stderr,"\tTask %d, sending vz[%d]=%f \n",ThisTask,NPartPrevSlice-1,PrevSliceVZ[NPartPrevSlice-1]);
  		#endif
  		MPI_Send( &PrevSliceVZ[0], NPartPrevSlice, MPI_FLOAT, PrevSlice, 2116, MPI_COMM_WORLD);
	}

	#ifdef VERB
		#ifndef DEBUG
		if (ThisTask==0)
		#endif
			fprintf(stderr,"\tTask %d received %d parts from task %d, and %d parts from task %d\n",ThisTask,NfromPrevSlice,PrevSlice,NfromNextSlice,NextSlice);
	#endif
		
      	P = (struct part_data *) realloc(P,sizeof(struct part_data) * (NumPart+ NfromPrevSlice+NfromNextSlice));

	if (P==NULL){
		fprintf(stderr,"ERROR reallocating memory for Particles in task %d \t", ThisTask);
		return -1;
	}
	#ifdef DEBUG
	fprintf(stderr,"\tTask %d, copying %d from Prev\n",ThisTask,NfromPrevSlice);
	#endif

	for (ii=0;ii<NfromPrevSlice;ii++){
		ipart =ii+NumPart;
		P[ipart].Pos[0]=FromPrevSliceX[ii];
		P[ipart].Pos[1]=FromPrevSliceY[ii];
		P[ipart].Pos[2]=FromPrevSliceZ[ii];
		P[ipart].Vel[0]=FromPrevSliceVX[ii];
		P[ipart].Vel[1]=FromPrevSliceVY[ii];
		P[ipart].Vel[2]=FromPrevSliceVZ[ii];
		
		i = (long) (invL * P[ipart].Pos[0]*Nlin);
		ilocal = (long) i - Nx*ThisTask;
		j = (long) (invL * P[ipart].Pos[1]*Nlin);
		k = (long) (invL * P[ipart].Pos[2]*Nlin);
		lin_ijk = k+j*Nlin+ilocal*Nlin*Nlin;
		
		if ((ilocal<0) || (ilocal>Nx) || (j<0) || (j>=Nlin) || (k<0) || (k>=Nlin)){
			fprintf(stderr,"\tTask %d ERROR: Transferred Particle %ld at [%f,%f,%f]->[%d,%d,%d] seems to be out of the right box interval [0,%f) ii=%d\n",ThisTask,ipart,P[ipart].Pos[0],P[ipart].Pos[2],P[ipart].Pos[2],i,j,k,Box,ii);
			return -1;
		}
		(*NPartPerCell)[lin_ijk]++;
	}
	#ifdef DEBUG	
	fprintf(stderr,"Task %d, parts from prev copied\n",ThisTask);
	fprintf(stderr,"Task %d, copying %d from Next\n",ThisTask,NfromNextSlice);
	#endif

	for (ii=0;ii<NfromNextSlice;ii++){
		ipart =ii+NumPart+NfromPrevSlice;
		P[ipart].Pos[0]=FromNextSliceX[ii];
		P[ipart].Pos[1]=FromNextSliceY[ii];
		P[ipart].Pos[2]=FromNextSliceZ[ii];
		P[ipart].Vel[0]=FromNextSliceVX[ii];
		P[ipart].Vel[1]=FromNextSliceVY[ii];
		P[ipart].Vel[2]=FromNextSliceVZ[ii];
		
		i = (long) (invL * P[ipart].Pos[0]*Nlin);
		ilocal = (long) i - Nx*ThisTask;
		j = (long) (invL * P[ipart].Pos[1]*Nlin);
		k = (long) (invL * P[ipart].Pos[2]*Nlin);
		lin_ijk = k+j*Nlin+ilocal*Nlin*Nlin;
		
		if ((ilocal<0) || (ilocal>Nx) || (j<0) || (j>=Nlin) || (k<0) || (k>=Nlin)){
			fprintf(stderr,"\tTask %d ERROR: Transferred Particle %ld at [%f,%f,%f]->[%d,%d,%d] seems to be out of the right box interval [0,%f) ii=%d\n",ThisTask,ipart,P[ipart].Pos[0],P[ipart].Pos[2],P[ipart].Pos[2],i,j,k,Box,ii);
		}
		(*NPartPerCell)[lin_ijk]++;
	}
	#ifdef DEBUG
	fprintf(stderr,"\tParticles Transferred Counted \n\n");
	#endif
	for (i=0;i<Nx;i++){
	for (j=0;j<Nlin;j++){
	for (k=0;k<Nlin;k++){
		lin_ijk = k + j*Nlin + i*Nlin*Nlin;
		(*ListOfPart)[lin_ijk] = (long *) calloc((*NPartPerCell)[lin_ijk],sizeof(long));
	}
	}
	}

	#ifdef DEBUG
	fprintf(stderr,"\t Grid allocated \n\n");
	#endif


	ipart=0;

	for (ipart=0;ipart<NumPart+NfromPrevSlice+NfromNextSlice;ipart++) {
		i = (long) (invL * P[ipart].Pos[0]*Nlin);
		ilocal = (long) i - Nx*ThisTask;
		j = (long) (invL * P[ipart].Pos[1]*Nlin);
		k = (long) (invL * P[ipart].Pos[2]*Nlin);
		if (ilocal>=0 && ilocal<Nx) {
			lin_ijk = k + j*Nlin + ilocal*Nlin*Nlin;

			(*ListOfPart)[lin_ijk][count[lin_ijk]] = ipart;
			count[lin_ijk]++;
		}
	}

	return 0;
}






//#else //ONLY_2LPT
void save_local_data(void);

void write_particle_data(void)
{
  int nprocgroup, groupTask, masterTask;

  if(ThisTask == 0)
    printf("\nwriting initial conditions... \n");


  if((NTask < NumFilesWrittenInParallel))
    {
      printf
	("Fatal error.\nNumber of processors must be a smaller or equal than `NumFilesWrittenInParallel'.\n");
      FatalError(24131);
    }


  nprocgroup = NTask / NumFilesWrittenInParallel;

  if((NTask % NumFilesWrittenInParallel))
    nprocgroup++;

  masterTask = (ThisTask / nprocgroup) * nprocgroup;


  for(groupTask = 0; groupTask < nprocgroup; groupTask++)
    {
      if(ThisTask == (masterTask + groupTask))	/* ok, it's this processor's turn */
	save_local_data();

      /* wait inside the group */
      MPI_Barrier(MPI_COMM_WORLD);
    }

  if(ThisTask == 0)
    printf("done with writing initial conditions.\n");
}

void save_local_data(void)
{
#define BUFFER 10
  size_t bytes;
  float *block;


  int blockmaxlen,  maxlongidlen;
  int4byte dummy;
  FILE *fd;
  char buf[300];
  int i, k, pc;



  if(NumPart == 0)
    return;

  if(NTaskWithN > 1)
    sprintf(buf, "%s/%s.%d", OutputDir, FileBase, ThisTask);
  else
    sprintf(buf, "%s/%s", OutputDir, FileBase);

  if(!(fd = fopen(buf, "w")))
    {
      printf("Error. Can't write in file '%s'\n", buf);
      FatalError(10);
    }

  for(i = 0; i < 6; i++)
    {
      header.npart[i] = 0;
      header.npartTotal[i] = 0;
      header.mass[i] = 0;
    }


#ifdef MULTICOMPONENTGLASSFILE
  qsort(P, NumPart, sizeof(struct part_data), compare_type);  // sort particles by type, because that's how they should be stored in a gadget binary file 

  for(i = 0; i < 3; i++)
    header.npartTotal[i] = header1.npartTotal[i + 1] * GlassTileFac * GlassTileFac * GlassTileFac;

  for(i = 0; i < NumPart; i++)
    header.npart[P[i].Type]++;

  if(header.npartTotal[0])
    header.mass[0] =
      (OmegaBaryon) * 3 * Hubble * Hubble / (8 * PI * G) * pow(Box, 3) / (header.npartTotal[0]);

  if(header.npartTotal[1])
    header.mass[1] =
      (Omega - OmegaBaryon - OmegaDM_2ndSpecies) * 3 * Hubble * Hubble / (8 * PI * G) * pow(Box,
											    3) /
      (header.npartTotal[1]);

  if(header.npartTotal[2])
    header.mass[2] =
      (OmegaDM_2ndSpecies) * 3 * Hubble * Hubble / (8 * PI * G) * pow(Box, 3) / (header.npartTotal[2]);


#else

  header.npart[1] = NumPart;
  header.npartTotal[1] = TotNumPart;
  header.npartTotal[2] = (TotNumPart >> 32);
  header.mass[1] = (Omega) * 3 * Hubble * Hubble / (8 * PI * G) * pow(Box, 3) / TotNumPart;

#ifdef  PRODUCEGAS
  header.npart[0] = NumPart;
  header.npartTotal[0] = TotNumPart;
  header.mass[0] = (OmegaBaryon) * 3 * Hubble * Hubble / (8 * PI * G) * pow(Box, 3) / TotNumPart;
  header.mass[1] = (Omega - OmegaBaryon) * 3 * Hubble * Hubble / (8 * PI * G) * pow(Box, 3) / TotNumPart;
#endif
#endif


  header.time = InitTime;
  header.redshift = 1.0 / InitTime - 1;

  header.flag_sfr = 0;
  header.flag_feedback = 0;
  header.flag_cooling = 0;
  header.flag_stellarage = 0;
  header.flag_metals = 0;

  header.num_files = NTaskWithN;

  header.BoxSize = Box;
  header.Omega0 = Omega;
  header.OmegaLambda = OmegaLambda;
  header.HubbleParam = HubbleParam;

  header.flag_stellarage = 0;
  header.flag_metals = 0;
  header.hashtabsize = 0;

  dummy = sizeof(header);
  my_fwrite(&dummy, sizeof(dummy), 1, fd);
  my_fwrite(&header, sizeof(header), 1, fd);
  my_fwrite(&dummy, sizeof(dummy), 1, fd);



#ifdef  PRODUCEGAS
  double meanspacing;
  double  shift_gas, shift_dm;
  meanspacing = Box / pow(TotNumPart, 1.0 / 3);
  shift_gas = -0.5 * (Omega - OmegaBaryon) / (Omega) * meanspacing;
  shift_dm = +0.5 * OmegaBaryon / (Omega) * meanspacing;
#endif



  if(!(block = malloc(bytes = BUFFER * 1024 * 1024)))
    {
      printf("failed to allocate memory for `block' (%g bytes).\n", (double)bytes);
      FatalError(24);
    }

  blockmaxlen = bytes / (3 * sizeof(float));

#ifdef NO64BITID
  int *blockid;
  blockid = (int *) block;
#else
  long long *blocklongid;
  blocklongid = (long long *) block;
#endif
  maxlongidlen = bytes / (sizeof(long long));

  // write coordinates 
  dummy = sizeof(float) * 3 * NumPart;
#ifdef  PRODUCEGAS
  dummy *= 2;
#endif
  my_fwrite(&dummy, sizeof(dummy), 1, fd);
  for(i = 0, pc = 0; i < NumPart; i++)
    {
      for(k = 0; k < 3; k++)
	{
	  block[3 * pc + k] = P[i].Pos[k];
#ifdef  PRODUCEGAS
	  block[3 * pc + k] = periodic_wrap(P[i].Pos[k] + shift_gas);
#endif
	}

      pc++;

      if(pc == blockmaxlen)
	{
	  my_fwrite(block, sizeof(float), 3 * pc, fd);
	  pc = 0;
	}
    }
  if(pc > 0)
    my_fwrite(block, sizeof(float), 3 * pc, fd);
#ifdef  PRODUCEGAS
  for(i = 0, pc = 0; i < NumPart; i++)
    {
      for(k = 0; k < 3; k++)
	{
	  block[3 * pc + k] = periodic_wrap(P[i].Pos[k] + shift_dm);
	}

      pc++;

      if(pc == blockmaxlen)
	{
	  my_fwrite(block, sizeof(float), 3 * pc, fd);
	  pc = 0;
	}
    }
  if(pc > 0)
    my_fwrite(block, sizeof(float), 3 * pc, fd);
#endif
  my_fwrite(&dummy, sizeof(dummy), 1, fd);



  // write velocities 
  dummy = sizeof(float) * 3 * NumPart;
#ifdef  PRODUCEGAS
  dummy *= 2;
#endif
  my_fwrite(&dummy, sizeof(dummy), 1, fd);
  for(i = 0, pc = 0; i < NumPart; i++)
    {
      for(k = 0; k < 3; k++)
	block[3 * pc + k] = P[i].Vel[k];

#ifdef MULTICOMPONENTGLASSFILE
      if(WDM_On == 1 && WDM_Vtherm_On == 1 && P[i].Type == 1)
	add_WDM_thermal_speeds(&block[3 * pc]);
#else
#ifndef PRODUCEGAS
      if(WDM_On == 1 && WDM_Vtherm_On == 1)
	add_WDM_thermal_speeds(&block[3 * pc]);
#endif
#endif

      pc++;

      if(pc == blockmaxlen)
	{
	  my_fwrite(block, sizeof(float), 3 * pc, fd);
	  pc = 0;
	}
    }
  if(pc > 0)
    my_fwrite(block, sizeof(float), 3 * pc, fd);
#ifdef PRODUCEGAS
  for(i = 0, pc = 0; i < NumPart; i++)
    {
      for(k = 0; k < 3; k++)
	block[3 * pc + k] = P[i].Vel[k];

      if(WDM_On == 1 && WDM_Vtherm_On == 1)
	add_WDM_thermal_speeds(&block[3 * pc]);

      pc++;

      if(pc == blockmaxlen)
	{
	  my_fwrite(block, sizeof(float), 3 * pc, fd);
	  pc = 0;
	}
    }
  if(pc > 0)
    my_fwrite(block, sizeof(float), 3 * pc, fd);
#endif
  my_fwrite(&dummy, sizeof(dummy), 1, fd);


  // write particle ID 
#ifdef NO64BITID
  dummy = sizeof(int) * NumPart;
#else
  dummy = sizeof(long long) * NumPart;
#endif
#ifdef  PRODUCEGAS
  dummy *= 2;
#endif
  my_fwrite(&dummy, sizeof(dummy), 1, fd);
  for(i = 0, pc = 0; i < NumPart; i++)
    {
#ifdef NO64BITID
      blockid[pc] = P[i].ID;
#else
      blocklongid[pc] = P[i].ID;
#endif

      pc++;

      if(pc == maxlongidlen)
	{
#ifdef NO64BITID
	  my_fwrite(blockid, sizeof(int), pc, fd);
#else
	  my_fwrite(blocklongid, sizeof(long long), pc, fd);
#endif
	  pc = 0;
	}
    }
  if(pc > 0)
    {
#ifdef NO64BITID
      my_fwrite(blockid, sizeof(int), pc, fd);
#else
      my_fwrite(blocklongid, sizeof(long long), pc, fd);
#endif
    }

#ifdef PRODUCEGAS
  for(i = 0, pc = 0; i < NumPart; i++)
    {
#ifdef NO64BITID
      blockid[pc] = P[i].ID + TotNumPart;
#else
      blocklongid[pc] = P[i].ID + TotNumPart;
#endif

      pc++;

      if(pc == maxlongidlen)
	{
#ifdef NO64BITID
	  my_fwrite(blockid, sizeof(int), pc, fd);
#else
	  my_fwrite(blocklongid, sizeof(long long), pc, fd);
#endif
	  pc = 0;
	}
    }
  if(pc > 0)
    {
#ifdef NO64BITID
      my_fwrite(blockid, sizeof(int), pc, fd);
#else
      my_fwrite(blocklongid, sizeof(long long), pc, fd);
#endif
    }
#endif

  my_fwrite(&dummy, sizeof(dummy), 1, fd);





  // write zero temperatures if needed 
#ifdef  PRODUCEGAS
  dummy = sizeof(float) * NumPart;
  my_fwrite(&dummy, sizeof(dummy), 1, fd);
  for(i = 0, pc = 0; i < NumPart; i++)
    {
      block[pc] = 0;

      pc++;

      if(pc == blockmaxlen)
	{
	  my_fwrite(block, sizeof(float), pc, fd);
	  pc = 0;
	}
    }
  if(pc > 0)
    my_fwrite(block, sizeof(float), pc, fd);
  my_fwrite(&dummy, sizeof(dummy), 1, fd);
#endif


  // write zero temperatures if needed 
#ifdef  MULTICOMPONENTGLASSFILE
  if(header.npart[0])
    {
      dummy = sizeof(float) * header.npart[0];
      my_fwrite(&dummy, sizeof(dummy), 1, fd);

      for(i = 0, pc = 0; i < header.npart[0]; i++)
	{
	  block[pc] = 0;

	  pc++;

	  if(pc == blockmaxlen)
	    {
	      my_fwrite(block, sizeof(float), pc, fd);
	      pc = 0;
	    }
	}
      if(pc > 0)
	my_fwrite(block, sizeof(float), pc, fd);
      my_fwrite(&dummy, sizeof(dummy), 1, fd);
    }
#endif



  free(block);

  fclose(fd);
}

//#endif //ONLY_2LPT



/* This catches I/O errors occuring for my_fwrite(). In this case we better stop.
 */
size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE * stream)
{
  size_t nwritten;

  if((nwritten = fwrite(ptr, size, nmemb, stream)) != nmemb)
    {
      printf("I/O error (fwrite) on task=%d has occured.\n", ThisTask);
      fflush(stdout);
      FatalError(777);
    }
  return nwritten;
}


/* This catches I/O errors occuring for fread(). In this case we better stop.
 */
size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE * stream)
{
  size_t nread;

  if((nread = fread(ptr, size, nmemb, stream)) != nmemb)
    {
      printf("I/O error (fread) on task=%d has occured.\n", ThisTask);
      fflush(stdout);
      FatalError(778);
    }
  return nread;
}


#ifdef MULTICOMPONENTGLASSFILE
int compare_type(const void *a, const void *b)
{
  if(((struct part_data *) a)->Type < (((struct part_data *) b)->Type))
    return -1;

  if(((struct part_data *) a)->Type > (((struct part_data *) b)->Type))
    return +1;

  return 0;
}
#endif
