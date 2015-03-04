#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>

#include "allvars.h"
#include "proto.h"





int distribute_part(int Nlin,int nthreads,float **Partx, float **Party, float **Partz, float **Partvx, float **Partvy, float **Partvz, long *NTotPart, float *Lbox, float *om_m, long ***ListOfPart,long **NPartPerCell)
{

  *om_m = Omega;
  *Lbox = Box;
  *NTotPart = (GlassTileFac * GlassTileFac * GlassTileFac);
  int nprocgroup, groupTask, masterTask;
  int i,j,k;
  int *PartPerFile;
  MPI_Status status;

  float **file_x, **file_y, **file_z, **file_vx,**file_vy,**file_vz;
  float *thisfile_x, *thisfile_y, *thisfile_z, *thisfile_vx,*thisfile_vy,*thisfile_vz;
  int Nthisfile;
  Nthisfile= (int) NumPart;

  
  if (ThisTask == 0){

      PartPerFile = (int*) calloc(NTask,sizeof(int));

      if(!((file_x)=(float **) malloc(NTask*sizeof(float *))))
      {
        fprintf(stderr,"\nfailed to allocate memory for GADGET data\n");
        exit(1);
      }
      if(!((file_y)=(float **) malloc(NTask*sizeof(float *))))
      {
        fprintf(stderr,"\nfailed to allocate memory for GADGET data\n");
        exit(1);
      }
      if(!((file_z)=(float **) malloc(NTask*sizeof(float *))))
      {
        fprintf(stderr,"\nfailed to allocate memory for GADGET data\n");
        exit(1);
      }
      if(!((file_vx)=(float **) malloc(NTask*sizeof(float *))))
      {
        fprintf(stderr,"\nfailed to allocate memory for GADGET data\n");
        exit(1);
      }
      if(!((file_vy)=(float **) malloc(NTask*sizeof(float *))))
      {
        fprintf(stderr,"\nfailed to allocate memory for GADGET data\n");
        exit(1);
      }
      if(!((file_vz)=(float **) malloc(NTask*sizeof(float *))))
      {
        fprintf(stderr,"\nfailed to allocate memory for GADGET data\n");
        exit(1);
      }


      if(!((file_x[0])=(float *) malloc(Nthisfile*sizeof(float))))
      {
        fprintf(stderr,"\nfailed to allocate memory for GADGET data\n");
        exit(1);
      }
      if(!((file_y[0])=(float *) malloc(Nthisfile*sizeof(float))))
      {
        fprintf(stderr,"\nfailed to allocate memory for GADGET data\n");
        exit(1);
      }
      if(!((file_z[0])=(float *) malloc(Nthisfile*sizeof(float))))
      {
        fprintf(stderr,"\nfailed to allocate memory for GADGET data\n");
        exit(1);
      }
      if(!((file_vx[0])=(float *) malloc(Nthisfile*sizeof(float))))
      {
        fprintf(stderr,"\nfailed to allocate memory for GADGET data\n");
        exit(1);
      }
      if(!((file_vy[0])=(float *) malloc(Nthisfile*sizeof(float))))
      {
        fprintf(stderr,"\nfailed to allocate memory for GADGET data\n");
        exit(1);
      }
      if(!((file_vz[0])=(float *) malloc(Nthisfile*sizeof(float))))
      {
        fprintf(stderr,"\nfailed to allocate memory for GADGET data\n");
        exit(1);
      } 
/*
     for (i=0;i<NumPart;i++){	
	//should be overcome
        file_x[0][i]=(float)partX[i];
        file_y[0][i]=(float)partY[i];
        file_z[0][i]=(float)partZ[i];
        file_vx[0][i]=(float)partVX[i];
        file_vy[0][i]=(float)partVY[i];
        file_vz[0][i]=(float)partVZ[i];
     }
*/
     memcpy(file_x[0],partX,Nthisfile*sizeof(float));
     memcpy(file_y[0],partY,Nthisfile*sizeof(float));
     memcpy(file_z[0],partZ,Nthisfile*sizeof(float));
     memcpy(file_vx[0],partVX,Nthisfile*sizeof(float));
     memcpy(file_vy[0],partVY,Nthisfile*sizeof(float));
     memcpy(file_vz[0],partVZ,Nthisfile*sizeof(float));

  }

//MPI_Barrier(MPI_COMM_WORLD); 
//if (ThisTask!=0) {
else {

  #ifdef DEBUG
  fprintf(stderr,"\tTask %d, sending %d particles\n",ThisTask,Nthisfile);
  #endif
  MPI_Send( &Nthisfile, 1, MPI_INT, 0, 1234, MPI_COMM_WORLD); 
  #ifdef DEBUG
  fprintf(stderr,"\tTask %d, sending x[%d]=%f \n",ThisTask,Nthisfile-1,partX[Nthisfile-1]);
  #endif
  MPI_Send( &partX[0], Nthisfile, MPI_FLOAT, 0, 1111, MPI_COMM_WORLD);
  #ifdef DEBUG
  fprintf(stderr,"\tTask %d, sending y[%d]=%f \n",ThisTask,Nthisfile-1,partY[Nthisfile-1]);
  #endif
  MPI_Send( &partY[0], Nthisfile, MPI_FLOAT, 0, 1112, MPI_COMM_WORLD);
  #ifdef DEBUG
  fprintf(stderr,"\tTask %d, sending z[%d]=%f \n",ThisTask,Nthisfile-1,partZ[Nthisfile-1]);
  #endif
  MPI_Send( &partZ[0], Nthisfile, MPI_FLOAT, 0, 1113, MPI_COMM_WORLD);
  #ifdef DEBUG
  fprintf(stderr,"\tTask %d, sending vz[%d]=%f \n",ThisTask,Nthisfile-1,partVZ[Nthisfile-1]);
  #endif
  MPI_Send( &partVZ[0], Nthisfile, MPI_FLOAT, 0, 1116, MPI_COMM_WORLD);
  #ifdef DEBUG
  fprintf(stderr,"\tTask %d, sending vx[%d]=%f \n",ThisTask,Nthisfile-1,partVX[Nthisfile-1]);
  #endif
  MPI_Send( &partVX[0], Nthisfile, MPI_FLOAT, 0, 1114, MPI_COMM_WORLD);
  #ifdef DEBUG
  fprintf(stderr,"\tTask %d, sending vy[%d]=%f \n",ThisTask,Nthisfile-1,partVY[Nthisfile-1]);
  #endif
  MPI_Send( &partVY[0], Nthisfile, MPI_FLOAT, 0, 1115, MPI_COMM_WORLD);

  

/*  free(thisfile_x);
  free(thisfile_y);
  free(thisfile_z);
  free(thisfile_vx);
  free(thisfile_vy);
  free(thisfile_vz);
*/
  //sleep(1000);
  #ifdef DEBUG
  fprintf(stderr,"\tDone with Task %d\n",ThisTask);
  #endif
}
//MPI_Barrier(MPI_COMM_WORLD);
  if (ThisTask == 0){
     int totN= NumPart;
     PartPerFile[groupTask]=NumPart;
     for(groupTask = 1; groupTask < NTask; groupTask++)
     {	
	#ifdef DEBUG
  	//sleep(2);
  	fprintf(stderr,"\tTask 0, starting contactwith task %d\n",groupTask);
  	//sleep(2);
	#endif
	MPI_Recv(&Nthisfile, 1, MPI_INT, groupTask, 1234, MPI_COMM_WORLD, &status);
	#ifdef DEBUG
  	//sleep(2);
  	fprintf(stderr,"\tTask 0, receiving %d particles from task %d\n",Nthisfile,groupTask);
  	//sleep(2);
	#endif
  	PartPerFile[groupTask]=Nthisfile;
     	totN+=Nthisfile;

      if(!((file_x[groupTask])=(float *) malloc(Nthisfile* sizeof(float))))
      {
        fprintf(stderr,"\nfailed to allocate memory for GADGET data\n");
        exit(1);
      }
      if(!((file_y[groupTask])=(float *) malloc(Nthisfile* sizeof(float))))
      {
        fprintf(stderr,"\nfailed to allocate memory for GADGET data\n");
        exit(1);
      }
      if(!((file_z[groupTask])=(float *) malloc(Nthisfile* sizeof(float))))
      {
        fprintf(stderr,"\nfailed to allocate memory for GADGET data\n");
        exit(1);
      }
      if(!((file_vx[groupTask])=(float *) malloc(Nthisfile* sizeof(float))))
      {
        fprintf(stderr,"\nfailed to allocate memory for GADGET data\n");
        exit(1);
      }
      if(!((file_vy[groupTask])=(float *) malloc(Nthisfile* sizeof(float))))
      {
        fprintf(stderr,"\nfailed to allocate memory for GADGET data\n");
        exit(1);
      }
      if(!((file_vz[groupTask])=(float *) malloc(Nthisfile* sizeof(float))))
      {
        fprintf(stderr,"\nfailed to allocate memory for GADGET data\n");
        exit(1);
      }

	/*
        thisfile_x=(float *) calloc(Nthisfile, sizeof(float));
        thisfile_y=(float *) calloc(Nthisfile, sizeof(float));
        thisfile_z=(float *) calloc(Nthisfile, sizeof(float));
        thisfile_vx=(float *) calloc(Nthisfile, sizeof(float));
        thisfile_vy=(float *) calloc(Nthisfile, sizeof(float));
        thisfile_vz=(float *) calloc(Nthisfile, sizeof(float));

	MPI_Recv( &thisfile_x[0], Nthisfile, MPI_FLOAT, groupTask, 1111, MPI_COMM_WORLD, &status);
	MPI_Recv( &thisfile_y[0], Nthisfile, MPI_FLOAT, groupTask, 1112, MPI_COMM_WORLD, &status);
	MPI_Recv( &thisfile_z[0], Nthisfile, MPI_FLOAT, groupTask, 1113, MPI_COMM_WORLD, &status);
	MPI_Recv( &thisfile_vx[0], Nthisfile, MPI_FLOAT, groupTask, 1114, MPI_COMM_WORLD, &status);
	MPI_Recv( &thisfile_vy[0], Nthisfile, MPI_FLOAT, groupTask, 1115, MPI_COMM_WORLD, &status);
	MPI_Recv( &thisfile_vz[0], Nthisfile, MPI_FLOAT, groupTask, 1116, MPI_COMM_WORLD, &status);
	*/

	MPI_Recv( &file_x[groupTask][0], Nthisfile, MPI_FLOAT, groupTask, 1111, MPI_COMM_WORLD, &status);
	#ifdef DEBUG
  	//sleep(2);
  	fprintf(stderr,"\tTask 0, received x[%d]=%f  from task %d\n",Nthisfile-1,file_x[groupTask][Nthisfile-1],groupTask);
  	//sleep(2);
	#endif
	MPI_Recv( &(file_y[groupTask][0]), Nthisfile, MPI_FLOAT, groupTask, 1112, MPI_COMM_WORLD, &status);
	#ifdef DEBUG
  	//sleep(2);
  	fprintf(stderr,"\tTask 0, received y[%d]=%f  from task %d\n",Nthisfile-1,file_y[groupTask][Nthisfile-1],groupTask);
  	//sleep(2);
	#endif
	MPI_Recv( &(file_z[groupTask][0]), Nthisfile, MPI_FLOAT, groupTask, 1113, MPI_COMM_WORLD, &status);
	#ifdef DEBUG
  	//sleep(2);
  	fprintf(stderr,"\tTask 0, received z[%d]=%f  from task %d\n",Nthisfile-1,file_z[groupTask][Nthisfile-1],groupTask);
  	//sleep(2);
	#endif
	MPI_Recv( &(file_vz[groupTask][0]), Nthisfile, MPI_FLOAT, groupTask, 1116, MPI_COMM_WORLD, &status);
	#ifdef DEBUG
  	//sleep(2);
  	fprintf(stderr,"\tTask 0, received Vz[%d]=%f  from task %d\n",Nthisfile-1,file_vz[groupTask][Nthisfile-1],groupTask);
  	//sleep(2);
	#endif
	MPI_Recv( &(file_vx[groupTask][0]), Nthisfile, MPI_FLOAT, groupTask, 1114, MPI_COMM_WORLD, &status);
	#ifdef DEBUG
  	//sleep(2);
  	fprintf(stderr,"\tTask 0, received Vx[%d]=%f  from task %d\n",Nthisfile-1,file_vx[groupTask][Nthisfile-1],groupTask);
  	//sleep(2);
	#endif
	MPI_Recv( &(file_vy[groupTask][0]), Nthisfile, MPI_FLOAT, groupTask, 1115, MPI_COMM_WORLD, &status);
	#ifdef DEBUG
  	//sleep(2);
  	fprintf(stderr,"\tTask 0, received Vy[%d]=%f  from task %d\n",Nthisfile-1,file_vy[groupTask][Nthisfile-1],groupTask);
  	//sleep(2);
	#endif

      /*
      for (i=0;i<NumPart;i++){
        file_x[groupTask][i]=thisfile_x[i];
        file_y[groupTask][i]=thisfile_y[i];
        file_z[groupTask][i]=thisfile_z[i];
        file_vx[groupTask][i]=thisfile_vx[i];
        file_vy[groupTask][i]=thisfile_vy[i];
        file_vz[groupTask][i]=thisfile_vz[i];
      }
	
      free(thisfile_x);
      free(thisfile_y);
      free(thisfile_z);
      free(thisfile_vx);
      free(thisfile_vy);
      free(thisfile_vz);
      */

     }

        //Distribute Particles
        (*NPartPerCell) = (long *) calloc(Nlin*Nlin*Nlin,sizeof(long));
        if( *NPartPerCell == NULL) {
                fprintf(stderr,"\tplace_halos(): could not allocate %d array for NPartPerCell[]\nABORTING",Nlin*Nlin*Nlin);
                exit(-1);
        }
        #ifdef DEBUG
        fprintf(stderr,"\t...NPartPerCell allocaled, continue...\n");
        #endif
        (*ListOfPart) = (long **) calloc(Nlin*Nlin*Nlin,sizeof(long *));
        if( *ListOfPart == NULL) {
                fprintf(stderr,"\tplace_halos(): could not allocate %d array for NPartPerCell[]\nABORTING",Nlin*Nlin*Nlin);
                exit(-1);
        }
        #ifdef DEBUG
        fprintf(stderr,"\t...ListOFPART allocaled, continue...\n");
        #endif
	long *count;
        count = (long *) calloc(Nlin*Nlin*Nlin,sizeof(long));
        if( count == NULL) {
                fprintf(stderr,"\tplace_halos(): could not allocate %d array for NPartPerCell[]\nABORTING",Nlin*Nlin*Nlin);
                exit(-1);
        }
        #ifdef DEBUG
        fprintf(stderr,"\t Counters allocated \n\n");
        #endif
        long ipart=0,ii,i_gadget_file,lin_ijk;
        double invL = 1.0/Box;

	#ifdef DEBUG
	fprintf(stderr,"\tL=%f, invL=%f, Nlin=%d\n",Box,invL,Nlin);
	#endif
        //Counting Particles in cells
        //t3=time(NULL);
        for(i_gadget_file=0; i_gadget_file<NTask; i_gadget_file++)
        for (ii=0;ii<PartPerFile[i_gadget_file];ii++) {
                if (file_x[i_gadget_file][ii]==Box)
                        file_x[i_gadget_file][ii]=0.;
                if (file_y[i_gadget_file][ii]==Box)
                        file_y[i_gadget_file][ii]=0.;
                if (file_z[i_gadget_file][ii]==Box)
                        file_z[i_gadget_file][ii]=0.;
                i = (long) (invL * file_x[i_gadget_file][ii]*Nlin);
                j = (long) (invL * file_y[i_gadget_file][ii]*Nlin);
                k = (long) (invL * file_z[i_gadget_file][ii]*Nlin);
                if ((i<0) || (i>=Nlin) || (j<0) || (j>=Nlin) || (k<0) || (k>=Nlin)){
                        fprintf(stderr,"\tERROR: Particle %ld at [%f,%f,%f]->[%d,%d,%d] seems to be out of the right box interval [0.,%f)\n",ipart,file_x[i_gadget_file][ii],file_y[i_gadget_file][ii],file_z[i_gadget_file][ii],i,j,k,Box);
			exit(0);
                }
                lin_ijk = k+j*Nlin+i*Nlin*Nlin;
                (*NPartPerCell)[lin_ijk]++;
                ipart++;
        }
        #ifdef DEBUG
         fprintf(stderr,"\tParticles Counted \n\n");
        #endif
        for (i=0;i<Nlin;i++){
        for (j=0;j<Nlin;j++){
        for (k=0;k<Nlin;k++){
                lin_ijk = k+j*Nlin+i*Nlin*Nlin;
                //fprintf(stderr,"lin_ijk=%d\n",lin_ijk);
                (*ListOfPart)[lin_ijk] = (long *) calloc((*NPartPerCell)[lin_ijk],sizeof(long));
        }
        }
        }
        #ifdef DEBUG
        fprintf(stderr,"\t Grid allocated \n\n");
        #endif
	(*NTotPart) = ipart;
        (*Partx) = (float *) calloc(ipart,sizeof(float));
        (*Party) = (float *) calloc(ipart,sizeof(float));
        (*Partz) = (float *) calloc(ipart,sizeof(float));

        (*Partvx) = (float *) calloc(ipart,sizeof(float));
        (*Partvy) = (float *) calloc(ipart,sizeof(float));
        (*Partvz) = (float *) calloc(ipart,sizeof(float));
        #ifdef DEBUG
        fprintf(stderr,"\t vector allocated \n\n");
         #endif
        ipart=0;
        //t5=time(NULL);
        //Distributing Particles
        //#pragma omp parallel for private(i_gadget_file,ipart,i,k,j,lin_ijk,ii) shared(PartPerFile,invL,Nlin,ListOfPart,count,out_x,out_z,out_y,out_vx,out_vz,out_vy,file_x,file_y,file_z,file_vx,file_vy,file_vz,stderr,gadget,no_gadget_files) default(none)
        for (i_gadget_file=0; i_gadget_file<NTask; i_gadget_file++)
        for (ii=0;ii<PartPerFile[i_gadget_file];ii++) {
                i = (long) (invL * file_x[i_gadget_file][ii]*Nlin);
                j = (long) (invL * file_y[i_gadget_file][ii]*Nlin);
                k = (long) (invL * file_z[i_gadget_file][ii]*Nlin);
                if (i<0 || i>=Nlin || j<0 || j>=Nlin || k<0 || k>=Nlin){
                        fprintf(stderr,"\tERROR: Particle %ld at [%f,%f,%f] seems to be out of the right box interval [0.,%f)",ipart,file_x[i_gadget_file][ii],file_y[i_gadget_file][ii],file_z[i_gadget_file][ii],Box);
                }
                (*Partx)[ipart] = file_x[i_gadget_file][ii];
                (*Party)[ipart] = file_y[i_gadget_file][ii];
                (*Partz)[ipart] = file_z[i_gadget_file][ii];

                (*Partvx)[ipart] = file_vx[i_gadget_file][ii];
                (*Partvy)[ipart] = file_vy[i_gadget_file][ii];
                (*Partvz)[ipart] = file_vz[i_gadget_file][ii];

                lin_ijk = k+j*Nlin+i*Nlin*Nlin;
                (*ListOfPart)[lin_ijk][count[lin_ijk]] = ipart;
                count[lin_ijk]++;
                ipart++;
        }

   }

MPI_Barrier(MPI_COMM_WORLD);

  return 0;
}

/*
void save_local_data(void)
{
#define BUFFER 10
  size_t bytes;
  float *block;
  int *blockid;
  long long *blocklongid;
  int blockmaxlen, maxidlen, maxlongidlen;
  int4byte dummy;
  FILE *fd;
  char buf[300];
  int i, k, pc;
  double meanspacing, shift_gas, shift_dm;


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


  meanspacing = Box / pow(TotNumPart, 1.0 / 3);
  shift_gas = -0.5 * (Omega - OmegaBaryon) / (Omega) * meanspacing;
  shift_dm = +0.5 * OmegaBaryon / (Omega) * meanspacing;


  if(!(block = malloc(bytes = BUFFER * 1024 * 1024)))
    {
      printf("failed to allocate memory for `block' (%g bytes).\n", (double)bytes);
      FatalError(24);
    }

  blockmaxlen = bytes / (3 * sizeof(float));

  blockid = (int *) block;
  blocklongid = (long long *) block;
  maxidlen = bytes / (sizeof(int));
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
*/

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
