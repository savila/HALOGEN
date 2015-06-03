
/*********************************************************************************** 
This code is a part of 2LPTic (http://arxiv.org/abs/astro-ph/9711187)
adapted to HALOGEN (http://arxiv.org/abs/1412.5228)
The original code 2LPTic was developped by 
Roman Scoccimarro and Sebastian Pueblas (http://cosmo.nyu.edu/roman/2LPT/)
************************************************************************************ */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "allvars.h"
#include "proto.h"


void read_parameterfile(char *fname)
{
#define FLOAT 1
#define STRING 2
#define INT 3
#define MAXTAGS 300

  FILE *fd;
  char buf[200], buf1[200], buf2[200], buf3[200];
  int i, j, nt;
  int id[MAXTAGS];
  void *addr[MAXTAGS];
  char tag[MAXTAGS][50];
  int errorFlag = 0;

  /* read parameter file on all processes for simplicty */

  nt = 0;

  strcpy(tag[nt], "Omega");
  addr[nt] = &Omega;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "OmegaLambda");
  addr[nt] = &OmegaLambda;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "OmegaBaryon");
  addr[nt] = &OmegaBaryon;
  id[nt++] = FLOAT;

#ifdef ONLY_2LPT
  strcpy(tag[nt], "OmegaDM_2ndSpecies");
  addr[nt] = &OmegaDM_2ndSpecies;
  id[nt++] = FLOAT;
#else
  OmegaDM_2ndSpecies=0.0;
#endif

  strcpy(tag[nt], "HubbleParam");
  addr[nt] = &HubbleParam;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "ShapeGamma");
  addr[nt] = &ShapeGamma;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "Sigma8");
  addr[nt] = &Sigma8;
  id[nt++] = FLOAT;

/*
  strcpy(tag[nt], "PrimordialIndex");
  addr[nt] = &PrimordialIndex;
  id[nt++] = FLOAT;
*/
//Assume Pk is already tilted
  PrimordialIndex=1;


  strcpy(tag[nt], "Box");
  addr[nt] = &Box;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "Redshift");
  addr[nt] = &Redshift;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "Nmesh");
  addr[nt] = &Nmesh;
  id[nt++] = INT;

  strcpy(tag[nt], "Nsample");
  addr[nt] = &Nsample;
  id[nt++] = INT;

  strcpy(tag[nt], "GlassFile");
  addr[nt] = GlassFile;
  id[nt++] = STRING;

  strcpy(tag[nt], "FileWithInputSpectrum");
  addr[nt] = FileWithInputSpectrum;
  id[nt++] = STRING;

  strcpy(tag[nt], "GlassTileFac");
  addr[nt] = &GlassTileFac;
  id[nt++] = INT;

  strcpy(tag[nt], "Seed");
  addr[nt] = &Seed;
  id[nt++] = INT;

/*
  strcpy(tag[nt], "SphereMode");
  addr[nt] = &SphereMode;
  id[nt++] = INT;
*/

  SphereMode = 0;


#ifdef ONLY_2LPT
  strcpy(tag[nt], "NumFilesWrittenInParallel");
  addr[nt] = &NumFilesWrittenInParallel;
  id[nt++] = INT;

  strcpy(tag[nt], "OutputDir");
  addr[nt] = OutputDir;
  id[nt++] = STRING;

  strcpy(tag[nt], "FileBase");
  addr[nt] = FileBase;
  id[nt++] = STRING;
#endif

  strcpy(tag[nt], "WhichSpectrum");
  addr[nt] = &WhichSpectrum;
  id[nt++] = INT;

  strcpy(tag[nt], "UnitVelocity_in_cm_per_s");
  addr[nt] = &UnitVelocity_in_cm_per_s;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "UnitLength_in_cm");
  addr[nt] = &UnitLength_in_cm;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "UnitMass_in_g");
  addr[nt] = &UnitMass_in_g;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "InputSpectrum_UnitLength_in_cm");
  addr[nt] = &InputSpectrum_UnitLength_in_cm;
  id[nt++] = FLOAT;


#ifdef ONLY_2LPT
  strcpy(tag[nt], "WDM_On");
  addr[nt] = &WDM_On;
  id[nt++] = INT;

  strcpy(tag[nt], "WDM_Vtherm_On");
  addr[nt] = &WDM_Vtherm_On;
  id[nt++] = INT;

  strcpy(tag[nt], "WDM_PartMass_in_kev");
  addr[nt] = &WDM_PartMass_in_kev;
  id[nt++] = FLOAT;
#endif

  WDM_On=0;
  
  
  if((fd = fopen(fname, "r")))
    {
      while(!feof(fd))
	{
	  buf[0] = 0;
	  fgets(buf, 200, fd);
	  

	  if(sscanf(buf, "%s%s%s", buf1, buf2, buf3) < 2)
	    continue;

	  if(buf1[0] == '%')
	    continue;

	  for(i = 0, j = -1; i < nt; i++)
	    if(strcmp(buf1, tag[i]) == 0)
	      {
		j = i;
		tag[i][0] = 0;
		break;
	      }

	  if(j >= 0)
	    {
	      switch (id[j])
		{
		case FLOAT:
		  *((double *) addr[j]) = atof(buf2);
		  break;
		case STRING:
		  strcpy(addr[j], buf2);
		  break;
		case INT:
		  *((int *) addr[j]) = atoi(buf2);
		  break;
		}
	    }
	  else
	    {
	      if(ThisTask == 0)
		fprintf(stderr, "\tError in file %s:   Tag '%s' not allowed or multiple defined.\n",fname,buf1);
	      	errorFlag = 1;
	    }
	}
      fclose(fd);

    }
  else
    {
      if(ThisTask == 0)
	fprintf(stderr, "\tParameter file %s not found.\n", fname);
      errorFlag = 1;
    }


  for(i = 0; i < nt; i++)
    {
      if(*tag[i])
	{
	  if(ThisTask == 0)
	    fprintf(stderr, "\tError. I miss a value for tag '%s' in parameter file '%s'.\n", tag[i], fname);
	  errorFlag = 1;
	}
    }

  if(errorFlag)
    {
      MPI_Finalize();
      exit(0);
    }


#undef FLOAT
#undef STRING
#undef INT
#undef MAXTAGS
}
