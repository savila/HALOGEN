///////////////////////////////////////////////////////////////////////
//                                                                   //
//   Copyright 2012 David Alonso                                     //
//                                                                   //
//                                                                   //
// This file is part of CUTE.                                        //
//                                                                   //
// CUTE is free software: you can redistribute it and/or modify it   //
// under the terms of the GNU General Public License as published by //
// the Free Software Foundation, either version 3 of the License, or //
// (at your option) any later version.                               //
//                                                                   //
// CUTE is distributed in the hope that it will be useful, but       //
// WITHOUT ANY WARRANTY; without even the implied warranty of        //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU //
// General Public License for more details.                          //
//                                                                   //
// You should have received a copy of the GNU General Public License //
// along with CUTE.  If not, see <http://www.gnu.org/licenses/>.     //
//                                                                   //
///////////////////////////////////////////////////////////////////////

/*********************************************************************/
//                               Main                                //
/*********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "cute.src/define.h"
#include "cute.src/common.h"
#include "cute.src/correlator.h"
#include "cute.src/io.h"
#include "cute.src/neighbors.h"


void run_monopole_corr_neighbors(Catalog cat_dat,double *corr,
								 double *ercorr,unsigned long long *DD,
								 int nb_r, double r_max, int nthreads)
{
  //////
  // Main routine for monopole using neighbor boxes
  
  // lint n_dat;
  int nside;
  // Catalog cat_dat;
 //fprintf(stderr,"a\n");
  NeighborBox *boxes;
 //fprintf(stderr,"a.1\n");

  nside=optimal_nside(l_box,r_max,cat_dat.np);

  boxes=catalog_to_boxes(nside,cat_dat);


 fprintf(stderr,"*** Correlating ***\n");
 // timer(0);
  corr_mono_box_neighbors(nside,boxes,cat_dat.np, cat_dat.pos,DD,nb_r,r_max,nthreads);
 fprintf(stderr,"***done***\n");
  
 // timer(1);

  make_CF(DD,cat_dat.np,corr,ercorr,nb_r,r_max);
  // write_CF(fnameOut,corr,ercorr,DD);


  free_catalog(cat_dat);
  free_boxes(nside,boxes);


 // timer(5);
}

int correlate(long N, float L, float *x, float *y, float *z,
			  double *corr, double *ercorr, unsigned long long *DD,
			  int nb_r, double r_max,int nthreads){
	
	Catalog cat = {0};
	long ii;
	
	// Allocate Catalog
	cat.np=N;
	#ifdef DEBUG
	fprintf(stderr,"\tAlloccing...\n");
	#endif
	cat.pos=(double *)malloc(3*N*sizeof(double));
	if (cat.pos==NULL){
		fprintf(stderr,"Couldn allocate cat.pos \n");
		exit(0);
	}
		
	fprintf(stderr,"...allocated\n");
	
	l_box = L;
	l_box_half = 0.5*L;
	// Copy data into catalog (seems like an unnecessary step...)
	fprintf(stderr,"Copying data...\n");
	for(ii=0;ii<N;ii++){
		cat.pos[3*ii] = (double) x[ii];
		cat.pos[3*ii+1] = (double) y[ii];
		cat.pos[3*ii+2] = (double) z[ii];
	}
	
	fprintf(stderr,"...data copied\n");
	
	// Run the 2PCF
	fprintf(stderr,"Runing monopole...\n");
	run_monopole_corr_neighbors(cat, corr, ercorr, DD,nb_r,r_max,nthreads);
	fprintf(stderr,"...Monopole done\n");
	
	return 0;
}

