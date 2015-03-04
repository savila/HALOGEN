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
//                      Correlators with OpenMP                      //
/*********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "define.h"
#include "common.h"

double I_DR;
double R2_MAX;

void corr_mono_box_neighbors(int nside,NeighborBox *boxes,
			     	 	 	 long np,double *pos,
			     	 	 	 unsigned long long *hh, int nb_r, double r_max,int nthreads)
{

  I_DR = nb_r/r_max;
  R2_MAX = r_max*r_max;

  //////
  // Correlator for monopole in the periodic-box case
  // using neighbor boxes  
  double agrid=l_box/nside;
  //double r_max=1/I_R_MAX;
  int index_max=(int)(r_max/agrid)+1;
  int i;

  fprintf(stderr,"  Boxes will be correlated up to %d box sizes \n",index_max);

  for(i=0;i<nb_r;i++)
    hh[i]=0; //Clear shared histogram

 #pragma omp parallel num_threads(nthreads) default(none) shared(index_max,nside,boxes,hh,l_box,np,pos,agrid,stderr,nb_r,I_DR,R2_MAX)
  {
    long ii;
    double a2grid=agrid*agrid;
    unsigned long long hthread[nb_r]; //Histogram filled by each thread

    for(ii=0;ii<nb_r;ii++)
      hthread[ii]=0; //Clear private histogram

 #pragma omp for  nowait schedule(dynamic)
    for(ii=0;ii<np;ii++) {
      int ix0,iy0,iz0;
      double x0,y0,z0;
      int idz;
      x0=pos[3*ii];
      y0=pos[3*ii+1];
      z0=pos[3*ii+2];
      ix0=(int)(x0/l_box*nside);
      iy0=(int)(y0/l_box*nside);
      iz0=(int)(z0/l_box*nside);

      for(idz=-index_max;idz<=index_max;idz++) {
	int idy,idz_dist2;
	int iwrapz=0;
	int iz1=iz0+idz;
	if(iz1<0) {
	  iz1+=nside;
	  iwrapz=1;
	}
	else if(iz1>=nside) {
	  iz1-=nside;
	  iwrapz=1;
	}
	idz_dist2=MAX(0,abs(idz)-1);
	idz_dist2=idz_dist2*idz_dist2;
	for(idy=-index_max;idy<=index_max;idy++) {
	  int idx,idy_dist2;
	  int iwrapy=0;
	  int iy1=iy0+idy;
	  if(iy1<0) {
	    iy1+=nside;
	    iwrapy=1;
	  }
	  else if(iy1>=nside) {
	    iy1-=nside;
	    iwrapy=1;
	  }
	  idy_dist2=MAX(0,abs(idy)-1);
	  idy_dist2=idy_dist2*idy_dist2;
	  for(idx=-index_max;idx<=index_max;idx++) {
	    int ibox,idx_dist;
	    int iwrapx=0;
	    int ix1=ix0+idx;
	    double d2max;
	    int jj;
	    if(ix1<0) {
	      ix1+=nside;
	      iwrapx=1;
	    }
	    else if(ix1>=nside) {
	      ix1-=nside;
	      iwrapx=1;
	    }
	    ibox=ix1+nside*(iy1+nside*iz1);
	    idx_dist=MAX(0,abs(idx)-1);
	    d2max=a2grid*(idx_dist*idx_dist+idy_dist2+idz_dist2);
	    if(d2max>R2_MAX) continue;
	    for(jj=0;jj<boxes[ibox].np;jj++) {
	      double xr[3];
	      double r2;
	      int ir;

	      xr[0]=fabs(x0-(boxes[ibox].pos)[3*jj]);
	      xr[1]=fabs(y0-(boxes[ibox].pos)[3*jj+1]);

	      xr[2]=fabs(z0-(boxes[ibox].pos)[3*jj+2]);


	      if(iwrapx) xr[0]=l_box-xr[0];

	      if(iwrapy) xr[1]=l_box-xr[1];

	      if(iwrapz) xr[2]=l_box-xr[2];

	      r2=xr[0]*xr[0]+xr[1]*xr[1]+xr[2]*xr[2];
	      if(r2>R2_MAX) continue;

//#ifdef _LOGBIN
//	      if(r2>0) {
//		ir=(int)(N_LOGINT*(0.5*log10(r2)-LOG_R_MAX)+NB_R);
//		if((ir<NB_R)&&(ir>=0))
//		  (hthread[ir])++;
//	      }
//#else //_LOGBIN
	      ir=(int)(sqrt(r2)*I_DR);
	      if(ir<nb_r) //Check bound
	    	  (hthread[ir])++;
//#endif //_LOGBIN
	    }
	  }
	}
      }
    } // end pragma omp for
    fprintf(stderr,"end of pragma omp for\n");
 #pragma omp critical
    {
      for(ii=0;ii<nb_r;ii++) //Check bound
	hh[ii]+=hthread[ii]; //Add private histograms to shared one
    } // end pragma omp critical
  } // end pragma omp parallel
}
