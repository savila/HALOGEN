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
//          Read-write routines for the periodic-box mode            //
/*********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "define.h"
#include "common.h"

void make_CF(unsigned long long *DD,int nD,
	     double *corr,double *ercorr,int nb_r,double rmax)
{
  //////
  // Creates correlation function and poisson errors
  // from pair counts DD, DR and RR
  double *edd;
  double rho_av=nD/(l_box*l_box*l_box);
  int ii;

  edd=(double *)malloc(sizeof(double)*nb_r);
  if(edd==NULL)
    error_mem_out();
  
#ifndef _LOGBIN
  DD[0]-=nD; //Substract diagonal
#endif //_LOGBIN
  for(ii=0;ii<nb_r;ii++)
    edd[ii]=1./sqrt((double)DD[ii]);

  for(ii=0;ii<nb_r;ii++) {
    if(DD[ii]==0) {
      corr[ii]=0;
      ercorr[ii]=0;
    }
    else {
      double r0,r1,vr,rho_r;
//#ifdef _LOGBIN
 //     r0=pow(10.,(((double)ii-NB_R)/N_LOGINT)+LOG_R_MAX);
  //    r1=pow(10.,(((double)ii+1-NB_R)/N_LOGINT)+LOG_R_MAX);
//#else //_LOGBIN
      r0=ii*rmax/(nb_r);
      r1=(ii+1)*rmax/(nb_r);
//#endif //_LOGBIN
      vr=4*M_PI*(r1*r1*r1-r0*r0*r0)/3;
      rho_r=DD[ii]/(nD*vr);
      corr[ii]=rho_r/rho_av-1;
      ercorr[ii]=(1+corr[ii])*edd[ii];
    }
  }

  free(edd);
}

