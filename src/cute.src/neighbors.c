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
//         Boxing routines for nearest-neighbor searching            //
/*********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "define.h"
#include "common.h"

void free_boxes(int nside,NeighborBox *boxes) 
{
  //////
  // Frees all memory associated with a box
  // set of size nside
  int ii;

  for(ii=0;ii<nside*nside*nside;ii++) {
    if(boxes[ii].np>0)
      free(boxes[ii].pos);
  }
  
  free(boxes);
}

NeighborBox *catalog_to_boxes(int n_box_side,Catalog cat)
{
  //////
  // Creates boxes for nearest-neighbor searching
  long ii;
  int nside;
  NeighborBox *boxes;


  fprintf(stderr,"*** Building neighbor boxes \n");
  nside=n_box_side;
  fprintf(stderr,"  There will be %d boxes per side with a size of %lf \n",
	 nside,l_box/nside);
  
  boxes=(NeighborBox *)malloc(nside*nside*nside*sizeof(NeighborBox));
  if(boxes==NULL) error_mem_out();
  fprintf(stderr,"\tboxes allocated\n");
  for(ii=0;ii<nside*nside*nside;ii++)
    boxes[ii].np=0;
  fprintf(stderr,"\tboxes init\n");

  fprintf(stderr,"nside=%d, tot=%d\n",nside,nside*nside*nside);

  for(ii=0;ii<cat.np;ii++) {
    int ix,iy,iz;
 
    ix=(int)(cat.pos[3*ii]/l_box*nside);
    iy=(int)(cat.pos[3*ii+1]/l_box*nside);
    iz=(int)(cat.pos[3*ii+2]/l_box*nside);
 
    if (ix==l_box)
	ix=0.;
    if (iy==l_box)
	iy=0.;
    if (iz==l_box)
	iz=0.;
    //fprintf(stderr,"ii=%ld,  [%d,%d,%d]->%d\n",ii,ix,iy,iz,ix+nside*(iy+nside*iz));

    (boxes[ix+nside*(iy+nside*iz)].np)++;
  }
  fprintf(stderr,"\t counted\n");

  for(ii=0;ii<nside*nside*nside;ii++) {
    int npar=boxes[ii].np;
    if(npar>0) {
      boxes[ii].pos=(double *)malloc(3*npar*sizeof(double));
      if(boxes[ii].pos==NULL) {
		error_mem_out();
      }
      boxes[ii].np=0;
    }
  }
  fprintf(stderr,"\t allocated. Distributing...\n");

  for(ii=0;ii<cat.np;ii++) {
    int ix,iy,iz,index,offset;
    
    //fprintf(stderr,"\tii = %ld ",ii);
    ix=(int)(cat.pos[3*ii]/l_box*nside);
    iy=(int)(cat.pos[3*ii+1]/l_box*nside);
    iz=(int)(cat.pos[3*ii+2]/l_box*nside);
    index=ix+nside*(iy+nside*iz);
    //fprintf(stderr,"\t[i,j,k] = [%ld,%ld,%ld] lin=%ld",ix,iy,iz,index);
    //fprintf(stderr,"\t[x,y,z] = [%f,%f,%f] ",cat.pos[3*ii],cat.pos[3*ii+1],cat.pos[3*ii+2]);
    offset=3*boxes[index].np;
  //  fprintf(stderr,"\toffset = %ld ",offset);
    (boxes[index].pos)[offset]=cat.pos[3*ii];
    (boxes[index].pos)[offset+1]=cat.pos[3*ii+1];
    (boxes[index].pos)[offset+2]=cat.pos[3*ii+2];
    (boxes[index].np)++;
//    fprintf(stderr,"\tnp = %ld \n",boxes[index].np);
  }
  fprintf(stderr,"\t ...distributed\n");

//  printf("\n");
  
  return boxes;
}
