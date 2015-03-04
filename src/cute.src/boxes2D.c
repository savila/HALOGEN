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
//                    Pixel functions and routines                   //
/*********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "define.h"
#include "common.h"

static double cth_min_bound;
static double cth_max_bound;
static double phi_min_bound;
static double phi_max_bound;

static int estimate_optimal_nside_radial(void)
{
  return 5*(int)(M_PI/aperture_los);
}

static int estimate_optimal_nside_angular(int np,double fsky)
{
  int n_side1=5*(int)(M_PI*I_THETA_MAX);
  int n_side2=(int)(sqrt(0.25*np/fsky));
  
  return MIN(n_side1,n_side2);
}

static int sph2pix(double cth,double phi)
{
  //////
  // Returns ipix to pixel with coordinates cth,phi

  int icth,iphi;
  if((cth<-1)||(cth>1)) {
    fprintf(stderr,"Wrong cos(theta) = %lf \n",cth);
    exit(1);
  }
  else if(cth==1)
    icth=n_side_cth-1;
  else
    icth=(int)(0.5*(1+cth)*n_side_cth);

  iphi=(int)(0.5*wrap_phi(phi)/M_PI*n_side_phi);
  
  return iphi+icth*n_side_phi;
}

//Cell2D
void free_Cells2D(int npix,Cell2D *cells)
{
  int ii;
  for(ii=0;ii<npix;ii++)
    if(cells[ii].np>0) free(cells[ii].ci);

  free(cells);
}

static Cell2D *init_Cells2D(int npix)
{
  int ii;
  Cell2D *cells=(Cell2D *)malloc(npix*sizeof(Cell2D));
  if(cells==NULL) error_mem_out();

  for(ii=0;ii<npix;ii++) {
    cells[ii].np=-1;
    cells[ii].ci=NULL;
  }

  return cells;
}

//Box2D
static void free_Box2DInfo(Box2DInfo *bi)
{
  free(bi->pos);
  free(bi);
}

void free_Boxes2D(int npix,Box2D *boxes)
{
  int ii;
  for(ii=0;ii<npix;ii++) {
    if(boxes[ii].np>0) 
      free_Box2DInfo(boxes[ii].bi);
  }

  free(boxes);
}

static Box2DInfo *init_Box2DInfo(int np)
{
  Box2DInfo *bi=(Box2DInfo *)malloc(sizeof(Box2DInfo));
  if(bi==NULL) error_mem_out();
  bi->pos=(double *)malloc(N_POS*np*sizeof(double));
  if(bi->pos==NULL) error_mem_out();

  return bi;
}

static Box2D *init_Boxes2D(int npix)
{
  int ii;
  Box2D *boxes=(Box2D *)malloc(npix*sizeof(Box2D));
  if(boxes==NULL) error_mem_out();

  for(ii=0;ii<npix;ii++) {
    boxes[ii].np=0;
    boxes[ii].bi=NULL;
  }

  return boxes;
}

//RadialPixels
static void free_RadialPixelInfo(RadialPixelInfo *pi)
{
  free(pi->pos);
  free(pi->redshifts);
  free(pi);
}

void free_RadialPixels(int npix,RadialPixel *pixrad)
{
  int ii;
  for(ii=0;ii<npix;ii++) {
    if(pixrad[ii].np>0) 
      free_RadialPixelInfo(pixrad[ii].pi);
  }

  free(pixrad);
}

static RadialPixelInfo *init_RadialPixelInfo(int np)
{
  RadialPixelInfo *pi=(RadialPixelInfo *)malloc(sizeof(RadialPixelInfo));
  if(pi==NULL) error_mem_out();
  pi->pos=(double *)malloc(N_POS*np*sizeof(double));
  if(pi->pos==NULL) error_mem_out();
  pi->redshifts=(double *)malloc(np*sizeof(double));
  if(pi->redshifts==NULL) error_mem_out();

  return pi;
}

static RadialPixel *init_RadialPixels(int npix)
{
  int ii;
  RadialPixel *pixrad=(RadialPixel *)malloc(npix*sizeof(RadialPixel));
  if(pixrad==NULL) error_mem_out();

  for(ii=0;ii<npix;ii++) {
    pixrad[ii].np=0;
    pixrad[ii].pi=NULL;
  }

  return pixrad;
}

void init_2D_params(Catalog cat_dat,Catalog cat_ran,int ctype)
{
  int ii;

  cth_min_bound=cat_dat.cth[0];
  cth_max_bound=cat_dat.cth[0];
  phi_min_bound=cat_dat.phi[0];
  phi_max_bound=cat_dat.phi[0];

  for(ii=0;ii<cat_dat.np;ii++) {
    double cth=cat_dat.cth[ii];
    double phi=cat_dat.phi[ii];

    if(cth<cth_min_bound) cth_min_bound=cth;
    if(phi<phi_min_bound) phi_min_bound=phi;
    if(cth>cth_max_bound) cth_max_bound=cth;
    if(phi>phi_max_bound) phi_max_bound=phi;
  }

  for(ii=0;ii<cat_ran.np;ii++) {
    double cth=cat_ran.cth[ii];
    double phi=cat_ran.phi[ii];

    if(cth<cth_min_bound) cth_min_bound=cth;
    if(phi<phi_min_bound) phi_min_bound=phi;
    if(cth>cth_max_bound) cth_max_bound=cth;
    if(phi>phi_max_bound) phi_max_bound=phi;
  }

  if(ctype==0) {
    n_side_cth=estimate_optimal_nside_radial();
    n_side_phi=2*n_side_cth;
  }
  else if(ctype==1) {
    if(!use_pm) {
      double fsky=(cth_max_bound-cth_min_bound)*
	(phi_max_bound-phi_min_bound)/(4*M_PI);
      n_side_cth=estimate_optimal_nside_angular(cat_dat.np,fsky);
      n_side_phi=2*n_side_cth;
    }
  }
  else if(ctype==5) {
    double fsky=(cth_max_bound-cth_min_bound)*
      (phi_max_bound-phi_min_bound)/(4*M_PI);
    n_side_cth=estimate_optimal_nside_angular(cat_dat.np,fsky);
    n_side_phi=2*n_side_cth;
  }
  else {
    fprintf(stderr,"WTF?? \n");
    exit(1);
  }

  n_boxes2D=n_side_phi*n_side_cth;

  double pixel_resolution=sqrt(4*M_PI/n_boxes2D)/DTORAD;
  printf("  There will be %d pixels in total\n",n_boxes2D);
  printf("  Pixel resolution is %.4lf deg \n",pixel_resolution);
}

static void get_pix_bounds(double alpha,int ipix,
			   int *icth_min,int *icth_max,
			   int *iphi_min,int *iphi_max)
{
  //////
  // Returns pixel bounds for all pixels within
  // theta_max=alpha
  int icth,iphi;
  double theta,th_hi,th_lo;
  double phi_hi,phi_lo;
  double cth_max,cth_min;

  icth=(int)(ipix/n_side_phi);
  iphi=(int)(ipix%n_side_phi);

  theta=acos(-1.0+2.0*((double)(icth+0.5))/n_side_cth);
  th_hi=acos(-1.0+2.0*((double)(icth+0.0))/n_side_cth);
  th_lo=acos(-1.0+2.0*((double)(icth+1.0))/n_side_cth);
  phi_hi=2*M_PI*((double)(iphi+1.0)/n_side_phi);
  phi_lo=2*M_PI*((double)(iphi+0.0)/n_side_phi);

  if(th_hi>M_PI-alpha) {
    cth_min=-1;
    cth_max=cos(th_lo-alpha);

    *iphi_min=0;
    *iphi_max=n_side_phi-1;
  }
  else if(th_lo<alpha) {
    cth_min=cos(th_hi+alpha);
    cth_max=1;

    *iphi_min=0;
    *iphi_max=n_side_phi-1;
  }
  else {
    double dphi;
    double calpha=cos(alpha);
    cth_min=cos(th_hi+alpha);
    cth_max=cos(th_lo-alpha);

    if(theta<0.5*M_PI) {
      double c_thlo=cos(th_lo);
      dphi=acos(sqrt((calpha*calpha-c_thlo*c_thlo)/
		     (1-c_thlo*c_thlo)));
    }
    else {
      double c_thhi=cos(th_hi);
      dphi=acos(sqrt((calpha*calpha-c_thhi*c_thhi)/
		     (1-c_thhi*c_thhi)));
    }

    if(dphi<M_PI) {
      double phi_max,phi_min;
      phi_min=phi_lo-dphi;
      phi_max=phi_hi+dphi;
      *iphi_min=(int)(0.5*phi_min/M_PI*n_side_phi);
      *iphi_max=(int)(0.5*phi_max/M_PI*n_side_phi);
    }
    else {
      *iphi_min=0;
      *iphi_max=n_side_phi-1;
    }
  }

  //Cut with mask
  cth_min=MAX((cth_min),(cth_min_bound));
  cth_max=MIN((cth_max),(cth_max_bound));

  *icth_min=(int)(0.5*(1+cth_min)*n_side_cth);
  *icth_max=(int)(0.5*(1+cth_max)*n_side_cth);
  if(*icth_max>=n_side_cth) *icth_max=n_side_cth-1;
  if(*icth_min<0) *icth_min=0;
}

Cell2D *mk_Cells2D_from_Catalog(Catalog cat,int **cell_indices,int *n_cell_full)
{
  int ii,nfull;
  Cell2D *cells;

  cells=init_Cells2D(n_boxes2D);
  
  nfull=0;
  for(ii=0;ii<cat.np;ii++) {
    double cth=cat.cth[ii];
    double phi=cat.phi[ii];
    int ipix=sph2pix(cth,phi);
    np_t np0=cells[ipix].np;
    if(np0<0) {
      nfull++;
      cells[ipix].np=0;
    }

#ifdef _WITH_WEIGHTS
    cells[ipix].np+=cat.weight[ii];
#else //_WITH_WEIGHTS
    cells[ipix].np++;
#endif //_WITH_WEIGHTS
  }
  
  *n_cell_full=nfull;
  printf("  There are objects in %d out of %d pixels\n",nfull,n_boxes2D);
  *cell_indices=(int *)malloc(nfull*sizeof(int));
  if(cell_indices==NULL) error_mem_out();

  nfull=0;
  for(ii=0;ii<n_boxes2D;ii++) {
    if(cells[ii].np>0) {
      int icth_min,icth_max;
      int iphi_min,iphi_max;
      int icth=(int)(ii/n_side_phi);
      int iphi=(int)(ii%n_side_phi);
      double cth=-1.0+2.0*((double)(icth+0.5))/n_side_cth;
      double phi=2*M_PI*((double)(iphi+0.5))/n_side_phi;
      double sth=sqrt(1-cth*cth);

      //Allocate cell info
      cells[ii].ci=(Cell2DInfo *)malloc(sizeof(Cell2DInfo));
      if(cells[ii].ci==NULL) error_mem_out();

      //Calculate cell bounds
      get_pix_bounds(1/I_THETA_MAX,ii,
		     &icth_min,&icth_max,&iphi_min,&iphi_max);
      (cells[ii].ci)->bounds[0]=icth_min;
      (cells[ii].ci)->bounds[1]=icth_max;
      (cells[ii].ci)->bounds[2]=iphi_min;
      (cells[ii].ci)->bounds[3]=iphi_max;

      //Calculate cell position
      (cells[ii].ci)->pos[0]=sth*cos(phi);
      (cells[ii].ci)->pos[1]=sth*sin(phi);
      (cells[ii].ci)->pos[2]=cth;

      //Get pixel index
      (*cell_indices)[nfull]=ii;
      nfull++;
    }
    else {
      cells[ii].np=0;
    }
  }

  return cells;
}

Box2D *mk_Boxes2D_from_Catalog(Catalog cat,int **box_indices,int *n_box_full)
{
  int ii,nfull;
  Box2D *boxes;

  boxes=init_Boxes2D(n_boxes2D);

  nfull=0;
  for(ii=0;ii<cat.np;ii++) {
    double cth=cat.cth[ii];
    double phi=cat.phi[ii];
    int ipix=sph2pix(cth,phi);
    int np0=boxes[ipix].np;
    if(np0==0) nfull++;
    boxes[ipix].np++;
  }

  *n_box_full=nfull;
  printf("  There are objects in %d out of %d pixels \n",nfull,n_boxes2D);
  *box_indices=(int *)malloc(nfull*sizeof(int));
  if(box_indices==NULL) error_mem_out();
  
  nfull=0;
  for(ii=0;ii<n_boxes2D;ii++) {
    if(boxes[ii].np>0) {
      int icth_min,icth_max,iphi_min,iphi_max;
      
      //Allocate box info
      boxes[ii].bi=init_Box2DInfo(boxes[ii].np);
      boxes[ii].np=0;

      //Calculate box bounds
      get_pix_bounds(1/I_THETA_MAX,ii,
		     &icth_min,&icth_max,&iphi_min,&iphi_max);
      (boxes[ii].bi)->bounds[0]=icth_min;
      (boxes[ii].bi)->bounds[1]=icth_max;
      (boxes[ii].bi)->bounds[2]=iphi_min;
      (boxes[ii].bi)->bounds[3]=iphi_max;

      //Get pixel index
      (*box_indices)[nfull]=ii;
      nfull++;
    }
  }

  for(ii=0;ii<cat.np;ii++) {
    double cth=cat.cth[ii];
    double phi=cat.phi[ii];
    double sth=sqrt(1-cth*cth);
    int ipix=sph2pix(cth,phi);
    int np0=boxes[ipix].np;
    (boxes[ipix].bi)->pos[N_POS*np0]=sth*cos(phi);
    (boxes[ipix].bi)->pos[N_POS*np0+1]=sth*sin(phi);
    (boxes[ipix].bi)->pos[N_POS*np0+2]=cth;
#ifdef _WITH_WEIGHTS
    (boxes[ipix].bi)->pos[N_POS*np0+3]=cat.weight[ii];
#endif //_WITH_WEIGHTS
    boxes[ipix].np++;
  }

  return boxes;
}

RadialPixel *mk_RadialPixels_from_Catalog(Catalog cat,int **pixrad_indices,
					  int *n_pixrad_full,int ctype)
{
  int ii,nfull;
  RadialPixel *pixrad;

  pixrad=init_RadialPixels(n_boxes2D);

  nfull=0;
  for(ii=0;ii<cat.np;ii++) {
    double cth=cat.cth[ii];
    double phi=cat.phi[ii];
    int ipix=sph2pix(cth,phi);
    int np0=pixrad[ipix].np;
    if(np0==0) nfull++;
    pixrad[ipix].np++;
  }

  *n_pixrad_full=nfull;
  printf("  There are objects in %d out of %d pixels \n",nfull,n_boxes2D);
  *pixrad_indices=(int *)malloc(nfull*sizeof(int));
  if(pixrad_indices==NULL) error_mem_out();
  
  nfull=0;
  double aperture;
  if(ctype==0) aperture=aperture_los;
  else if(ctype==5) aperture=1./I_THETA_MAX;
  else {
    fprintf(stderr,"WTF??\n");
    exit(1);
  }
  for(ii=0;ii<n_boxes2D;ii++) {
    if(pixrad[ii].np>0) {
      int icth_min,icth_max,iphi_min,iphi_max;
      
      //Allocate box info
      pixrad[ii].pi=init_RadialPixelInfo(pixrad[ii].np);
      pixrad[ii].np=0;

      //Calculate box bounds
      get_pix_bounds(aperture,ii,
		     &icth_min,&icth_max,&iphi_min,&iphi_max);
      (pixrad[ii].pi)->bounds[0]=icth_min;
      (pixrad[ii].pi)->bounds[1]=icth_max;
      (pixrad[ii].pi)->bounds[2]=iphi_min;
      (pixrad[ii].pi)->bounds[3]=iphi_max;

      //Get pixel index
      (*pixrad_indices)[nfull]=ii;
      nfull++;
    }
  }

  for(ii=0;ii<cat.np;ii++) {
    double zz=cat.red[ii];
    double cth=cat.cth[ii];
    double phi=cat.phi[ii];
    double sth=sqrt(1-cth*cth);
    int ipix=sph2pix(cth,phi);
    int np0=pixrad[ipix].np;
    (pixrad[ipix].pi)->pos[N_POS*np0]=sth*cos(phi);
    (pixrad[ipix].pi)->pos[N_POS*np0+1]=sth*sin(phi);
    (pixrad[ipix].pi)->pos[N_POS*np0+2]=cth;
#ifdef _WITH_WEIGHTS
    (pixrad[ipix].pi)->pos[N_POS*np0+3]=cat.weight[ii];
#endif //_WITH_WEIGHTS
    (pixrad[ipix].pi)->redshifts[np0]=zz;
    pixrad[ipix].np++;
  }

  return pixrad;
}
