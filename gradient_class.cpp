/*
 *Implements the necessary class for gradient optimaztion in register3d
 *Contains coordinate of control points, gradient values, coefficients, etc.
 *Gustavo Rohde, fall 2000
 */


#include <stdlib.h>
#include <math.h>
#include "gradient_class.h"
#include "register3d.h"
#include <stdio.h>

void gradient_class::init(int lx,int rx,int ty,int by,int fz,int bz,int nptsx, int nptsy, int nptsz,
						  int sizex, int sizey, int sizez, data_type FLEX){
	
  int i,j,x,y,z;
  int vx_start,vx_end,vy_start,vy_end,vz_start,vz_end;
  //	data_type spcx,spcy,spcz;
  
  nx = nptsx;
  ny = nptsy;
  nz = nptsz;
  num_points = nptsx*nptsy*nptsz;
  
  //allocate space
  control_points = (int *)malloc(num_points*3*sizeof(int));	
  gradient = (data_type *)malloc(num_points*3*sizeof(data_type));
  coefficients = (data_type *)malloc(num_points*3*sizeof(data_type));
  
  
  if ((control_points || gradient || coefficients)==NULL){
    printf("not enough memory for gradient structure.\n");
    exit(1);
  }
  
  
  //calculate support
  spcx = (data_type)(myabs(((data_type)(rx-lx))/((data_type)(nptsx-1))));
  spcy = (data_type)(myabs(((data_type)(ty-by))/((data_type)(nptsy-1))));
  spcz = (data_type)(myabs(((data_type)(bz-fz))/((data_type)(nptsz-1))));
  
  
  supp = (int)mymax(spcx,spcy);
  supp = (int)ceil(FLEX*supp);
  suppz = (int)ceil(FLEX*spcz);
  //printf("Inside grd.init. Supp = %d , Suppz = %d, spcz = %f, spcx = %f, spcy = %f \n", supp, suppz,spcz,spcx,spcy);

  
  //supp = 40;
  //initialize structures
  /*i=0;
    for(x=lx;x<=rx;x=x+(int)spcx){
    for(y=ty;y<=by;y=y+(int)spcy){
    for(z=fz;z<=bz;z=z+(int)spcz){
    set(i,0,(data_type)floor(x),'p');
    set(i,1,(data_type)floor(y),'p');
    set(i,2,(data_type)floor(z),'p');
    i++;
    }
    }
    }*/
  i=0;
  for(x=0;x<nptsx;x++){
    for(y=0;y<nptsy;y++){
      for(z=0;z<nptsz;z++){
	
	set(i,0,(data_type)floor(lx+x*spcx),'p');
	set(i,1,(data_type)floor(ty+y*spcy),'p');
	set(i,2,(data_type)floor(fz+z*spcz),'p');
	i++;
      }
    }
  }
  
  for(i=0;i<num_points;i++){
    for(j=0;j<3;j++){
      set(i,j,0,'c');
      set(i,j,0,'g');
    }
  }
  
  
  minx = sizex-1;
  maxx = 0;
  miny = sizey-1;
  maxy = 0;
  minz = sizez-1;
  maxz = 0;
  
  
  for(i=0;i<num_points;i++){
    
    //find appropriate indexes for acessing
    vx_start = (int)mymax((get(i,0,'p')-supp),0);
    vx_end = (int)mymin((get(i,0,'p')+supp+1),(sizex-1));
    vy_start = (int)mymax((get(i,1,'p')-supp),0);
    vy_end = (int)mymin((get(i,1,'p')+supp+1),(sizey-1));
    vz_start = (int)mymax((get(i,2,'p')-suppz),0);
    vz_end = (int)mymin((get(i,2,'p')+suppz+1),(sizez-1));
    
    if (vx_start<minx)
      minx = vx_start;
    if (vx_end>maxx)
      maxx = vx_end;
    if (vy_start<miny)
      miny = vy_start;
    if (vy_end>maxy)
      maxy = vy_end;
    if (vz_start<minz)
      minz = vz_start;
    if (vz_end>maxz)
      maxz = vz_end;
    
  }
  
  
}

/**
 *
 * Secondary initialization function. Meant to be used after hot spots in deformation field have been identified.
 * Takes in the centre coordinate of the hot spot and the spacings with which it will build the new gradient structure.
 * The new gradient structure is allways a cube ( 8 points )
 *
 */

void gradient_class::init2(int x, int y, int z, int sx, int sy, int sz, int sp, int spz, int sizex, int sizey, int sizez){

  int vx_start,vx_end,vy_start,vy_end,vz_start,vz_end;
  int lx,rx,ty,by,fz,bz;

  supp = sp;
  suppz = spz;

  nx = 2;
  ny = 2;
  nz = 2;
  num_points = 8; //cube around location x,y,z

  //allocate space
  control_points = (int *)malloc(num_points*3*sizeof(int));	
  gradient = (data_type *)malloc(num_points*3*sizeof(data_type));
  coefficients = (data_type *)malloc(num_points*3*sizeof(data_type));

  if ((control_points || gradient || coefficients)==NULL){
    printf("Not enough memory for gradient structure.\n");
    exit(1);
  }

  spcx = sx;
  spcy = sy;
  spcz = sz;

  /*
  set(0,0,(data_type)max(0,floor(x-(data_type)sx/2)),'p');
  set(0,1,(data_type)max(0,floor(y-(data_type)sy/2)),'p');
  set(0,2,(data_type)max(0,floor(z-(data_type)sz/2)),'p');

  set(1,0,(data_type)min((sizex-1),floor(x+(data_type)sx/2)),'p');
  set(1,1,(data_type)max(0,floor(y-(data_type)sy/2)),'p');
  set(1,2,(data_type)max(0,floor(z-(data_type)sz/2)),'p');
 
  set(2,0,(data_type)max(0,floor(x-(data_type)sx/2)),'p');
  set(2,1,(data_type)min((sizey-1),floor(y+(data_type)sy/2)),'p');
  set(2,2,(data_type)max(0,floor(z-(data_type)sz/2)),'p');

  set(3,0,(data_type)min((sizex-1),floor(x+(data_type)sx/2)),'p');
  set(3,1,(data_type)min((sizey-1),floor(y+(data_type)sy/2)),'p');
  set(3,2,(data_type)max(0,floor(z-(data_type)sz/2)),'p');

  set(4,0,(data_type)max(0,floor(x-(data_type)sx/2)),'p');
  set(4,1,(data_type)max(0,floor(y-(data_type)sy/2)),'p');
  set(4,2,(data_type)min((sizez-1),floor(z+(data_type)sz/2)),'p');

  set(5,0,(data_type)min((sizex-1),floor(x+(data_type)sx/2)),'p');
  set(5,1,(data_type)max(0,floor(y-(data_type)sy/2)),'p');
  set(5,2,(data_type)min((sizez-1),floor(z+(data_type)sz/2)),'p');

  set(6,0,(data_type)max(0,floor(x-(data_type)sx/2)),'p');
  set(6,1,(data_type)min((sizey-1),floor(y+(data_type)sy/2)),'p');
  set(6,2,(data_type)min((sizez-1),floor(z+(data_type)sz/2)),'p');

  set(7,0,(data_type)min((sizex-1),floor(x+(data_type)sx/2)),'p');
  set(7,1,(data_type)min((sizey-1),floor(y+(data_type)sy/2)),'p');
  set(7,2,(data_type)min((sizez-1),floor(z+(data_type)sz/2)),'p');
  */

  int i=0;
  
  lx  = (int)mymax(0,floor(x-(data_type)sx/2.0));
  rx =  (int)mymin((sizex-1),floor(x+(data_type)sx/2.0));
  ty = (int)mymax(0,floor(y-(data_type)sy/2.0));
  by = (int)mymin((sizey-1),floor(y+(data_type)sy/2.0));
  fz = (int)mymax(0,floor(z-(data_type)sz/2.0));
  bz = (int)mymin((sizez-1),floor(z+(data_type)sz/2.0));
  
  int xx,yy,zz;
  for(xx=0;xx<nx;xx++){
    for(yy=0;yy<ny;yy++){
      for(zz=0;zz<nz;zz++){
	
	set(i,0,(data_type)floor(lx+xx*spcx),'p');
	set(i,1,(data_type)floor(ty+yy*spcy),'p');
	set(i,2,(data_type)floor(fz+zz*spcz),'p');
	i++;
      }
    }
  }
 
  
  for( i=0; i<num_points; i++){
    for(int j=0; j<3; j++){
      set(i,j,0,'c');
      set(i,j,0,'g');
    }
  }


  minx = sizex-1;
  maxx = 0;
  miny = sizey-1;
  maxy = 0;
  minz = sizez-1;
  maxz = 0;

  
  for(i=0;i<num_points;i++){
    
    //find appropriate indexes for acessing
    vx_start = (int)mymax((get(i,0,'p')-supp),0);
    vx_end = (int)mymin((get(i,0,'p')+supp+1),(sizex-1));
    vy_start = (int)mymax((get(i,1,'p')-supp),0);
    vy_end = (int)mymin((get(i,1,'p')+supp+1),(sizey-1));
    vz_start = (int)mymax((get(i,2,'p')-suppz),0);
    vz_end = (int)mymin((get(i,2,'p')+suppz+1),(sizez-1));
    
    if (vx_start<minx)
      minx = vx_start;
    if (vx_end>maxx)
      maxx = vx_end;
    if (vy_start<miny)
      miny = vy_start;
    if (vy_end>maxy)
      maxy = vy_end;
    if (vz_start<minz)
      minz = vz_start;
    if (vz_end>maxz)
      maxz = vz_end;
    
  }

  //printf("\n");
  //printf("   Inside init2.\n");
  //printf("   Centre point %d, %d, %d\n", x, y, z);
  //printf("   Spacings %d, %d, %d\n", sx, sy, sz);
  //printf("   Active zone: %d, %d, %d, %d, %d, %d \n", minx, maxx, miny, maxy, minz, maxz);
  //for (i=0;i<num_points;i++){
  //  printf("      %d, %d, %d \n", (int)get(i,0,'p'), (int)get(i,1,'p'), (int)get(i,2,'p') );
  //}
  

}//of function


/*
 *updates the coefficients field of this class
 */


void gradient_class::coeff_update(data_type konst){

	int i,j;
	data_type new_val;

	for(i=0;i<num_points;i++){
		for(j=0;j<3;j++){
			
			new_val = get(i,j,'c') - konst*get(i,j,'g');
			set(i,j,new_val,'c');

		}
	}

}

void gradient_class::normalize(void){

  int i,j;
  data_type normmm;
  normmm = norm();
  for (i=0;i<num_points;i++){
    for (j=0;j<3;j++){
      set(i,j,get(i,j,'g')/normmm,'g');
    }
  }

}


void gradient_class::max_d_compute(data_type trsh){

  long i,j;
  data_type ck,ckp1,ckm1,gk,gkp1,gkm1;
  data_type d1,d2, K;
  data_type a1,a2;
  data_type a,b;
  
  //printf("Threshold recieved %f \n", trsh);

  min_a = 5000;
  //trsh = 1.0/3.0;

  //loop through all values
  for (i=0;i<num_points;i++){
    
    //loop through X,Y,Z components
    for (j=0;j<3;j++){
      
      //x
      gk = get(i,j,'g');
      ck = get(i,j,'c');
     
      if ( (i+nz*ny)<num_points ){
	ckp1 = get(i+nz*ny,j,'c');
	gkp1 = get(i+nz*ny,j,'g');
      }else{
	ckp1 = 0;
	gkp1 = 0;
      }

      if ( (i-nz*ny)>0 ){
	ckm1 = get(i-nz*ny,j,'c');
	gkm1 = get(i-nz*ny,j,'g');
      }else{
	ckm1 = 0;
	gkm1 = 0;
      }
      
      a = (ck-ckm1);
      b = (gkm1-gk);
      K = supp*trsh/8.0;
      if (b<0)
	K = -K;
      a1 = (K - a)/b;

      a = (ckp1-ck);
      b = (gk-gkp1);
      K = supp*trsh/8.0;
      if (b<0)
	K = -K;
      a2 = (K - a)/b;

      
      min_a=keep_min(min_a,a1,a2);
      //keep minimum alpha instead

      //y
     
      if ( (i+nz)<num_points ){
	ckp1 = get(i+nz,j,'c');
	gkp1 = get(i+nz,j,'g');
      }else{
	ckp1 = 0;
	gkp1 = 0;
      }
      
      if ( (i-nz)>0 ){
	ckm1 = get(i-nz,j,'c');
	gkm1 = get(i-nz,j,'g');
      }else{
	ckm1 = 0;
	gkm1 = 0;
      }
      
      a = (ck-ckm1);
      b = (gkm1-gk);
      K = supp*trsh/8.0;
      if (b<0)
	K = -K;
      a1 = (K - a)/b;

      a = (ckp1-ck);
      b = (gk-gkp1);
      K = supp*trsh/8.0;
      if (b<0)
	K = -K;
      a2 = (K - a)/b;

      
      min_a=keep_min(min_a,a1,a2);
      
      //z
     
      if ( (i+1)<num_points ){
	ckp1 = get(i+1,j,'c');
	gkp1 = get(i+1,j,'g');
      }else{
	ckp1 = 0;
	gkp1 = 0;
      }
      
      if ( (i-1)>0 ){
	ckm1 = get(i-1,j,'c');
	gkm1 = get(i-1,j,'g');
      }else{
	ckm1 = 0;
	gkm1 = 0;
      }
      
      a = (ck-ckm1);
      b = (gkm1-gk);
      K = suppz*trsh/8.0;
      if (b<0)
	K = -K;
      a1 = (K - a)/b;

      a = (ckp1-ck);
      b = (gk-gkp1);
      K = supp*trsh/8.0;
      if (b<0)
	K = -K;
      a2 = (K - a)/b;

      
      min_a=keep_min(min_a,a1,a2);

    }


  }

  
  
  //printf("Maximum alpha allowed: %f\n", min_a);

}

data_type gradient_class::keep_min(data_type a1, data_type a2, data_type a3){

  data_type result;

  result = a1;
  if (a2<result)
    result = a2;
  if (a3<result)
    result = a3;

  return result;

}

data_type gradient_class::keep_max(data_type m1, data_type m2, data_type m3){

  data_type result;

  if (m1<0)
    m1 = -m1;
  if (m2<0)
    m2 = -m2;
  if (m3<0)
    m3 = -m3;

  result = m1;
  if (m2>result)
    result = m2;
  if (m3>result)
    result = m3;

  return result;


}



data_type gradient_class::norm(void){

  int i,j;
  data_type res;

  res = 0;
  for (i=0;i<num_points;i++){
    for (j=0;j<3;j++){
      res = res + get(i,j,'g')*get(i,j,'g');
    }
  }

  res = res/(num_points*3);
  return sqrt(res);


}

data_type gradient_class::get(int x, int y, char w){

	long coord;
	data_type result;

	coord = x*3 + y;
	

	switch (w){

	case 'p':
		result = (data_type)control_points[coord];
		break;
	case 'c':
		result = coefficients[coord];
		break;
	case 'g':
		result = gradient[coord];
		break;

	}

	return result;

}

void gradient_class::set(int x, int y, data_type val, char w){

	long coord;

	coord = x*3 + y;
		
	switch (w){

	case 'p':
		control_points[coord]=(int)val;
		break;
	case 'c':
		coefficients[coord]=(data_type)val;
		break;
	case 'g':
		gradient[coord]=(data_type)val;
		break;

	}

}

void gradient_class::destroy(){

	free(coefficients);
	free(control_points);
	free(gradient);
	num_points = 0;

}


void gradient_class::reset(){

	free(coefficients);
	free(control_points);
	free(gradient);
	num_points = 0;

}
