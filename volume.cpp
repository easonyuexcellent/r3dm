/*file: volume.cpp
Implements data structure and functions to hold a volume of data.
Uses matrix3d.cpp
Gustavo Rohde
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include "register3d.h"
#include "volume.h"


/*
 *allocates memory
 */
void volume::init(int x, int y,int z, int hs, char *flname){

  int i,j;//,k;

  /*ny= y;
  nx= x;
  nz= z;
  nslice = nx*ny;
  ntot=(long)(nslice*nz);*/

  header_size = hs;
  nyR= y;
  nxR= x;
  nzR= z;
  nsliceR = nxR*nyR;
  ntotR=(long)(nsliceR*nzR);
  
  filename = flname;


  //NOTE: allow for out of memory errors!
  //m =(input_vol_type ***)malloc(nx*sizeof(input_vol_type**));
  Rm =(input_vol_type ***)malloc(nxR*sizeof(input_vol_type**));

  // printf("test 1 ok\n");

  for(i=0;i<nxR;i++){
    
    //m[i] = (input_vol_type **)malloc(ny*sizeof(input_vol_type*));
    Rm[i] = (input_vol_type **)malloc(nyR*sizeof(input_vol_type*));

  }

  // printf("test 2 ok\n");

  for (i=0;i<nxR;i++) {

    for (j=0;j<nyR;j++) {

      //m[i][j] = (input_vol_type *)malloc(nz*sizeof(input_vol_type));
      Rm[i][j] = (input_vol_type *)malloc(nzR*sizeof(input_vol_type));

    }

  }

  // printf("test 3 ok\n");
  
  //old unidimensional access implrementation
  //m = (input_vol_type *)malloc(ntot*sizeof(input_vol_type));
  //if (m==NULL){
  //  printf("Not enough memory to hold volume\n");
  //  exit(1);
  //}else{
  //  setzero();
  //}

}



/**
 *
 *Given the coordinates of a bounding box, the source and target volumes, builds a binary volume 
 * containing the union of the binary source and target volumes.
 *
 *
 */

void volume::init2(int resol, int x, int y, int z, int x1, int x2, int y1, int y2, int z1, int z2, volume &source, volume &target, char *flname){

  int i,j,k;
	
  data_type trsh = 0.25; //in terms of the max of the volume
  input_vol_type mas, mat; //max of volume
  data_type csts, cstt;
  int res;

  res = (int)pow(float(2),resol);

  nyR= y/2;
  nxR= x/2;
  nzR= z/2;
  nsliceR = nxR*nyR;
  ntotR=(long)(nsliceR*nzR);

  x1=x1/res;
  x2=x2/res;
  y1=y1/res;
  y2=y2/res;
  z1=z1/res;
  z2=z2/res;

  filename = flname;
  
  mas = source.max_vol();
  csts = mas*trsh;
  mat = target.max_vol();
  cstt = mat*trsh;

  //printf("inside init2(), %f, %f \n", csts, cstt);

  Rm =(input_vol_type ***)malloc(nxR*sizeof(input_vol_type**));

  for(i=0;i<nxR;i++){
    Rm[i] = (input_vol_type **)malloc(nyR*sizeof(input_vol_type*));
  }

  for (i=0;i<nxR;i++) {
    for (j=0;j<nyR;j++) {
      Rm[i][j] = (input_vol_type *)malloc(nzR*sizeof(input_vol_type));
    }
  }

  //set all to zero
  for (k=0;k<nzR;k++){
    for (j=0;j<nyR;j++){
      for (i=0;i<nxR;i++){
	Rm[i][j][k]=0;
      }
    }
  }

 //set voxels within bounding box to one
  for (k=z1;k<z2;k++){
    for (j=y1;j<y2;j++){
      for (i=x1;i<x2;i++){
	Rm[i][j][k]=1;
      }
    }
  }
  
  //set voxels within bounding box and in max of the images to one
  // start from the left first
  for (k=z1;k<nzR;k++){
    for (j=y1;j<nyR;j++){
      for (i=x1;i<nzR;i++){
	
	if ( (source.get(i,j,k)<csts) && (target.get(i,j,k)<cstt) ){
	  Rm[i][j][k]=0;
	}else{
	  break;
	}
	
      }
    }
  }
  
  //from the right now
  
 for (k=z1;k<nzR;k++){
   for (j=y1;j<nyR;j++){
     for (i=x2;i>-1;i--){
       
       if ( (source.get(i,j,k)<csts) && (target.get(i,j,k)<cstt) ){
	 Rm[i][j][k]=0;
       }else{
	 break;
       }
       
     }
   }
 }



}//of init2 


/**
 *
 *Performs binary dilation.
 *Meant to be called imediately after init2 only.
 *
 */

void volume::dilate(void){

  int i,j,k,nnx,nny,nnz;
  int ii,jj,kk;

  nnx = nxR;
  nny = nyR;
  nnz = nzR;

  //allocate memory
  m = (input_vol_type ***)malloc(nnx*sizeof(input_vol_type**));
  if(m==NULL){
    printf("Not enough memory to build binary image.\n");
    exit(1);
  }

  for(i=0;i<nnx;i++){
    m[i] = (input_vol_type **)malloc(nny*sizeof(input_vol_type*));
    //Rm[i] = (input_vol_type **)malloc(nyR*sizeof(input_vol_type*));
    if(m[i]==NULL){
      printf("Not enough memory to build binary image, 2.\n");
      exit(1);
    }
  }

  
  for (i=0;i<nnx;i++) {
    for (j=0;j<nny;j++) {
      m[i][j] = (input_vol_type *)malloc(nnz*sizeof(input_vol_type));
      //Rm[i][j] = (input_vol_type *)malloc(nzR*sizeof(input_vol_type));
      if(m[i][j]==NULL){
	printf("Not enough memory to build binary image., 3\n");
	exit(1);
      }
    }
  }


  //perform dilation
  for (k=0;k<nzR;k++){
    for (j=0;j<nyR;j++){
      for (i=0;i<nxR;i++){
	m[i][j][k]=0;
	if(Rm[i][j][k]>0){
	 
      for(kk=mymax(0,k-2);kk<mymin(nzR-1,k+2);kk++){
        for(jj=mymax(0,j-2);jj<mymin(nyR-1,j+2);jj++){
          for(ii=mymax(0,i-2);ii<mymin(nxR-1,i+2);ii++){
		m[ii][jj][kk]=1;
	      }
	    }
	  }

	}
	
      }
    }
  }


  //copy m into Rm
  for (k=0;k<nzR;k++){
    for (j=0;j<nyR;j++){
      for (i=0;i<nxR;i++){
	Rm[i][j][k]=m[i][j][k];
      }
    }
  }

 
  //dealocate m
 
  for (i=0; i<nnx; i++) {
    for (j=0; j<nny; j++){
      free(m[i][j]);
    }
  }
  
  for (i=0;i<nnx;i++){
    free(m[i]);
  }
  
  free(m);
  
  

}//of dilate


/**
 *
 *Given an integer n by which to subsample, computes the nth lower resolution.
 *n typically goes from high to low.
 *NOTE: The image at full resolution must be divisible by n!
 */

void volume::build_resolution(int n){

  int i,j,k;

  data_type jmp;

  //printf("Subsampling volume.\n");

  jmp = pow(float(2),n);

  nx = (int) (((data_type)nxR)/((data_type)jmp));
  ny = (int) (((data_type)nyR)/((data_type)jmp));
  nz = (int) (((data_type)nzR)/((data_type)jmp));
  //nz = (int) (((data_type)nzR)/((1)));

  nslice = nx*ny;
  ntot=(long)(nslice*nz);

  //clean up first if m already used
  //if(m!=NULL){
  //printf("M was not null!\n");
  //destroy();
  // }

  //allocate necessary memory for mx
  m = (input_vol_type ***)malloc(nx*sizeof(input_vol_type**));
  if(m==NULL){
    printf("Not enough memory to build resolution.\n");
    exit(1);
  }
    
  for(i=0;i<nx;i++){
    m[i] = (input_vol_type **)malloc(ny*sizeof(input_vol_type*));
    //Rm[i] = (input_vol_type **)malloc(nyR*sizeof(input_vol_type*));
    if(m[i]==NULL){
      printf("Not enough memory to build resolution, 2.\n");
      exit(1);
    }
  }

  
  for (i=0;i<nx;i++) {
    for (j=0;j<ny;j++) {
      m[i][j] = (input_vol_type *)malloc(nz*sizeof(input_vol_type));
      //Rm[i][j] = (input_vol_type *)malloc(nzR*sizeof(input_vol_type));
      if(m[i][j]==NULL){
	printf("Not enough memory to build resolution., 3\n");
	exit(1);
      }
    }
  }
  
  setzero();//for safety

  

  //now do subsampling

  if (n==0){
    int coordx,coordy,coordz;
    for (k=0;k<nz;k++){
      for (j=0;j<ny;j++){
	for (i=0;i<nx;i++){
	  coordx = (int)(i*jmp);
	  coordy = (int)(j*jmp);
	  coordz = (int)(k*jmp);
	  set(i,j,k,Rm[coordx][coordy][coordz]);
	}
      }
    }
  }else{

    int coordx,coordy,coordz;
    data_type new_val=0;
    data_type fctr = (jmp+1)*(jmp+1)*(jmp+1);
    int ii,jj,kk;
    
    for (k=0;k<nz;k++){
      for (j=0;j<ny;j++){
	for (i=0;i<nx;i++){
	  coordx = (int)(i*jmp);
	  coordy = (int)(j*jmp);
	  coordz = (int)(k*jmp);
	  new_val = 0;
	  for(ii=(coordx-(int)jmp/2);ii<(coordx+(int)jmp/2+1);ii++){
	    for(jj=(coordy-(int)jmp/2);jj<(coordy+(int)jmp/2+1);jj++){
	      for(kk=(coordz-(int)jmp/2);kk<(coordz+(int)jmp/2+1);kk++){
		new_val = new_val+(data_type)get_Rm(ii,jj,kk);
	      }
	    }
	  }
	  new_val = new_val/fctr; 
	  set(i,j,k,new_val);
	  
	}
      }
    }

  }
  //printf("Done building res.\n");

}


/*
 * interpolates
 */

data_type volume::get_intrpl(float x,float y,float z){

	int x1,x2,y1,y2,z1,z2;
	data_type ur,ul,lr,ll,l,r;

	if ( (x<nx-1) & (y<ny-1) & (z<nz-1) & (x>=0) & (y>=0) & (z>=0) ){
		//if coordinates within image
			
		/*trilinear interpolation (spline of order 1 in 3D)*/
		
		/*determine octant*/
		x1=(int)floor(x);
		x2=(int)ceil(x);
		y1=(int)floor(y);
		y2=(int)ceil(y);
		z1=(int)floor(z);
		z2=(int)ceil(z);

		/*do interpolation as separable b-splines of order 1*/
		ur = (z-z1)*(get(x2,y1,z2)-get(x2,y1,z1))+get(x2,y1,z1);
		ul = (z-z1)*(get(x1,y1,z2)-get(x1,y1,z1))+get(x1,y1,z1);
		lr = (z-z1)*(get(x2,y2,z2)-get(x2,y2,z1))+get(x2,y2,z1);
		ll = (z-z1)*(get(x1,y2,z2)-get(x1,y2,z1))+get(x1,y2,z1);

		r = (y-y1)*(lr-ur) + ur;
		l = (y-y1)*(ll-ul) + ul;

		return ((x-x1)*(r-l)+l);


	}else{
		return 0;
	}
}

/*
 * reads in a volume
 */
int volume::vol_read(void){

  FILE *fp;
  input_vol_type *temp;
  temp=(input_vol_type *)malloc(nxR*sizeof(input_vol_type));
  int i,j,k;
  char buff;
  
  if ((fp=fopen(filename,"r+b"))==NULL) {
    printf("Could not open input file.\n");
    exit(1);
  }else{
    
    //skip header_bytes
    for(i=0;i<header_size;i++){
      fread(&buff,sizeof(char),1,fp);
    }
      

    for (k=0;k<nzR;k++){
      for (j=0;j<nyR;j++){
		  fread(temp,sizeof(input_vol_type),nxR,fp);
	for (i=0;i<nxR;i++){
	  set_Rm(i,j,k,temp[i]);
	  //	printf("inside read ok\n");
	}
      }
    }
    //printf("Input file read.\n");
    fclose(fp);
    return 1;
  }

}

/*
 * writes a volume to a file
 */
void volume::vol_write(void){

	FILE *fp;
	input_vol_type *temp;
	temp=(input_vol_type *)malloc(nxR*sizeof(input_vol_type));
	int i,j,k;

	if ((fp=fopen(filename,"wb"))==NULL) {
		printf("Could not open output file.\n");
		free (temp);
		exit(1);
	}else{
	  

	  //	for (i=0;i<ntot;i++){
	  //  temp = (input_vol_type)m[i];
	  //  fwrite(&temp,sizeof(input_vol_type),1,fp);
	  //}
	  printf("Writing output to :");
	  printf(filename);
	  printf("\n");
	  printf("Image of size %d %d %d\n", nx,ny,nz);

	  for (k=0;k<nz;k++){
	    // printf("%d ",k);
	    for (j=0;j<ny;j++){
	      for (i=0;i<nx;i++){
			  temp[i] = get(i,j,k);
			  //m[i]=(input_vol_type)temp;
	      }
		fwrite(temp,sizeof(input_vol_type),nxR,fp);
	    }
	  }

	  printf("File written.\n");
		fclose(fp);
	
	}
		
	
}

/*
 *computes the boinding box of the volume
 *accorind to a threshold herein specified
 */

void volume::compt_bbox(void){

	int x,y,z;
	
	data_type trsh = data_type(0.15); //in terms of the max of the volume
	input_vol_type ma; //max of volume
	data_type cst;

	lx = nx;
	rx = 0;
	ty = ny;
	by = 0;
	fz = nz;
	bz = 0;

	ma = max_vol();
	cst = ma*trsh;

	cut_trsh = (int)cst;

	for(x=0;x<nx;x++){
		for(y=0;y<ny;y++){
			for(z=0;z<nz;z++){
				if(get(x,y,z)>cst){
			
					if(x<lx)
						lx=x;
					if(y<ty)
						ty=y;
					if(z<fz)
						fz=z;

					if(x>rx)
						rx=x;
					if(y>by)
						by=y;
					if(z>bz)
						bz=z;

				}
			}
		}
	}


	//cout << lx << " " << rx<<"\n";
	//cout << ty << " " << by<<"\n";
	//cout << fz << " " << bz<<"\n";
}

/*
 *computes the maximum value of the volume
 */

input_vol_type volume::max_vol(){

	int x,y,z;
	input_vol_type currmax;
	input_vol_type val;

	currmax = get(0,0,0);
	for(x=0;x<nx;x++){
		for(y=0;y<ny;y++){
			for(z=0;z<nz;z++){
				val = get(x,y,z);
				if(val>currmax)
					currmax = val;
			}
		}
	}

	return currmax;

}

/*
 * Destroy 
 */
void volume::destroy()
{
  // free(m);
 int i,j;
 
 for (i=0; i<nx; i++) {
   for (j=0; j<ny; j++){
     free(m[i][j]);
   }
 }

 for (i=0;i<nx;i++){
   free(m[i]);
 }

 free(m);


}


/*
 * Destroy_Rm 
 */
void volume::destroy_Rm()
{
  // free(m);
 int i,j;
 
 for (i=0; i<nxR; i++) {
   for (j=0; j<nyR; j++){
     free(Rm[i][j]);
   }
 }

 for (i=0;i<nxR;i++){
   free(Rm[i]);
 }

 free(Rm);


}
/*
 * Set all matrix elements to zero.
 */
void volume::setzero()
{
  //long i;

  // for(i=0;i<ntot;i++)
  //	  m[i]=0;
  
  int i,j,k;

  for(i=0;i<nx;i++){
    for(j=0;j<ny;j++){
      for(k=0;k<nz;k++){
	set(i,j,k,(data_type)0);
      }
    }
  }

}


/*
 * returns the value stored in indexes i, j, s
 */
input_vol_type volume::get(int i, int j, int s){

  //long kk;

  //	kk = s*nslice + j*nx + i;
	
  //if ( (i>nx-1) || (j>ny-1) || (s>nz-1))
  //printf("inside volume get!\n");

  return m[i][j][s];

}


/*
 * set the value in indexes i, j to val
 */
void volume::set(int i, int j, int s, data_type val){

  //long kk;

  //	kk = s*nslice + j*nx + i;
  //m[kk]=(input_vol_type)val;
	
  //return 1;

  m[i][j][s]=(input_vol_type)val;

}

/*
 * returns the value stored in indexes i, j, s
 */
input_vol_type volume::get_Rm(int i, int j, int s){

  //long kk;
  
  //	kk = s*nslice + j*nx + i;
  
  //if ( (i>nx-1) || (j>ny-1) || (s>nz-1))
  //printf("inside volume get!\n");
  if ( (i<nxR) && (j<nyR) && (s<nzR) && (i>=0) && (j>=0) && (s>=0) ){
   
    return Rm[i][j][s];
  }else{
    //printf("flag");
    return 0;
  }

}



/*
 * set the value in indexes i, j in Rm to val
 */
void volume::set_Rm(int i, int j, int s, data_type val){

  //long kk;

  //	kk = s*nslice + j*nx + i;
  //m[kk]=(input_vol_type)val;
	
  //return 1;

  Rm[i][j][s]=(input_vol_type)val;

}
