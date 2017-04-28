/*
Defines a compactly supported radial basis function.
Gustavo Rohde, fall 2000
*/


#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "basisfctn.h"


void basis_function::init(int support,int supportz, int num_init){

	
  int i,j,k;

        //Allocates memory for array, and builds it
	
	//NOTE: suppz must be smaller than support!

  //	supp=support; //actually 1/2 of the real support
  //suppz=supportz;
  //nslice = (2*supp+1)*(2*supp+1);
  //ntot = nslice*(2*supp+1);
  //nx = (2*supp+1);
  //ny = nx;
  //nz = nx;

	if(num_init<1){

	  // printf("first init\n");
	  supp=support; //actually 1/2 of the real support
	  suppz=supportz;
	  nslice = (2*supp+1)*(2*supp+1);
	  ntot = nslice*(2*supp+1);
	  nx = (2*supp+1);
	  ny = nx;
	  nz = nx;

	  //m = (data_type *)malloc(ntot*sizeof(data_type));
	  
	  //NOTE: allow for out of memory errors!
	  m =(data_type ***)malloc(nx*sizeof(data_type**));

	  for(i=0;i<nx;i++){
	    m[i] = (data_type **)malloc(ny*sizeof(data_type*));
	  }

	  for (i=0;i<nx;i++) {
	    for (j=0;j<ny;j++) {
	      m[i][j] = (data_type *)malloc(nz*sizeof(data_type));
	    }
	  }

	}else{

	  destroy();
	  supp=support; //actually 1/2 of the real support
	  suppz=supportz;
	  nslice = (2*supp+1)*(2*supp+1);
	  ntot = nslice*(2*supp+1);
	  nx = (2*supp+1);
	  ny = nx;
	  nz = nx;


	  //m = (data_type *)realloc(m,ntot*sizeof(data_type));
	  
	  //NOTE: allow for out of memory errors!
	  m =(data_type ***)malloc(nx*sizeof(data_type**));

	  for(i=0;i<nx;i++){
	    m[i] = (data_type **)malloc(ny*sizeof(data_type*));
	  }

	  for (i=0;i<nx;i++) {
	    for (j=0;j<ny;j++) {
	      m[i][j] = (data_type *)malloc(nz*sizeof(data_type));
	    }
	  }


	}

	if (m==NULL){
	  printf("Not enough memory for RBF.\n");
	  exit(1);
	}else{
	  build();
	}

}

void basis_function::build(void){
	
	//builds basis function into array "m"

	int x,y,z;
	data_type r; //distance from center of array
	data_type dx,dy,dz,ddz,vz;

	for(z=0;z<(nz);z++){

	  ddz = ( (data_type)abs(z - supp) ) / ( (data_type)suppz );
	  vz = WU2_np(ddz);

		for(y=0;y<(ny);y++){
			for(x=0;x<(nx);x++){

				dx = (data_type)( x - supp );
				dy = (data_type)( y - supp );
				dz = (data_type)( z - supp );

				r = (data_type)sqrt( dx*dx + dy*dy + dz*dz );
				set(x,y,z,vz*(WU2(r)));

			}
		}
	}


}

data_type basis_function::WU2(data_type r){

	//implements Wu's radial positive definite B-spline phi(2,1) widh C^2 smoothness in 3 or less dimensions
	//max hight is 4 at radius = 0;

	data_type temp;
	 
	r = r/((data_type)supp);

	temp = max(1-r,0); 
	return (data_type)(pow(temp,4)*(4 + 16*r + 12*pow(r,2) + 3*pow(r,3)));
	
}

data_type basis_function::WU2_np(data_type r){

	//implements Wu's radial positive definite B-spline phi(2,1) widh C^2 smoothness in 3 or less dimensions
	//max hight is 4 at radius = 0;
	//does not use supp implicitly

	data_type temp;
	 
	//r = r/((data_type)suppz);

	temp = max(1-r,0); 
	return (data_type)(pow(temp,4)*(4 + 16*r + 12*pow(r,2) + 3*pow(r,3)));
	
}


void basis_function::set(int x, int y, int z, data_type val){

  //long coord;

  //	coord = z*nslice + y*(2*supp+1) + x;
	
  //	m[coord] = val;

  //if( (x>nx-1)||(y>ny-1)||(z>nz-1))
  //printf("inside set!\n");
  m[x][y][z] = val;
  

}

data_type basis_function::get(int x, int y, int z){

  //long coord;

  //coord = z*nslice + y*(2*supp+1) + x;
	
  //return m[coord];
  //if ((x>=nx)||(y>=ny)||(z>=nz))
  //  printf("inside Get!\n");

  return m[x][y][z];

}

void basis_function::destroy(void){

 int i,j;

	if(m==NULL){
	}else{
	  //free(m);

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

}


void basis_function::reset(void){

	//if(m==NULL){
	//}else{
	//free(m);
	//}
	supp = 0;
	nslice = 0;
	ntot=0;
	nx = 0;
	ny = 0;
	nz = 0;

}


void basis_function::write(void){

	FILE *fp;
	data_type temp;
	long i;

	//printf("in write\n");
	//printf("%d\n",sizeof(data_type));

	if ((fp=fopen("rbf.vol","wb"))==NULL) {
		printf("Could not open output file.\n");
		exit(1);
	}else{
		for (i=0;i<ntot;i++){
		  //temp = m[i];
		  // fwrite(&temp,sizeof(data_type),1,fp);
		}
		fclose(fp);
	
	}

}
