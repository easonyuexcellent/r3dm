/* FILE: splcoord.cpp
Implements a class to hold a 3D array of sampling coordinates (deformation field).
Gustavo Rohde, fall 2000
*/


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "register3d.h"
#include "splcoord.h"
#include "called_crd_upsp.h"

void sampling_coord::init(int x, int y, int z){

    //WARNING: calls build_identity automatically
    long ntot;
    int i,j;//,k;

    ny= y;
    nx= x;
    nz= z;

    nslice = x*y;
    ntot=nslice*z;

    //X = (data_type *)malloc(ntot*sizeof(data_type));
    X = (data_type ***)malloc(nx*sizeof(data_type**));
    for(i=0;i<nx;i++){
        X[i] = (data_type **)malloc(ny*sizeof(data_type*));
    }

    for (i=0;i<nx;i++) {
        for (j=0;j<ny;j++) {
            X[i][j] = (data_type *)malloc(nz*sizeof(data_type));
        }
    }

    if (X==NULL){
        printf("Not enough memory for deformation field\n");
        exit(1);
    }

    // Y = (data_type *)malloc(ntot*sizeof(data_type));

    Y = (data_type ***)malloc(nx*sizeof(data_type**));
    for(i=0;i<nx;i++){
        Y[i] = (data_type **)malloc(ny*sizeof(data_type*));
    }

    for (i=0;i<nx;i++) {
        for (j=0;j<ny;j++) {
            Y[i][j] = (data_type *)malloc(nz*sizeof(data_type));
        }
    }
    if (Y==NULL){
        printf("Not enough memory for deformation field\n");
        exit(1);
    }


    // Z = (data_type *)malloc(ntot*sizeof(data_type));
    Z = (data_type ***)malloc(nx*sizeof(data_type**));
    for(i=0;i<nx;i++){
        Z[i] = (data_type **)malloc(ny*sizeof(data_type*));
    }

    for (i=0;i<nx;i++) {
        for (j=0;j<ny;j++) {
            Z[i][j] = (data_type *)malloc(nz*sizeof(data_type));
        }
    }
    if (Z==NULL){
        printf("Not enough memory for deformation field\n");
        exit(1);
    }

    //printf("Def_field allocated\n");
    build_identity();


}

/**
 *
 * Upsample the three arrays of the deformation field by using interpolation
 *
 *
 */

void sampling_coord::upsample(){

    int nx_new, ny_new, nz_new; //new dimensions
    int i,j,k;

    //printf("Upsampling deformation field.\n");

    nx_new = 2*nx;
    ny_new = 2*ny;
    nz_new = 2*nz;

    //allocate arrays to hold X,Y,Z temporarily
    // data_type ***XX, ***YY, ***ZZ;

    XX = (data_type ***)malloc(nx*sizeof(data_type**));
    YY = (data_type ***)malloc(nx*sizeof(data_type**));
    ZZ = (data_type ***)malloc(nx*sizeof(data_type**));
    for(i=0;i<nx;i++){
        XX[i] = (data_type **)malloc(ny*sizeof(data_type*));
        YY[i] = (data_type **)malloc(ny*sizeof(data_type*));
        ZZ[i] = (data_type **)malloc(ny*sizeof(data_type*));
    }

    for (i=0;i<nx;i++) {
        for (j=0;j<ny;j++) {
            XX[i][j] = (data_type *)malloc(nz*sizeof(data_type));
            YY[i][j] = (data_type *)malloc(nz*sizeof(data_type));
            ZZ[i][j] = (data_type *)malloc(nz*sizeof(data_type));
        }
    }

    if ((XX==NULL)||(YY==NULL)||(ZZ==NULL)){
        printf("Not enough memory for deformation field copy in sampling_coord::upsample.\n");
        exit(1);
    }

    //printf("Allocated memory ok.\n");

    //copy existing X, Y, Z into XX, YY, and ZZ arrays
    for(i=0;i<nx;i++){
        for(j=0;j<ny;j++){
            for(k=0;k<nz;k++){
                XX[i][j][k]=X[i][j][k];
                YY[i][j][k]=Y[i][j][k];
                ZZ[i][j][k]=Z[i][j][k];
            }
        }
    }

    //printf("Copied fine.\n");

    //now free X, Y, and Z
    destroy();

    //printf("Destroyed old X,Y,Z ok.\n");

    //allocate X,Y, and Z of new sizes
    X = (data_type ***)malloc(nx_new*sizeof(data_type**));
    Y = (data_type ***)malloc(nx_new*sizeof(data_type**));
    Z = (data_type ***)malloc(nx_new*sizeof(data_type**));

    for(i=0;i<nx_new;i++){
        X[i] = (data_type **)malloc(ny_new*sizeof(data_type*));
        Y[i] = (data_type **)malloc(ny_new*sizeof(data_type*));
        Z[i] = (data_type **)malloc(ny_new*sizeof(data_type*));
    }

    for (i=0;i<nx_new;i++) {
        for (j=0;j<ny_new;j++) {
            X[i][j] = (data_type *)malloc(nz_new*sizeof(data_type));
            Y[i][j] = (data_type *)malloc(nz_new*sizeof(data_type));
            Z[i][j] = (data_type *)malloc(nz_new*sizeof(data_type));
        }
    }

    if ((X==NULL)||(Y==NULL)||(Z==NULL)){
        printf("Not enough memory for deformation field copy in sampling_coord::upsample.\n");
        exit(1);
    }

    //printf("Allocated new X,Y,Z fine.\n");

    //now populate new X,Y, and Z
    data_type coordx,coordy,coordz;

    for(i=0;i<nx_new;i++){
        for(j=0;j<ny_new;j++){
            for(k=0;k<nz_new;k++){

                coordx = ((data_type)i)/2.0;
                coordy = ((data_type)j)/2.0;
                coordz = ((data_type)k)/2.0;

                X[i][j][k] = 2.0*get_interp_o(coordx,coordy,coordz,'x');
                Y[i][j][k] = 2.0*get_interp_o(coordx,coordy,coordz,'y');
                Z[i][j][k] = 2.0*get_interp_o(coordx,coordy,coordz,'z');

            }
        }
    }

    //printf("Upsampled ok.\n");
    //now clean up XX,YY,ZZ
    for (i=0; i<nx; i++) {
        for (j=0; j<ny; j++){
            free(XX[i][j]);
            free(YY[i][j]);
            free(ZZ[i][j]);
        }
    }

    for (i=0;i<nx;i++){
        free(XX[i]);
        free(YY[i]);
        free(ZZ[i]);
    }

    free(XX);
    free(YY);
    free(ZZ);

    //printf("Freed memory ok.\n");
    //update fields nx, ny, nz
    nx = nx_new;
    ny = ny_new;
    nz = nz_new;

    //now correct for errors in the get_interp_o function

    //X
    /* for (k=0;k<nz;k++){
    for (j=0;j<ny;j++){
      if ( X[nx-1][j][k] <= X[nx-2][j][k] )
    //X[nx-1][j][k] = X[nx-2][j][k]+1;
    nx-1;
    }
    }*/

    //Y
    /*
  for (k=0;k<nz;k++){
    for (i=0;i<nx;i++){
      if ( Y[i][ny-1][k] <= Y[i][ny-2][k] )
    //Y[i][ny-1][k] = Y[i][ny-2][k]+1;
    ny-1;
    }
    }*/

    //Z
    /*
  for (i=0;i<nx;i++){
    for (j=0;j<ny;j++){
      if( Z[i][j][nz-1] <= Z[i][j][nz-2] )
    //Z[i][j][nz-1] = Z[i][j][nz-2] +1;
    nz-1;
    }
    }*/

    //printf("Done Resampling.\n");

}//of upsample


void sampling_coord::upsample_multi(){

    int nx_new, ny_new, nz_new; //new dimensions
    int i,j,k;

    //printf("Upsampling deformation field.\n");

    nx_new = 2*nx;
    ny_new = 2*ny;
    nz_new = 2*nz;

    //allocate arrays to hold X,Y,Z temporarily
    // data_type ***XX, ***YY, ***ZZ;

    XX = (data_type ***)malloc(nx*sizeof(data_type**));
    YY = (data_type ***)malloc(nx*sizeof(data_type**));
    ZZ = (data_type ***)malloc(nx*sizeof(data_type**));
    for(i=0;i<nx;i++){
        XX[i] = (data_type **)malloc(ny*sizeof(data_type*));
        YY[i] = (data_type **)malloc(ny*sizeof(data_type*));
        ZZ[i] = (data_type **)malloc(ny*sizeof(data_type*));
    }

    for (i=0;i<nx;i++) {
        for (j=0;j<ny;j++) {
            XX[i][j] = (data_type *)malloc(nz*sizeof(data_type));
            YY[i][j] = (data_type *)malloc(nz*sizeof(data_type));
            ZZ[i][j] = (data_type *)malloc(nz*sizeof(data_type));
        }
    }

    if ((XX==NULL)||(YY==NULL)||(ZZ==NULL)){
        printf("Not enough memory for deformation field copy in sampling_coord::upsample.\n");
        exit(1);
    }

    //printf("Allocated memory ok.\n");

    //copy existing X, Y, Z into XX, YY, and ZZ arrays
    for(i=0;i<nx;i++){
        for(j=0;j<ny;j++){
            for(k=0;k<nz;k++){
                XX[i][j][k]=X[i][j][k];
                YY[i][j][k]=Y[i][j][k];
                ZZ[i][j][k]=Z[i][j][k];
            }
        }
    }

    //printf("Copied fine.\n");

    //now free X, Y, and Z
    destroy();

    //printf("Destroyed old X,Y,Z ok.\n");

    //allocate X,Y, and Z of new sizes
    X = (data_type ***)malloc(nx_new*sizeof(data_type**));
    Y = (data_type ***)malloc(nx_new*sizeof(data_type**));
    Z = (data_type ***)malloc(nx_new*sizeof(data_type**));

    for(i=0;i<nx_new;i++){
        X[i] = (data_type **)malloc(ny_new*sizeof(data_type*));
        Y[i] = (data_type **)malloc(ny_new*sizeof(data_type*));
        Z[i] = (data_type **)malloc(ny_new*sizeof(data_type*));
    }

    for (i=0;i<nx_new;i++) {
        for (j=0;j<ny_new;j++) {
            X[i][j] = (data_type *)malloc(nz_new*sizeof(data_type));
            Y[i][j] = (data_type *)malloc(nz_new*sizeof(data_type));
            Z[i][j] = (data_type *)malloc(nz_new*sizeof(data_type));
        }
    }

    if ((X==NULL)||(Y==NULL)||(Z==NULL)){
        printf("Not enough memory for deformation field copy in sampling_coord::upsample.\n");
        exit(1);
    }

    //printf("Allocated new X,Y,Z fine.\n");

    //now populate new X,Y, and Z
    data_type coordx,coordy,coordz;

    called_crd_upsp **p;
    p=(called_crd_upsp **)malloc(nx_new*sizeof(called_crd_upsp *));
    for(i=0;i<nx_new;i++){
        p[i]=new(called_crd_upsp);
        p[i]->set(this,i);
        p[i]->start();
        if ((i+1)%core_number==0)
            p[i]->wait();
    }
    for(i=0;i<nx_new;i++){
        p[i]->wait();
    }
    for(i=0;i<nx_new;i++){
        delete(p[i]);
    }

    //printf("Upsampled ok.\n");
    //now clean up XX,YY,ZZ
    for (i=0; i<nx; i++) {
        for (j=0; j<ny; j++){
            free(XX[i][j]);
            free(YY[i][j]);
            free(ZZ[i][j]);
        }
    }

    for (i=0;i<nx;i++){
        free(XX[i]);
        free(YY[i]);
        free(ZZ[i]);
    }

    free(XX);
    free(YY);
    free(ZZ);

    //printf("Freed memory ok.\n");
    //update fields nx, ny, nz
    nx = nx_new;
    ny = ny_new;
    nz = nz_new;

    //now correct for errors in the get_interp_o function

    //X
    /* for (k=0;k<nz;k++){
    for (j=0;j<ny;j++){
      if ( X[nx-1][j][k] <= X[nx-2][j][k] )
    //X[nx-1][j][k] = X[nx-2][j][k]+1;
    nx-1;
    }
    }*/

    //Y
    /*
  for (k=0;k<nz;k++){
    for (i=0;i<nx;i++){
      if ( Y[i][ny-1][k] <= Y[i][ny-2][k] )
    //Y[i][ny-1][k] = Y[i][ny-2][k]+1;
    ny-1;
    }
    }*/

    //Z
    /*
  for (i=0;i<nx;i++){
    for (j=0;j<ny;j++){
      if( Z[i][j][nz-1] <= Z[i][j][nz-2] )
    //Z[i][j][nz-1] = Z[i][j][nz-2] +1;
    nz-1;
    }
    }*/

    //printf("Done Resampling.\n");

}//of upsample

/**
 *
 * Meant to be used within upsample() function only;
 * Given a coordinate, not necessarily on a grid, it gives the interpolated value
 *
 */

data_type sampling_coord::get_interp_o(data_type x, data_type y, data_type z, char w){


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
        ur = (z-z1)*(get_o(x2,y1,z2,w)-get_o(x2,y1,z1,w))+get_o(x2,y1,z1,w);
        ul = (z-z1)*(get_o(x1,y1,z2,w)-get_o(x1,y1,z1,w))+get_o(x1,y1,z1,w);
        lr = (z-z1)*(get_o(x2,y2,z2,w)-get_o(x2,y2,z1,w))+get_o(x2,y2,z1,w);
        ll = (z-z1)*(get_o(x1,y2,z2,w)-get_o(x1,y2,z1,w))+get_o(x1,y2,z1,w);

        r = (y-y1)*(lr-ur) + ur;
        l = (y-y1)*(ll-ul) + ul;

        return ((x-x1)*(r-l)+l);


    }else{
        //return 0;

        if (x>nx-1){
            x1 = (int)floor(x);
            x2 = x1;
        }else{
            x1=(int)floor(x);
            x2=(int)ceil(x);
        }

        if (y>ny-1){
            y1 = (int)floor(y);
            y2 = y1;
        }else{
            y1=(int)floor(y);
            y2=(int)ceil(y);
        }

        if (z>nz-1){
            z1 = (int)floor(z);
            z2 = z1;
        }else{
            z1=(int)floor(z);
            z2=(int)ceil(z);
        }



        /*do interpolation as separable b-splines of order 1*/
        ur = (z-z1)*(get_o(x2,y1,z2,w)-get_o(x2,y1,z1,w))+get_o(x2,y1,z1,w);
        ul = (z-z1)*(get_o(x1,y1,z2,w)-get_o(x1,y1,z1,w))+get_o(x1,y1,z1,w);
        lr = (z-z1)*(get_o(x2,y2,z2,w)-get_o(x2,y2,z1,w))+get_o(x2,y2,z1,w);
        ll = (z-z1)*(get_o(x1,y2,z2,w)-get_o(x1,y2,z1,w))+get_o(x1,y2,z1,w);

        r = (y-y1)*(lr-ur) + ur;
        l = (y-y1)*(ll-ul) + ul;

        return ((x-x1)*(r-l)+l);

        //if (x>nx-1){


        //return 0;



    }


}//of get_interp_o


/**
 *
 * Used to save the final deformation field.
 * Saves three volumes, each with the displacement in each dimension.
 * Adds suffix dx, dy, and dz to the filename
 *
 */

void sampling_coord::save(char *fln){

    int pp;
    FILE *fp;
    data_type temp;
    int i,j,k;
    char sdx[200];
    char sdy[200];
    char sdz[200];

    pp = strcspn(fln,".");

    strncpy(sdx,fln,pp);
    strncpy(&sdx[pp],"\0",2);

    strncpy(sdy,fln,pp);
    strncpy(&sdy[pp],"\0",2);

    strncpy(sdz,fln,pp);
    strncpy(&sdz[pp],"\0",2);

    strncat(sdx,".dx",strlen(".dx"));
    strncat(sdy,".dy",strlen(".dy"));
    strncat(sdz,".dz",strlen(".dz"));

    //save x deformation field
    if ((fp=fopen(sdx,"wb"))==NULL) {
        printf("Could not open output file for x deformation field.\n");
        exit(1);
    }else{

        for (k=0;k<nz;k++){
            for (j=0;j<ny;j++){
                for (i=0;i<nx;i++){
                    temp = get(i,j,k,'x');
                    fwrite(&temp,sizeof(data_type),1,fp);
                }
            }
        }

        fclose(fp);

    }

    //save y deformatino field
    if ((fp=fopen(sdy,"wb"))==NULL) {
        printf("Could not open output file for y deformatino field.\n");
        exit(1);
    }else{

        for (k=0;k<nz;k++){
            for (j=0;j<ny;j++){
                for (i=0;i<nx;i++){
                    temp = get(i,j,k,'y');
                    fwrite(&temp,sizeof(data_type),1,fp);
                }
            }
        }

        fclose(fp);

    }

    //save z deformation field
    if ((fp=fopen(sdz,"wb"))==NULL) {
        printf("Could not open output file for z deformation field.\n");
        exit(1);
    }else{

        for (k=0;k<nz;k++){
            for (j=0;j<ny;j++){
                for (i=0;i<nx;i++){
                    temp = get(i,j,k,'z');
                    fwrite(&temp,sizeof(data_type),1,fp);
                }
            }
        }

        fclose(fp);

    }



    //printf("%d\n",strlen(s2));
    //printf(s2);
    //printf("\n");


}


/**
 * Used by get_interp_o only
 * Returns the value at x,y,z of XX, or YY, or ZZ, depending on character w.
 *
 */

data_type sampling_coord::get_o(int x, int y, int z, char w){

    data_type result;

    switch (w){
    
    case 'x':

        result=XX[x][y][z];
        break;

    case 'y':

        result=YY[x][y][z];
        break;

    case 'z':

        result=ZZ[x][y][z];
        break;

    }

    return result;

}



void sampling_coord::set(int x, int y, int z, char w, data_type val){

    //long coord;

    //coord = z*nslice + y*nx + x;

    switch (w){
    
    case 'x':
        //X[coord]=val;
        X[x][y][z]=val;
        break;
    case 'y':
        //Y[coord]=val;
        Y[x][y][z]=val;
        break;
    case 'z':
        //Z[coord]=val;
        Z[x][y][z]=val;
        break;

    }

}


data_type sampling_coord::get(int x, int y, int z, char w){

    //long coord;
    data_type result;

    //coord = z*nslice + y*nx + x;

    switch (w){

    case 'x':
        //result=X[coord];
        result=X[x][y][z];
        break;
    case 'y':
        //result=Y[coord];
        result=Y[x][y][z];
        break;
    case 'z':
        //result=Z[coord];
        result=Z[x][y][z];
        break;

    }

    return result;
}

void sampling_coord::build_identity(void){

    int x,y,z;

    //X
    for(z=0;z<nz;z++){
        for(y=0;y<ny;y++){
            for(x=0;x<nx;x++)
                set(x,y,z,'x',(data_type)x);
        }
    }

    //Y
    for(z=0;z<nz;z++){
        for(x=0;x<nx;x++){
            for(y=0;y<ny;y++)
                set(x,y,z,'y',(data_type)y);
        }
    }

    //Z
    for(x=0;x<nx;x++){
        for(y=0;y<ny;y++){
            for(z=0;z<nz;z++)
                set(x,y,z,'z',(data_type)z);
        }
    }


}

void sampling_coord::destroy(void){

    //free(X);
    //free(Y);
    //free(Z);


    int i,j;

    for (i=0; i<nx; i++) {
        for (j=0; j<ny; j++){
            free(X[i][j]);
            free(Y[i][j]);
            free(Z[i][j]);
        }
    }

    for (i=0;i<nx;i++){
        free(X[i]);
        free(Y[i]);
        free(Z[i]);
    }

    free(X);
    free(Y);
    free(Z);

    //nx=0;
    //ny=0;
    //nz=0;
    //nslice=0;
    //ntot=0;

}
