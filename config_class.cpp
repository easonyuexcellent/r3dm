/*
FILE: config_class.c
Part of the register3d.c program.
Gustavo Rohde, winter 2001

The purpose of this class is to provide storage for all the registration parameters.
Methods that read these parameters from an ASCII input file and store them in the appropriate variables are also provided.

*/

#include "config_class.h"
#include "register3d.h"
#include <stdio.h>
#include <stdlib.h>

//memory dealocation
void config_class::init(char *cfg_file){

  FILE *fp;
  
  if(( fp=fopen(cfg_file,"r")) == NULL ){
    printf("Cannot open input configuration file.\n");
    exit(1);
  }
  
  fscanf(fp,"%s", junk);
  fscanf(fp,"%s", flname_src);
  fscanf(fp,"%d", &header_src);
  fscanf(fp,"%s", junk);
  fscanf(fp,"%s", flname_trg);
  fscanf(fp,"%d", &header_trg);

  fscanf(fp,"%d", &nx);
  fscanf(fp,"%d", &ny);
  fscanf(fp,"%d", &nz);

  fscanf(fp,"%s", junk);
  fscanf(fp,"%d", &num_levels);
  fscanf(fp,"%d", &num_resolutions);
  fscanf(fp,"%d", &max_res);
  ctrpts_array=new int[num_levels];
  ctrpts_arrayZ=new int[num_levels];
  res_array = new int[num_resolutions];

  int i;
  for(i=0;i<num_levels;i++){
    fscanf(fp,"%d", &ctrpts_array[i]);
  }


  for(i=0;i<num_levels;i++){
    fscanf(fp,"%d", &ctrpts_arrayZ[i]);
  }


  for(i=0;i<max_res;i++){
    fscanf(fp,"%d", &res_array[i]);
  }
  
  fscanf(fp,"%d", &mmax);
  fscanf(fp,"%d", &mmin);
  fscanf(fp,"%d", &nbins_src);
  fscanf(fp,"%d", &nbins_trg);
  fscanf(fp,"%f", &s_trsh);
  //printf("Threshold read %f\n",s_trsh);

  fscanf(fp,"%s", junk);
  fscanf(fp,"%d", &OP_MODE);
  fscanf(fp,"%d", &B_SEG);
  fscanf(fp, "%f", &FLEX);
  fscanf(fp, "%f", &br_EPS);
  fscanf(fp, "%f", &min_EPS);
  fscanf(fp, "%f", &br_STEP);
  fscanf(fp, "%f", &grd_STEP);
  fscanf(fp, "%f", &pick_const);
  fscanf(fp, "%d", &X_opt);
  fscanf(fp, "%d", &Y_opt);
  fscanf(fp, "%d", &Z_opt);
  fscanf(fp, "%d", &SAVE_DEF_FIELD);

  fscanf(fp,"%s", junk);
  fscanf(fp,"%s", out);
  //printf("here5: %d \n", nx);
  fclose(fp);
  //printf("here6: %d \n", nx);
}


//memory dealocation
void config_class::destroy(void){

  delete [] ctrpts_array;
  delete [] ctrpts_arrayZ;
  delete [] res_array;

}
