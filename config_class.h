/*
FILE: config_class.h
Part of register3d.c program.
The purpose of this class is to provide storage for all the registration parameters.
Methods that read these parameters from an ASCII input file and store them in the appropriate variables are also provided.

Gustavo Rohde
*/


#ifndef __CONFIG_CLASS_H__
#define __CONFIG_CLASS_H__

#include "register3d.h"

class config_class {

 public:

  //used for initialization
  void init(char *cfg_file);

  //used for memory dealocation
  void destroy(void);

  //used for reading junk
  char junk[200];

  //these are the image variables
  char flname_src[200];
  int header_src;
  char flname_trg[200];
  int header_trg;
  char out[200];
  
  int nx,ny,nz;

  //these are the parameters
  int num_levels;
  int num_resolutions;
  int max_res;
  int *res_array;
  int *ctrpts_array;
  int *ctrpts_arrayZ;
  int mmax;
  int mmin;
  int nbins_src;
  int nbins_trg;
  data_type s_trsh;

  

  //these are the configuration variables
  int OP_MODE;
  int B_SEG;
  data_type FLEX;
  data_type br_EPS;
  data_type min_EPS;
  data_type br_STEP;
  data_type grd_STEP;
  data_type pick_const;
  int X_opt, Y_opt, Z_opt;
  int SAVE_DEF_FIELD;

};


#endif
