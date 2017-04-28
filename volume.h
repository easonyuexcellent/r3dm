/*
Implements data structure and functions to hold a volume of data.
Uses matrix3d.cpp
Gustavo Rohde
*/

#ifndef __VOLUME_H__
#define __VOLUME_H__

#include "register3d.h"
//#include "matrix3d.h"

class volume {

private:

  int header_size; //in bytes
  long nslice; //number of pixels in a slice
  long ntot;
  
  long nsliceR;
  long ntotR;
  
  input_vol_type ***m; //holds pixel data at current resolution
  
  input_vol_type ***Rm; //holds the pixel data at full resolution
  
  
  
 public:
  
  //constructor and destructor
  //volume(void); 
  //~volume(void);
  
  //initialization
  void init(int,int,int,int,char*);
  void init2(int,int,int,int,int,int,int,int,int,int,volume &, volume &, char *); //initializes a binary volume 
  void setzero(void);//initializes the array to zero
  void destroy(void);//deallocates the memory for m
  void destroy_Rm(void);//deallocates the memory for Rm
  
  //utility routines
  void build_resolution(int res); //resamples the data in Rm 
  void compt_bbox(void);
  input_vol_type max_vol(void);
  void dilate(); //performs binary dilation, meant to be called imediatly after init2(...) only
  
  //Access routines
  input_vol_type get(int, int, int);
  void set(int, int, int, data_type);
  void set_Rm(int, int, int, data_type); //modifies data in Rm
  input_vol_type get_Rm(int, int, int); //access to data in Rm
  data_type get_intrpl(float,float,float);
  
  //read and write routines
  int vol_read(void);
  void vol_write(void);
  int lx,rx,ty,by,fz,bz; //coordinates for the bounding box
  
  //useful variables
  int nx; //number of columns
  int ny; // of rows
  int nz; //of slices
  
  int cut_trsh;
  
  int nxR;
  int nyR;
  int nzR;
  
  char *filename;
};

#endif

