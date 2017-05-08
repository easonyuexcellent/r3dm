/* FILE: splcoord.h
Implements a class to hold a 3D array of sampling coordinates (deformation field).
Gustavo Rohde, fall 2000
*/

#ifndef __SAMPLING_COORD_H__
#define __SAMPLING_COORD_H__

#include "register3d.h"
#include <string.h>

class sampling_coord {

private:

  data_type ***X;
  data_type ***Y;
  data_type ***Z;
  data_type ***XX, ***YY, ***ZZ;//temporary arrays used in upsampling
  
  long nslice; //number of pixels in a slice
  long ntot;
  
 public:
  
  int nx;
  int ny;
  int nz;
  
  //initialization
  void init(int,int,int);
  void build_identity(void);
  
  //upsampling
  void upsample(void);
  void upsample_multi(void);
  data_type get_interp_o(data_type, data_type, data_type, char);
  data_type get_o(int, int, int, char);
  
  //data acess
  data_type get(int,int,int,char);
  void set(int,int,int,char, data_type);
  
  //deallocation
  void destroy(void);
  
  //output
  void save(char *fln);
};




#endif
