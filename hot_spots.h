/**
 *
 * FILE: hot_spots.h
 * Gustavo Rohde, fall 2000;
 * Class definition for structure which identifies the hot_spots of a regular grid optimizaion and
 *			provides a place to hold them.
 *
 */

#ifndef __HOT_SPOTS_H__
#define __HOT_SPOTS_H__

#include "register3d.h"
#include "gradient_class.h"
#include "volume.h"

class hot_spots {
  
 public:
  
  //number of hot regions
  int n_pts;
  
  //access to hot points
  int get(int i, int j);
  void set(int i, int j, int val);
  
  //based on konst threshold alone
  void prune1(gradient_class &gradient, data_type konst, int num_init, volume &source, volume &target);

  void prune12(gradient_class &gradient, data_type konst, int num_init, volume &source, volume &target);
  
  //destructor
  void destroy(void);
  
  //Contains the number of hot spots found
  int num_hp;
  
 private:
  
  //A 3Xn array containing coordinates of hot spots;
  int *hp;
  
  //find maximum value
  void find_max(gradient_class &grd,int *index, data_type *max);


  void find_max_multi(gradient_class &grd,int *index, data_type *max);
  
  //check to see if point is already included (or close to) in one of the points in hp
  int check_distance(gradient_class &grd, int index);


  int check_distance_multi(gradient_class &grd, int index);
  
  

};


#endif
