/*FILE: gradient_class.h
 *Implements the necessary class for gradient optimaztion in register3d
 *Contains coordinate of control points, gradient values, coefficients, etc.
 *Gustavo Rohde, fall 2000
 */

#ifndef __GRADIENT_CLASS_H__
#define __GRADIENT_CLASS_H__

#include "register3d.h"

class gradient_class{

private:

  //all these vectors are of size num_points by three,
  //use get(point,coordinate) to access them
  int *control_points; //control point coordinate
  data_type *gradient; //gradient of the cost function with respect to the ith coefficient
  data_type *coefficients;//coefficient values

public:
	
  int minx,maxx,miny,maxy,minz,maxz; //bounding box of current active region
  data_type spcx,spcy,spcz; //the spacing between the knots in each dimmension
  data_type min_a;//this is the maximum value of the constant by which the gradient can be multiplied without infringing the
  //specified threshold for the jacobian
  int nx,ny,nz;//number of control points in the x, y, and z directions
	
  //main initialization function, initializes a regular grid of control points
  void init(int lx,int rx,int ty,int by,int fz,int bz,int nptsx, int nptsy, int nptsz, int sizex, int sizey, int sizez, data_type FLEX);

  //secondary initialization function, initializes a cube around a centre location
  void init2(int x, int y, int z, int sx, int sy, int sz, int sp, int spz, int sizex, int sizey, int sizez);

  void destroy(void);//deallocation procedure
  void reset(void); //zeroes everything
  void set(int i,int j,data_type val,char w); //here i refers to the point, and j to which coordinate
  void coeff_update(data_type konst);//updates the values of the coefficients accordint to the gradient descent formula
  data_type get(int,int,char);//here i refers to the point, and j to which coordinate
  int num_points; //number of control points
  int supp; //support of the basis, NOTE: the support actually is 2*supp+1;
  int suppz; //must be smaller than supp
  
  void normalize(void); //sets the L2 norm of the gradient to 1
  data_type norm(void); //computes the L2 norm of the gradient
  void max_d_compute(data_type trsh); //computes the maximum derivative (absolute value); Also determines what maximum coefficient can be
  data_type keep_max(data_type,data_type,data_type);//given three values, return the maximum of these
  data_type keep_min(data_type,data_type,data_type);//given three values, return the minimum of these

};




#endif 
