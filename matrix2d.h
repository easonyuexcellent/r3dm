/*
File: MATRIX2D.H
intended to be used to store a 2d joint histogram
Gustavo Rohde
*/


#ifndef __MATRIX_2D_H__
#define __MATRIX_2D_H__

#include "register3d.h"

typedef long int JHtype;

class matrix_2d {

public:
  //constructor and destructor
  //	matrix_2d(void); 
  //	~matrix_2d(void);
  
  //initialization/destruction routines
  void init(int, int);//allocates memory for the array
  void setzero(void);//initializes the array to zero
  void destroy(void);//deallocates the memory
  
  //Access routines
  JHtype get(int, int);
  void set(int, int, JHtype);//access to values in x, and y
  void increment(int,int);//adds 1 to coordinate x, y
  void subtract(int,int,JHtype);//subtracts 1 from coordinate x, y
  void add(int,int,JHtype);//adds JHtype to the current coordinate value
  data_type2 mi(void);//computes the mutual information for current array
  JHtype sum(void);//computes the total of the array
  //void trg_en_cpt(void);//computes the entropy of the target image
  int nx; //number of histogram bins in source image
  int ny; //number of histogram bins in target image
  
 private:
  
  JHtype ntot; //total number of elements in Joint histogram matrix
  data_type nl2; //stores natural log of 2, fbor speed
  JHtype **m; //stores the actual Joinh Histogram data
  JHtype *tmp_src_hist; //array to contain the temporary histogram of the current source image
  JHtype *tmp_trg_hist; //array to contain histogram of the target image
  //data_type En_Trg; //Entropy of target image, needs to be computed only once.
} ;

#endif /* __MATRIX_2D_H__ */

