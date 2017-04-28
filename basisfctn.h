/*
Defines a compactly supported radial basis function.
Gustavo Rohde, fall 2000
*/


#ifndef __BASIS_FUNCTION_H__
#define __BASIS_FUNCTION_H__

#include "register3d.h"

class basis_function{

private:

  //	data_type *m;//basis function 3D array
        data_type ***m;
	//	int supp; //support must be an integer 
	//int suppz;
	long ntot; //number of total elements in three dimensional array
	long nslice; //number of elements per slice
	

public:
	int supp; //support must be an integer 
	int suppz;

	int nx;
	int ny;
	int nz;

	//initialization
	void init(int,int,int); 
	void build(void);

	//access routines
	void set(int,int,int,data_type);
	data_type get(int,int,int);

	//radial basis
	data_type WU2(data_type r);
	data_type WU2_np(data_type r);
	
	//deallocation routines
	void reset(void);
	void destroy(void);


	void write(void);
};



#endif
