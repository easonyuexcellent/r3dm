/*type definitions for register3d.c
Gustavo Rohde
*/



#ifndef __REGISTER_3D_H__
#define __REGISTER_3D_H__

typedef struct {
	int lx;
	int rx;
	int ty;
	int by;
	int fz;
	int bz;
} bounding_box;


typedef short int input_vol_type;
typedef float data_type;
typedef double data_type2;
#define min_data_type -32768
#define max_data_type  32767
#define myabs(b) ((b>0)? (b):(-b))
#define max(a,b) ((a>b)? (a):(b))
#define min(a,b) ((a>b)? (b):(a))
#define sign(a) ((a>0)? (1):(0))

/*
data_type const FLEX = (data_type)1.6;
data_type const br_EPS = (data_type)0.0001;
data_type const min_EPS = (data_type)0.005;
data_type const br_STEP = (data_type)0.15;
*/

#endif
