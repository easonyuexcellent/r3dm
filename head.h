#ifndef HEAD_H
#define HEAD_H

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
#define core_number 4
#define switch_number 32
#define myabs(b) ((b>0)? (b):(-b))
#define mymax(a,b) ((a>b)? (a):(b))
#define mymin(a,b) ((a>b)? (b):(a))
#define sign(a) ((a>0)? (1):(0))

#endif // HEAD_H
