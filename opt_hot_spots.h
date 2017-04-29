#ifndef MULTITHREAD_HOT_SPOTS_H
#define MULTITHREAD_HOT_SPOTS_H

#include <head.h>
#include "register3d.h"
#include "volume.h"
#include "splcoord.h"
#include "matrix2d.h"
#include "basisfctn.h"
#include "gradient_class.h"
#include "config_class.h"
#include <stdlib.h>
#include <time.h>
#include "multithread_hot_spots.h"
#include "hot_spots.h"
#include <stdio.h>
#include <string.h>
#include <QMetaType>
#include <QThread>

class multithread_hot_spots : public QThread
{
public:
    //functions
    multithread_hot_spots(){endflag=0;runflag=0;}
    void init(hot_spots *, matrix_2d *,volume *, volume *,sampling_coord *,basis_function *,input_vol_type,input_vol_type,int,int );
    void grd(void);
    void run();
    data_type2 line_minimize();
    data_type2 f_eval(data_type konst, matrix_2d &JH_temp, matrix_2d &JH_temp_copy,
             sampling_coord &def_field_temp, int lb, int rb, int up, int lob,
             int fb, int bb);
    void JH_inst_up(matrix_2d &JH_temp,
            sampling_coord &def_field_temp, int lb, int rb, int up, int lob,
            int fb, int bb);
    //data
    data_type grd_STEP;
    data_type min_EPS;
    data_type br_STEP;
    int X_opt, Y_opt, Z_opt;
    data_type s_trsh;
    gradient_class gradient; //the gradient structure
    int endflag;
    int runflag;

private:  char std_msg[256];
    hot_spots *htspts;
    //Algorithm variables
    matrix_2d * JH; //the instantaneous (current) joint histogram
    volume * source; //source volume
    volume * target; //target volume
    sampling_coord * def_field; //the deformation field
    basis_function * RBF; //the radial basis function for current level
    input_vol_type mmax,mmin; //mymax and mymin values to consider in 2d histogram
    int nbins_src, nbins_trg; //number of bins in the histogram

};

#endif // MULTITHREAD_HOT_SPOTS_H
