#ifndef CALLED_GRD_BUILD_H
#define CALLED_GRD_BUILD_H

#include "head.h"
#include "register3d.h"
#include "volume.h"
#include "splcoord.h"
#include "matrix2d.h"
#include "basisfctn.h"
#include "gradient_class.h"
#include "config_class.h"
#include <stdlib.h>
#include <time.h>
#include "hot_spots.h"
#include <stdio.h>
#include <string.h>
#include <QThread>

class called_grd_build : public QThread
{
public:
    void init(gradient_class *_grd, int _left, int _right,gradient_class *_gradient,volume *_source,volume *_target,sampling_coord *_def_field,basis_function *_RBF,matrix_2d *_JH);
    called_grd_build();
    void run();
    int left;
    int right;
    int X_opt, Y_opt, Z_opt;
    int nbins_src, nbins_trg;
    input_vol_type mmax,mmin; //max and min values to consider in 2d histogram
    data_type grd_STEP;
    data_type fpp;

private:
    int i;
    gradient_class *gradient;
    volume *source; //source volume
    volume *target; //target volume
    sampling_coord *def_field; //the deformation field
    basis_function *RBF; //the radial basis function for current level
    matrix_2d *JH; //the instantaneous (current) joint histogram
};

#endif // CALLED_GRD_BUILD_H
