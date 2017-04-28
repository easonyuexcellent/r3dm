/*FILE: OPTIMIZATION.h
 *class responsible for implementing the 3D optimization algorithm
 *Gustavo Rohde, fall 2000

 This class provides the variable as well as methods to perform spline based nonrigid registration.
 Usage of this class by the outside world should follow:

 -init()
 -run()
 -output()
 -destroy()

 */


#ifndef __OPTIMIZATION_H__
#define __OPTIMIZATION_H__

#include "register3d.h"
#include "volume.h"
#include "splcoord.h"
#include "matrix2d.h"
#include "basisfctn.h"
#include "gradient_class.h"
#include "config_class.h"
#include <stdlib.h>
#include <time.h>

class OPTIMIZATION {

 public:
  
  /**Initialization method. Given a config_class object it will read input images and
     initialize all appropriate variables
  */
  void init(config_class &cfg); //initialization routine
  
  /**Memory cleanup method
   */
  void destroy(void); //desctructor
  
  /**Methods that run the algorithms here supported
   */
  void run(void); //runs the old algorithm
  void run_gus(void); //runs the new one
  void run_multi(void); //runs the multi-thread one
  
  /**Given a filename for output, it will save the result image after registration. If so desired it will
     also save the deformation fields
  */
  void output(char *fn); //outputs the result
  
  
  //these are the configuration variables
  data_type FLEX;
  data_type br_EPS;
  data_type min_EPS;
  data_type br_STEP;
  data_type grd_STEP;
  data_type pick_const;
  int X_opt, Y_opt, Z_opt;
  int SAVE_DEF_FIELD;
  data_type s_trsh;
  
 private:
  
  //Algorithm variables
  matrix_2d JH; //the instantaneous (current) joint histogram
  volume source; //source volume
  volume target; //target volume
  volume bin_vol; //binary volume containing the union of the source and target volumes (not used)
  sampling_coord def_field; //the deformation field
  basis_function RBF; //the radial basis function for current level
  gradient_class gradient; //the gradient structure
  bounding_box bbox; // the bounding box coordinates structure (bounding box of union of source and target objects)
  input_vol_type mmax,mmin; //max and min values to consider in 2d histogram
  int nbins_src, nbins_trg; //number of bins in the histogram 
  int npoints_x; //number of control points in the x dimension
  int npoints_y; //number of control points in the y demension
  int npoints_z; //number of control points in the z dimension
  int num_levels;//total number of levels 
  int curr_level;// the current level;
  int num_rbf_init;//number of times the RBF has been initialized

  int num_resolutions;//number of resolutions to be used
  int curr_resolution;//index indicating current resolution
  int curr_res_index;//yet another index indicating current resolution
  int max_res_index;//maximum resolution index allowed
  int *res_array;//array that contains the sequence of resolutions
  int *ctrpts_array;//array that contains the sequence of control points in the x and y directions
  int *ctrptsZ_array;//array that contains the sequence of control points in the x and y directions
  time_t time_init;//initial time
  time_t time_now;//current time

  /**Method that evaluates the gradient of the costfunction with respect to the current set of basis
     functions
   */
  void grd(void); //given that the gradient has been initialized, computes the gradient

  /**Method that updates the Joint Histogram Structure
   */
  void JH_inst(); //computes the whole 2d Joint histogram for the current def_field

  /**Method that updates the Joint Histogram Structure (no longer used)
   */
  void JH_inst1(); //computes the whole 2d Joint histogram with the current def_field and coefficients

  /**Intended to aid line_minimize, takes in a constant and outputs the mutual information based
     cost function.
  */
  void JH_inst_up(matrix_2d &JH_temp,
		  sampling_coord &def_field_temp, int lb, int rb, int up, int lob,
		  int fb, int bb);

  /**Procedure that performs line monimization along the direction of the negative gradient.
   */
  data_type2 line_minimize(void);

  /**Intended to aid line_minimize, takes in a constant and outputs the mutual information.
   */
  data_type2 f_eval(data_type konst, matrix_2d &JH_temp, matrix_2d &JH_temp_copy,
		   sampling_coord &def_field_temp, int lb, int rb, int up, int lob,
		   int fb, int bb); 

  /**Method that uses the current gradient structure to update and actually change the current deformation
     field.
   */
  void def_field_update(void);
  void def_field_update2(gradient_class *);

  /**Method that applyes the current deformation field to the source image and stores the result in
     the source image itself
  */
  void generate_result(void); 
	
};



#endif
