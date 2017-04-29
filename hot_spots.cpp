/**
 *
 * FILE: hot_spots.cpp
 * Gustavo Rohde, fall 2000
 *
 * Implements structure that identifies the hot regions of a regular grid set of
 *		radial basis function coefficients. 
 *
 */


#include "hot_spots.h"
#include "register3d.h"
#include <math.h>
#include <stdlib.h>
#include "called_grd_max.h"
#include "called_grd_diff.h"

/**
 *
 * The constructor for the class.
 * Its tasks are to, upon recieving a gradient and a threshold constant identify the hot spots of
 *			the gradient. 
 *
 * Does so by using konst as a hard (explicit) threshold
 *
 */


void hot_spots::prune1(gradient_class &gradient, data_type konst, int num_init, volume &source, volume &target){
  
  int i;
  int index;
  int hp_in=0;
  data_type coeff_mag;
  //int res;
  
  num_hp = 0;	
  
  //    res = (int)pow(2,resol);
  
  //allocate memory to hold the hot spots (NOTE: the right way to do this would be with lists)
  if (num_init<1){
    hp = (int *)malloc(gradient.num_points*3*sizeof(int));
  } else {
    hp = (int *)realloc(hp,gradient.num_points*3*sizeof(int));
  }
  find_max(gradient,&index,&coeff_mag);	
  //	printf("inside prune, mag is %f, and index is %d \n", coeff_mag, index);
  //	printf("that is right 2\n");
  //initialize to zero
  for (long ii=0;ii<gradient.num_points*3;ii++)
    hp[ii]=0;
  
  //go through each point and compare it to konst
  coeff_mag = data_type(1.05432);
  while (coeff_mag>0){
    
    //find maximum magnitude of coefficients in gradient class
    find_max(gradient,&index,&coeff_mag);
    
    //	printf("inside prune, mag is %f, and index is %d \n", coeff_mag, index);
    //	printf("that is right \n");
    
    if ((coeff_mag>=konst)&&(source.get(gradient.get(index,0,'p'),gradient.get(index,1,'p'),gradient.get(index,2,'p'))>0.5*source.cut_trsh)&&(target.get(gradient.get(index,0,'p'),gradient.get(index,1,'p'),gradient.get(index,2,'p'))>0.5*target.cut_trsh)){
      //if (coeff_mag>=konst){
      //check to see if it doesn't neighbor any other already in hp
      if (check_distance(gradient,index)==0) {
	//add point to hp
	set(num_hp,0,(int)gradient.get(index,0,'p'));
	set(num_hp,1,(int)gradient.get(index,1,'p'));
	set(num_hp,2,(int)gradient.get(index,2,'p'));
	
	num_hp++;	
	
      }
      
      //set to zero
      gradient.set(index,0,0,'c');
      gradient.set(index,1,0,'c');
      gradient.set(index,2,0,'c');
      
      
      
    } else {
      //set all to zero
      gradient.set(index,0,0,'c');
      gradient.set(index,1,0,'c');
      gradient.set(index,2,0,'c');
    }
    
  }
  
  //printf("\n");
  //printf("Inside prune1\n");
  //for (i=0;i<num_hp;i++){
  // printf("      %d, %d, %d \n", (int)get(i,0), (int)get(i,1), (int)get(i,2) );
  //}
  //printf("\n");
  
}



void hot_spots::prune12(gradient_class &gradient, data_type konst, int num_init, volume &source, volume &target){

  int i;
  int index;
  int hp_in=0;
  data_type coeff_mag;
  //int res;

  num_hp = 0;

  //    res = (int)pow(2,resol);

  //allocate memory to hold the hot spots (NOTE: the right way to do this would be with lists)
  if (num_init<1){
    hp = (int *)malloc(gradient.num_points*3*sizeof(int));
  } else {
    hp = (int *)realloc(hp,gradient.num_points*3*sizeof(int));
  }
  find_max_multi(gradient,&index,&coeff_mag);
  /*
  if (gradient.num_points>switch_number){
    find_max_multi(gradient,&index,&coeff_mag);
  }else{
    find_max(gradient,&index,&coeff_mag);
  }*/
  //	printf("inside prune, mag is %f, and index is %d \n", coeff_mag, index);
  //	printf("that is right 2\n");
  //initialize to zero
  for (long ii=0;ii<gradient.num_points*3;ii++)
    hp[ii]=0;

  //go through each point and compare it to konst
  coeff_mag = data_type(1.05432);
  while (coeff_mag>0){

    //find maximum magnitude of coefficients in gradient class
	find_max_multi(gradient,&index,&coeff_mag);
    //if (gradient.num_points>switch_number){
    //  find_max_multi(gradient,&index,&coeff_mag);
    //}else{
    //  find_max(gradient,&index,&coeff_mag);
    //}

    //	printf("inside prune, mag is %f, and index is %d \n", coeff_mag, index);
    //	printf("that is right \n");

	if ((coeff_mag>=konst)&&(source.get(gradient.get(index,0,'p'),gradient.get(index,1,'p'),gradient.get(index,2,'p'))>0.5*source.cut_trsh)&&(target.get(gradient.get(index,0,'p'),gradient.get(index,1,'p'),gradient.get(index,2,'p'))>0.5*target.cut_trsh)){
		//if (coeff_mag>=konst){
		//check to see if it doesn't neighbor any other already in hp
		if (check_distance_multi(gradient,index)==0) {
			//add point to hp
			set(num_hp,0,(int)gradient.get(index,0,'p'));
			set(num_hp,1,(int)gradient.get(index,1,'p'));
			set(num_hp,2,(int)gradient.get(index,2,'p'));
			
			num_hp++;

		}

		//set to zero
		gradient.set(index,0,0,'c');
		gradient.set(index,1,0,'c');
		gradient.set(index,2,0,'c');
		
		

	} else {
		//set all to zero
		gradient.set(index,0,0,'c');
		gradient.set(index,1,0,'c');
		gradient.set(index,2,0,'c');
	}
/*
    if (gradient.num_points>switch_number){
        if ((coeff_mag>=konst)&&(source.get(gradient.get(index,0,'p'),gradient.get(index,1,'p'),gradient.get(index,2,'p'))>0.5*source.cut_trsh)&&(target.get(gradient.get(index,0,'p'),gradient.get(index,1,'p'),gradient.get(index,2,'p'))>0.5*target.cut_trsh)){
          //if (coeff_mag>=konst){
          //check to see if it doesn't neighbor any other already in hp
          if (check_distance_multi(gradient,index)==0) {
        //add point to hp
        set(num_hp,0,(int)gradient.get(index,0,'p'));
        set(num_hp,1,(int)gradient.get(index,1,'p'));
        set(num_hp,2,(int)gradient.get(index,2,'p'));

        num_hp++;

          }

          //set to zero
          gradient.set(index,0,0,'c');
          gradient.set(index,1,0,'c');
          gradient.set(index,2,0,'c');



        } else {
          //set all to zero
          gradient.set(index,0,0,'c');
          gradient.set(index,1,0,'c');
          gradient.set(index,2,0,'c');
        }
    }else{
        if ((coeff_mag>=konst)&&(source.get(gradient.get(index,0,'p'),gradient.get(index,1,'p'),gradient.get(index,2,'p'))>0.5*source.cut_trsh)&&(target.get(gradient.get(index,0,'p'),gradient.get(index,1,'p'),gradient.get(index,2,'p'))>0.5*target.cut_trsh)){
          //if (coeff_mag>=konst){
          //check to see if it doesn't neighbor any other already in hp
          if (check_distance(gradient,index)==0) {
        //add point to hp
        set(num_hp,0,(int)gradient.get(index,0,'p'));
        set(num_hp,1,(int)gradient.get(index,1,'p'));
        set(num_hp,2,(int)gradient.get(index,2,'p'));

        num_hp++;

          }

          //set to zero
          gradient.set(index,0,0,'c');
          gradient.set(index,1,0,'c');
          gradient.set(index,2,0,'c');



        } else {
          //set all to zero
          gradient.set(index,0,0,'c');
          gradient.set(index,1,0,'c');
          gradient.set(index,2,0,'c');
        }
    }*/

  }

  //printf("\n");
  //printf("Inside prune1\n");
  //for (i=0;i<num_hp;i++){
  // printf("      %d, %d, %d \n", (int)get(i,0), (int)get(i,1), (int)get(i,2) );
  //}
  //printf("\n");

}



/**
 *
 * Finds the maximum value of the coefficients withing the gradient class.
 *
 */

void hot_spots::find_max(gradient_class &grd,int *index, data_type *max_val){


  int i;
  data_type cx,cy,cz,curr_val;

  *index = 0;
  cx = grd.get(0,0,'c');
  cy = grd.get(0,1,'c');
  cz = grd.get(0,2,'c');
  *max_val = (data_type)sqrt(cx*cx + cy*cy + cz*cz);
  //	printf("in find_max, max_val is %f \n", *max_val);
  //	printf("in find_max, grd.numl is %d \n", grd.num_points);

  for (i=0;i<grd.num_points;i++){

    cx = grd.get(i,0,'c');
    cy = grd.get(i,1,'c');
    cz = grd.get(i,2,'c');
    curr_val = (data_type)sqrt(cx*cx + cy*cy + cz*cz);
    //	printf("i : %d \n", i);
    //	printf("curr_val is %f \n", curr_val);

    if (curr_val>(*max_val)){
      *max_val = curr_val;
      *index = i;
    }

  }
}



void hot_spots::find_max_multi(gradient_class &grd,int *_index, data_type *_max_val){
    int number=grd.num_points;
    int core=core_number;
    called_grd_max **p;
    p =(called_grd_max **)malloc(number*sizeof(called_grd_max *));
    int gap=(number+core-1)/core;
    int left=0;
    int right=gap;
    for (int i=0;i<core;i++){
        p[i]=new called_grd_max;
        p[i]->init(&grd,left,right);
        p[i]->start();
        left+=gap;
        right+=gap;
        right=mymin(right,number);
    }
    p[0]->wait();
    int index=p[0]->index;
    data_type max_val=p[0]->max_val;
    for (int i=1;i<core;i++){
        p[i]->wait();
        if (p[i]->max_val>max_val){
            index=p[i]->index;
            max_val=p[i]->max_val;
        }
    }
    *_index=index;
    *_max_val=max_val;
}


/**
 *
 * Given the index of a point in the gradient structure, determine if that point
 *       already exists in the hp list of hot spots, or is too close to any of the
 *       points in hp.
 *
 * Returns zero or one. Zero indicates the point specified by index is not close to 
 *       any other already in hp. One indicates it is.
 *
 *
 */


int hot_spots::check_distance(gradient_class &grd, int index){

  int i;
  int spcx, spcy, spcz;
  int diffx, diffy, diffz;

  //spc = grd.supp;
  // spcz = grd.suppz;
  spcx = (int)ceil(grd.spcx);
  spcy = (int)ceil(grd.spcy);
  spcz = (int)ceil(grd.spcz);

  //printf("index in check_distance : %d \n", index);

  for (i=0;i<num_hp;i++){

    diffx = (int)(grd.get(index,0,'p')-get(i,0));
    if(diffx<0) diffx=-diffx;

    diffy = (int)(grd.get(index,1,'p')-get(i,1));
    if(diffy<0) diffy=-diffy;

    diffz = (int)(grd.get(index,2,'p')-get(i,2));
    if(diffz<0) diffz=-diffz;

    // printf(" ---------------------------> check_distance: %f , %f , %f \n", grd.get(index,0,'p'), grd.get(index,1,'p'), grd.get(index,2,'p'));
    // printf(" --------------------------->               : %d , %d , %d \n", get(i,0), get(i,1) , get(i,2));
    // printf("\n");

    if ( (diffx <= spcx) && (diffy <= spcy) && (diffz<=spcz) ){

      return 1;

    }

  }

  return 0;
}



int hot_spots::check_distance_multi(gradient_class &grd, int index){
    int number=grd.num_points;
    int core=core_number;
    called_grd_diff **p;
    p =(called_grd_diff **)malloc(number*sizeof(called_grd_diff *));
    int gap=(number+core-1)/core;
    int left=0;
    int right=gap;
    for (int i=0;i<core;i++){
        p[i]=new called_grd_diff;
        p[i]->init(&grd,this,index,left,right);
        p[i]->start();
        left+=gap;
        right+=gap;
        right=mymin(right,number);
    }
    while (1){
        int count = 0;
        for (int i=0;i<core;i++){
            if (p[i]->flag==1){
                for (int j=0;j<core;j++){
                    p[j]->quit();
                }
                return 1;
            }
            count += p[i]->flag;
        }
        if (count == 0)
            return 0;

    }
}


/**
 *
 * Allows acces to modifying hp.
 *
 */

void hot_spots::set(int x, int y, int val){
	
   long coord;

   coord = x*3 + y;
		
   hp[coord]=(int)val;

}

/**
 *
 * Allows access to coordinate values hp
 *
 */

int hot_spots::get(int x, int y) {

  long coord;

  coord = x*3 + y;
		
  return hp[coord];

}


/**
 *
 * Destructor
 *
 */

void hot_spots::destroy(void){

  free(hp);

}
