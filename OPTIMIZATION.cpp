/**
 *
 * FILE: OPTIMIZATION.h
 * Gustavo Rohde, fall 2000
 * Class responsible for running Non Rigid Mutual Information based registration algorithm as
 * described in Rohde Et all.
 * Instantiates all data, and implements all necessary functions.
 * The propper order of call is:
 *									 init();
 *									 run();// or run_gus
 *									 output("outputfile");
 *									 destroy();
 */
//#define UNIX
#define PC

#ifdef UNIX
#include <stream.h>
#else
#include <iostream>
#endif


#include "register3d.h"
#include "OPTIMIZATION.h"
#include "matrix2d.h"
#include "hot_spots.h"
#include "config_class.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <string.h>

#include "opt_hot_spots.h"
#include "called_grd_build.h"




/**
 * Initializes necessary data for optimization algorithm.
 * Pre=requisite for the run() method.
 * Tasks:
 *			Load volumes, which will be known hereafter as "source" and "tartet".
 *			Load initial settings:
 *									Starting number of control points in each dimmension.
 *									Maximum and minimum values to be considered in volumes.
 *									Number of levels allowed.
 *									Number of histogram bins for both volumes.
 *
 *			Allocate and initialize the deformation field known hereafter as "def_field".
 *			Evaluate the Joint Histogram to be optimized, known hereafter as "JH".
 *			Compute the bounding box of the union of the source and target values.
 *
 *
 */

void OPTIMIZATION::init(config_class &cfg){



    int size_x,size_y,size_z;

    FLEX = cfg.FLEX;
    br_EPS = cfg.br_EPS;
    min_EPS = cfg.min_EPS;
    br_STEP = cfg.br_STEP;
    grd_STEP = cfg.grd_STEP;
    pick_const = cfg.pick_const;
    X_opt = cfg.X_opt;
    Y_opt = cfg.Y_opt;
    Z_opt = cfg.Z_opt;

    size_x = cfg.nx;
    size_y = cfg.ny;
    size_z = cfg.nz;
    SAVE_DEF_FIELD = cfg.SAVE_DEF_FIELD;
    s_trsh = cfg.s_trsh;

    printf("\n");
    printf("       Nonrigid Image Registration with Mutual Information (r3d).\n");
    printf("                        by Gustavo Rohde, 2001.\n");

    num_levels = cfg.num_levels;
    //printf("Number of levels requested: %d \n",num_levels);


    printf("\n");
    printf("\n");
    printf("Registering \n");
    printf(cfg.flname_src);
    printf("\n");
    printf(" to ");
    printf("\n");
    printf(cfg.flname_trg);
    printf(".\n");
    printf("Output in \n");
    printf(cfg.out);
    printf(".\n");


    printf("Number of levels requested = %d\n",cfg.num_levels);
    int i;
    //for(i=0;i<cfg.num_levels;i++){
    //  printf(" %d", cfg.ctrpts_array[i]);
    //}
    //printf("\n");

    //for(i=0;i<cfg.num_levels;i++){
    //  printf(" %d", cfg.ctrpts_arrayZ[i]);
    // }
    //printf("\n");

    printf("Number of resolutions = %d\n",cfg.num_resolutions);
    printf("Maximum resolution at = %d\n",cfg.max_res);

    //for(i=0;i<cfg.max_res;i++){
    // printf(" %d", cfg.res_array[i]);
    //}
    //printf("\n");
    //printf("%d\n",cfg_obj.header_src);
    //printf(cfg_obj.flname_trg);
    //printf("\n");
    //printf("%d\n",cfg_obj.header_trg);
    //printf("%d\n",cfg_obj.nx);
    //printf("%d\n",cfg_obj.ny);
    //printf("%d\n",cfg_obj.nz);


    //printf("\n");
    printf("Maximum data value to include in cost function computation: %d \n", cfg.mmax);
    printf("Minimum data value to include in cost function computation: %d \n", cfg.mmin);

    printf("Number of histogram bins in Source image: %d \n", cfg.nbins_src);
    printf("Number of histogram bins in Target image: %d \n", cfg.nbins_trg);

    if (cfg.OP_MODE==2)
        printf("Using multi Adaptive Grid Algorithm.\n");
    else if (cfg.OP_MODE==1)
        printf("Using Adaptive Grid Algorithm.\n");
    else
        printf("Using Regular Grid Approach.\n");

    if (cfg.B_SEG==1)
        printf("Optimizing region within bounding box of Source image only.\n");

    //printf("\n");
    //printf("OP_MODE: %d \n", cfg_obj.OP_MODE);
    //printf("FLEX: %f \n", cfg_obj.FLEX);
    //printf("br_EPS: %f \n", cfg_obj.br_EPS);
    //printf("min_EPS: %f \n", cfg_obj.min_EPS);
    //printf("br_EPS: %f \n", cfg_obj.br_STEP);
    //printf("grd_STEP: %f \n", cfg_obj.grd_STEP);
    //printf("pick_const: %f \n", cfg_obj.pick_const);
    //printf("opt: %d %d %d \n", cfg_obj.X_opt, cfg_obj.Y_opt, cfg_obj.Z_opt);


    num_resolutions = cfg.num_resolutions;
    curr_resolution = num_resolutions;

    curr_res_index = 0;
    max_res_index = cfg.max_res;

    res_array = new int[num_resolutions];
    ctrpts_array = new int[num_levels];
    ctrptsZ_array = new int[num_levels];

    //int i;
    for(i=0;i<cfg.num_levels;i++){
        ctrpts_array[i]=cfg.ctrpts_array[i];
    }

    for(i=0;i<cfg.num_levels;i++){
        ctrptsZ_array[i]=cfg.ctrpts_arrayZ[i];
    }

    for(i=0;i<cfg.max_res;i++){
        res_array[i]=cfg.res_array[i];
    }

    npoints_x = ctrpts_array[0];
    npoints_y = ctrpts_array[0];
    npoints_z = ctrptsZ_array[0];

    num_rbf_init = 0;

    mmax = cfg.mmax;
    mmin = cfg.mmin;

    nbins_src = cfg.nbins_src;
    nbins_trg = cfg.nbins_trg;

    JH.init(nbins_src,nbins_trg);
    JH.setzero();

    //read in volumes

    source.init(size_x,size_y,size_z,cfg.header_src,cfg.flname_src);
    target.init(size_x,size_y,size_z,cfg.header_trg,cfg.flname_trg);

    printf("Reading images ... \n");
    source.vol_read();
    target.vol_read();
    printf("Images read.\n");
    printf("\n");
    printf("initializing ...\n");

    source.build_resolution(num_resolutions);
    target.build_resolution(num_resolutions);
    //source.build_resolution(0);
    //target.build_resolution(0);

    //	printf("Building resolution finished.Size is %d,%d,%d\n",source.nx,source.ny,source.nz);

    //initialize deformation field
    //def_field.init(size_x,size_y,size_z);
    def_field.init(source.nx,source.ny,source.nz);
    //printf("Deformation field initialized.\n");

    time(&time_init);
    //build 2d joint histogram
    JH_inst();
    //printf("JH_inst computed.\n");


    //establish bounding box of source and target combined
    source.compt_bbox();
    target.compt_bbox();


    if (cfg.B_SEG==1){

        bbox.lx = mymax(0,target.lx-5);
        bbox.rx = mymin(target.rx+5,target.nx-1);
        bbox.ty = mymax(0,target.ty-5);
        bbox.by = mymin(target.by+5,target.ny-1);
        bbox.fz = mymax(0,target.fz-5);
        bbox.bz = mymin(target.bz+5,target.nz-1);

        printf("The bounding box is: \n");
        printf("%d,%d\n",bbox.lx,bbox.rx);
        printf("%d,%d\n",bbox.ty,bbox.by);
        printf("%d,%d\n",bbox.fz,bbox.bz);
        printf("\n");

    }else{

        //x
        if (source.lx < target.lx)
            bbox.lx = source.lx;
        else
            bbox.lx = target.lx;

        if (source.rx > target.rx)
            bbox.rx = source.rx;
        else
            bbox.rx = target.rx;

        //y
        if (source.ty < target.ty)
            bbox.ty = source.ty;
        else
            bbox.ty = target.ty;

        if (source.by > target.by)
            bbox.by = source.by;
        else
            bbox.by = target.by;

        //z
        if (source.fz < target.fz)
            bbox.fz = source.fz;
        else
            bbox.fz = target.fz;

        if (source.bz > target.bz)
            bbox.bz = source.bz;
        else
            bbox.bz = target.bz;
    }

    //printf("The bounding box is: \n");
    //printf("%d,%d\n",bbox.lx,bbox.rx);
    //printf("%d,%d\n",bbox.ty,bbox.by);
    //printf("%d,%d\n",bbox.fz,bbox.bz);
    //printf("\n");

    //bin_vol.init2(num_resolutions,size_x, size_y, size_z, bbox.lx,bbox.rx,bbox.ty,bbox.by,bbox.fz,bbox.bz, source, target, "bin_im2.vol");
    //bin_vol.dilate();
    //printf("initialized ok\n");
    //bin_vol.build_resolution(0);
    //printf("build resolution\n");
    //bin_vol.vol_write();
    //printf("wrote volume\n");

    //printf("Initialization finished.\n");
    //printf("\n");



}




/**
 *
 * Updates the JH (Joint Histogram) for the current deformation field.
 * Tasks:
 *			Go through entire domain of deformation field regardless, and recompute JH in its entirety.
 *
 *
 */

void OPTIMIZATION::JH_inst(){

    //NOTE: target image histogram does not Change, it is stored as the second argument
    int x,y,z,jhx,jhy;
    data_type del;
    data_type src_val, trg_val;

    JH.setzero();
    del = (data_type)(mmax-mmin);

    //do simple count
    for(x=0;x<source.nx;x++){
        for(y=0;y<source.ny;y++){
            for(z=0;z<source.nz;z++){

                src_val = source.get_intrpl(def_field.get(x,y,z,'x'),def_field.get(x,y,z,'y'),def_field.get(x,y,z,'z'));
                trg_val = (data_type)target.get(x,y,z);

                if ((src_val<mmin)||(trg_val<mmin)||(src_val>mmax)||(trg_val>mmax)) {

                    //do nothing

                }else{

                    jhx = (int)((((data_type)(src_val-mmin))/del)*(JH.nx-1));
                    jhy = (int)((((data_type)(trg_val-mmin))/del)*(JH.ny-1));
                    JH.increment(jhx,jhy);
                }

            }
        }
    }

    //printf("True current Mutual Information: %f\n",JH.mi());
    //time(&time_now);
    //printf("%d %f\n",(time_now-time_init),JH.mi());


}




/**
 *
 * NOTE: This method is entirely outdated and should not be used anywhere in the code.
 *			It is kept here for syntax purposes. JH_inst_update is its replacement.
 *
 */

void OPTIMIZATION::JH_inst1(){

    //NOTE: Please note that efficiency could be highly improved here by taking advantage of locality
    int x,y,z,jhx,jhy;
    data_type del;
    data_type src_val, trg_val;
    sampling_coord def_field_temp;

    int i;
    int ii,jj,kk,bx,by,bz;
    int vx_start,vx_end,vy_start,vy_end,vz_start,vz_end;
    int bx_start,bx_end,by_start,by_end,bz_start,bz_end;
    data_type new_valx,new_valy,new_valz;


    //first create temporary deformation field
    def_field_temp.init(def_field.nx,def_field.ny,def_field.nz);
    for(x=0;x<source.nx;x++){
        for(y=0;y<source.ny;y++){
            for(z=0;z<source.nz;z++){


                def_field_temp.set(x,y,z,'x',def_field.get(x,y,z,'x'));
                def_field_temp.set(x,y,z,'y',def_field.get(x,y,z,'y'));
                def_field_temp.set(x,y,z,'z',def_field.get(x,y,z,'z'));

            }
        }
    }


    //def_field_temp = def_field_temp + coeff*RBF
    for(i=0;i<gradient.num_points;i++){


        //find appropriate indexes for acessing
        vx_start = (int)mymax((gradient.get(i,0,'p')-gradient.supp),0);
        vx_end = (int)mymin((gradient.get(i,0,'p')+gradient.supp),(source.nx));
        vy_start = (int)mymax((gradient.get(i,1,'p')-gradient.supp),0);
        vy_end = (int)mymin((gradient.get(i,1,'p')+gradient.supp),(source.ny));
        vz_start = (int)mymax((gradient.get(i,2,'p')-gradient.supp),0);
        vz_end = (int)mymin((gradient.get(i,2,'p')+gradient.supp),(source.nz));

        //find appropriate indexes for RBF
        bx_start = (int)(gradient.supp - (gradient.get(i,0,'p') - vx_start));
        bx_end = (int)(gradient.supp + (vx_end - gradient.get(i,0,'p')));
        by_start = (int)(gradient.supp - (gradient.get(i,1,'p') - vy_start));
        by_end = (int)(gradient.supp + (vy_end - gradient.get(i,1,'p')));
        bz_start = (int)(gradient.supp - (gradient.get(i,2,'p') - vz_start));
        bz_end = (int)(gradient.supp + (vz_end - gradient.get(i,2,'p')));


        for(ii=vx_start,bx=bx_start;ii<vx_end,bx<bx_end;ii++,bx++){
            for(jj=vy_start,by=by_start;jj<vy_end,by<by_end;jj++,by++){
                for(kk=vz_start,bz=bz_start;kk<vz_end,bz<bz_end;kk++,bz++){

                    new_valx = def_field_temp.get(ii,jj,kk,'x') + gradient.get(i,0,'c')*RBF.get(bx,by,bz);
                    new_valy = def_field_temp.get(ii,jj,kk,'y') + gradient.get(i,1,'c')*RBF.get(bx,by,bz);
                    new_valz = def_field_temp.get(ii,jj,kk,'z') + gradient.get(i,2,'c')*RBF.get(bx,by,bz);

                    def_field_temp.set(ii,jj,kk,'x',new_valx);
                    def_field_temp.set(ii,jj,kk,'y',new_valy);
                    def_field_temp.set(ii,jj,kk,'z',new_valz);



                }
            }
        }


    }


    JH.setzero();
    del = (data_type)(mmax-mmin);

    //do simple count
    for(x=0;x<source.nx;x++){
        for(y=0;y<source.ny;y++){
            for(z=0;z<source.nz;z++){

                src_val = source.get_intrpl(def_field_temp.get(x,y,z,'x'),def_field_temp.get(x,y,z,'y'),def_field_temp.get(x,y,z,'z'));
                trg_val = (data_type)target.get(x,y,z);

                if ((src_val<mmin)||(trg_val<mmin)||(src_val>mmax)||(trg_val>mmax)) {

                    //do nothing

                }else{

                    jhx = (int)((((data_type)(src_val-mmin))/del)*(JH.nx-1));
                    jhy = (int)((((data_type)(trg_val-mmin))/del)*(JH.ny-1));
                    JH.increment(jhx,jhy);
                }

            }
        }
    }

    def_field_temp.destroy();
    //printf("Current Mutual Information: %f\n",JH.mi());

}



void OPTIMIZATION::run_multi(){

    int i;
    int STOP = 0;
    //data_type2 last_min,curr_min,diff;
    curr_level = 0;
    hot_spots htspts;

    printf("Registration has started.\n");
    time(&time_now);
    time_t time_last=time_now;
    printf("Time: %ld\n",(time_now-time_init));

    //loop on all levels of algorithm
    for(i=0;i<num_levels;i++){

        npoints_x = ctrpts_array[i];
        npoints_y = ctrpts_array[i];
        npoints_z = ctrptsZ_array[i];

        printf("LEVEL %d out of %d.\n",i+1,num_levels);


        //first establish regular grid based on current level
        gradient.init(bbox.lx,bbox.rx,bbox.ty,bbox.by,bbox.fz,bbox.bz,npoints_x,npoints_y,npoints_z,source.nx,source.ny,source.nz,FLEX);

        //create appropriate radial basis function
        RBF.init(gradient.supp,gradient.suppz,num_rbf_init);
        printf("     Current support of RBF %d,%d.\n",gradient.supp,gradient.suppz);

        grd_multi();

        time_last=time_now;
        time(&time_now);
        printf("Time %ld cost in grd_multi...\n",(time_now-time_last));

        line_minimize();

        time_last=time_now;
        time(&time_now);
        printf("Time %ld cost in line_minimize...\n",(time_now-time_last));



        /*
    last_min = JH.mi();
    STOP=0;

    while(STOP<1){ //evaluate G once, allow for choosing more than once

      //sprintf(std_msg,"     Computing gradient of regular grid ... ");
      grd();
      //sprintf(std_msg,"     Doing line minimization of regular grid ... ");
      curr_min = line_minimize();


      diff = last_min-curr_min;
      if(diff<0) diff=-diff;
      if (diff<min_EPS) {
    //sprintf(std_msg,"should have stoped");
    STOP=1;
      } else {
    last_min = curr_min;
      }

      STOP = 1;

    }*/


        //hold temporarily
        int tx, ty, tz, support, supportz;
        tx = (int)ceil(gradient.spcx);
        ty = (int)ceil(gradient.spcy);
        tz = (int)ceil(gradient.spcz);
        support = gradient.supp;
        supportz = gradient.suppz;

        def_field_update();

        JH_inst();

        time_last=time_now;
        time(&time_now);
        printf("Time %ld cost in hold temporarity...\n",(time_now-time_last));

        //now identify "hot spots"
        data_type2 konst = pick_const;
        htspts.prune3(gradient, konst, num_rbf_init, source, target);
        printf("     Number of hot spots is %d \n", htspts.num_hp);

        gradient.destroy();

        time_last=time_now;
        time(&time_now);
        printf("Time %ld cost in identify hotspots...\n",(time_now-time_last));

        printf("Registration has started.\n");
        time(&time_now);
        printf("Time: %ld\n",(time_now-time_init));


        //multi_processing "hot spots"
        int j=0;
        if (htspts.num_hp>0){
            opt_hot_spots **p;
            p =(opt_hot_spots **)malloc(htspts.num_hp*sizeof(opt_hot_spots *));
            for (j=0;j<htspts.num_hp;j++){
                p[j]=new opt_hot_spots;
                p[j]->init( &htspts,&JH,&source,&target,&def_field,&RBF,mmax,mmin,nbins_src,nbins_trg);
                p[j]->grd_STEP=grd_STEP;
                p[j]->br_STEP=br_STEP;
                p[j]->X_opt=X_opt;
                p[j]->Y_opt=Y_opt;
                p[j]->Z_opt=Z_opt;
                p[j]->min_EPS=min_EPS;
                p[j]->s_trsh=s_trsh;
                p[j]->gradient.init2(htspts.get(j,0),htspts.get(j,1),htspts.get(j,2),tx,ty,tz,support,supportz,source.nx,source.ny,source.nz);
                p[j]->start();
                //if ((j+1)%core_number==0)
                //    p[j]->wait();
            }
            for (j=0;j<htspts.num_hp;j++){
                p[j]->wait();
            }
            for (j=0;j<htspts.num_hp;j++){
                def_field_update2(&(p[j]->gradient));
                p[j]->gradient.destroy();
            }
        }

        time_t time_last=time_now;
        time(&time_now);
        printf("Time %ld cost in processing hot spots...\n",(time_now-time_last));

        //reset
        STOP=0;
        num_rbf_init++;
        //npoints_x = npoints_x + (npoints_x - 1);
        //npoints_y = npoints_y + (npoints_y - 1);
        //npoints_z = npoints_z + (npoints_z - 1);
        RBF.reset();
        curr_level++;

        //sprintf(std_msg," The new test: res_array[curr_res_index] = %d, npoints_x = %d.", res_array[curr_res_index], npoints_x);

        time_last=time_now;
        time(&time_now);
        printf("Time %ld cost in reset RBF...\n",(time_now-time_last));

        //upsample resolution if necessary
        if ( (res_array[curr_res_index]<=npoints_x)&&(max_res_index>curr_res_index) ){
            printf("\n");

            printf("Changing resolutions %d.\n", curr_resolution-1);
            curr_resolution--;

            source.destroy();//destroys only current resolution
            target.destroy();//destroys only current resolution

            source.build_resolution(curr_resolution);
            target.build_resolution(curr_resolution);


            def_field.upsample_multi();

            bbox.lx = bbox.lx*2;
            bbox.rx = bbox.rx*2;
            bbox.by = bbox.by*2;
            bbox.ty = bbox.ty*2;
            bbox.fz = bbox.fz*2;
            bbox.bz = bbox.bz*2;

            JH_inst();

            time(&time_now);
            printf("Time,  Cost Function\n");
            printf("%d %f\n",(time_now-time_init),JH.mi());


            curr_res_index++;

        }

        time_last=time_now;
        time(&time_now);
        printf("Time %ld cost in upsample...\n",(time_now-time_last));


        //end
    }

    //sprintf(std_msg,"Generating result.");
    if (curr_resolution!=0) {

        source.destroy();//destroys only current resolution
        target.destroy();//destroys only current resolution
        source.build_resolution(0);
        target.build_resolution(0);

        for (i=0;i<curr_resolution;i++){
            def_field.upsample_multi();
        }
        JH_inst();

    }

    printf("Optimization completed.\n");

    time(&time_now);
    printf("Time: %ld seconds.\n",(time_now-time_init));

    printf("Final value for cost function %f .\n", JH.mi());

    generate_result();

}



/**
 *
 * This is the method that actually control all instantiated classes and optimizes the deformation_field.
 * Tasks:
 *	      	Loop through all levels of algorithm and optimize "def_field" in each.
 *		For each level do:
 *		      		Place control points in a regular grid.
 *     				Build appropriate Radial Basis Function "RBF".
 *		      		Evaluate gradient of cost function with respect to grid. (call "grd()" )
 *		       		Do line minimization once on entire grid. (call "line_minimize()" )
 *
 *					 Loop through gradient and choose "hot spots".
 *					      For each hot spot do:
 *						Build appropriate RBF.
 *					       	Optimize hot spot (evaluate gradient, line minimize ...)
 *
 *               Generate final result by resampling the source volume with optimized "def_field". Store the result in "source".
 *
 *
 */

void OPTIMIZATION::run_gus(){

    int i;
    int STOP = 0;
    data_type2 last_min,curr_min,diff;
    curr_level = 0;
    hot_spots htspts;

    printf("Registration has started.\n");

    //loop on all levels of algorithm
    for(i=0;i<num_levels;i++){

        npoints_x = ctrpts_array[i];
        npoints_y = ctrpts_array[i];
        npoints_z = ctrptsZ_array[i];

        printf("LEVEL %d out of %d.\n",i+1,num_levels);

        //first establish regular grid based on current level
        gradient.init(bbox.lx,bbox.rx,bbox.ty,bbox.by,bbox.fz,bbox.bz,npoints_x,npoints_y,npoints_z,source.nx,source.ny,source.nz,FLEX);

        //create appropriate radial basis function
        RBF.init(gradient.supp,gradient.suppz,num_rbf_init);
        printf("     Current support of RBF %d,%d.\n",gradient.supp,gradient.suppz);

        last_min = JH.mi();
        STOP=0;

        while(STOP<1){ //evaluate G once, allow for choosing more than once

            //printf("     Computing gradient of regular grid ... \n");
            grd();
            //printf("     Doing line minimization of regular grid ... \n");
            curr_min = line_minimize();


            diff = last_min-curr_min;
            if(diff<0) diff=-diff;
            if (diff<min_EPS) {
                //printf("should have stoped");
                STOP=1;
            } else {
                last_min = curr_min;
            }

            STOP = 1;

        }


        //hold temporarily
        int tx, ty, tz, support, supportz;
        tx = (int)ceil(gradient.spcx);
        ty = (int)ceil(gradient.spcy);
        tz = (int)ceil(gradient.spcz);
        support = gradient.supp;
        supportz = gradient.suppz;

        def_field_update();

        JH_inst();


        //now identify "hot spots"
        data_type2 konst = pick_const;
        htspts.prune1(gradient, konst, num_rbf_init, source, target);
        printf("     Number of hot spots is %d \n", htspts.num_hp);

        gradient.destroy();

        //loop through "hot spots"
        int j;
        for (j=0;j<htspts.num_hp;j++){

            gradient.init2(htspts.get(j,0),htspts.get(j,1),htspts.get(j,2),tx,ty,tz,support,supportz,source.nx,source.ny,source.nz);
            last_min = JH.mi();

            printf("     Point %d out of %d. Point coordinate: %d, %d, %d.\n", j+1,htspts.num_hp,htspts.get(j,0),htspts.get(j,1),htspts.get(j,2));
            STOP=0;
            while(STOP==0){


                grd();
                //printf("          Inside line minimize ... \n");
                curr_min = line_minimize();


                diff = last_min-curr_min;
                if(diff<0) diff=-diff;

                if ( (diff<0.001) || (curr_min>last_min)) {

                    //printf("should have stoped");
                    STOP=2;

                } else {
                    last_min = curr_min;
                }
                //STOP++;

            }

            def_field_update();
            //JH_inst();
            gradient.destroy();

        }

        //reset
        STOP=0;
        num_rbf_init++;
        //npoints_x = npoints_x + (npoints_x - 1);
        //npoints_y = npoints_y + (npoints_y - 1);
        //npoints_z = npoints_z + (npoints_z - 1);
        RBF.reset();
        curr_level++;

        //printf(" The new test: res_array[curr_res_index] = %d, npoints_x = %d.\n", res_array[curr_res_index], npoints_x);

        //upsample resolution if necessary
        if ( (res_array[curr_res_index]<=npoints_x)&&(max_res_index>curr_res_index) ){
            printf("\n");
            printf("Changing resolutions %d.\n", curr_resolution-1);
            curr_resolution--;

            source.destroy();//destroys only current resolution
            target.destroy();//destroys only current resolution

            source.build_resolution(curr_resolution);
            target.build_resolution(curr_resolution);


            def_field.upsample();

            bbox.lx = bbox.lx*2;
            bbox.rx = bbox.rx*2;
            bbox.by = bbox.by*2;
            bbox.ty = bbox.ty*2;
            bbox.fz = bbox.fz*2;
            bbox.bz = bbox.bz*2;

            JH_inst();

            time(&time_now);
            printf("Time,  Cost Function\n");
            printf("%ld %f\n",(time_now-time_init),JH.mi());

            curr_res_index++;

        }



        //end
    }

    //printf("Generating result.\n");
    if (curr_resolution!=0) {

        source.destroy();//destroys only current resolution
        target.destroy();//destroys only current resolution
        source.build_resolution(0);
        target.build_resolution(0);

        for (i=0;i<curr_resolution;i++){
            def_field.upsample();
        }
        JH_inst();

    }

    printf("Optimization completed.\n");
    time(&time_now);
    //printf("%d %f\n",(time_now-time_init),JH.mi());
    printf("Time: %ld seconds.\n",(time_now-time_init));
    printf("Final value for cost function %f .\n", JH.mi());
    generate_result();

}




/**
 *
 * Runs the old version of the algoritm. Basic regular grid optimization.
 *
 */

void OPTIMIZATION::run(){

    int i;
    int STOP = 0;
    data_type2 last_min,curr_min,diff;
    curr_level = 0;

    printf("Registration has started.\n");
    //loop on all levels of algorithm
    //while (curr_level <= num_levels )
    for(i=0;i<num_levels;i++){

        printf("LEVEL %d out of %d.\n",i+1,num_levels);

        npoints_x = ctrpts_array[i];
        npoints_y = ctrpts_array[i];
        npoints_z = ctrptsZ_array[i];

        gradient.init(bbox.lx,bbox.rx,bbox.ty,bbox.by,bbox.fz,bbox.bz,npoints_x,npoints_y,npoints_z,source.nx,source.ny,source.nz,FLEX);
        RBF.init(gradient.supp,gradient.suppz,num_rbf_init);
        printf("     Current support of RBF %d,%d.\n",gradient.supp,gradient.suppz);

        num_rbf_init++;
        last_min = JH.mi();
        while(STOP==0){

            //printf("     Computing gradient of regular grid ... \n");
            grd();
            //printf("     Doing line minimization of regular grid ... \n");

            curr_min = line_minimize();

            diff = last_min-curr_min;
            if(diff<0) diff=-diff;

            if ( (diff<min_EPS) || (curr_min>last_min)) {
                //printf("should have stoped");
                STOP=1;
            } else {
                //STOP=1;
                last_min = curr_min;
            }
            //STOP=1;

        }

        STOP=0;
        def_field_update();
        JH_inst();
        gradient.destroy();
        //npoints_x = npoints_x + (npoints_x - 1);
        //npoints_y = npoints_y + (npoints_y - 1);
        //npoints_z = npoints_z + (npoints_z - 1);
        RBF.reset();
        curr_level++;

        //upsample resolution if necessary
        if ( (res_array[curr_res_index]<=npoints_x)&&(max_res_index>curr_res_index) ){

            printf("Changing resolutions %d.\n", curr_resolution-1);

            curr_resolution--;
            source.destroy();//destroys only current resolution
            target.destroy();//destroys only current resolution

            source.build_resolution(curr_resolution);
            target.build_resolution(curr_resolution);

            def_field.upsample();

            bbox.lx = bbox.lx*2;
            bbox.rx = bbox.rx*2;
            bbox.by = bbox.by*2;
            bbox.ty = bbox.ty*2;
            bbox.fz = bbox.fz*2;
            bbox.bz = bbox.bz*2;

            JH_inst();

            time(&time_now);
            printf("Time,  Cost Function\n");
            printf("%ld %f\n",(time_now-time_init),JH.mi());

            curr_res_index++;

        }
        //end
    }


    //printf("Generating result.\n");
    if (curr_resolution!=0) {

        source.destroy();//destroys only current resolution
        target.destroy();//destroys only current resolution
        source.build_resolution(0);
        target.build_resolution(0);

        for (i=0;i<curr_resolution;i++){
            def_field.upsample();
        }
        JH_inst();

    }

    printf("Optimization completed.\n");
    time(&time_now);
    //printf("%d %f\n",(time_now-time_init),JH.mi());
    printf("Time: %ld seconds.\n",(time_now-time_init));
    printf("Final value for cost function &f .\n", JH.mi());
    generate_result();



}


/**
 *
 * Function called from "line_minimize()".
 * Updates the JH (joint histogram) based on the newly computed coefficients in "line_minimize()".
 * Tasks:
 *			Update JH based on new coefficients without actually changing "def_field."
 *
 *
 */

void OPTIMIZATION::JH_inst_up(matrix_2d &JH_temp,
                              sampling_coord &def_field_temp, int minx, int maxx, int miny, int maxy,
                              int minz, int maxz){



    int i;
    int ii,jj,kk;
    int vx_start,vx_end,vy_start,vy_end,vz_start,vz_end;
    int bx_start,bx_end,by_start,by_end,bz_start,bz_end;
    int bx, by, bz;
    data_type new_valx,new_valy,new_valz;
    data_type cx,cy,cz;
    int x,y,z,jhx,jhy;
    data_type del;
    data_type src_val, trg_val;

    del = (data_type)(mmax-mmin);

    //now do def_field_temp = def_field (only for the current region) again
    for(x=0;x<def_field_temp.nx;x++){
        for(y=0;y<def_field_temp.ny;y++){
            for(z=0;z<def_field_temp.nz;z++){

                def_field_temp.set(x,y,z,'x',def_field.get(x+minx,y+miny,z+minz,'x'));
                def_field_temp.set(x,y,z,'y',def_field.get(x+minx,y+miny,z+minz,'y'));
                def_field_temp.set(x,y,z,'z',def_field.get(x+minx,y+miny,z+minz,'z'));

            }
        }
    }


    //redo def_field_temp
    for(i=0;i<gradient.num_points;i++){


        //find appropriate indexes for acessing
        //		vx_start = (int)max((gradient.get(i,0,'p')-gradient.supp),0);
        //vx_end = (int)min((gradient.get(i,0,'p')+gradient.supp+1),(def_field_temp.nx));
        //vy_start = (int)max((gradient.get(i,1,'p')-gradient.supp),0);
        //vy_end = (int)min((gradient.get(i,1,'p')+gradient.supp+1),(def_field_temp.ny));
        //vz_start = (int)max((gradient.get(i,2,'p')-gradient.suppz),0);
        //vz_end = (int)min((gradient.get(i,2,'p')+gradient.suppz+1),(def_field_temp.nz));
        vx_start = (int)mymax((gradient.get(i,0,'p')-gradient.supp-minx),0);
        vx_end = (int)mymin((gradient.get(i,0,'p')+gradient.supp+1-minx),(def_field_temp.nx));
        vy_start = (int)mymax((gradient.get(i,1,'p')-gradient.supp-miny),0);
        vy_end = (int)mymin((gradient.get(i,1,'p')+gradient.supp+1-miny),(def_field_temp.ny));
        vz_start = (int)mymax((gradient.get(i,2,'p')-gradient.suppz-minz),0);
        vz_end = (int)mymin((gradient.get(i,2,'p')+gradient.suppz+1-minz),(def_field_temp.nz));

        //find appropriate indexes for RBF
        //bx_start = (int)(gradient.supp - (gradient.get(i,0,'p') - vx_start));
        //	bx_end = (int)(gradient.supp + (vx_end - gradient.get(i,0,'p')));
        //by_start = (int)(gradient.supp - (gradient.get(i,1,'p') - vy_start));
        //by_end = (int)(gradient.supp + (vy_end - gradient.get(i,1,'p')));
        //bz_start = (int)(gradient.supp - (gradient.get(i,2,'p') - vz_start));
        //bz_end = (int)(gradient.supp + (vz_end - gradient.get(i,2,'p')));

        //find appropriate indexes for RBF
        bx_start = (int)(gradient.supp - ((gradient.get(i,0,'p')-minx) - vx_start));
        bx_end = (int)(gradient.supp + (vx_end - (gradient.get(i,0,'p')-minx)));
        by_start = (int)(gradient.supp - ((gradient.get(i,1,'p')-miny) - vy_start));
        by_end = (int)(gradient.supp + (vy_end - (gradient.get(i,1,'p')-miny)));
        bz_start = (int)(gradient.supp - ((gradient.get(i,2,'p')-minz) - vz_start));
        bz_end = (int)(gradient.supp + (vz_end - (gradient.get(i,2,'p')-minz)));

        cx=gradient.get(i,0,'c');
        cy=gradient.get(i,1,'c');
        cz=gradient.get(i,2,'c');

        for(ii=vx_start,bx=bx_start;ii<vx_end,bx<bx_end;ii++,bx++){
            for(jj=vy_start,by=by_start;jj<vy_end,by<by_end;jj++,by++){
                for(kk=vz_start,bz=bz_start;kk<vz_end,bz<bz_end;kk++,bz++){

                    new_valx = def_field_temp.get(ii,jj,kk,'x') + ( cx )*RBF.get(bx,by,bz);
                    new_valy = def_field_temp.get(ii,jj,kk,'y') + ( cy )*RBF.get(bx,by,bz);
                    new_valz = def_field_temp.get(ii,jj,kk,'z') + ( cz )*RBF.get(bx,by,bz);

                    //def_field_temp.set(ii-minx,jj-miny,kk-minz,'x',new_valx);
                    //def_field_temp.set(ii-minx,jj-miny,kk-minz,'y',new_valy);
                    //def_field_temp.set(ii-minx,jj-miny,kk-minz,'z',new_valz);

                    def_field_temp.set(ii,jj,kk,'x',new_valx);
                    def_field_temp.set(ii,jj,kk,'y',new_valy);
                    def_field_temp.set(ii,jj,kk,'z',new_valz);
                }
            }
        }


    }


    //printf("In JH_inst_up, the JH of images - current region is %f \n", JH_temp.mi());

    //now calculate JH_temp
    for(x=0;x<def_field_temp.nx;x++){
        for(y=0;y<def_field_temp.ny;y++){
            for(z=0;z<def_field_temp.nz;z++){

                src_val = source.get_intrpl(def_field_temp.get(x,y,z,'x'),def_field_temp.get(x,y,z,'y'),def_field_temp.get(x,y,z,'z'));
                trg_val = (data_type)target.get(x+minx,y+miny,z+minz);

                if ((src_val<mmin)||(trg_val<mmin)||(src_val>mmax)||(trg_val>mmax)) {

                    //do nothing

                }else{

                    jhx = (int)((((data_type)(src_val-mmin))/del)*(JH.nx-1));
                    jhy = (int)((((data_type)(trg_val-mmin))/del)*(JH.ny-1));
                    JH_temp.increment(jhx,jhy);
                }

            }
        }
    }


    //and finally JH = JH_temp
    for(ii=0;ii<JH_temp.nx;ii++){
        for(jj=0;jj<JH_temp.ny;jj++){

            JH.set(ii,jj,JH_temp.get(ii,jj));

        }
    }

    //printf("New mi is %f \n", JH.mi());

}





/**
 *
 * Built in aid to the function "line_minimize()";
 * Given a value "konst" it evaluates what the cost function (MI) value would be if
 *			a new deformation field was created by using the current coefficients plus
 *			the "konst" times the gradients.
 *
 * Takes in many argument to avoid recomputation and reallocation, etc.
 *
 *
 */

data_type2 OPTIMIZATION::f_eval(data_type konst, matrix_2d &JH_temp, matrix_2d &JH_temp_copy,
                                sampling_coord &def_field_temp, int minx, int maxx, int miny, int maxy,
                                int minz, int maxz){

    /* Inputs:

     1)konst = the constant by which to multiply the current gradient
     2)JH_temp = the current JH minus the joint histogram of the current region being worked on
     3)JH_temp_copy = just a space for storing JH_temp_coopy = JH_temp in the first few lines
     meant so that this function, which is called quite often doesn't have to allocate new
     memory every time.
     4)def_field_temp = a copy of the current def_field but only for the rigion which is currently
     being worked on.
     5,6 ...) the bounds (in coordinates) of the current region
     
     
     
  */


    // printf("                           Inside f_eval %d , %d , %d , %d , %d , %d \n", minx, maxx, miny, maxy, minz, maxz);


    int i;
    int ii,jj,kk;
    int vx_start,vx_end,vy_start,vy_end,vz_start,vz_end;
    int bx_start,bx_end,by_start,by_end,bz_start,bz_end;
    int bx, by, bz;
    data_type new_valx,new_valy,new_valz;
    data_type cx,cy,cz;
    int x,y,z,jhx,jhy;
    data_type del;
    data_type src_val, trg_val;

    del = (data_type)(mmax-mmin);

    //copy JH_temp into JH_Temp_copy
    for(ii=0;ii<JH_temp.nx;ii++){
        for(jj=0;jj<JH_temp.ny;jj++){

            JH_temp_copy.set(ii,jj,JH_temp.get(ii,jj));

        }
    }

    //printf("In f_eval, the JH of images - current region is %f \n", JH_temp.mi());
    //printf("In f_eval, the JH of images - current region is %f \n", JH_temp_copy.mi());
    //now do def_field_temp = def_field (only for the current region) again
    for(x=0;x<def_field_temp.nx;x++){
        for(y=0;y<def_field_temp.ny;y++){
            for(z=0;z<def_field_temp.nz;z++){

                def_field_temp.set(x,y,z,'x',def_field.get(x+minx,y+miny,z+minz,'x'));
                def_field_temp.set(x,y,z,'y',def_field.get(x+minx,y+miny,z+minz,'y'));
                def_field_temp.set(x,y,z,'z',def_field.get(x+minx,y+miny,z+minz,'z'));

            }
        }
    }

    //printf("              Inside f_eval checkpoint 1. \n");

    //redo def_field_temp
    for(i=0;i<gradient.num_points;i++){

        //printf("     point %d out of %d\n", i, gradient.num_points);

        //find appropriate indexes for acessing
        //vx_start = (int)max((gradient.get(i,0,'p')-gradient.supp),0);
        //vx_end = (int)min((gradient.get(i,0,'p')+gradient.supp+1),(def_field_temp.nx));
        //vy_start = (int)max((gradient.get(i,1,'p')-gradient.supp),0);
        //vy_end = (int)min((gradient.get(i,1,'p')+gradient.supp+1),(def_field_temp.ny));
        //vz_start = (int)max((gradient.get(i,2,'p')-gradient.suppz),0);
        //vz_end = (int)min((gradient.get(i,2,'p')+gradient.suppz+1),(def_field_temp.nz));
        vx_start = (int)mymax((gradient.get(i,0,'p')-gradient.supp-minx),0);
        vx_end = (int)mymin((gradient.get(i,0,'p')+gradient.supp+1-minx),(def_field_temp.nx));
        vy_start = (int)mymax((gradient.get(i,1,'p')-gradient.supp-miny),0);
        vy_end = (int)mymin((gradient.get(i,1,'p')+gradient.supp+1-miny),(def_field_temp.ny));
        vz_start = (int)mymax((gradient.get(i,2,'p')-gradient.suppz-minz),0);
        vz_end = (int)mymin((gradient.get(i,2,'p')+gradient.suppz+1-minz),(def_field_temp.nz));

        //find appropriate indexes for RBF
        bx_start = (int)(gradient.supp - ((gradient.get(i,0,'p')-minx) - vx_start));
        bx_end = (int)(gradient.supp + (vx_end - (gradient.get(i,0,'p')-minx)));
        by_start = (int)(gradient.supp - ((gradient.get(i,1,'p')-miny) - vy_start));
        by_end = (int)(gradient.supp + (vy_end - (gradient.get(i,1,'p')-miny)));
        bz_start = (int)(gradient.supp - ((gradient.get(i,2,'p')-minz) - vz_start));
        bz_end = (int)(gradient.supp + (vz_end - (gradient.get(i,2,'p')-minz)));

        cx=gradient.get(i,0,'c')-konst*gradient.get(i,0,'g');
        cy=gradient.get(i,1,'c')-konst*gradient.get(i,1,'g');
        cz=gradient.get(i,2,'c')-konst*gradient.get(i,2,'g');

        //printf("\n");
        //printf("in f_eval() V, %d, %d, %d, %d, %d, %d\n",vx_start, vx_end, vy_start, vy_end, vz_start, vz_end);
        //printf("in f_eval() B, %d, %d, %d, %d, %d, %d\n",bx_start, bx_end, by_start, by_end, bz_start, bz_end);
        //printf("in f_eval() def_field_temp size is: %d, %d, %d\n", def_field_temp.nx, def_field_temp.ny, def_field_temp.nz);
        //printf("in f_eval() RBF size is: %d, %d, %d\n", RBF.nx, RBF.ny, RBF.nz);
        //printf("\n");

        for(ii=vx_start,bx=bx_start;ii<vx_end,bx<bx_end;ii++,bx++){
            for(jj=vy_start,by=by_start;jj<vy_end,by<by_end;jj++,by++){
                for(kk=vz_start,bz=bz_start;kk<vz_end,bz<bz_end;kk++,bz++){


                    //  printf("------------------------>>>>>>>>>>>> %d,%d %d,%d %d,%d\n", ii,bx,jj,by,kk,bz);
                    new_valx = def_field_temp.get(ii,jj,kk,'x') + ( cx )*RBF.get(bx,by,bz);
                    new_valy = def_field_temp.get(ii,jj,kk,'y') + ( cy )*RBF.get(bx,by,bz);
                    new_valz = def_field_temp.get(ii,jj,kk,'z') + ( cz )*RBF.get(bx,by,bz);

                    //new_valx = def_field_temp.get(ii,jj,kk,'x') + ( -konst*gradient.get(i,0,'g') )*RBF.get(bx,by,bz);
                    //new_valy = def_field_temp.get(ii,jj,kk,'y') + ( -konst*gradient.get(i,1,'g') )*RBF.get(bx,by,bz);
                    //new_valz = def_field_temp.get(ii,jj,kk,'z') + ( -konst*gradient.get(i,2,'g') )*RBF.get(bx,by,bz);

                    //def_field_temp.set(ii-minx,jj-miny,kk-minz,'x',new_valx);
                    //def_field_temp.set(ii-minx,jj-miny,kk-minz,'y',new_valy);
                    //def_field_temp.set(ii-minx,jj-miny,kk-minz,'z',new_valz);

                    def_field_temp.set(ii,jj,kk,'x',new_valx);
                    def_field_temp.set(ii,jj,kk,'y',new_valy);
                    def_field_temp.set(ii,jj,kk,'z',new_valz);
                }
            }
        }


    }

    //printf("              Inside f_eval checkpoint 2. \n");

    //JH_temp_copy.setzero();
    //now calculate JH_temp_copy
    for(x=0;x<def_field_temp.nx;x++){
        for(y=0;y<def_field_temp.ny;y++){
            for(z=0;z<def_field_temp.nz;z++){

                src_val = source.get_intrpl(def_field_temp.get(x,y,z,'x'),def_field_temp.get(x,y,z,'y'),def_field_temp.get(x,y,z,'z'));
                trg_val = (data_type)target.get(x+minx,y+miny,z+minz);

                if ((src_val<mmin)||(trg_val<mmin)||(src_val>mmax)||(trg_val>mmax)) {

                    //do nothing

                }else{

                    jhx = (int)((((data_type)(src_val-mmin))/del)*(JH.nx-1));
                    jhy = (int)((((data_type)(trg_val-mmin))/del)*(JH.ny-1));
                    JH_temp_copy.increment(jhx,jhy);
                }

            }
        }
    }

    //	printf("              Inside f_eval checkpoint 3. \n");


    return JH_temp_copy.mi();
}




/**
 *
 * Called from within "run()" method.
 * Does line minimization along the computed gradient.
 * Tasks:
 *			Do preliminary calculations to be able to call "f_eval(...)" efficiently.
 *			Bracket the interval, use something similar to Brent's method.
 *			Do quadratic model based line minimization with four point bracketing update method.
 *			Once the "result" minimizer has been found update coefficients.
 *			Update JH by calling JH_inst_up(...) method.
 *			Return new minimum found.
 *
 *
 */

data_type2 OPTIMIZATION::line_minimize(void){

    int i;
    int ii,jj,kk;
    int vx_start,vx_end,vy_start,vy_end,vz_start,vz_end;
    int bx_start,bx_end,by_start,by_end,bz_start,bz_end;
    int bx, by, bz;
    int minx, maxx, miny, maxy, minz, maxz;
    matrix_2d JH_temp;  //the local joint histogram
    matrix_2d JH_temp_copy;
    sampling_coord def_field_temp;
    data_type new_valx,new_valy,new_valz;
    int x,y,z,jhx,jhy;
    data_type del;
    data_type src_val, trg_val;


    //first check that gradient is not zero
    data_type ssum=0;
    for(i=0;i<gradient.num_points;i++){
        ssum = ssum + myabs(gradient.get(i,0,'g')) + myabs(gradient.get(i,1,'g')) + myabs(gradient.get(i,2,'g'));
    }

    if (ssum==0)
        return 0;


    //----------------------------initial settings-----------------------------------

    JH_temp.init(nbins_src,nbins_trg);
    JH_temp_copy.init(nbins_src,nbins_trg);

    JH_temp.setzero();
    JH_temp_copy.setzero();

    //first get the boundaries of the region currently being worked on
    minx = gradient.minx;
    maxx = gradient.maxx;
    miny = gradient.miny;
    maxy = gradient.maxy;
    minz = gradient.minz;
    maxz = gradient.maxz;

    //create a temporary deformation field
    def_field_temp.init(maxx-minx+1,maxy-miny+1,maxz-minz+1);


    //now do def_field_temp = def_field (only for the current region)
    for(x=0;x<def_field_temp.nx;x++){
        for(y=0;y<def_field_temp.ny;y++){
            for(z=0;z<def_field_temp.nz;z++){

                //printf(" -------> %d, %d, %d \n", x+minx,y+miny,z+minz);
                def_field_temp.set(x,y,z,'x',def_field.get(x+minx,y+miny,z+minz,'x'));
                def_field_temp.set(x,y,z,'y',def_field.get(x+minx,y+miny,z+minz,'y'));
                def_field_temp.set(x,y,z,'z',def_field.get(x+minx,y+miny,z+minz,'z'));

            }
        }
    }

    //printf("--------------------> Inside line minimize, check point 1\n");

    //def_field_temp = def_field_temp + coeff*RBF
    for(i=0;i<gradient.num_points;i++){


        //find appropriate indexes for acessing
        vx_start = (int)mymax((gradient.get(i,0,'p')-gradient.supp),0);
        vx_end = (int)mymin((gradient.get(i,0,'p')+gradient.supp+1),(source.nx));
        vy_start = (int)mymax((gradient.get(i,1,'p')-gradient.supp),0);
        vy_end = (int)mymin((gradient.get(i,1,'p')+gradient.supp+1),(source.ny));
        vz_start = (int)mymax((gradient.get(i,2,'p')-gradient.suppz),0);
        vz_end = (int)mymin((gradient.get(i,2,'p')+gradient.suppz+1),(source.nz));

        //find appropriate indexes for RBF
        bx_start = (int)(gradient.supp - (gradient.get(i,0,'p') - vx_start));
        bx_end = (int)(gradient.supp + (vx_end - gradient.get(i,0,'p')));
        by_start = (int)(gradient.supp - (gradient.get(i,1,'p') - vy_start));
        by_end = (int)(gradient.supp + (vy_end - gradient.get(i,1,'p')));
        bz_start = (int)(gradient.supp - (gradient.get(i,2,'p') - vz_start));
        bz_end = (int)(gradient.supp + (vz_end - gradient.get(i,2,'p')));


        for(ii=vx_start,bx=bx_start;ii<vx_end,bx<bx_end;ii++,bx++){
            for(jj=vy_start,by=by_start;jj<vy_end,by<by_end;jj++,by++){
                for(kk=vz_start,bz=bz_start;kk<vz_end,bz<bz_end;kk++,bz++){

                    new_valx = def_field_temp.get(ii-minx,jj-miny,kk-minz,'x') + gradient.get(i,0,'c')*RBF.get(bx,by,bz);
                    new_valy = def_field_temp.get(ii-minx,jj-miny,kk-minz,'y') + gradient.get(i,1,'c')*RBF.get(bx,by,bz);
                    new_valz = def_field_temp.get(ii-minx,jj-miny,kk-minz,'z') + gradient.get(i,2,'c')*RBF.get(bx,by,bz);

                    def_field_temp.set(ii-minx,jj-miny,kk-minz,'x',new_valx);
                    def_field_temp.set(ii-minx,jj-miny,kk-minz,'y',new_valy);
                    def_field_temp.set(ii-minx,jj-miny,kk-minz,'z',new_valz);



                }
            }
        }
    }

    //now do JH_temp = JH - the histogram from the current region

    del = (data_type)(mmax-mmin);

    //first built JH_temp to contain the histogram of the current region
    for(x=0;x<def_field_temp.nx;x++){
        for(y=0;y<def_field_temp.ny;y++){
            for(z=0;z<def_field_temp.nz;z++){

                src_val = source.get_intrpl(def_field_temp.get(x,y,z,'x'),def_field_temp.get(x,y,z,'y'),def_field_temp.get(x,y,z,'z'));
                trg_val = (data_type)target.get(x+minx,y+miny,z+minz);

                if ((src_val<mmin)||(trg_val<mmin)||(src_val>mmax)||(trg_val>mmax)) {

                    //do nothing

                }else{

                    jhx = (int)((((data_type)(src_val-mmin))/del)*(JH.nx-1));
                    jhy = (int)((((data_type)(trg_val-mmin))/del)*(JH.ny-1));
                    JH_temp.increment(jhx,jhy);
                }

            }
        }
    }



    //printf("--------------------> Inside line minimize, check point 2\n");

    //now JH_temp = JH - JH_temp
    for(ii=0;ii<JH.nx;ii++){
        for(jj=0;jj<JH.ny;jj++){

            JH_temp.set(ii,jj,(JH.get(ii,jj)-JH_temp.get(ii,jj)));
            //if(JH_temp.get(ii,jj)>0)
            // 	printf("jh %d\n",JH_temp.get(ii,jj));

        }
    }

    //printf("In line minimize, the JH of images - current region is %f \n", JH_temp.mi());

    //now do def_field_temp = def_field (only for the current region) again
    for(x=0;x<def_field_temp.nx;x++){
        for(y=0;y<def_field_temp.ny;y++){
            for(z=0;z<def_field_temp.nz;z++){


                def_field_temp.set(x,y,z,'x',def_field.get(x+minx,y+miny,z+minz,'x'));
                def_field_temp.set(x,y,z,'y',def_field.get(x+minx,y+miny,z+minz,'y'));
                def_field_temp.set(x,y,z,'z',def_field.get(x+minx,y+miny,z+minz,'z'));

            }
        }
    }

    //printf("--------------------> Inside line minimize, check point 3\n");


    //-----------------------bracket the interval---------------------------------


    int STOP=0;
    int count=1;
    //data_type br_step=br_STEP*(curr_level/5.0+1);
    data_type br_step=br_STEP;
    data_type konst=0;
    data_type konst_min = 0;
    data_type2 mi_min;
    data_type konst_ini = 0;
    data_type2 mi_ini;
    data_type konst_last=0;
    data_type2 mi_last;
    data_type konst_curr=0;
    data_type2 mi_curr;
    data_type VIOL=0;

    mi_ini = f_eval(0,JH_temp,JH_temp_copy,def_field_temp,minx,maxx,miny,maxy,minz,maxz);

    //printf("--------------------> Inside line minimize, check point 4\n");

    mi_last = mi_ini;
    mi_min = mi_ini;

    //for(i=0;i<gradient.num_points;i++){
    //	printf("    %f,%f,%f,\n",gradient.get(i,0,'g'),gradient.get(i,1,'g'),gradient.get(i,2,'g'));
    //}

    while(STOP==0){

        konst_curr = count*count*br_step;
        mi_curr = f_eval(konst_curr,JH_temp,JH_temp_copy,def_field_temp,minx,maxx,miny,maxy,minz,maxz);

        //printf("                 Bracketing print outs : : %f,%f\n",konst_curr,mi_curr);
        if(mi_curr<mi_min){
            konst_min=konst_curr;
            mi_min=mi_curr;
        }

        if((mi_curr-mi_min)>0)//used to be br_EPS
            STOP=1;

        //if(count>20)
        // STOP=1;

        if(konst_curr >= gradient.min_a){
            VIOL=1;
            //printf("     JACOBIAN THRESHOLD VIOLATION STOPPED\n");
            STOP=1;
        }

        mi_last=mi_curr;
        konst_last=konst_curr;


        count++;
    }

    //correct for a particular situation from the above loop
    if(konst_min==konst_ini){
        konst_min = (konst_last-konst_ini)/((data_type)2.0);
        mi_min = f_eval(konst_min,JH_temp,JH_temp_copy,def_field_temp,minx,maxx,miny,maxy,minz,maxz);
    }


    //test printing
    //printf("\n");
    //printf("                    LINE MINIMIZATION        \n");
    //printf("                 Bracketing print outs : : %f,%f\n",konst_curr,mi_curr);
    //printf("                       %f,%f  |  %f,%f  |  %f,%f\n",konst_ini,mi_ini,konst_min,mi_min,konst_last,mi_last);
    //printf("                       %f,%f\n",konst_min,mi_min);
    //printf("                       %f,%f\n",konst_last,mi_last);
    //printf("\n");

    //----------------------minimize within interval------------------------------
    //use quadratic model allyed with four point bracketing update

    data_type2 result, mi_result;
    data_type2 a,b,c,fa,fb,fc,d,fd,num,den;
    //result = konst_min;

    //first check that function has been bracketed correctly and is decreasing
    if ( (mi_min<mi_ini) && (mi_min<mi_last) && (VIOL<1) ) {
        //do line minimization
        int count=0;
        STOP=0;
        while (STOP==0) {

            a = konst_ini;
            b = konst_last;
            c = konst_min;
            fa = mi_ini;
            fb = mi_last;
            fc = mi_min;
            result = (data_type)c;
            mi_result = (data_type)fc;

            //the quadratic model minimization
            num = fa*(b*b-c*c) - fc*(b*b-a*a) + fb*(c*c-a*a);
            den = fa*(b-c) - fc*(b-a) + fb*(c-a);
            d = (data_type)0.5*(num/den);
            fd = f_eval((data_type)d,JH_temp,JH_temp_copy,def_field_temp,minx,maxx,miny,maxy,minz,maxz);



            //four point bracketing update
            if ( (fd>fa)||(fd>fb) ) {
                break;
            } else {

                if (fc<fd) {
                    konst_ini = a;
                    mi_ini = fa;
                    konst_last = d;
                    mi_last = fd;
                    konst_min = c;
                    mi_min = fc;
                } else if (fc>=fd) {
                    konst_ini = c;
                    mi_ini = fc;
                    konst_last = b;
                    mi_last = fb;
                    konst_min = d;
                    mi_min = fd;
                }
                //printf("                       %f,%f  |  %f,%f  |  %f,%f\n",konst_ini,mi_ini,konst_min,mi_min,konst_last,mi_last);

                if ( (myabs(mi_min-mi_result)<min_EPS)||count>10 )
                    STOP=1;

            }

            count++;
        }//of while loop

    } else {
        //something went bad
        result = konst_min;
        mi_result = mi_min;
    }

    //make updating calls
    gradient.coeff_update(result);
    //printf("updating the JH ... \n");
    //JH_inst1();
    JH_inst_up(JH_temp,def_field_temp,minx,maxx,miny,maxy,minz,maxz);

    //deallocation

    //JH_temp.setzero();
    //JH_temp_copy.setzero();

    JH_temp.destroy();
    JH_temp_copy.destroy();
    def_field_temp.destroy();


    return mi_result;

}




/**
 *
 * Called from within "run()" method.
 * Updates the deformation field based on the new coefficients found for current control points
 *
 *
 */

void OPTIMIZATION::def_field_update(void){

    int i;
    int ii,jj,kk;
    int vx_start,vx_end,vy_start,vy_end,vz_start,vz_end;
    int bx_start,bx_end,by_start,by_end,bz_start,bz_end;

    int bx, by, bz;
    data_type new_valx, new_valy, new_valz;

    //printf("Updating deformation field ... \n");

    int vxmin,vxmax,vymin,vymax,vzmin,vzmax;
    vxmin = source.nx-1;
    vxmax = 0;
    vymin = source.ny-1;
    vymax = 0;
    vzmin = source.nz-1;
    vzmax = 0;

    //loop through each point in instance gradient
    for (i=0;i<gradient.num_points;i++){

        //find appropriate indexes for acessing
        vx_start = (int)mymax((gradient.get(i,0,'p')-gradient.supp),0);
        vx_end = (int)mymin((gradient.get(i,0,'p')+gradient.supp+1),(source.nx));
        vy_start = (int)mymax((gradient.get(i,1,'p')-gradient.supp),0);
        vy_end = (int)mymin((gradient.get(i,1,'p')+gradient.supp+1),(source.ny));
        vz_start = (int)mymax((gradient.get(i,2,'p')-gradient.suppz),0);
        vz_end = (int)mymin((gradient.get(i,2,'p')+gradient.suppz+1),(source.nz));

        if(vx_start<vxmin)
            vxmin = vx_start;
        if(vx_end>vxmax)
            vxmax = vx_end;

        if(vy_start<vymin)
            vymin = vy_start;
        if(vy_end>vymax)
            vymax = vy_end;

        if(vz_start<vzmin)
            vzmin = vz_start;
        if(vz_end>vzmax)
            vzmax = vz_end;

        //find appropriate indexes for RBF
        bx_start = (int)(gradient.supp - (gradient.get(i,0,'p') - vx_start));
        bx_end = (int)(gradient.supp + (vx_end - gradient.get(i,0,'p')));
        by_start = (int)(gradient.supp - (gradient.get(i,1,'p') - vy_start));
        by_end = (int)(gradient.supp + (vy_end - gradient.get(i,1,'p')));
        bz_start = (int)(gradient.supp - (gradient.get(i,2,'p') - vz_start));
        bz_end = (int)(gradient.supp + (vz_end - gradient.get(i,2,'p')));


        for(ii=vx_start,bx=bx_start;ii<vx_end,bx<bx_end;ii++,bx++){
            for(jj=vy_start,by=by_start;jj<vy_end,by<by_end;jj++,by++){
                for(kk=vz_start,bz=bz_start;kk<vz_end,bz<bz_end;kk++,bz++){

                    new_valx = def_field.get(ii,jj,kk,'x') + gradient.get(i,0,'c')*RBF.get(bx,by,bz);
                    new_valy = def_field.get(ii,jj,kk,'y') + gradient.get(i,1,'c')*RBF.get(bx,by,bz);
                    new_valz = def_field.get(ii,jj,kk,'z') + gradient.get(i,2,'c')*RBF.get(bx,by,bz);

                    def_field.set(ii,jj,kk,'x',new_valx);
                    def_field.set(ii,jj,kk,'y',new_valy);
                    def_field.set(ii,jj,kk,'z',new_valz);

                }
            }
        }


    }

    //printf(" Zone %d,%d %d,%d %d,%d \n", vxmin,vxmax,vymin,vymax,vzmin,vzmax);

}



void OPTIMIZATION::def_field_update2(gradient_class *p){

    int i;
    int ii,jj,kk;
    int vx_start,vx_end,vy_start,vy_end,vz_start,vz_end;
    int bx_start,bx_end,by_start,by_end,bz_start,bz_end;

    int bx, by, bz;
    data_type new_valx, new_valy, new_valz;

    //printf("Updating deformation field ... \n");

    int vxmin,vxmax,vymin,vymax,vzmin,vzmax;
    vxmin = source.nx-1;
    vxmax = 0;
    vymin = source.ny-1;
    vymax = 0;
    vzmin = source.nz-1;
    vzmax = 0;

    //loop through each point in instance gradient
    for (i=0;i<p->num_points;i++){

        //find appropriate indexes for acessing
        vx_start = (int)mymax((p->get(i,0,'p')-p->supp),0);
        vx_end = (int)mymin((p->get(i,0,'p')+p->supp+1),(source.nx));
        vy_start = (int)mymax((p->get(i,1,'p')-p->supp),0);
        vy_end = (int)mymin((p->get(i,1,'p')+p->supp+1),(source.ny));
        vz_start = (int)mymax((p->get(i,2,'p')-p->suppz),0);
        vz_end = (int)mymin((p->get(i,2,'p')+p->suppz+1),(source.nz));

        if(vx_start<vxmin)
            vxmin = vx_start;
        if(vx_end>vxmax)
            vxmax = vx_end;

        if(vy_start<vymin)
            vymin = vy_start;
        if(vy_end>vymax)
            vymax = vy_end;

        if(vz_start<vzmin)
            vzmin = vz_start;
        if(vz_end>vzmax)
            vzmax = vz_end;

        //find appropriate indexes for RBF
        bx_start = (int)(p->supp - (p->get(i,0,'p') - vx_start));
        bx_end = (int)(p->supp + (vx_end - p->get(i,0,'p')));
        by_start = (int)(p->supp - (p->get(i,1,'p') - vy_start));
        by_end = (int)(p->supp + (vy_end - p->get(i,1,'p')));
        bz_start = (int)(p->supp - (p->get(i,2,'p') - vz_start));
        bz_end = (int)(p->supp + (vz_end - p->get(i,2,'p')));


        for(ii=vx_start,bx=bx_start;ii<vx_end,bx<bx_end;ii++,bx++){
            for(jj=vy_start,by=by_start;jj<vy_end,by<by_end;jj++,by++){
                for(kk=vz_start,bz=bz_start;kk<vz_end,bz<bz_end;kk++,bz++){

                    new_valx = def_field.get(ii,jj,kk,'x') + p->get(i,0,'c')*RBF.get(bx,by,bz);
                    new_valy = def_field.get(ii,jj,kk,'y') + p->get(i,1,'c')*RBF.get(bx,by,bz);
                    new_valz = def_field.get(ii,jj,kk,'z') + p->get(i,2,'c')*RBF.get(bx,by,bz);

                    def_field.set(ii,jj,kk,'x',new_valx);
                    def_field.set(ii,jj,kk,'y',new_valy);
                    def_field.set(ii,jj,kk,'z',new_valz);

                }
            }
        }


    }

    //printf(" Zone %d,%d %d,%d %d,%d \n", vxmin,vxmax,vymin,vymax,vzmin,vzmax);

}



void OPTIMIZATION::grd_multi(void){
    int number=gradient.num_points;
    int core=core_number;
    called_grd_build **p;
    p =(called_grd_build **)malloc(number*sizeof(called_grd_build *));
    int gap=(number+core-1)/core;
    int left=0;
    int right=gap;
    data_type fpp=JH.mi();
    for (int i=0;i<core;i++){
        p[i]=new called_grd_build;
        p[i]->init(&gradient,left,right,&gradient,&source,&target,&def_field,&RBF,&JH);
        p[i]->X_opt=X_opt;
        p[i]->Y_opt=Y_opt;
        p[i]->Z_opt=Z_opt;
        p[i]->nbins_src=nbins_src;
        p[i]->nbins_trg=nbins_trg;
        p[i]->mmax=mmax;
        p[i]->mmin=mmin;
        p[i]->grd_STEP=grd_STEP;
        p[i]->fpp=fpp;
        p[i]->start();
        left+=gap;
        right+=gap;
        right=mymin(right,number);
    }
    for (int i=0;i<core;i++){
        p[i]->wait();
    }

    gradient.normalize();
    gradient.max_d_compute(s_trsh);
}



/**
 *
 * Called from within "run()" method.
 * Computes the current gradient (for the current control points);
 * Long and complicated function. Tries to maximize use of locality through the
 *			compactly supported radial basis funcion.
 *
 * See code for actual detail.
 *
 */

void OPTIMIZATION::grd(void){

    data_type STEP = gradient.supp*grd_STEP/8.0; //this is the h in (f(x+h)-f(x-h))/(2*h)
    data_type mic;
    int i;
    int ii,jj,kk,iii,jjj;
    int jhx, jhy;
    data_type del;
    data_type fmh, fph, fpp;
    data_type src_val, trg_val;
    matrix_2d JH_temp;  //the local joint histogram
    matrix_2d JH_i;  //JH- local joint histogram
    int vx_start,vx_end,vy_start,vy_end,vz_start,vz_end;
    int bx_start,bx_end,by_start,by_end,bz_start,bz_end;
    int bx, by, bz;

    JH_temp.init(nbins_src,nbins_trg);
    JH_i.init(nbins_src,nbins_trg);

    fpp = JH.mi();

    JH_temp.setzero();
    JH_i.setzero();


    //printf("Computing gradient ... \n");


    //printf( "number of points %d", gradient.num_points);
    //printf("\n");

    del = (data_type)(mmax-mmin);

    //loop through each point in instance gradient
    for (i=0;i<gradient.num_points;i++){

        //printf("Point %d out of %d .\n",i,gradient.num_points);
        //printf("                  In grd(), point %d out of %d.\n",i, gradient.num_points);

        JH_temp.setzero();
        JH_i.setzero();



        //subtract the current subportion of the volume from the current histogram

        //for that build one first

        //find appropriate indexes for acessing
        vx_start = (int)mymax((gradient.get(i,0,'p')-gradient.supp),0);
        vx_end = (int)mymin((gradient.get(i,0,'p')+gradient.supp),(source.nx));
        vy_start = (int)mymax((gradient.get(i,1,'p')-gradient.supp),0);
        vy_end = (int)mymin((gradient.get(i,1,'p')+gradient.supp),(source.ny));
        vz_start = (int)mymax((gradient.get(i,2,'p')-gradient.suppz),0);
        vz_end = (int)mymin((gradient.get(i,2,'p')+gradient.suppz),(source.nz));

        //find appropriate indexes for RBF
        bx_start = (int)(gradient.supp - (gradient.get(i,0,'p') - vx_start));
        bx_end = (int)(gradient.supp + (vx_end - gradient.get(i,0,'p')));

        by_start = (int)(gradient.supp - (gradient.get(i,1,'p') - vy_start));
        by_end = (int)(gradient.supp + (vy_end - gradient.get(i,1,'p')));

        bz_start = (int)(gradient.supp - (gradient.get(i,2,'p') - vz_start));
        bz_end = (int)(gradient.supp + (vz_end - gradient.get(i,2,'p')));

        bx = bx_start;
        by = by_start;
        bz = bz_start;




        for(ii=vx_start,bx=bx_start;ii<vx_end,bx<bx_end;ii++,bx++){
            for(jj=vy_start,by=by_start;jj<vy_end,by<by_end;jj++,by++){
                for(kk=vz_start,bz=bz_start;kk<vz_end,bz<bz_end;kk++,bz++){

                    src_val = source.get_intrpl((def_field.get(ii,jj,kk,'x')+gradient.get(i,0,'c')*RBF.get(bx,by,bz)),(def_field.get(ii,jj,kk,'y')+gradient.get(i,1,'c')*RBF.get(bx,by,bz)),(def_field.get(ii,jj,kk,'z')+gradient.get(i,2,'c')*RBF.get(bx,by,bz)));
                    trg_val = (data_type)target.get(ii,jj,kk);

                    if ((src_val<mmin)||(trg_val<mmin)||(src_val>mmax)||(trg_val>mmax)) {
                        //do nothing
                    }else{

                        jhx = (int)((((data_type)(src_val-mmin))/del)*(JH.nx-1));
                        jhy = (int)((((data_type)(trg_val-mmin))/del)*(JH.ny-1));

                        JH_temp.increment(jhx,jhy);


                    }

                }
            }
        }


        //do subtraction
        for(ii=0;ii<JH.nx;ii++){
            for(jj=0;jj<JH.ny;jj++){

                JH_i.set(ii,jj,JH.get(ii,jj)-JH_temp.get(ii,jj));

                //put JH_i into JH_temp for later reuse
                //JH_temp.set(ii,jj,(JH_i.get(ii,jj)));
            }
        }


        if (X_opt==1){
            //-------------- X direction gradient component ---------------//
            //printf("......   x\n");

            //printf("xdirection gradient: \n");
            //reset JH_temp to equal JH_inst
            for(iii=0;iii<JH.nx;iii++){
                for(jjj=0;jjj<JH.ny;jjj++){
                    JH_temp.set(iii,jjj,(JH_i.get(iii,jjj)));
                }
            }

            //build new "TEST" histogram (reuse JH_temp)

            for(ii=vx_start,bx=bx_start;ii<vx_end,bx<bx_end;ii++,bx++){
                for(jj=vy_start,by=by_start;jj<vy_end,by<by_end;jj++,by++){
                    for(kk=vz_start,bz=bz_start;kk<vz_end,bz<bz_end;kk++,bz++){

                        //would need to generate the new local test deformation field based on the test coefficients
                        //then apply it to the image
                        //then build the new local 2d joint histogram
                        //add joint histogram to the image
                        //compute the mutual information

                        src_val = source.get_intrpl(  ((gradient.get(i,0,'c')+STEP)*RBF.get(bx,by,bz) + def_field.get(ii,jj,kk,'x')) , ((gradient.get(i,1,'c'))*RBF.get(bx,by,bz) + def_field.get(ii,jj,kk,'y')) , ((gradient.get(i,2,'c'))*RBF.get(bx,by,bz) + def_field.get(ii,jj,kk,'z')));
                        trg_val = (data_type)target.get(ii,jj,kk);

                        if ((src_val<mmin)||(trg_val<mmin)||(src_val>mmax)||(trg_val>mmax)) {
                            //do nothing
                        }else{
                            jhx = (int)((((data_type)(src_val-mmin))/del)*(JH.nx-1));
                            jhy = (int)((((data_type)(trg_val-mmin))/del)*(JH.ny-1));
                            JH_temp.increment(jhx,jhy);

                        }


                    }
                }
            }


            fph = JH_temp.mi();


            //repeat procedure for f(x-h)

            //reset JH_temp to equal JH_inst
            for(iii=0;iii<JH.nx;iii++){
                for(jjj=0;jjj<JH.ny;jjj++){
                    JH_temp.set(iii,jjj,(JH_i.get(iii,jjj)));
                }
            }


            for(ii=vx_start,bx=bx_start;ii<vx_end,bx<bx_end;ii++,bx++){
                for(jj=vy_start,by=by_start;jj<vy_end,by<by_end;jj++,by++){
                    for(kk=vz_start,bz=bz_start;kk<vz_end,bz<bz_end;kk++,bz++){


                        src_val = source.get_intrpl( ((gradient.get(i,0,'c')-STEP)*RBF.get(bx,by,bz) + def_field.get(ii,jj,kk,'x')) , ((gradient.get(i,1,'c'))*RBF.get(bx,by,bz) + def_field.get(ii,jj,kk,'y')) , ((gradient.get(i,2,'c'))*RBF.get(bx,by,bz) + def_field.get(ii,jj,kk,'z')));
                        trg_val = (data_type)target.get(ii,jj,kk);
                        if ((src_val<mmin)||(trg_val<mmin)||(src_val>mmax)||(trg_val>mmax)) {
                            //do nothing
                        }else{
                            jhx = (int)((((data_type)(src_val-mmin))/del)*(JH.nx-1));
                            jhy = (int)((((data_type)(trg_val-mmin))/del)*(JH.ny-1));

                            JH_temp.increment(jhx,jhy);
                        }


                    }
                }
            }

            fmh = JH_temp.mi();
            /*if( sign(fph-fpp)==sign(fpp-fmh) ){
     mic = (fph-fmh)/(2*STEP);
     }else{
     mic = 0;
     }*/
            mic = (fph-fmh)/(2*STEP);
            //mic = (fph-fpp)/(STEP);
        }else{
            mic = 0;
        }
        gradient.set(i,0,mic,'g');



        //-------------- Y direction gradient component ---------------//
        if (Y_opt==1){
            //printf("......  y\n");

            //printf("ydirection gradient: \n");
            //reset JH_temp to equal JH_inst
            for(iii=0;iii<JH.nx;iii++){
                for(jjj=0;jjj<JH.ny;jjj++){
                    JH_temp.set(iii,jjj,(JH_i.get(iii,jjj)));
                }
            }

            //build new "TEST" histogram (reuse JH_temp)

            for(ii=vx_start,bx=bx_start;ii<vx_end,bx<bx_end;ii++,bx++){
                for(jj=vy_start,by=by_start;jj<vy_end,by<by_end;jj++,by++){
                    for(kk=vz_start,bz=bz_start;kk<vz_end,bz<bz_end;kk++,bz++){

                        //would need to generate the new local test deformation field based on the test coefficients
                        //then apply it to the image
                        //then build the new local 2d joint histogram
                        //add joint histogram to the image
                        //compute the mutual information

                        src_val = source.get_intrpl(  ((gradient.get(i,0,'c'))*RBF.get(bx,by,bz) + def_field.get(ii,jj,kk,'x')) , ((gradient.get(i,1,'c')+STEP)*RBF.get(bx,by,bz) + def_field.get(ii,jj,kk,'y')) , ((gradient.get(i,2,'c'))*RBF.get(bx,by,bz) + def_field.get(ii,jj,kk,'z')));
                        trg_val = (data_type)target.get(ii,jj,kk);

                        if ((src_val<mmin)||(trg_val<mmin)||(src_val>mmax)||(trg_val>mmax)) {
                            //do nothing
                        }else{
                            jhx = (int)((((data_type)(src_val-mmin))/del)*(JH.nx-1));
                            jhy = (int)((((data_type)(trg_val-mmin))/del)*(JH.ny-1));

                            JH_temp.increment(jhx,jhy);
                        }

                    }
                }
            }


            fph = JH_temp.mi();


            //repeat procedure for f(x-h)

            //reset JH_temp to equal JH_inst
            for(iii=0;iii<JH.nx;iii++){
                for(jjj=0;jjj<JH.ny;jjj++){
                    JH_temp.set(iii,jjj,(JH_i.get(iii,jjj)));
                }
            }


            for(ii=vx_start,bx=bx_start;ii<vx_end,bx<bx_end;ii++,bx++){
                for(jj=vy_start,by=by_start;jj<vy_end,by<by_end;jj++,by++){
                    for(kk=vz_start,bz=bz_start;kk<vz_end,bz<bz_end;kk++,bz++){


                        src_val = source.get_intrpl( ((gradient.get(i,0,'c'))*RBF.get(bx,by,bz) + def_field.get(ii,jj,kk,'x')) , ((gradient.get(i,1,'c')-STEP)*RBF.get(bx,by,bz) + def_field.get(ii,jj,kk,'y')) , ((gradient.get(i,2,'c'))*RBF.get(bx,by,bz) + def_field.get(ii,jj,kk,'z')));
                        trg_val = (data_type)target.get(ii,jj,kk);
                        if ((src_val<mmin)||(trg_val<mmin)||(src_val>mmax)||(trg_val>mmax)) {
                            //do nothing
                        }else{
                            jhx = (int)((((data_type)(src_val-mmin))/del)*(JH.nx-1));
                            jhy = (int)((((data_type)(trg_val-mmin))/del)*(JH.ny-1));

                            JH_temp.increment(jhx,jhy);
                        }


                    }
                }
            }

            fmh = JH_temp.mi();

            /*if( sign(fph-fpp)==sign(fpp-fmh) ){
    mic = (fph-fmh)/(2*STEP);
    }else{
    mic = 0;
    }*/

            mic = (fph-fmh)/(2*STEP);
            //mic = (fph-fpp)/(STEP);
        }else{
            mic = 0;
        }
        gradient.set(i,1,mic,'g');



        //-------------- Z direction gradient component ---------------//
        if (Z_opt==1){
            //printf("......  z \n");

            //printf("zdirection gradient: \n");
            //reset JH_temp to equal JH_inst
            for(iii=0;iii<JH.nx;iii++){
                for(jjj=0;jjj<JH.ny;jjj++){
                    JH_temp.set(iii,jjj,(JH_i.get(iii,jjj)));
                }
            }

            //build new "TEST" histogram (reuse JH_temp)

            for(ii=vx_start,bx=bx_start;ii<vx_end,bx<bx_end;ii++,bx++){
                for(jj=vy_start,by=by_start;jj<vy_end,by<by_end;jj++,by++){
                    for(kk=vz_start,bz=bz_start;kk<vz_end,bz<bz_end;kk++,bz++){

                        //would need to generate the new local test deformation field based on the test coefficients
                        //then apply it to the image
                        //then build the new local 2d joint histogram
                        //add joint histogram to the image
                        //compute the mutual information

                        src_val = source.get_intrpl(  ((gradient.get(i,0,'c'))*RBF.get(bx,by,bz) + def_field.get(ii,jj,kk,'x')) , ((gradient.get(i,1,'c'))*RBF.get(bx,by,bz) + def_field.get(ii,jj,kk,'y')) , ((gradient.get(i,2,'c')+STEP)*RBF.get(bx,by,bz) + def_field.get(ii,jj,kk,'z')));
                        trg_val = (data_type)target.get(ii,jj,kk);

                        if ((src_val<mmin)||(trg_val<mmin)||(src_val>mmax)||(trg_val>mmax)) {
                            //do nothing
                        }else{
                            jhx = (int)((((data_type)(src_val-mmin))/del)*(JH.nx-1));
                            jhy = (int)((((data_type)(trg_val-mmin))/del)*(JH.ny-1));

                            JH_temp.increment(jhx,jhy);
                        }

                    }
                }
            }


            fph = JH_temp.mi();


            //repeat procedure for f(x-h)

            //reset JH_temp to equal JH_inst
            for(iii=0;iii<JH.nx;iii++){
                for(jjj=0;jjj<JH.ny;jjj++){
                    JH_temp.set(iii,jjj,(JH_i.get(iii,jjj)));
                }
            }


            for(ii=vx_start,bx=bx_start;ii<vx_end,bx<bx_end;ii++,bx++){
                for(jj=vy_start,by=by_start;jj<vy_end,by<by_end;jj++,by++){
                    for(kk=vz_start,bz=bz_start;kk<vz_end,bz<bz_end;kk++,bz++){


                        src_val = source.get_intrpl( ((gradient.get(i,0,'c'))*RBF.get(bx,by,bz) + def_field.get(ii,jj,kk,'x')) , ((gradient.get(i,1,'c'))*RBF.get(bx,by,bz) + def_field.get(ii,jj,kk,'y')) , ((gradient.get(i,2,'c')-STEP)*RBF.get(bx,by,bz) + def_field.get(ii,jj,kk,'z')));
                        trg_val = (data_type)target.get(ii,jj,kk);
                        if ((src_val<mmin)||(trg_val<mmin)||(src_val>mmax)||(trg_val>mmax)) {
                            //do nothing
                        }else{
                            jhx = (int)((((data_type)(src_val-mmin))/del)*(JH.nx-1));
                            jhy = (int)((((data_type)(trg_val-mmin))/del)*(JH.ny-1));

                            JH_temp.increment(jhx,jhy);
                        }


                    }
                }
            }

            fmh = JH_temp.mi();

            /*if( sign(fph-fpp)==sign(fpp-fmh) ){
    mic = (fph-fmh)/(2*STEP);
    }else{
    mic = 0;
    }*/

            mic = (fph-fmh)/(2*STEP);

            //mic = (fph-fpp)/(STEP);
        }else{
            mic = 0;
        }
        gradient.set(i,2,mic,'g');

        //done with all three components of gradient for this control point
    }

    gradient.normalize();
    gradient.max_d_compute(s_trsh);
    JH_temp.destroy();
    JH_i.destroy();
}




/*
 *
 * Called from within "run()" method.
 * Not a complicated function.
 * Applyes the current deformation field to the source image and stores the result in the
 *			source image itself.
 *
 */

void OPTIMIZATION::generate_result(void){

    int ii,jj,kk;
    data_type new_val;
    volume outvol;

    char str_result[] = "result.vol";
    outvol.init(source.nx,source.ny,source.nz,0,str_result);

    for(ii=0;ii<source.nx;ii++){
        for(jj=0;jj<source.ny;jj++){
            for(kk=0;kk<source.nz;kk++){

                new_val = source.get_intrpl(def_field.get(ii,jj,kk,'x'),def_field.get(ii,jj,kk,'y'),def_field.get(ii,jj,kk,'z'));
                outvol.set_Rm(ii,jj,kk,new_val);

            }
        }
    }

    for(ii=0;ii<source.nx;ii++){
        for(jj=0;jj<source.ny;jj++){
            for(kk=0;kk<source.nz;kk++){

                source.set(ii,jj,kk,(data_type)outvol.get_Rm(ii,jj,kk));

            }
        }
    }
    outvol.destroy_Rm();

}




/**
 *
 * One of the public methods.
 * Called to output results in a file.
 *
 */

void OPTIMIZATION::output(char *fn){

    source.filename = fn;
    source.vol_write();

    if (SAVE_DEF_FIELD==1)
        def_field.save(fn);


}




/**
 * One of the public methods
 * Called to perform clean up
 */

void OPTIMIZATION::destroy(){

    JH.destroy();
    source.destroy();
    source.destroy_Rm();
    target.destroy();
    target.destroy_Rm();
    //bin_vol.destroy();
    //bin_vol.destroy_Rm();
    //gradient.destroy();
    RBF.destroy();

    //
    def_field.destroy();
    //

    delete [] ctrpts_array;
    delete [] ctrptsZ_array;
    delete [] res_array;

}
