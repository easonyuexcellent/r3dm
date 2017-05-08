#include "called_grd_build.h"

void called_grd_build::init(gradient_class *_grd, int _left, int _right,gradient_class *_gradient,volume *_source,volume *_target,sampling_coord *_def_field,basis_function *_RBF,matrix_2d *_JH){
    gradient=_grd;
    left=_left;
    right=_right;
    gradient=_gradient;
    source=_source;
    target=_target;
    def_field=_def_field;
    RBF=_RBF;
    JH=_JH;
}

called_grd_build::called_grd_build()
{

}

void called_grd_build::run(){


    data_type STEP = gradient->supp*grd_STEP/8.0; //this is the h in (f(x+h)-f(x-h))/(2*h)
    data_type mic;
    int i;
    int ii,jj,kk,iii,jjj;
    int jhx, jhy;
    data_type del;
    data_type fmh, fph;
    data_type src_val, trg_val;
    matrix_2d JH_temp;  //the local joint histogram
    matrix_2d JH_i;  //JH- local joint histogram
    int vx_start,vx_end,vy_start,vy_end,vz_start,vz_end;
    int bx_start,bx_end,by_start,by_end,bz_start,bz_end;
    int bx, by, bz;

    JH_temp.init(nbins_src,nbins_trg);
    JH_i.init(nbins_src,nbins_trg);

    JH_temp.setzero();
    JH_i.setzero();
    del = (data_type)(mmax-mmin);

    //loop through each point in instance gradient
    for (i=left;i<right;i++){

        //printf("Point %d out of %d .\n",i,gradient->num_points);
        //printf("                  In grd(), point %d out of %d.\n",i, gradient->num_points);

        JH_temp.setzero();
        JH_i.setzero();



        //subtract the current subportion of the volume from the current histogram

        //for that build one first

        //find appropriate indexes for acessing
        vx_start = (int)mymax((gradient->get(i,0,'p')-gradient->supp),0);
        vx_end = (int)mymin((gradient->get(i,0,'p')+gradient->supp),(source->nx));
        vy_start = (int)mymax((gradient->get(i,1,'p')-gradient->supp),0);
        vy_end = (int)mymin((gradient->get(i,1,'p')+gradient->supp),(source->ny));
        vz_start = (int)mymax((gradient->get(i,2,'p')-gradient->suppz),0);
        vz_end = (int)mymin((gradient->get(i,2,'p')+gradient->suppz),(source->nz));

        //find appropriate indexes for RBF
        bx_start = (int)(gradient->supp - (gradient->get(i,0,'p') - vx_start));
        bx_end = (int)(gradient->supp + (vx_end - gradient->get(i,0,'p')));

        by_start = (int)(gradient->supp - (gradient->get(i,1,'p') - vy_start));
        by_end = (int)(gradient->supp + (vy_end - gradient->get(i,1,'p')));

        bz_start = (int)(gradient->supp - (gradient->get(i,2,'p') - vz_start));
        bz_end = (int)(gradient->supp + (vz_end - gradient->get(i,2,'p')));

        bx = bx_start;
        by = by_start;
        bz = bz_start;




        for(ii=vx_start,bx=bx_start;ii<vx_end,bx<bx_end;ii++,bx++){
            for(jj=vy_start,by=by_start;jj<vy_end,by<by_end;jj++,by++){
                for(kk=vz_start,bz=bz_start;kk<vz_end,bz<bz_end;kk++,bz++){

                    src_val = source->get_intrpl((def_field->get(ii,jj,kk,'x')+gradient->get(i,0,'c')*RBF->get(bx,by,bz)),(def_field->get(ii,jj,kk,'y')+gradient->get(i,1,'c')*RBF->get(bx,by,bz)),(def_field->get(ii,jj,kk,'z')+gradient->get(i,2,'c')*RBF->get(bx,by,bz)));
                    trg_val = (data_type)target->get(ii,jj,kk);

                    if ((src_val<mmin)||(trg_val<mmin)||(src_val>mmax)||(trg_val>mmax)) {
                        //do nothing
                    }else{
                        jhx = (int)((((data_type)(src_val-mmin))/del)*(JH->nx-1));
                        jhy = (int)((((data_type)(trg_val-mmin))/del)*(JH->ny-1));

                        JH_temp.increment(jhx,jhy);

                    }

                }
            }
        }


        //do subtraction
        for(ii=0;ii<JH->nx;ii++){
            for(jj=0;jj<JH->ny;jj++){

                JH_i.set(ii,jj,JH->get(ii,jj)-JH_temp.get(ii,jj));

                //put JH_i into JH_temp for later reuse
                //JH_temp.set(ii,jj,(JH_i.get(ii,jj)));
            }
        }


        if (X_opt==1){
            //-------------- X direction gradient component ---------------//
            //printf("......   x\n");

            //printf("xdirection gradient: \n");
            //reset JH_temp to equal JH_inst
            for(iii=0;iii<JH->nx;iii++){
                for(jjj=0;jjj<JH->ny;jjj++){
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

                        src_val = source->get_intrpl(  ((gradient->get(i,0,'c')+STEP)*RBF->get(bx,by,bz) + def_field->get(ii,jj,kk,'x')) , ((gradient->get(i,1,'c'))*RBF->get(bx,by,bz) + def_field->get(ii,jj,kk,'y')) , ((gradient->get(i,2,'c'))*RBF->get(bx,by,bz) + def_field->get(ii,jj,kk,'z')));
                        trg_val = (data_type)target->get(ii,jj,kk);

                        if ((src_val<mmin)||(trg_val<mmin)||(src_val>mmax)||(trg_val>mmax)) {
                            //do nothing
                        }else{
                            jhx = (int)((((data_type)(src_val-mmin))/del)*(JH->nx-1));
                            jhy = (int)((((data_type)(trg_val-mmin))/del)*(JH->ny-1));
                            JH_temp.increment(jhx,jhy);

                        }


                    }
                }
            }


            fph = JH_temp.mi();


            //repeat procedure for f(x-h)

            //reset JH_temp to equal JH_inst
            for(iii=0;iii<JH->nx;iii++){
                for(jjj=0;jjj<JH->ny;jjj++){
                    JH_temp.set(iii,jjj,(JH_i.get(iii,jjj)));
                }
            }


            for(ii=vx_start,bx=bx_start;ii<vx_end,bx<bx_end;ii++,bx++){
                for(jj=vy_start,by=by_start;jj<vy_end,by<by_end;jj++,by++){
                    for(kk=vz_start,bz=bz_start;kk<vz_end,bz<bz_end;kk++,bz++){


                        src_val = source->get_intrpl( ((gradient->get(i,0,'c')-STEP)*RBF->get(bx,by,bz) + def_field->get(ii,jj,kk,'x')) , ((gradient->get(i,1,'c'))*RBF->get(bx,by,bz) + def_field->get(ii,jj,kk,'y')) , ((gradient->get(i,2,'c'))*RBF->get(bx,by,bz) + def_field->get(ii,jj,kk,'z')));
                        trg_val = (data_type)target->get(ii,jj,kk);
                        if ((src_val<mmin)||(trg_val<mmin)||(src_val>mmax)||(trg_val>mmax)) {
                            //do nothing
                        }else{
                            jhx = (int)((((data_type)(src_val-mmin))/del)*(JH->nx-1));
                            jhy = (int)((((data_type)(trg_val-mmin))/del)*(JH->ny-1));

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
        gradient->set(i,0,mic,'g');



        //-------------- Y direction gradient component ---------------//
        if (Y_opt==1){
            //printf("......  y\n");

            //printf("ydirection gradient: \n");
            //reset JH_temp to equal JH_inst
            for(iii=0;iii<JH->nx;iii++){
                for(jjj=0;jjj<JH->ny;jjj++){
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

                        src_val = source->get_intrpl(  ((gradient->get(i,0,'c'))*RBF->get(bx,by,bz) + def_field->get(ii,jj,kk,'x')) , ((gradient->get(i,1,'c')+STEP)*RBF->get(bx,by,bz) + def_field->get(ii,jj,kk,'y')) , ((gradient->get(i,2,'c'))*RBF->get(bx,by,bz) + def_field->get(ii,jj,kk,'z')));
                        trg_val = (data_type)target->get(ii,jj,kk);

                        if ((src_val<mmin)||(trg_val<mmin)||(src_val>mmax)||(trg_val>mmax)) {
                            //do nothing
                        }else{
                            jhx = (int)((((data_type)(src_val-mmin))/del)*(JH->nx-1));
                            jhy = (int)((((data_type)(trg_val-mmin))/del)*(JH->ny-1));

                            JH_temp.increment(jhx,jhy);
                        }

                    }
                }
            }


            fph = JH_temp.mi();


            //repeat procedure for f(x-h)

            //reset JH_temp to equal JH_inst
            for(iii=0;iii<JH->nx;iii++){
                for(jjj=0;jjj<JH->ny;jjj++){
                    JH_temp.set(iii,jjj,(JH_i.get(iii,jjj)));
                }
            }


            for(ii=vx_start,bx=bx_start;ii<vx_end,bx<bx_end;ii++,bx++){
                for(jj=vy_start,by=by_start;jj<vy_end,by<by_end;jj++,by++){
                    for(kk=vz_start,bz=bz_start;kk<vz_end,bz<bz_end;kk++,bz++){


                        src_val = source->get_intrpl( ((gradient->get(i,0,'c'))*RBF->get(bx,by,bz) + def_field->get(ii,jj,kk,'x')) , ((gradient->get(i,1,'c')-STEP)*RBF->get(bx,by,bz) + def_field->get(ii,jj,kk,'y')) , ((gradient->get(i,2,'c'))*RBF->get(bx,by,bz) + def_field->get(ii,jj,kk,'z')));
                        trg_val = (data_type)target->get(ii,jj,kk);
                        if ((src_val<mmin)||(trg_val<mmin)||(src_val>mmax)||(trg_val>mmax)) {
                            //do nothing
                        }else{
                            jhx = (int)((((data_type)(src_val-mmin))/del)*(JH->nx-1));
                            jhy = (int)((((data_type)(trg_val-mmin))/del)*(JH->ny-1));

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
        gradient->set(i,1,mic,'g');



        //-------------- Z direction gradient component ---------------//
        if (Z_opt==1){
            //printf("......  z \n");

            //printf("zdirection gradient: \n");
            //reset JH_temp to equal JH_inst
            for(iii=0;iii<JH->nx;iii++){
                for(jjj=0;jjj<JH->ny;jjj++){
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

                        src_val = source->get_intrpl(  ((gradient->get(i,0,'c'))*RBF->get(bx,by,bz) + def_field->get(ii,jj,kk,'x')) , ((gradient->get(i,1,'c'))*RBF->get(bx,by,bz) + def_field->get(ii,jj,kk,'y')) , ((gradient->get(i,2,'c')+STEP)*RBF->get(bx,by,bz) + def_field->get(ii,jj,kk,'z')));
                        trg_val = (data_type)target->get(ii,jj,kk);

                        if ((src_val<mmin)||(trg_val<mmin)||(src_val>mmax)||(trg_val>mmax)) {
                            //do nothing
                        }else{
                            jhx = (int)((((data_type)(src_val-mmin))/del)*(JH->nx-1));
                            jhy = (int)((((data_type)(trg_val-mmin))/del)*(JH->ny-1));

                            JH_temp.increment(jhx,jhy);
                        }

                    }
                }
            }


            fph = JH_temp.mi();


            //repeat procedure for f(x-h)

            //reset JH_temp to equal JH_inst
            for(iii=0;iii<JH->nx;iii++){
                for(jjj=0;jjj<JH->ny;jjj++){
                    JH_temp.set(iii,jjj,(JH_i.get(iii,jjj)));
                }
            }


            for(ii=vx_start,bx=bx_start;ii<vx_end,bx<bx_end;ii++,bx++){
                for(jj=vy_start,by=by_start;jj<vy_end,by<by_end;jj++,by++){
                    for(kk=vz_start,bz=bz_start;kk<vz_end,bz<bz_end;kk++,bz++){


                        src_val = source->get_intrpl( ((gradient->get(i,0,'c'))*RBF->get(bx,by,bz) + def_field->get(ii,jj,kk,'x')) , ((gradient->get(i,1,'c'))*RBF->get(bx,by,bz) + def_field->get(ii,jj,kk,'y')) , ((gradient->get(i,2,'c')-STEP)*RBF->get(bx,by,bz) + def_field->get(ii,jj,kk,'z')));
                        trg_val = (data_type)target->get(ii,jj,kk);
                        if ((src_val<mmin)||(trg_val<mmin)||(src_val>mmax)||(trg_val>mmax)) {
                            //do nothing
                        }else{
                            jhx = (int)((((data_type)(src_val-mmin))/del)*(JH->nx-1));
                            jhy = (int)((((data_type)(trg_val-mmin))/del)*(JH->ny-1));

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
        gradient->set(i,2,mic,'g');

        //done with all three components of gradient for this control point
    }
    JH_temp.destroy();
    JH_i.destroy();
}
