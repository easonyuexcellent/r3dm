#include "multithread_hot_spots.h"

void multithread_hot_spots::init(hot_spots * _htspts, matrix_2d * _JH,volume * _source, volume * _target,sampling_coord * _def_field,basis_function * _RBF,input_vol_type _mmax, input_vol_type _mmin,int _nbins_src,int _nbins_trg)
{
    htspts=_htspts;
    JH=_JH;
    source=_source;
    target=_target;
    def_field=_def_field;
    RBF=_RBF;
    mmax=_mmax;
    mmin=_mmin;
    nbins_src=_nbins_src;
    nbins_trg=_nbins_trg;
}

void multithread_hot_spots::run(){
    runflag=1;

    data_type2 last_min,curr_min,diff;
    int STOP=0;
    while(STOP==0){

        grd();
        //sprintf(std_msg,"          Inside line minimize ... ");
        curr_min = line_minimize();


        diff = last_min-curr_min;
        if(diff<0) diff=-diff;

        if ( (diff<0.001) || (curr_min>last_min)) {

            //sprintf(std_msg,"should have stoped");
            STOP=2;
            runflag=0;
            endflag=1;
            return;

        } else {
            last_min = curr_min;
        }
    }
}

data_type2 multithread_hot_spots::line_minimize(){

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

                //sprintf(std_msg," -------> %d, %d, %d ", x+minx,y+miny,z+minz);
                def_field_temp.set(x,y,z,'x',def_field->get(x+minx,y+miny,z+minz,'x'));
                def_field_temp.set(x,y,z,'y',def_field->get(x+minx,y+miny,z+minz,'y'));
                def_field_temp.set(x,y,z,'z',def_field->get(x+minx,y+miny,z+minz,'z'));

            }
        }
    }

    //sprintf(std_msg,"--------------------> Inside line minimize, check point 1");

    //def_field_temp = def_field_temp + coeff*RBF
    for(i=0;i<gradient.num_points;i++){


        //find appropriate indexes for acessing
        vx_start = (int)mymax((gradient.get(i,0,'p')-gradient.supp),0);
        vx_end = (int)mymin((gradient.get(i,0,'p')+gradient.supp+1),(source->nx));
        vy_start = (int)mymax((gradient.get(i,1,'p')-gradient.supp),0);
        vy_end = (int)mymin((gradient.get(i,1,'p')+gradient.supp+1),(source->ny));
        vz_start = (int)mymax((gradient.get(i,2,'p')-gradient.suppz),0);
        vz_end = (int)mymin((gradient.get(i,2,'p')+gradient.suppz+1),(source->nz));

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

                    new_valx = def_field_temp.get(ii-minx,jj-miny,kk-minz,'x') + gradient.get(i,0,'c')*RBF->get(bx,by,bz);
                    new_valy = def_field_temp.get(ii-minx,jj-miny,kk-minz,'y') + gradient.get(i,1,'c')*RBF->get(bx,by,bz);
                    new_valz = def_field_temp.get(ii-minx,jj-miny,kk-minz,'z') + gradient.get(i,2,'c')*RBF->get(bx,by,bz);

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

                src_val = source->get_intrpl(def_field_temp.get(x,y,z,'x'),def_field_temp.get(x,y,z,'y'),def_field_temp.get(x,y,z,'z'));
                trg_val = (data_type)target->get(x+minx,y+miny,z+minz);

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


    //sprintf(std_msg,"--------------------> Inside line minimize, check point 2");
    //now JH_temp = JH - JH_temp
    for(ii=0;ii<JH->nx;ii++){
        for(jj=0;jj<JH->ny;jj++){

            JH_temp.set(ii,jj,(JH->get(ii,jj)-JH_temp.get(ii,jj)));
            //if(JH_temp.get(ii,jj)>0)
            // 	sprintf(std_msg,"jh %d",JH_temp.get(ii,jj));

        }
    }

    //sprintf(std_msg,"In line minimize, the JH of images - current region is %f ", JH_temp.mi());

    //now do def_field_temp = def_field (only for the current region) again
    for(x=0;x<def_field_temp.nx;x++){
        for(y=0;y<def_field_temp.ny;y++){
            for(z=0;z<def_field_temp.nz;z++){


                def_field_temp.set(x,y,z,'x',def_field->get(x+minx,y+miny,z+minz,'x'));
                def_field_temp.set(x,y,z,'y',def_field->get(x+minx,y+miny,z+minz,'y'));
                def_field_temp.set(x,y,z,'z',def_field->get(x+minx,y+miny,z+minz,'z'));

            }
        }
    }

    //sprintf(std_msg,"--------------------> Inside line minimize, check point 3");


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

    //sprintf(std_msg,"--------------------> Inside line minimize, check point 4");
    mi_last = mi_ini;
    mi_min = mi_ini;

    //for(i=0;i<gradient.num_points;i++){
    //	sprintf(std_msg,"    %f,%f,%f,",gradient.get(i,0,'g'),gradient.get(i,1,'g'),gradient.get(i,2,'g'));
    //}

    while(STOP==0){

        konst_curr = count*count*br_step;
        mi_curr = f_eval(konst_curr,JH_temp,JH_temp_copy,def_field_temp,minx,maxx,miny,maxy,minz,maxz);

        //sprintf(std_msg,"                 Bracketing print outs : : %f,%f",konst_curr,mi_curr);
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
            //sprintf(std_msg,"     JACOBIAN THRESHOLD VIOLATION STOPPED");
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
    //sprintf(std_msg,"");
    //sprintf(std_msg,"                    LINE MINIMIZATION        ");
    //sprintf(std_msg,"                 Bracketing print outs : : %f,%f",konst_curr,mi_curr);
    //sprintf(std_msg,"                       %f,%f  |  %f,%f  |  %f,%f",konst_ini,mi_ini,konst_min,mi_min,konst_last,mi_last);
    //sprintf(std_msg,"                       %f,%f",konst_min,mi_min);
    //sprintf(std_msg,"                       %f,%f",konst_last,mi_last);
    //sprintf(std_msg,"");

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
                //sprintf(std_msg,"                       %f,%f  |  %f,%f  |  %f,%f",konst_ini,mi_ini,konst_min,mi_min,konst_last,mi_last);

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
    //sprintf(std_msg,"updating the JH ... ");
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

void multithread_hot_spots::JH_inst_up(matrix_2d &JH_temp,
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

    def_field_temp.set(x,y,z,'x',def_field->get(x+minx,y+miny,z+minz,'x'));
    def_field_temp.set(x,y,z,'y',def_field->get(x+minx,y+miny,z+minz,'y'));
    def_field_temp.set(x,y,z,'z',def_field->get(x+minx,y+miny,z+minz,'z'));

      }
    }
  }


  //redo def_field_temp
  for(i=0;i<gradient.num_points;i++){


    //find appropriate indexes for acessing
    //		vx_start = (int)mymax((gradient.get(i,0,'p')-gradient.supp),0);
    //vx_end = (int)mymin((gradient.get(i,0,'p')+gradient.supp+1),(def_field_temp.nx));
    //vy_start = (int)mymax((gradient.get(i,1,'p')-gradient.supp),0);
    //vy_end = (int)mymin((gradient.get(i,1,'p')+gradient.supp+1),(def_field_temp.ny));
    //vz_start = (int)mymax((gradient.get(i,2,'p')-gradient.suppz),0);
    //vz_end = (int)mymin((gradient.get(i,2,'p')+gradient.suppz+1),(def_field_temp.nz));
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

      new_valx = def_field_temp.get(ii,jj,kk,'x') + ( cx )*RBF->get(bx,by,bz);
      new_valy = def_field_temp.get(ii,jj,kk,'y') + ( cy )*RBF->get(bx,by,bz);
      new_valz = def_field_temp.get(ii,jj,kk,'z') + ( cz )*RBF->get(bx,by,bz);

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


  //sprintf(std_msg,"In JH_inst_up, the JH of images - current region is %f ", JH_temp.mi());

  //now calculate JH_temp
  for(x=0;x<def_field_temp.nx;x++){
    for(y=0;y<def_field_temp.ny;y++){
      for(z=0;z<def_field_temp.nz;z++){

    src_val = source->get_intrpl(def_field_temp.get(x,y,z,'x'),def_field_temp.get(x,y,z,'y'),def_field_temp.get(x,y,z,'z'));
    trg_val = (data_type)target->get(x+minx,y+miny,z+minz);

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


  //and finally JH = JH_temp
  for(ii=0;ii<JH_temp.nx;ii++){
    for(jj=0;jj<JH_temp.ny;jj++){

      JH->set(ii,jj,JH_temp.get(ii,jj));

    }
  }

  //sprintf(std_msg,"New mi is %f ", JH->mi());

}

data_type2 multithread_hot_spots::f_eval(data_type konst, matrix_2d &JH_temp, matrix_2d &JH_temp_copy,
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


  // sprintf(std_msg,"                           Inside f_eval %d , %d , %d , %d , %d , %d ", minx, maxx, miny, maxy, minz, maxz);


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

  //sprintf(std_msg,"In f_eval, the JH of images - current region is %f ", JH_temp.mi());
  //sprintf(std_msg,"In f_eval, the JH of images - current region is %f ", JH_temp_copy.mi());
  //now do def_field_temp = def_field (only for the current region) again
  for(x=0;x<def_field_temp.nx;x++){
    for(y=0;y<def_field_temp.ny;y++){
      for(z=0;z<def_field_temp.nz;z++){

    def_field_temp.set(x,y,z,'x',def_field->get(x+minx,y+miny,z+minz,'x'));
    def_field_temp.set(x,y,z,'y',def_field->get(x+minx,y+miny,z+minz,'y'));
    def_field_temp.set(x,y,z,'z',def_field->get(x+minx,y+miny,z+minz,'z'));

      }
    }
  }

  //sprintf(std_msg,"              Inside f_eval checkpoint 1. ");

  //redo def_field_temp
  for(i=0;i<gradient.num_points;i++){

    //sprintf(std_msg,"     point %d out of %d", i, gradient.num_points);

    //find appropriate indexes for acessing
    //vx_start = (int)mymax((gradient.get(i,0,'p')-gradient.supp),0);
    //vx_end = (int)mymin((gradient.get(i,0,'p')+gradient.supp+1),(def_field_temp.nx));
    //vy_start = (int)mymax((gradient.get(i,1,'p')-gradient.supp),0);
    //vy_end = (int)mymin((gradient.get(i,1,'p')+gradient.supp+1),(def_field_temp.ny));
    //vz_start = (int)mymax((gradient.get(i,2,'p')-gradient.suppz),0);
    //vz_end = (int)mymin((gradient.get(i,2,'p')+gradient.suppz+1),(def_field_temp.nz));
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

    //sprintf(std_msg,"");
    //sprintf(std_msg,"in f_eval() V, %d, %d, %d, %d, %d, %d",vx_start, vx_end, vy_start, vy_end, vz_start, vz_end);
    //sprintf(std_msg,"in f_eval() B, %d, %d, %d, %d, %d, %d",bx_start, bx_end, by_start, by_end, bz_start, bz_end);
    //sprintf(std_msg,"in f_eval() def_field_temp size is: %d, %d, %d", def_field_temp.nx, def_field_temp.ny, def_field_temp.nz);
    //sprintf(std_msg,"in f_eval() RBF size is: %d, %d, %d", RBF->nx, RBF->ny, RBF->nz);
    //sprintf(std_msg,"");

    for(ii=vx_start,bx=bx_start;ii<vx_end,bx<bx_end;ii++,bx++){
      for(jj=vy_start,by=by_start;jj<vy_end,by<by_end;jj++,by++){
    for(kk=vz_start,bz=bz_start;kk<vz_end,bz<bz_end;kk++,bz++){


      //  sprintf(std_msg,"------------------------>>>>>>>>>>>> %d,%d %d,%d %d,%d", ii,bx,jj,by,kk,bz);
      new_valx = def_field_temp.get(ii,jj,kk,'x') + ( cx )*RBF->get(bx,by,bz);
      new_valy = def_field_temp.get(ii,jj,kk,'y') + ( cy )*RBF->get(bx,by,bz);
      new_valz = def_field_temp.get(ii,jj,kk,'z') + ( cz )*RBF->get(bx,by,bz);

      //new_valx = def_field_temp.get(ii,jj,kk,'x') + ( -konst*gradient.get(i,0,'g') )*RBF->get(bx,by,bz);
      //new_valy = def_field_temp.get(ii,jj,kk,'y') + ( -konst*gradient.get(i,1,'g') )*RBF->get(bx,by,bz);
      //new_valz = def_field_temp.get(ii,jj,kk,'z') + ( -konst*gradient.get(i,2,'g') )*RBF->get(bx,by,bz);

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

  //sprintf(std_msg,"              Inside f_eval checkpoint 2. ");

  //JH_temp_copy.setzero();
  //now calculate JH_temp_copy
  for(x=0;x<def_field_temp.nx;x++){
    for(y=0;y<def_field_temp.ny;y++){
      for(z=0;z<def_field_temp.nz;z++){

    src_val = source->get_intrpl(def_field_temp.get(x,y,z,'x'),def_field_temp.get(x,y,z,'y'),def_field_temp.get(x,y,z,'z'));
    trg_val = (data_type)target->get(x+minx,y+miny,z+minz);

    if ((src_val<mmin)||(trg_val<mmin)||(src_val>mmax)||(trg_val>mmax)) {

      //do nothing

    }else{

      jhx = (int)((((data_type)(src_val-mmin))/del)*(JH->nx-1));
      jhy = (int)((((data_type)(trg_val-mmin))/del)*(JH->ny-1));
      JH_temp_copy.increment(jhx,jhy);
    }

      }
    }
  }

  //	sprintf(std_msg,"              Inside f_eval checkpoint 3. ");


  return JH_temp_copy.mi();
}

void multithread_hot_spots::grd(void){

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

  fpp = JH->mi();

  JH_temp.setzero();
  JH_i.setzero();


  //sprintf(std_msg,"Computing gradient ... ");


  //sprintf(std_msg, "number of points %d", gradient.num_points);
  //sprintf(std_msg,"");

  del = (data_type)(mmax-mmin);

  //loop through each point in instance gradient
  for (i=0;i<gradient.num_points;i++){

    //sprintf(std_msg,"Point %d out of %d .",i,gradient.num_points);
    //sprintf(std_msg,"                  In grd(), point %d out of %d.",i, gradient.num_points);

    JH_temp.setzero();
    JH_i.setzero();



    //subtract the current subportion of the volume from the current histogram

    //for that build one first

    //find appropriate indexes for acessing
    vx_start = (int)mymax((gradient.get(i,0,'p')-gradient.supp),0);
    vx_end = (int)mymin((gradient.get(i,0,'p')+gradient.supp),(source->nx));
    vy_start = (int)mymax((gradient.get(i,1,'p')-gradient.supp),0);
    vy_end = (int)mymin((gradient.get(i,1,'p')+gradient.supp),(source->ny));
    vz_start = (int)mymax((gradient.get(i,2,'p')-gradient.suppz),0);
    vz_end = (int)mymin((gradient.get(i,2,'p')+gradient.suppz),(source->nz));

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

      src_val = source->get_intrpl((def_field->get(ii,jj,kk,'x')+gradient.get(i,0,'c')*RBF->get(bx,by,bz)),(def_field->get(ii,jj,kk,'y')+gradient.get(i,1,'c')*RBF->get(bx,by,bz)),(def_field->get(ii,jj,kk,'z')+gradient.get(i,2,'c')*RBF->get(bx,by,bz)));
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
      //sprintf(std_msg,"......   x");

      //sprintf(std_msg,"xdirection gradient: ");
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

        src_val = source->get_intrpl(  ((gradient.get(i,0,'c')+STEP)*RBF->get(bx,by,bz) + def_field->get(ii,jj,kk,'x')) , ((gradient.get(i,1,'c'))*RBF->get(bx,by,bz) + def_field->get(ii,jj,kk,'y')) , ((gradient.get(i,2,'c'))*RBF->get(bx,by,bz) + def_field->get(ii,jj,kk,'z')));
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


        src_val = source->get_intrpl( ((gradient.get(i,0,'c')-STEP)*RBF->get(bx,by,bz) + def_field->get(ii,jj,kk,'x')) , ((gradient.get(i,1,'c'))*RBF->get(bx,by,bz) + def_field->get(ii,jj,kk,'y')) , ((gradient.get(i,2,'c'))*RBF->get(bx,by,bz) + def_field->get(ii,jj,kk,'z')));
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
    gradient.set(i,0,mic,'g');



    //-------------- Y direction gradient component ---------------//
    if (Y_opt==1){
      //sprintf(std_msg,"......  y");

       //sprintf(std_msg,"ydirection gradient: ");
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

        src_val = source->get_intrpl(  ((gradient.get(i,0,'c'))*RBF->get(bx,by,bz) + def_field->get(ii,jj,kk,'x')) , ((gradient.get(i,1,'c')+STEP)*RBF->get(bx,by,bz) + def_field->get(ii,jj,kk,'y')) , ((gradient.get(i,2,'c'))*RBF->get(bx,by,bz) + def_field->get(ii,jj,kk,'z')));
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


        src_val = source->get_intrpl( ((gradient.get(i,0,'c'))*RBF->get(bx,by,bz) + def_field->get(ii,jj,kk,'x')) , ((gradient.get(i,1,'c')-STEP)*RBF->get(bx,by,bz) + def_field->get(ii,jj,kk,'y')) , ((gradient.get(i,2,'c'))*RBF->get(bx,by,bz) + def_field->get(ii,jj,kk,'z')));
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
    gradient.set(i,1,mic,'g');



    //-------------- Z direction gradient component ---------------//
    if (Z_opt==1){
      //sprintf(std_msg,"......  z ");

      //sprintf(std_msg,"zdirection gradient: ");
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

        src_val = source->get_intrpl(  ((gradient.get(i,0,'c'))*RBF->get(bx,by,bz) + def_field->get(ii,jj,kk,'x')) , ((gradient.get(i,1,'c'))*RBF->get(bx,by,bz) + def_field->get(ii,jj,kk,'y')) , ((gradient.get(i,2,'c')+STEP)*RBF->get(bx,by,bz) + def_field->get(ii,jj,kk,'z')));
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


        src_val = source->get_intrpl( ((gradient.get(i,0,'c'))*RBF->get(bx,by,bz) + def_field->get(ii,jj,kk,'x')) , ((gradient.get(i,1,'c'))*RBF->get(bx,by,bz) + def_field->get(ii,jj,kk,'y')) , ((gradient.get(i,2,'c')-STEP)*RBF->get(bx,by,bz) + def_field->get(ii,jj,kk,'z')));
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
    gradient.set(i,2,mic,'g');

    //done with all three components of gradient for this control point
  }

  gradient.normalize();
  gradient.max_d_compute(s_trsh);
  JH_temp.destroy();
  JH_i.destroy();
}
