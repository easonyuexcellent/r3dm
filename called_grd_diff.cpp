#include "called_grd_diff.h"
#include "math.h"
called_grd_diff::called_grd_diff()
{
}


void called_grd_diff::init(gradient_class *_grd,hot_spots *_htspts,int _index,int _left,int _right){
    grd=_grd;
    htspts=_htspts;
    index=_index;
    left=_left;
    right=_right;
    flag=2;
}

void called_grd_diff::run(){

    int i;
    int spcx, spcy, spcz;
    int diffx, diffy, diffz;

    //spc = grd.supp;
    // spcz = grd.suppz;
    spcx = (int)ceil(grd->spcx);
    spcy = (int)ceil(grd->spcy);
    spcz = (int)ceil(grd->spcz);

    //printf("index in check_distance : %d \n", index);

    for (i=left;i<right;i++){

      diffx = (int)(grd->get(index,0,'p')-htspts->get(i,0));
      if(diffx<0) diffx=-diffx;

      diffy = (int)(grd->get(index,1,'p')-htspts->get(i,1));
      if(diffy<0) diffy=-diffy;

      diffz = (int)(grd->get(index,2,'p')-htspts->get(i,2));
      if(diffz<0) diffz=-diffz;

      // printf(" ---------------------------> check_distance: %f , %f , %f \n", grd.get(index,0,'p'), grd.get(index,1,'p'), grd.get(index,2,'p'));
      // printf(" --------------------------->               : %d , %d , %d \n", get(i,0), get(i,1) , get(i,2));
      // printf("\n");

      if ( (diffx <= spcx) && (diffy <= spcy) && (diffz<=spcz) ){

        flag = 1;
        return;

      }

    }

    flag =0;
    return;
}
