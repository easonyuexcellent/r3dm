#include "called_grd_max.h"
#include "math.h"

called_grd_max::called_grd_max()
{

}

void called_grd_max::init(gradient_class *_grd, int _lflag, int _rflag){
    grd=_grd;
    lflag=_lflag;
    rflag=_rflag;
}

void called_grd_max::get(int *_index, data_type *_max_val){
    *_index=index;
    *_max_val=max_val;
}

void called_grd_max::run(){
    data_type curr_val;
    index = lflag;
    data_type cx,cy,cz;
    cx = grd->get(lflag,0,'c');
    cy = grd->get(lflag,1,'c');
    cz = grd->get(lflag,2,'c');
    max_val = (data_type)sqrt(cx*cx + cy*cy + cz*cz);
    //	printf("in find_max, max_val is %f \n", *max_val);
    //	printf("in find_max, grd.numl is %d \n", grd.num_points);

    for (int i=lflag+1;i<rflag;i++){

        cx = grd->get(i,0,'c');
        cy = grd->get(i,1,'c');
        cz = grd->get(i,2,'c');
        curr_val = (data_type)sqrt(cx*cx + cy*cy + cz*cz);
        //	printf("i : %d \n", i);
        //	printf("curr_val is %f \n", curr_val);

        if (curr_val>(max_val)){
            max_val = curr_val;
            index = i;
        }
    }
}
