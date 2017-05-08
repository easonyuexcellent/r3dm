#include "called_crd_upsp.h"

called_crd_upsp::called_crd_upsp()
{

}

void called_crd_upsp::run()
{
    int j,k;
    data_type coordx,coordy,coordz;

    for(j=0;j<2*p->ny;j++){
        for(k=0;k<2*p->nz;k++){

            coordx = ((data_type)i)/2.0;
            coordy = ((data_type)j)/2.0;
            coordz = ((data_type)k)/2.0;

            p->set(i,j,k,'x',2.0*p->get_interp_o(coordx,coordy,coordz,'x'));
            p->set(i,j,k,'y',2.0*p->get_interp_o(coordx,coordy,coordz,'y'));
            p->set(i,j,k,'z',2.0*p->get_interp_o(coordx,coordy,coordz,'z'));

        }
    }
}


void called_crd_upsp::set(sampling_coord *_p, int _i)
{
    p=_p;
    i=_i;
}
