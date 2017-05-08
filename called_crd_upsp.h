#ifndef CALLED_CRD_UPSP_H
#define CALLED_CRD_UPSP_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "register3d.h"
#include "splcoord.h"
#include <QThread>

class called_crd_upsp : public QThread
{
public:
    called_crd_upsp();
    void run();
    void set(sampling_coord *_p,int _i);
private:
    int i;
    sampling_coord *p;
};

#endif // CALLED_CRD_UPSP_H
