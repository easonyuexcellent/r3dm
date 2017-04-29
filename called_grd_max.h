#ifndef GRD_DIST_H
#define GRD_DIST_H
#include <QThread>
#include "gradient_class.h"
class called_grd_max : public QThread
{
public:
    called_grd_max();
    void run();
    void get(int *,data_type *);
    data_type max_val;
    int index;
    void init(gradient_class *_grd,int _lflag,int _rflag);
private:
    int lflag;
    int rflag;
    gradient_class *grd;
};

#endif // GRD_DIST_H
