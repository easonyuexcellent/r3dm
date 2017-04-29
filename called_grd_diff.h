#ifndef CALLED_GRD_DIFF_H
#define CALLED_GRD_DIFF_H
#include <QThread>
#include "gradient_class.h"
#include "hot_spots.h"

class called_grd_diff : public QThread
{
public:
    called_grd_diff();
    void run();
    int flag;
    void init(gradient_class *_grd,hot_spots *_htspts,int _index,int _left, int _right);
 private:
    int index;
    int left;
    int right;
    hot_spots *htspts;
    gradient_class *grd;
};
#endif // CALLED_GRD_DIFF_H
