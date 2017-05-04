#ifndef NODE_H
#define NODE_H
#include "head.h"

class node{
public:
    int index;
    data_type val;
    node (int a,data_type b):index(a),val(b){;}
    bool operator < (const node &m) const{
        return val>m.val;
    }
};

#endif // NODE_H
