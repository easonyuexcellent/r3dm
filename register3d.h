/*type definitions for register3d.c
Gustavo Rohde
*/



#ifndef __REGISTER_3D_H__
#define __REGISTER_3D_H__
#include "head.h"
#include <QThread>
#include <QMetaType>
#include <QString>

class Register3d : public QThread
{
public:
    Register3d(int, char **);

protected:
    void run();
    void printMessage();

private:
    int argc; //argument count
    char ** argv; //argument value
};

#endif
