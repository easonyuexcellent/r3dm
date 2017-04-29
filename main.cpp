

#include "register3d.h"
int main(int argc, char *argv[])
{
    Register3d r3d(argc,argv);
    r3d.start();
    r3d.wait();
}
