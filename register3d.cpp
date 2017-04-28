#define PC

#ifdef UNIX
#include <stream.h>
#else
#include <iostream>
#endif

#include <stdlib.h>
#include <stdio.h>
#include "register3d.h"
#include "volume.h"
#include "basisfctn.h"
#include "splcoord.h"
#include "OPTIMIZATION.h"
#include "config_class.h"
#include <QMetaType>
#include <QString>

Register3d::Register3d(int _a, char ** _b)
{
    argc=_a;
    argv=_b;
}

void Register3d::run()
{

    fprintf(stderr," Entering program.\n");

    config_class cfg_obj;

    if(argc!=2){
      printf("Please enter a valid configuration file.\n");
      exit(1);
    }

    cfg_obj.init(argv[1]);

    OPTIMIZATION program;
    program.init(cfg_obj);
    if (cfg_obj.OP_MODE==2){
        program.run_multi();
    }else if (cfg_obj.OP_MODE==1){
        program.run_gus();
    }else{
        program.run();
    }

    program.output(cfg_obj.out);
    program.destroy();
    cfg_obj.destroy();
    printf("Run procedure ended.\n");
}
