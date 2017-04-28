/*
The Adaptive Grid Registration Algorithm: A New Spline Modeling Aproach for Nonrigid Image Registration.
Program to perform nonrigid data registration in 3d.
Adapted to suit medical images.
By Gustavo Kunde Rohde, fall&winter 2000,2001.

This is the "main" function.
Its purpose is to initialize the optimization class given the input file.
Then it performs registration by calling the "run" method in the optimization class.
Memory cleanup is performed afterwards


*/

//#define UNIX
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

int main(int argc, char *argv[]){

  fprintf(stderr," Entering program.\n");

  config_class cfg_obj;

  if(argc!=2){
    printf("Please enter a valid configuration file.\n");
    exit(1);
  }

  cfg_obj.init(argv[1]);


  OPTIMIZATION program;

  program.init(cfg_obj);
  fprintf(stderr," Calling run procedure.\n");

  if (cfg_obj.OP_MODE==0){
    program.run();
  }else{
    program.run_gus();
  }

  program.output(cfg_obj.out);
  program.destroy();
  cfg_obj.destroy();

  return 1;
}
