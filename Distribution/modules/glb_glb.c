#include <stdio.h>
#include <globes/globes.h>
#include "source/glb_wrapper.h"

void glbInit(char *name)
{
  
  fprintf(stderr,"\n************************************\n");
  fprintf(stderr,"* This is Super Init               *\n");
  fprintf(stderr,"************************************\n\n");
  glb_init(name);
  glbLoadPrior("glb_prior_module"); 
}
