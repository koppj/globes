#include <stdio.h>
#include <globes/globes.h>


void glbSuperInit(char *name)
{
  
  fprintf(stderr,"\n************************************\n");
  fprintf(stderr,"* This is Super Init               *\n");
  fprintf(stderr,"************************************\n\n");
  glbInit(name);
  glbLoadPrior("glb_prior_module"); 
}
