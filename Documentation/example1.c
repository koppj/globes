/* Example: Correlation between s22th13 and deltacp
 * Copyright 2004 GLoBES collaboration
 * Compile with ``make example1''
 */ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <globes/globes.h>   /* GLoBES library */
#include "myio.h"             /* my input-output routines */

/* If filename given, write to file; for empty filename write to screen */
char MYFILE[]="test1.dat";

int main(int argc, char *argv[])
{ 
  /* Initialize libglobes */
  GLBInit(argv[0]); 

  /* Initialize experiment NuFact.glb */
  GLBInitExperiment("NuFact.glb",&ExpList[0],&numofexps); 
 
  /* Intitialize output */
  InitOutput(MYFILE,"Format: Log(10,s22th13)   deltacp   chi^2 \n"); 
  
  /* Define standard oscillation parameters */
  double theta12 = asin(sqrt(0.8))/2;
  double theta13 = asin(sqrt(0.001))/2;
  double theta23 = M_PI/4;
  double deltacp = M_PI/2;
  double sdm = 7e-5;
  double ldm = 2e-3;
  
	/* Initialize parameter vector(s) */
  glb_params true_values = alloc_params();
  glb_params test_values = alloc_params();

  true_values = 
    GLBDefineParams(true_values,theta12,theta13,theta23,deltacp,sdm,ldm);
  test_values =  
    GLBDefineParams(test_values,theta12,theta13,theta23,deltacp,sdm,ldm);  

  /* The simulated data are computed */
  GLBSetOscillationParameters(true_values);
  MSetRates();

  /* Iteration over all values to be computed */
  double thetheta13,x,y,res;    
    
  for(x=-4.0;x<-2.0+0.01;x=x+2.0/50)
  for(y=0.0;y<200.0+0.01;y=y+200.0/50)
  {
      /* Set vector of test=fit values */
      thetheta13=asin(sqrt(pow(10,x)))/2;
      test_values=set_osc_params(test_values,thetheta13,1);
      test_values=set_osc_params(test_values,y*M_PI/180.0,3);
  
      /* Compute Chi^2 for all loaded experiments and all rules */
      res=GLBChiSys(test_values,GLB_ALL,GLB_ALL);
      AddToOutput(x,y,res);
  }
   
  /* Destroy parameter vector(s) */
  free_params(true_values);
  free_params(test_values); 
  
  exit(0);
}
