/* Example: Projection of two- and n-dimensional manifold onto stheta-axis
 * Copyright 2004 GLoBES collaboration
 * Compile with ``make example2''
 */ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <globes/globes.h>   /* GLoBES library */
#include "myio.h"             /* my input-output routines */

/* If filename given, write to file; for empty filename write to screen */
char MYFILE[]="test2.dat";

int main(int argc, char *argv[])
{ 
  /* Initialize libglobes */
  GLBInit(argv[0]); 

  /* Initialize experiment NuFact.glb */
  GLBInitExperiment("NuFact.glb",&ExpList[0],&numofexps); 

  /* Intitialize output */
  InitOutput(MYFILE,"Format: Log(10,s22th13)   chi^2 one param   chi^2 all params \n"); 

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
  glb_params input_errors = alloc_params();

  true_values = 
    GLBDefineParams(true_values,theta12,theta13,theta23,deltacp,sdm,ldm);
  test_values =  
    GLBDefineParams(test_values,theta12,theta13,theta23,deltacp,sdm,ldm);  

  /* The simulated data are computed */
  GLBSetOscillationParameters(true_values);
  MSetRates();

  /* Set starting values and input errors for all projections */  
  input_errors =  
    GLBDefineParams(input_errors,theta12*0.1,10,10,10,sdm*0.1,ldm);  
  input_errors = set_density_params(input_errors,0.05,GLB_ALL);
  GLBSetStartingValues(true_values);
  GLBSetInputErrors(input_errors);

  /* Iteration over all values to be computed */
  double thetheta13,x,res1,res2;    
  int i,projection[GLB_OSCP+32]; 
  for(x=-4;x<-2.0+0.001;x=x+2.0/50)
  {
      /* Set vector of test=fit values */
      thetheta13=asin(sqrt(pow(10,x)))/2;
      test_values=set_osc_params(test_values,thetheta13,1);
     
      /* Guess fit value for deltacp in order to safely find minimum */
      test_values=set_osc_params(test_values,200.0/2*(x+4)*M_PI/180,3);
 
      /* Compute Chi^2 for two-parameter correlation: minimize over deltacp only */
      /* Note that this has to be done every time a different GLB.. is used! */
      for(i=0;i<GLB_OSCP+numofexps;i++) projection[i]=GLB_FIXED;
      projection[3]=GLB_FREE;
      SelectProjection(&projection[0]);        
      res1=GLBChiNP(test_values,NULL,GLB_ALL);
      
      /* Compute Chi^2 for full correlation: minimize over all but theta13 */
      /* Calls SelectProjection automatically! */
      res2=GLBChiTheta(test_values,NULL,GLB_ALL);
  
      AddToOutput(x,res1,res2);
  }
  
  /* Destroy parameter vector(s) */
  free_params(true_values);
  free_params(test_values); 
  free_params(input_errors); 
  
  exit(0);
}
