/* Example: Finding the sgn(dm31^2)-degeneracy
 * Copyright 2004 GLoBES collaboration
 * Compile with ``make example3''
 */ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <globes/globes.h>   /* GLoBES library */
#include "myio.h"             /* my input-output routines */

/* If filename given, write to file; for empty filename write to screen */
char MYFILE[]="test3.dat";

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
  glb_params starting_values = alloc_params();
  glb_params input_errors = alloc_params();
  glb_params deg_pos = alloc_params();

  /* The simulated data are computed */
  true_values = 
    GLBDefineParams(true_values,theta12,theta13,theta23,deltacp,sdm,ldm);
  GLBSetOscillationParameters(true_values);
  MSetRates();

  /* Find sgn-degeneracy */  
  starting_values=
    GLBDefineParams(starting_values,theta12,theta13,theta23,deltacp,sdm,-ldm);  
  input_errors =  
    GLBDefineParams(input_errors,theta12*0.1,10,10,10,sdm*0.1,ldm/3);  
  input_errors = set_density_params(input_errors,0.05,GLB_ALL);
  GLBSetStartingValues(starting_values);
  GLBSetInputErrors(input_errors);
  /* BUGFIX: later the following will be replaced by 
    double CL=GLBChiALL(starting_values,deg_pos,GLB_ALL); */ 
  int i,projection[GLB_OSCP+32]; 
  for(i=0;i<GLB_OSCP+numofexps;i++) projection[i]=GLB_FREE;
  SelectProjection(&projection[0]);        
  double CL=GLBChiNP(starting_values,deg_pos,GLB_ALL);
   
  printf("Position of degeneracy: s22th13=%g, deltacp=%g; Confidence level: %g \n",
    get_osc_params(deg_pos,1),get_osc_params(deg_pos,3),CL);
  
  /* If degeneracy at low enough confidence level: compute section */
  if(CL<9.0)
  {
    double thetheta13,x,y,res;    
    for(x=-4.0;x<-2.0+0.01;x=x+2.0/50)
    for(y=0.0;y<200.0+0.01;y=y+200.0/50)
    {
        /* Set vector of test=fit values */
        thetheta13=asin(sqrt(pow(10,x)))/2;
        deg_pos=set_osc_params(deg_pos,thetheta13,1);
        deg_pos=set_osc_params(deg_pos,y*M_PI/180.0,3);
    
        /* Compute Chi^2 for all loaded experiments and all rules */
        res=GLBChiSys(deg_pos,GLB_ALL,GLB_ALL);
        AddToOutput(x,y,res);
    }
  }
   
  /* Destroy parameter vector(s) */
  free_params(true_values);
  free_params(starting_values); 
  free_params(input_errors); 
  free_params(deg_pos); 
  
  exit(0);
}
