/* Example: Creating bar charts
 * Copyright 2004 GLoBES collaboration
 * Compile with ``make example4''
 */ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <globes/globes.h>   /* GLoBES library */
#include "myio.h"             /* my input-output routines */

/* If filename given, write to file; for empty filename write to screen */
char MYFILE1[]="test4a.dat";
char MYFILE2[]="test4b.dat";
char MYFILE3[]="test4c.dat";
char MYFILE4[]="test4d.dat";

double theta12,theta13,theta23,deltacp,sdm,ldm;

double Sign(double n)
{
 if(n>=0) return +1.0;
 else return -1.0;   
}

void SetOscParams(double theldm)
{
 glb_params true_values = alloc_params();
 true_values = 
   GLBDefineParams(true_values,theta12,theta13,theta23,deltacp,sdm,theldm);
 GLBSetOscillationParameters(true_values);
 MSetRates();
 free_params(true_values);  
}

double CalcSystematics(double theldm,double thex)
{
  SetOscParams(theldm);
  double thetheta13=asin(sqrt(pow(10,thex)))/2;
  glb_params test_values = alloc_params();
  test_values=
    GLBDefineParams(test_values,theta12,thetheta13,theta23,0.0,sdm,theldm);
  double res=GLBChiSys(test_values,GLB_ALL,GLB_ALL);
  free_params(test_values);  
  return res;
}

double CalcNoSystematics(double theldm,double thex)
{
 int rules=GLBGetNumberOfRules(0);
 int j;
  /* Later with systematics off/on */
 for(j=0;j<rules;j++) MSetErrorDim(2,0,j);
 double res=CalcSystematics(theldm,thex);
 for(j=0;j<rules;j++) MSetErrorDim(0,0,j);
 return res; 
}

double CalcProjection(double theldm,double thex,glb_params start_vector)
{
  SetOscParams(theldm);

  double thetheta13=asin(sqrt(pow(10,thex)))/2;
  double thesign=Sign(get_osc_params(start_vector,5));

  glb_params input_errors = alloc_params();
  glb_params starting_values = alloc_params();
  glb_params minimum = alloc_params();
  input_errors =  
    GLBDefineParams(input_errors,theta12*0.1,10,10,10,sdm*0.1,theldm/3);  
  input_errors = set_density_params(input_errors,0.05,GLB_ALL);
  /* Set starting values to respective +-dm31^2 to avoid falling into unwanted solution;
     Note that the error on dm31^2 should not be too small in order to avoid large
     contributions from the priors */  
  starting_values =GLBDefineParams(starting_values,theta12,theta13,
                                   theta23,deltacp,sdm,theldm*thesign);
  GLBSetStartingValues(starting_values);
  GLBSetInputErrors(input_errors);

  start_vector=set_osc_params(start_vector,thetheta13,1);
  double res=GLBChiTheta(start_vector,minimum,GLB_ALL);
  
  /* Trick for next run: Find position of minimum better by using deltacp of prior 
     calculation; note that this modifies the org/deg-vectors in the main routine! */
  start_vector=set_osc_params(start_vector,get_osc_params(minimum,3),3);
 
  free_params(input_errors);
  free_params(starting_values);
  free_params(minimum);
  return res; 
}

double FindDeg(glb_params deg_pos,double theldm)
{
  SetOscParams(theldm);

  glb_params input_errors = alloc_params();
  glb_params starting_values = alloc_params();

  starting_values=
    GLBDefineParams(starting_values,theta12,theta13,theta23,deltacp,sdm,-theldm);  
  input_errors =  
    GLBDefineParams(input_errors,theta12*0.1,10,10,10,sdm*0.1,theldm/3);  
  input_errors = set_density_params(input_errors,0.05,GLB_ALL);
  GLBSetStartingValues(starting_values);
  GLBSetInputErrors(input_errors);
  /* BUGFIX: later the following will be replaced by 
    double CL=GLBChiALL(starting_values,deg_pos,GLB_ALL); */ 
  int i,projection[GLB_OSCP+32]; 
  for(i=0;i<GLB_OSCP+numofexps;i++) projection[i]=GLB_FREE;
  SelectProjection(&projection[0]);        
  double CL=GLBChiNP(starting_values,deg_pos,GLB_ALL);

  free_params(input_errors);
  free_params(starting_values);
  return CL;
}

int main(int argc, char *argv[])
{ 
  /* Initialize libglobes */
  GLBInit(argv[0]); 

  /* Initialize experiment NuFact.glb */
  GLBInitExperiment("NuFact.glb",&ExpList[0],&numofexps); 

  /* Set standard oscillation parameters */
  theta12 = asin(sqrt(0.8))/2;
  theta13 = 0;
  theta23 = M_PI/4;
  deltacp = 0;
  sdm = 7e-5;
  ldm = 2e-3;

  glb_params org1 = alloc_params();
  glb_params org2 = alloc_params();
  glb_params deg1 = alloc_params();
  glb_params deg2 = alloc_params();
   
  double thetheta13,x,res1,res2;
  
  /* Compute 1st edges of bars: systematics off, fit value of deltacp=0 */
  InitOutput(MYFILE1,"Format: Log(10,s22th13)  chi^2:ldm=2e-3  chi^2:ldm=3.0e-3 \n"); 
  for(x=-7.0;x<-2.0+0.01;x=x+5.0/50)
  {
     res1=CalcNoSystematics(2.0e-3,x);
     res2=CalcNoSystematics(3.0e-3,x);
     AddToOutput(x,res1,res2);
  }
  
  /* Compute 2nd edges of bars: systematics on, fit value of deltacp=0 */
  InitOutput(MYFILE2,"Format: Log(10,s22th13)  chi^2:ldm=2e-3  chi^2:ldm=3.0e-3 \n"); 
  for(x=-7.0;x<-2.0+0.01;x=x+5.0/50)
  {
     res1=CalcSystematics(2.0e-3,x);
     res2=CalcSystematics(3.0e-3,x);
     AddToOutput(x,res1,res2);
  }
 
  /* Compute 3rd edges of bars: systematics on, projection onto s22th13 axis */
  InitOutput(MYFILE3,"Format: Log(10,s22th13)  chi^2:ldm=2e-3  chi^2:ldm=3.0e-3 \n"); 
  org1 =GLBDefineParams(org1,theta12,theta13,
                                   theta23,deltacp,sdm,2.0e-3);
  org2 =GLBDefineParams(org2,theta12,theta13,
                                   theta23,deltacp,sdm,3.0e-3);
  for(x=-7.0;x<-2.0+0.01;x=x+5.0/50)
  {
     res1=CalcProjection(2.0e-3,x,org1);
     res2=CalcProjection(3.0e-3,x,org2);
     AddToOutput(x,res1,res2);
  }

  /* Find sgn-degeneracies */
  FindDeg(deg1,2.0e-3);
  FindDeg(deg2,3.0e-3);
 
  /* Compute 4th edges of bars: systematics on, degeneracy, projection onto s22th13 axis */
  InitOutput(MYFILE4,"Format: Log(10,s22th13)  chi^2:ldm=2e-3  chi^2:ldm=3.0e-3 \n"); 
  for(x=-7.0;x<-2.0+0.01;x=x+5.0/30)
  {
     res1=CalcProjection(2.0e-3,x,deg1);
     res2=CalcProjection(3.0e-3,x,deg2);
     AddToOutput(x,res1,res2);
  }

  free_params(org1);
  free_params(org2);
  free_params(deg1);
  free_params(deg2);
  exit(0);
}
