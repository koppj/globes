/* GLoBES -- General LOng Baseline Experiment Simulator
 * (C) 2002 - 2004,  The GLoBES Team
 *
 * GLoBES is mainly intended for academic purposes. Proper
 * credit must be given if you use GLoBES or parts of it. Please
 * read the section 'Credit' in the README file.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */
 
/* This example illustrates the use of modules with globes 
 *
 * Compile with `make -f MakefileExample'
 */ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <globes/globes.h>   /* GLoBES library */



int main(int argc, char *argv[])
{ 
  double (*new_prior)(const glb_params)=NULL;
  /* char* MYFILE=""; */ 
  char* MYFILE=""; /* if empty, write to screen */
  FILE* stream;
 double chi2;

  /* Define my standard oscillation parameters */
  double theta12 = asin(sqrt(0.8))/2;
  double theta13 = asin(sqrt(0.001))/2;
  double theta23 = M_PI/4;
  double deltacp = M_PI/2;
  double sdm = 7e-5;
  double ldm = 2.0e-3;
  
  glb_params true_values;
  glb_params fit_values;
  glb_params starting_values;
  glb_params input_errors;
  glb_params minimum;

  if(strlen(MYFILE)>0) stream=fopen(MYFILE, "w");
  else stream = stdout;
  
  /* Initialize libglobes */
  glbInit(argv[0]); 

 
 
  /* Initialize one experiment NuFact.glb */
  glbInitExperiment("NuFact1.glb",&glb_experiment_list[0],&glb_num_of_exps); 
  glbInitExperiment("NuFact1.glb",&glb_experiment_list[0],&glb_num_of_exps); 
   
  /* Initialize a number of parameter vector(s) */
  true_values = glbAllocParams();
  fit_values = glbAllocParams();
  starting_values = glbAllocParams();
  input_errors = glbAllocParams();
  minimum = glbAllocParams();
    
  /* Assign standard oscillation parameters */
  glbDefineParams(true_values,theta12,theta13,theta23,deltacp,sdm,ldm);

  /* The simulated data are computed with "true_values" */
  glbSetOscillationParameters(true_values);
  glbSetRates();
     
  
  /* Set test/fit values slightly off true values at s22th13=0.0015 */
  glbCopyParams(true_values,fit_values);
  glbSetOscParams(fit_values,asin(sqrt(0.0015))/2,GLB_THETA_13); 
  
  //glbLoadPrior("prior-template");
 

  /* Prepare minimizors: Set errors for external parameters: */
  /* 10% for each of the solar parameters, 5% for the matter density */
  glbDefineParams(input_errors,theta12*0.1,0,0,0,sdm*0.1,0);
  glbSetDensityParams(input_errors,0.05,GLB_ALL);
  glbSetStartingValues(true_values);
  glbSetInputErrors(input_errors);
  
 

  
  
  
  /* Check what the difference is if we keep in addition deltacp fixed. */
  /* This corresponds to projection onto theta13-deltacp-plane */
  fprintf(stream,"Single function ...\n");
  chi2 = glbChiThetaDelta(fit_values,minimum,1);
  fprintf(stream,"chi2 with correlations other than with deltacp: %g \n\n",
	  chi2);
  glbPrintParams(stream,minimum);
  fprintf(stream,"ALL function ...\n");
  chi2 = glbChiThetaDelta(fit_values,minimum,GLB_ALL);
  fprintf(stream,"chi2 with correlations other than with deltacp: %g \n\n",
	  chi2);
  glbPrintParams(stream,minimum);
  
  if(strlen(MYFILE)>0) fclose(stream);

  exit(0);
}
