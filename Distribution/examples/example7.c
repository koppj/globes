/* GLoBES -- General LOng Baseline Experiment Simulator
 * (C) 2002 - 2007,  The GLoBES Team
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

/* 
 * Example: Atmospheric neutrinos
 * Compile with ``make example7''
 */ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <globes/globes.h>   /* GLoBES library */
#include "myio.h"             /* my input-output routines */

/* If filename given, write to file; for empty filename write to screen */
char MYFILE[]="";

int main(int argc, char *argv[])
{ 
  /* Initialize libglobes */
  glbInit(argv[0]); 

  /* Initialize experiment NFstandard.glb */
  glbInitExperiment("pingu-mockup.glb",&glb_experiment_list[0],&glb_num_of_exps); 

  /* Intitialize output */
  InitOutput(MYFILE,"Format: Log(10,s22th13)   chi^2 one param   chi^2 all params \n"); 

  /* Define standard oscillation parameters */
  double theta12 = asin(sqrt(0.8))/2;
  double theta13 = asin(sqrt(0.1))/2;
  double theta23 = M_PI/4;
  double deltacp = M_PI/2;
  double sdm = 7e-5;
  double ldm = 2e-3;
  
  /* Initialize parameter and projection vector(s) */
  glb_params true_values = glbAllocParams();
  glb_params test_values = glbAllocParams();
  glb_params input_errors = glbAllocParams();
  glb_projection th23dm31_projection = glbAllocProjection();

  glbDefineParams(true_values,theta12,theta13,theta23,deltacp,sdm,ldm);
  glbSetDensityParams(true_values,1.0,GLB_ALL);
  glbDefineParams(test_values,theta12,theta13,theta23,deltacp,sdm,ldm);  
  glbSetDensityParams(test_values,1.0,GLB_ALL);

  /* The simulated data are computed */
  glbSetOscillationParameters(true_values);
  glbSetRates();

  /* Set starting values and input errors for all projections */  
  glbDefineParams(input_errors,theta12*0.1,0,0,0,sdm*0.1,0);  
  glbSetDensityParams(input_errors,0.05,GLB_ALL);
  glbSetCentralValues(true_values);
  glbSetInputErrors(input_errors);

  /* Set two-parameter projection onto th23-dm31 plane: leave th13 and deltacp free. */
//  glbDefineProjection(th23dm31_projection,GLB_FIXED,GLB_FREE,GLB_FIXED,
//    GLB_FREE,GLB_FIXED,GLB_FIXED);
  glbDefineProjection(th23dm31_projection,GLB_FIXED,GLB_FIXED,GLB_FIXED,
    GLB_FIXED,GLB_FIXED,GLB_FIXED);//FIXME
  glbSetDensityProjectionFlag(th23dm31_projection, GLB_FIXED, GLB_ALL);
  glbSetProjection(th23dm31_projection); 

  /* Test interpolation of atmospheric fluxes */
//  double cza, logE;
//  for (cza=-1.0; cza <= 1.0+0.00001; cza += 0.1)
//    for (logE=-1; logE <= 4.+0.00001; logE += 0.05)
//    {
//      double E = exp(logE*log(10.));
//      printf("%10.5g %10.5g   %10.5g %10.5g %10.5g   %10.5g %10.5g %10.5g\n", cza, E,
//          glbFlux(0, 0, E, glbCosThetaToL(cza), GLB_NU_E,   +1),
//          glbFlux(0, 0, E, glbCosThetaToL(cza), GLB_NU_MU,  +1),
//          glbFlux(0, 0, E, glbCosThetaToL(cza), GLB_NU_TAU, +1),
//          glbFlux(0, 0, E, glbCosThetaToL(cza), GLB_NU_E,   -1),
//          glbFlux(0, 0, E, glbCosThetaToL(cza), GLB_NU_MU,  -1),
//          glbFlux(0, 0, E, glbCosThetaToL(cza), GLB_NU_TAU, -1));
//    }
//  getchar();

  /* Iteration over all values to be computed */
  double s2th23, dm31;
  for (s2th23=0.0; s2th23 < 1.0+0.0001; s2th23 += 0.02)
  {
    for (dm31=2e-3; dm31 < 3e-3+0.0001; dm31 += 0.05e-3)
    {
      glbSetOscParams(test_values, asin(fmin(1.0, sqrt(s2th23))), GLB_THETA_23);
      glbSetOscParams(test_values, dm31, GLB_DM_31);

      /* Compute Chi^2 for two-parameter correlation: minimize over deltacp only */
      AddToOutput(s2th23, dm31, glbChiNP(test_values,NULL,GLB_ALL));
    }
  }
  
  /* Destroy parameter and projection vector(s) */
  glbFreeParams(true_values);
  glbFreeParams(test_values); 
  glbFreeParams(input_errors); 

  exit(0);
}
