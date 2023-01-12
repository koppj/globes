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
 * Example: using the SMEFT plugin
 * Compile with ``make example-smeft''
 */ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <globes/globes.h>   /* GLoBES library */
#include "myio.h"             /* my input-output routines */

/* If filename given, write to file; for empty filename write to screen */
char MYFILE[]="test-smeft.dat";

int main(int argc, char *argv[])
{ 
  /* Initialize libglobes */
  glbInit(argv[0]);

  /* Initialize experiment NFstandard.glb */
  glbInitExperiment("NF-smeft.glb",&glb_experiment_list[0],&glb_num_of_exps); 
 
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
  glb_params true_values = glbAllocParams();
  glb_params test_values = glbAllocParams();

  glbDefineParams(true_values,theta12,theta13,theta23,deltacp,sdm,ldm);
  glbSetDensityParams(true_values,1.0,GLB_ALL);
  glbDefineParams(test_values,theta12,theta13,theta23,deltacp,sdm,ldm);  
  glbSetDensityParams(test_values,1.0,GLB_ALL);

  /* The simulated data are computed */
  glbSetOscillationParameters(true_values);
  glbSetRates();

  /* verify that production/detection coefficients have been
   * loaded correctly */
  double E = 1.5;
  printf("E=%g, prod. coeff=(%g %g %g %g %g)\n",
         E,
         glbEFTFluxCoeff(0, 0, GLB_EFT_L, GLB_EFT_L, 0, E),
         glbEFTFluxCoeff(0, 0, GLB_EFT_L, GLB_EFT_R, 0, E),
         glbEFTFluxCoeff(0, 0, GLB_EFT_S, GLB_EFT_P, 0, E),
         glbEFTFluxCoeff(0, 0, GLB_EFT_P, GLB_EFT_T, 0, E),
         glbEFTFluxCoeff(0, 0, GLB_EFT_T, GLB_EFT_T, 2, E));
  printf("E=%g, det. coeff=(%g %g %g %g %g)\n",
         E,
         glbEFTXSecCoeff(0, 0, GLB_EFT_L, GLB_EFT_L, 0, E),
         glbEFTXSecCoeff(0, 0, GLB_EFT_L, GLB_EFT_R, 0, E),
         glbEFTXSecCoeff(0, 0, GLB_EFT_S, GLB_EFT_T, 0, E),
         glbEFTXSecCoeff(0, 0, GLB_EFT_P, GLB_EFT_T, 0, E),
         glbEFTXSecCoeff(0, 0, GLB_EFT_T, GLB_EFT_T, 2, E));
  int *q = glbEFTFluxQuarkFlavors(0, 0);
  printf("quark flavors flux #0: %d %d\n", q[0], q[1]);
  q = glbEFTFluxQuarkFlavors(0, 1);
  printf("quark flavors flux #1: %d %d\n", q[0], q[1]);
  q = glbEFTXSecQuarkFlavors(0, 0);
  printf("quark flavors xsec #0: %d %d\n", q[0], q[1]);
  q = glbEFTXSecQuarkFlavors(0, 1);
  printf("quark flavors xsec #1: %d %d\n", q[0], q[1]);
  getchar();


  /* Iteration over all values to be computed */
  double thetheta13,x,y,res;    
    
  for(x=-4.0;x<-2.0+0.01;x=x+2.0/50)
  for(y=0.0;y<200.0+0.01;y=y+200.0/50)
  {
      /* Set vector of test values */
      thetheta13=asin(sqrt(pow(10,x)))/2;
      glbSetOscParams(test_values,thetheta13,GLB_THETA_13);
      glbSetOscParams(test_values,y*M_PI/180.0,GLB_DELTA_CP);
  
      /* Compute Chi^2 for all loaded experiments and all rules */
      res=glbChiSys(test_values,GLB_ALL,GLB_ALL);
      AddToOutput(x,y,res);

      int i, j;
      double *p;
//      for (j=0; j < 4; j++)
//      {
//        p = glbGetSignalFitRatePtr(0,j);
//        printf("Sig %1d:  ", j);
//        for (i=0; i < glbGetNumberOfBins(0); i++)
//          printf("%7.4g ", p[i]);
//        printf("\n");
//
//        p = glbGetBGFitRatePtr(0,j);
//        printf("BG %1d:   ", j);
//        for (i=0; i < glbGetNumberOfBins(0); i++)
//          printf("%7.4g ", p[i]);
//        printf("\n");
//      }
//      for (j=0; j < 6; j++)
//      {
//        p = glbGetChannelFitRatePtr(0,j,GLB_POST);
//        printf("Ch %1d:  ", j);
//        for (i=0; i < glbGetNumberOfBins(0); i++)
//          printf("%7.4g ", p[i]);
//        printf("\n");
//      }
//      getchar();
  }
   
  /* Destroy parameter vector(s) */
  glbFreeParams(true_values);
  glbFreeParams(test_values); 
  
  exit(0);
}
