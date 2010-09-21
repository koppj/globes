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

/***************************************************************************
 *            G L o B E S   -   M I N I - W O R K S H O P                  *
 *                                                                         *
 *                     May 13, 2010, Fermilab                              *
 *                                                                         *
 * Hands-on session on LBNE: Computing chi^2 for the wrong hierarchy       *
 * as a function of th13                                                   *
 ***************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <globes/globes.h>            /* GLoBES library */

/* Output file, empty string = stdout */
char MYFILE[]="th13chi2.dat";
FILE *outfile = NULL;

int main(int argc, char *argv[])
{
  /* Initialize libglobes */
  glbInit(argv[0]);

  /* New minimization algorithm -> significant speed-up */
  glbSelectMinimizer(GLB_MIN_POWELL);

  /* Load experiment definition */
  glbInitExperiment("NFstandard.glb",&glb_experiment_list[0],&glb_num_of_exps);

  /* Intitialize output */
  if (strlen(MYFILE) > 0)
  {
    outfile = fopen(MYFILE, "w");
    if (outfile == NULL)
    {
      printf("Error opening output file.\n");
      return -1;
    }
  }
  else
    outfile = stdout;  

  /* ``True'' oscillation parameters */
  double theta12 = asin(sqrt(0.3));
  double theta13 = asin(sqrt(0.001))/2.0;
  double theta23 = 45.0 * M_PI/180.0;
  double deltacp = M_PI/4.0;
  double sdm = 8.0e-5;
  double ldm = 2.5e-3;               /* Note: True hierarchy normal */
  
  /* Allocate oscillation parameter vectors */
  glb_params test_values = glbAllocParams();
  glb_params true_values = glbAllocParams();
  glb_params input_errors = glbAllocParams();
  
  /* Define oscillation parameter vector */
  glbDefineParams(true_values,theta12,theta13,theta23,deltacp,sdm,ldm);
  glbSetDensityParams(true_values,1.0,GLB_ALL);                  

  /* Define external prior vector */
  glbDefineParams(input_errors,theta12*0.1,0,0,0,sdm*0.1,0);/* 10% external error on solar params */
  glbSetDensityParams(input_errors,0.05,GLB_ALL); /* 5% matter density uncertainty */

  /* Tell GLoBES about oscillation parameters and priors */
  glbSetCentralValues(true_values);
  glbSetInputErrors(input_errors);
  glbSetOscillationParameters(true_values);

  /* Compute ``true'' event rates */
  glbSetRates();
 
  /* Scan parameter space as function of x = log sin^2 2 theta_13 */
  double x_lower = -4.0;
  double x_upper = -2.0;
  double delta_x =  0.05;
  double x;
  double chi2;
  for(x=x_lower; x <= x_upper+0.00001; x+=delta_x)
  {
    /* Define test values of oscillation parameters for IH fit */
    glbCopyParams(true_values, test_values);
    glbSetOscParams(test_values, -ldm, GLB_DM_31);/* Wrong hierarchy fit */
    glbSetOscParams(test_values, asin(sqrt(pow(10,x)))/2.0, GLB_THETA_13);

    chi2=glbChiTheta13(test_values,NULL,GLB_ALL);
    fprintf(outfile, " %15.10g %15.10g\n", pow(10.0,x), chi2);
  }
 
  fclose(outfile);
  glbFreeParams(test_values);                     /* Release allocated memory */
  glbFreeParams(true_values);
  glbFreeParams(input_errors);
  return 0;
}

