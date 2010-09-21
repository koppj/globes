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
 *                                                                         *
 * Solution to problem 6: Usage of degfinder                               *
 ***************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <globes/globes.h>            /* GLoBES library */
#include "degfinder.h"

/* Output file, empty string = stdout */
char MYFILE[]="th13chi2.dat";
FILE *outfile = NULL;

int main(int argc, char *argv[])
{
  int i;

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

  /* Data structures for degfinder -- prescan over delta */
  const int n_prescan_params = 1;                 /* # of params for prescan        */
  const int prescan_params[] = { GLB_DELTA_CP };  /* List of params to scan over    */
  const double prescan_min[] = { 0.0          };  /* Lower edge of scan range       */
  const double prescan_max[] = { 2*M_PI       };  /* Upper edge of scan range       */
  const int prescan_steps[]  = { 16           };  /* Number of steps                */
  const int max_deg = 32;                         /* Max # of degeneracies expected */
  int n_deg;                                      /* # of degeneracies found        */
  double deg_chi2[max_deg];                       /* chi^2 values at degeneracies   */
  glb_params deg_pos[max_deg];                    /* Locations of degeneracies      */
  for (i=0; i < max_deg; i++)
    deg_pos[i] = glbAllocParams();

  /* Define projection for prescan -- keep all parameters fixed */
  /* (only scan over delta_{CP} -> fastest                      */
  glb_projection prescan_proj = glbAllocProjection();
  glbDefineProjection(prescan_proj, GLB_FIXED, GLB_FIXED, GLB_FIXED,
                                    GLB_FIXED, GLB_FIXED, GLB_FIXED);
  glbSetDensityProjectionFlag(prescan_proj, GLB_FIXED, GLB_ALL);  
  
  /* Define projection for production run -- marginalize over all */
  /* oscillation parameters (except theta-130                    */
  glb_projection proj = glbAllocProjection();
  glbDefineProjection(proj, GLB_FREE, GLB_FIXED, GLB_FREE,
                            GLB_FREE, GLB_FREE, GLB_FREE);
  glbSetDensityProjectionFlag(proj, GLB_FREE, GLB_ALL);  

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

    /* Invoke degfinder. This will first perform ther prescan over delta_{CP}
     * (switching off systamtical errors and minimization over other oscillation
     * parameters) to get the rough locations of the degenerate solutions. These
     * are then refined using a full-fledged marginalization.
     * We use the flag DEG_NO_NH to tell degfinder that we are interested only
     * in inverted hierarchy solutions */
    fprintf(stderr, "Running degfinder at sin^2 2th13 = %g ... ", pow(10.0,x));
    degfinder(test_values, n_prescan_params, prescan_params, prescan_min, prescan_max,
              prescan_steps, prescan_proj, proj, &n_deg, deg_pos, deg_chi2,
              DEG_NO_NH);
    fprintf(stderr, "found %d local minima.\n", n_deg);

    /* Analyze degeneracies */
    double chi2 = 1.0e10;
    for (i=0; i < n_deg; i++)
    {
      fprintf(stderr, "  delta_{CP} = %10.7g degrees, chi^2 = %10.7g\n",
              glbGetOscParams(deg_pos[i], GLB_DELTA_CP) * 180.0/M_PI, deg_chi2[i]);
      if (glbGetOscParams(deg_pos[i], GLB_DM_31) < 0.0)
      {
        if (deg_chi2[i] < chi2)
          chi2 = deg_chi2[i];
      }
      else
        fprintf(stderr, "Error: Minimizer converged into NH minimum."
                        "This should not have happened.\n");
    }

    fprintf(outfile, " %15.10g %15.10g\n", pow(10.0,x), chi2);
  }
 
  fclose(outfile);
  for (i=0; i < max_deg; i++)
    if (deg_pos[i] != NULL)  { glbFreeParams(deg_pos[i]); deg_pos[i] = NULL; }
  glbFreeProjection(proj);
  glbFreeProjection(prescan_proj);

  glbFreeParams(test_values);                     /* Release allocated memory */
  glbFreeParams(true_values);
  glbFreeParams(input_errors);
  return 0;
}

