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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <globes/globes.h>
#include "source/glb_wrapper.h"
#include "source/glb_error.h"
#include "glb_prior.h"
#include "glb_probability.h"
#include "glb_rate_engine.h"
#include "glb_sys.h"


/* Moved from glb_wrapper.c, in order to have the module support stuff
 * only used in this file.
 */ 
static void final_clean() 
{
  glb_builtin_prior_clean();
  glb_clean_up(); 
  return; 
}

void glbInit(char *name)
{
  glb_dlhandle prior;

  atexit(final_clean);
  glb_init(name);

  glb_builtin_prior_init();

  /* Initialize some global variables */
  glb_sys_list = NULL;

  /* Setup default probability engine */
  glbRegisterProbabilityEngine(-1, NULL, NULL, NULL);
  glbRegisterProbabilityEngine(6, &glb_probability_matrix,
                                  &glb_set_oscillation_parameters,
                                  &glb_get_oscillation_parameters);

  /* Setup default priors */  
  glbRegisterPriorFunction( glb_builtin_prior_prior,
			    glb_builtin_prior_starting_values,
			    glb_builtin_prior_input_errors);
  
  /* Register built-in chi^2 functions */
  /* When making any changes here, don't forget to update the array
   * glb_2011_compatible_chi_functions in glb_sys.c and the functions
   * glbConvertErrorDim and glbGetErrorDim */
  glbDefineChiFunction(&glb_chi_sys_w_bg,       4, "sysSpectrumTilt");
  glbDefineChiFunction(&glb_chi_no_sys,         0, "sysNoSysSpectrum");
  glbDefineChiFunction(&glb_chi_sys_w_bgtot,    4, "sysTotalRatesTilt");
  glbDefineChiFunction(&glb_chi_spec,           1, "sysSpectrumOnly");
  glbDefineChiFunction(&glb_chi_no_sys_tot,     0, "sysNoSysTotalRates");
  glbDefineChiFunction(&glb_chi_sys_w_bg_calib, 4, "sysSpectrumCalib");
}



