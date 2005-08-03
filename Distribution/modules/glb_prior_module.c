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
 *
 *
 *
 * This is a template for prior modules 
 *
 * In order to avoid conflicting symbol names, replace `prior_template'
 * in the macro definition below by the module name which should match
 * the name of the executable (without postfix). Use the GLB_MOD_NSP 
 * macro to define *any* exported symbols like this
 *
 *  #define my_exported_symbol GLB_MOD_NSP(my_exported_symbol)
 *
 * Using this macro ensures that the namespace resolution works
 * independently of the use of libtool.
 *
 * Modules for GLoBES are most easily compiled and linked by using
 * the `globes-module' script (which is in the the bin directory of
 * your GLoBES installation).
 *
 *  globes-module my-module.c
 *
 * In order for the module loader of GLoBES to find my-module, you
 * have to add the directory where my-module can be found to the
 * GLB_MODULES environment variable.
 *
 * Prior modules are loaded by invoking 
 *
 *  glbLoadPrior("my-module")
 *
 * which sets the prior function (ant its relatives) to the functions
 * provided by my-module. glbLoadPrior returns the module handle for
 * further access to the module. Note the module name given to glbLoadPrior
 * is without extension, since this function tries all possible extensions
 * by itself.
 */

#define GLB_MODULE_NAME glb_prior_module

/* This tells GLoBES that this module provides a prior function */
#define GLB_PRIOR_MODULE

/* Headers needed by this module */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <globes/globes.h>


#ifndef GLB_WO_MODULES

#endif /* !GLB_WO_MODULES */

/* The module header has to be inlcuded after GLB_MODULE_NAME and
 * GLB_XXX_MODULE have been defined.
 */
#include <globes/glb-modules.h> /* Module header */ 



int glb_module_id = GLB_PRIOR_MODULE_ID;



/* File-scope static variables to store starting values and input errors */
static glb_params sv;
static glb_params er;

/* glb_module_init is mandatory. Any module must have glb_module_init, 
 * which returns 0 on succes or -1 on failure. This function is called 
 * when the module is loaded.
 *
 * In our case it allocates memory for sv and er and sets sv/er to any 
 * input errors starting values which have already been defined 
 */   
int glb_module_init()
{
  sv=glbAllocParams();
  er=glbAllocParams(); 
  glbGetStartingValues(sv);
  glbGetInputErrors(er);


  return 0;
}

/* glb_module_clean is optional. It should free all memory reserved by
 * glb_module_init and also reverse all side effects of glb_module_init.
 *
 * Here we just free the memory for sv and er.
 */
void glb_module_clean()
{
  glbFreeParams(er);
  glbFreeParams(sv);
}

/* This is a local function which implements the prior calculation
 * and ensure that the builtin prior behaves as described in the 
 * documentation, i.e. it returns 0 when the error is below 1E-12
 */
static double sprior(double x, double center, double sigma)
{
  if(fabs(sigma-0)<1E-12) return 0;
  return (x-center)*(x-center)/sigma/sigma;
}  

/* glb_module_prior is mandatory (for prior modules).
 *
 * This is the default prior, i.e. priors are only added for free
 * parameters, and the prior accepts starting values as stored in sv
 * and input errors as stored in er. 
 */
double glb_module_prior(const glb_params in)
{
  int i;
  double res=0;
  glb_projection pro;
  pro=glbAllocProjection();
  glbGetProjection(pro);
  for(i=0;i<GLB_OSCP;i++)
    {   
    if(glbGetProjectionFlag(pro,i)==GLB_FREE) res += 
						sprior(glbGetOscParams(in,i),
						       glbGetOscParams(sv,i),
						       glbGetOscParams(er,i)
						       );
    }
  
  for(i=0;i<glb_num_of_exps;i++)
    {
      
      if(glbGetDensityProjectionFlag(pro,i)==GLB_FREE) 
	res += 
	  sprior(glbGetDensityParams(in,i),
		 glbGetDensityParams(sv,i),
		 glbGetDensityParams(er,i)
		 );
    }

  glbFreeProjection(pro);
 
  return res;
}

/* glb_module_starting_values is optional. It is used to synchronize
 * the values set by glbSetStartingValues with the values used in here.
 * This function (if it exists) is called by glbSetStartingValues
 * with the same argument.
 *
 * Here we just copy the argument to our local buffer.
 */ 
int glb_module_starting_values(const glb_params in)
{
  if(glbCopyParams(in,sv)!=NULL) return 0;
  return -1;
}

/* glb_module_input_errors is optional. It is used to synchronize
 * the values set by glbSetInputErrors with the values used in here.
 * This function (if it exists) is called by glbSetInputErrors
 * with the same argument.
 *
 * Here we just copy the argument to our local buffer.
 */ 
int glb_module_input_errors(const glb_params in)
{ 
  if(glbCopyParams(in,er)!=NULL) return 0;
  return -1;
}
