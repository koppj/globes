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

  /*setup default priors */  
  glbRegisterPriorFunction( glb_builtin_prior_prior,
			    glb_builtin_prior_starting_values,
			    glb_builtin_prior_input_errors);
  

}

