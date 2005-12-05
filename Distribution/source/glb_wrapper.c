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




#include <stdio.h>
#include <string.h>
#include <obstack.h>
#include <globes/globes.h>


#include "glb_probability.h"
#include "glb_fluxes.h"
#include "glb_rate_engine.h"
#include "glb_minimize.h"
#include "glb_types.h"
#include "glb_multiex.h"
#include "glb_version.h"
#include "glb_error.h"
#include "glb_smear.h"
#include "glb_lexer.h"
#include "glb_parser_addons.h"
#include "glb_path.h"


#include "glb_wrapper.h"

/* The global variables */
int glb_num_of_exps;
glb_exp glb_experiment_list[32];
int glb_rule_number;



#define obstack_chunk_alloc glb_malloc
#define obstack_chunk_free glb_free

/* This  is a bunch of function in order to deal with the
 * glb_params_type type
 */

glb_osc_type *glb_alloc_osc_type()
{
  glb_osc_type *temp;
  temp=(glb_osc_type *) glb_malloc(sizeof(glb_osc_type));
  temp->osc_params=(double *) glb_malloc(sizeof(double) *  GLB_OSCP); 
  temp->length = (size_t) GLB_OSCP;
  return temp;
}

void glb_free_osc_type(glb_osc_type *stale)
{
  if(stale!=NULL) glb_free(stale->osc_params);
  glb_free(stale);
}

glb_density_type *glb_alloc_density_type()
{
  int i;
  glb_density_type *temp;
  temp=(glb_density_type *) glb_malloc(sizeof(glb_density_type));
  temp->density_params=(double *) glb_malloc(sizeof(double)*glb_num_of_exps);
  temp->length=(size_t) glb_num_of_exps;
  for(i=0;i<glb_num_of_exps;i++) temp->density_params[i]=1.0;
  return temp;
}

void glb_free_density_type(glb_density_type *stale)
{
  if(stale!=NULL) glb_free(stale->density_params);
  glb_free(stale);
}

glb_params glbAllocParams()
{
  glb_params_type *temp;
  temp=(glb_params_type *) glb_malloc(sizeof(glb_params_type));
  temp->osc=glb_alloc_osc_type();
  temp->density=glb_alloc_density_type();
  temp->iterations=0;
  return (glb_params) temp;
}

void glbFreeParams(glb_params stale)
{
  glb_free_osc_type(((glb_params_type *) stale)->osc);
  glb_free_density_type(((glb_params_type *) stale)->density);
  glb_free(stale);
  stale = NULL;
}

glb_params glbSetIteration(glb_params in, int iter)
{
  glb_params_type *s;
  s=(glb_params_type *) in;
  s->iterations=iter;
  return (glb_params) s;
}

int glbGetIteration(const glb_params in)
{
  int iter;
  glb_params_type *s;
  s=(glb_params_type *) in;
  iter = s->iterations;
  return iter;
}

glb_params glbSetDensityParams(glb_params in, 
				 double dens, int which)
{
  int i;
  glb_params_type *s;
  s=(glb_params_type *) in;
  if(which==GLB_ALL)
    {
      for(i=0;i<(s->density)->length;i++) (s->density)->density_params[i]=dens;
    }
  else if(0 <= which&&which < (s->density)->length )
    {
      (s->density)->density_params[which]=dens;
    }
  else 
    {
      glb_error("Density list length mismatch");
      return NULL;
    }
  return (glb_params) s;
}

double glbGetDensityParams(const glb_params in, int which)
{
  double out;
  glb_params_type *s;
  s=(glb_params_type *) in;
  if(0 <= which&&which < (s->density)->length )
    {
      out=(s->density)->density_params[which];
    }
  else 
    {
      glb_error("Density list length mismatch");
      return 0;
    }
  return out;
}

glb_params glbSetOscParams(glb_params in,
				 double osc, int which)
{

  glb_params_type *s;
  s=(glb_params_type *) in;
   if(0 <= which&&which < GLB_OSCP )
    {
      (s->osc)->osc_params[which]=osc;
    }
  else 
    {
      glb_error("Oscillation list length mismatch");
      return NULL;
    }
  return (glb_params) s;
}

double glbGetOscParams(const glb_params in, int which)
{
  double out;
  glb_params_type *s;
  s=(glb_params_type *) in;
 if(0 <= which&&which < GLB_OSCP )
    {
      out=(s->osc)->osc_params[which];
    }
  else 
    {
      glb_error("Oscillation list length mismatch");
      return 0;
    }
  return out;
}

glb_params 
glbDefineParams(glb_params in,double theta12, double theta13, double theta23,
	      double delta, double dms, double dma)
{
  in=glbSetOscParams(in,theta12,0);
  in=glbSetOscParams(in,theta13,1);
  in=glbSetOscParams(in,theta23,2);
  in=glbSetOscParams(in,delta,3);
  in=glbSetOscParams(in,dms,4);
  in=glbSetOscParams(in,dma,5);
  return in;
}


glb_params 
glbCopyParams(const glb_params source, glb_params dest)
{
  int i;
  const glb_params_type *in;
  glb_params_type *out;
  in=(glb_params_type *) source;
  out=(glb_params_type *) dest;
  
  
  if(out==NULL||in==NULL) return NULL;
  if(out==in) return (glb_params) out;
  out->iterations=in->iterations;
  
  /* BUG #13 -- The use of realloc should solve a memory leak and makes this
     function pretty featureful, since it supports copying parameters
     even if the number of parameters and/or experiments has changed.
  */

  if((in->osc)->length>0&&
     (in->osc)->length!=(out->osc)->length)
    {
      (out->osc)->osc_params=glb_realloc((out->osc)->osc_params,(in->osc)->length * sizeof(double));
      (out->osc)->length=(in->osc)->length;
    }

  if((in->density)->length>0&&
     (in->density)->length!=(out->density)->length)
    {
      (out->density)->density_params=
	glb_realloc((out->density)->density_params,(in->density)->length * sizeof(double));
      (out->density)->length=(in->density)->length;
    }
  
  for(i=0;i<(out->osc)->length;i++) 
    (out->osc)->osc_params[i]=(in->osc)->osc_params[i];

  for(i=0;i<(out->density)->length;i++) 
    (out->density)->density_params[i]=(in->density)->density_params[i];

  return (glb_params) out;
  
}

void glbPrintParams(FILE *stream, const glb_params in)
{
  size_t k;
  const glb_params_type *source;
  source=(glb_params_type *) in;
  for(k=0;k<(source->osc)->length;k++) 
	fprintf(stream,"%g ",(source->osc)->osc_params[k]);   
  fprintf(stream,"\n");
  for(k=0;k<(source->density)->length;k++) 
	fprintf(stream,"%g ",(source->density)->density_params[k]);    
  fprintf(stream,"\n");
  fprintf(stream,"Iterations: %d\n",source->iterations);
}

/* This  is a bunch of function in order to deal with the
 * glb_projection_type type
 */

glb_osc_proj_type *glb_alloc_osc_proj_type()
{
  int i;
  glb_osc_proj_type *temp;
  temp=(glb_osc_proj_type *) glb_malloc(sizeof(glb_osc_proj_type));
  temp->osc_params=(int *) glb_malloc(sizeof(int) *  GLB_OSCP); 
  temp->length = (size_t) GLB_OSCP;
  for(i=0;i<GLB_OSCP;i++) temp->osc_params[i]=GLB_FREE;
  return temp;
}

void glb_free_osc_proj_type(glb_osc_proj_type *stale)
{
  if(stale!=NULL) glb_free(stale->osc_params);
  glb_free(stale);
}

glb_density_proj_type *glb_alloc_density_proj_type()
{
  int i;
  glb_density_proj_type *temp;
  temp=(glb_density_proj_type *) glb_malloc(sizeof(glb_density_proj_type));
  temp->density_params=(int *) glb_malloc(sizeof(int)*glb_num_of_exps);
  temp->length=(size_t) glb_num_of_exps;
  for(i=0;i<glb_num_of_exps;i++) temp->density_params[i]=GLB_FREE;
  return temp;
}

void glb_free_density_proj_type(glb_density_proj_type *stale)
{
  if(stale!=NULL) glb_free(stale->density_params);
  glb_free(stale);
}

glb_projection glbAllocProjection()
{
  glb_projection_type *temp;
  temp=(glb_projection_type *) glb_malloc(sizeof(glb_projection_type));
  temp->osc=glb_alloc_osc_proj_type();
  temp->density=glb_alloc_density_proj_type();
  return (glb_projection) temp;
}

void glbFreeProjection(glb_projection stale)
{
  glb_free_osc_proj_type(((glb_projection_type *) stale)->osc);
  glb_free_density_proj_type(((glb_projection_type *) stale)->density);
  /* BUG #14 -- This line solves a serious leak in connection with the
   * prior-module.
   */
  glb_free(stale);
  stale = NULL;
}


glb_projection glbSetDensityProjectionFlag(glb_projection in, 
				       int flag, int which)
{
  int i;
  glb_projection_type *s;
  if((flag!=GLB_FREE)&&(flag!=GLB_FIXED)) {
    glb_error("Projection flag must be either GLB_FREE or GLB_FIXED");
    return NULL;
  }

  s=(glb_projection_type *) in;
  if(which==GLB_ALL)
    {
      for(i=0;i<glb_num_of_exps;i++) (s->density)->density_params[i]=flag;
    }
  else if(0 <= which&&which < (s->osc)->length )
    {
      (s->density)->density_params[which]=flag;
    }
  else 
    {
      glb_error("Density list length mismatch");
      return NULL;
    }
  return (glb_projection) s;
}

int glbGetDensityProjectionFlag(const glb_projection in, int which)
{
  int out;
  glb_projection_type *s;
  s=(glb_projection_type *) in;
  if(0 <= which&&which < (s->density)->length )
    {
      out=(s->density)->density_params[which];
    }
  else 
    {
      glb_error("Density list length mismatch");
      return -1;
    }
  return out;
}

glb_projection glbSetProjectionFlag(glb_projection in, 
				 int flag, int which)
{

  glb_projection_type *s;
  if((flag!=GLB_FREE)&&(flag!=GLB_FIXED)) {
    glb_error("Projection flag must be either GLB_FREE or GLB_FIXED");
    return NULL;
  }

  s=(glb_projection_type *) in;
   if(0 <= which&&which < GLB_OSCP )
    {
      (s->osc)->osc_params[which]=flag;
    }
  else 
    {
      glb_error("Oscillation list length mismatch");
      return NULL;
    }
  return (glb_projection) s;
}

int glbGetProjectionFlag(const glb_projection in, int which)
{
  int out;
  glb_projection_type *s;
  s=(glb_projection_type *) in;
 if(0 <= which&&which < GLB_OSCP )
    {
      out=(s->osc)->osc_params[which];
    }
  else 
    {
      glb_error("Oscillation list length mismatch");
      return -1;
    }
  return out;
}


glb_projection 
glbDefineProjection(glb_projection in,int theta12, int theta13, int theta23,
	      int delta, int dms, int dma)
{
  in=glbSetProjectionFlag(in,theta12,0);
  in=glbSetProjectionFlag(in,theta13,1);
  in=glbSetProjectionFlag(in,theta23,2);
  in=glbSetProjectionFlag(in,delta,3);
  in=glbSetProjectionFlag(in,dms,4);
  in=glbSetProjectionFlag(in,dma,5);
  return in;
}


glb_projection
glbCopyProjection(const glb_projection source, glb_projection dest)
{
  int i;
  const glb_projection_type *in;
  glb_projection_type *out;
  in=(glb_projection_type *) source;
  out=(glb_projection_type *) dest;

  if(out==NULL||in==NULL) return NULL;
  if(out==in) return (glb_projection) out;

  /* BUG #13 -- The use of realloc should solve a memory leak and makes this
     function pretty featureful, since it supports copying parameters
     even if the number of parameters and/or experiments has changed.
  */

  if((in->osc)->length>0&&
     (in->osc)->length!=(out->osc)->length)
    {
      (out->osc)->osc_params=glb_realloc((out->osc)->osc_params,(in->osc)->length * sizeof(int));
      (out->osc)->length=(in->osc)->length;
    }

  if((in->density)->length>0&&
     (in->density)->length!=(out->density)->length)
    {
      (out->density)->density_params=
	glb_realloc((out->density)->density_params,(in->density)->length * sizeof(int));
      (out->density)->length=(in->density)->length;
    }

  for(i=0;i<(in->osc)->length;i++) 
    (out->osc)->osc_params[i]=(in->osc)->osc_params[i];

  for(i=0;i<(in->density)->length;i++) 
    (out->density)->density_params[i]=(in->density)->density_params[i];

  return (glb_projection) out;
}


void glbPrintProjection(FILE *stream, const glb_projection in)
{
  size_t k;
  const glb_projection_type *source;
  source=(glb_projection_type *) in;
  for(k=0;k<(source->osc)->length;k++) 
    {
      if((source->osc)->osc_params[k]==GLB_FREE) 
	fprintf(stream,"free  ");
      if((source->osc)->osc_params[k]==GLB_FIXED) 
	fprintf(stream,"fixed ");
      if(((source->osc)->osc_params[k]!=GLB_FIXED) &&
	 ((source->osc)->osc_params[k]!=GLB_FREE) )
	fprintf(stream,"undef ");
    }
  fprintf(stream,"\n");
  for(k=0;k<(source->density)->length;k++) 
    {
      if((source->density)->density_params[k]==GLB_FREE) 
	fprintf(stream,"free  ");
      if((source->density)->density_params[k]==GLB_FIXED) 
	fprintf(stream,"fixed ");
      if(((source->density)->density_params[k]!=GLB_FIXED) &&
	 ((source->density)->density_params[k]!=GLB_FREE) )
	fprintf(stream,"undef ");
    }
  fprintf(stream,"\n");
}


/* Wrappers for the old functions dealing with setting the oscillation 
 * parameters
 */


int 
glbSetOscillationParameters(const glb_params in)
{
  int i;
  double nsp[GLB_OSCP-6+1];
  if(in==NULL) return -1;
  glb_set_c_vacuum_parameters(glbGetOscParams(in,0),glbGetOscParams(in,1)
		       ,glbGetOscParams(in,2),glbGetOscParams(in,3));
  glb_set_c_squared_masses(0,glbGetOscParams(in,4),glbGetOscParams(in,5));
  if(GLB_OSCP>6)
  {
	for(i=0;i<GLB_OSCP-6;i++) nsp[i]=glbGetOscParams(in,6+i);
	glb_set_c_ns_params(nsp);
  }
  for(i=0;i<glb_num_of_exps;i++) glb_set_profile_scaling(glbGetDensityParams(in,i),i);
  return 0;

}

int
glbSetStartingValues(const glb_params in)
{
  int i,s=0;
  double nsp[GLB_OSCP-6+1];
  if(in==NULL) return -1;

  glb_set_solar_starting_values(glbGetOscParams(in,0));
  glb_set_starting_values(glbGetOscParams(in,1),glbGetOscParams(in,2),
		    glbGetOscParams(in,3),glbGetOscParams(in,4),
		    glbGetOscParams(in,5),0);
  if(GLB_OSCP>6)
  {
	for(i=0;i<GLB_OSCP-6;i++) nsp[i]=glbGetOscParams(in,6+i);
	glb_set_ns_starting_values(nsp);
  }

  for(i=0;i<glb_num_of_exps;i++) 
    glbSetDensityStartingValue(glbGetDensityParams(in,i),i);
  s=glb_user_defined_starting_values(in);
  return s;

}

int
glbSetInputErrors(const glb_params in)
{
  int i,s=0;
  double nsp[GLB_OSCP-6+1];
  if(in==NULL) return -1;

  glb_set_solar_input_errors(glbGetOscParams(in,0));
  glb_set_input_errors(glbGetOscParams(in,1),glbGetOscParams(in,2),
		    glbGetOscParams(in,3),glbGetOscParams(in,4),
		    glbGetOscParams(in,5),0);
  if(GLB_OSCP>6)
  {
	for(i=0;i<GLB_OSCP-6;i++) nsp[i]=glbGetOscParams(in,6+i);
	glb_set_ns_input_errors(nsp);
  }

  /* FIXME - the density glb_error.has the same size than the
   * density_center, if the user does not use glbSetDensityParams.
   * Furthermore the setting in the gls-file is ignored...!
   */
  for(i=0;i<glb_num_of_exps;i++) 
    glbSetDensityInputError(glbGetDensityParams(in,i),i);  
  s=glb_user_defined_input_errors(in);
  return s;

}


int 
glbGetOscillationParameters(glb_params in)
{
  int i;
  double *t;
  if(in==NULL) return -1;
  t=glb_get_vacuum_parameters();
  for(i=0;i<4;i++) glbSetOscParams(in,t[i],i);
  glb_free(t);
  t=glb_get_squared_masses();
  for(i=0;i<2;i++) glbSetOscParams(in,t[i],i+4);
  glb_free(t);
  if(GLB_OSCP>6)
  {  
	  t=glb_get_ns_params();
      for(i=0;i<GLB_OSCP-6;i++) glbSetOscParams(in,t[i],i+6);
      glb_free(t);
  }
  return 0;
}


int
glbGetStartingValues(glb_params in)
{
  int i;
  double *t;
  if(in==NULL) return -1;
  t=glb_return_input_values();
  for(i=0;i<GLB_OSCP;i++) glbSetOscParams(in,t[i],i);
  for(i=0;i<glb_num_of_exps;i++) glbSetDensityParams(in,t[i+GLB_OSCP],i);
  glb_free(t);
  return 0;

}

int
glbGetInputErrors(glb_params in)
{
  int i;
  double *t;
  if(in==NULL) return -1;
  t=glb_return_input_errors();
  for(i=0;i<GLB_OSCP;i++) glbSetOscParams(in,t[i],i);
  for(i=0;i<glb_num_of_exps;i++) glbSetDensityParams(in,t[i+GLB_OSCP],i);
  glb_free(t);
  return 0;

}




/* Here we have some "memory managament" functions, ie. things
 * to be called at start-up or exit...
 */

static struct obstack glb_rate_stack;

void
glbClearExperimentList()
{
  int i;
  for(i=0;i<32;i++) glbFreeExp(glb_experiment_list[i]);
  for(i=0;i<32;i++) glb_experiment_list[i]=glbAllocExp();
  glb_num_of_exps=0;
  glbResetCounters();
}

int
glbGetChannelRates(double **data, size_t *length,
		   int exp, int channel, int smearing)
{
  size_t l;
  int i;
  double *temp;
  if(exp>=glb_num_of_exps) 
    {
      glb_error("Error in glbGetChannelRates: 4th argument is larger than"
		" glb_num_of_exps");
      *data=NULL;
      *length=0;
      return -1;
    }
  if(channel>= ((struct glb_experiment *) glb_experiment_list[exp])->numofchannels) 
    {
      glb_error("Error in glbGetChannelRates: 5th argument is larger than"
		" numofchannels");
      *data=NULL;
      *length=0;
      return -1;
    }

  if(exp<0) 
    {
      glb_error("Error in glbGetChannelRates: 4th argument must be larger than"
		" 0");
      *data=NULL;
      *length=0;
      return -1;
    }
 if(channel<0) 
    {
      glb_error("Error in glbGetChannelRates: 5th argument must be larger than"
		" 0");
      *data=NULL;
      *length=0;
      return -1;
    }
 if((smearing!=GLB_PRE) && (smearing!=GLB_POST))
   {
     glb_error("Error in glbGetChannelRates: 3rd argument must be either"
	       " GLB_PRE or GLB_POST");
     *data=NULL;
     *length=0;
     return -1;
   }


  if(smearing==GLB_PRE) l=(size_t)
			  ((struct glb_experiment *) glb_experiment_list[exp])->simbins;
  else l=(size_t) 
	 ((struct glb_experiment *) glb_experiment_list[exp])->numofbins;
  if(l<=0)
    {
      glb_error("Error in glbGetChannelRates: Could not determine rate vector"
		" length");
      *data=NULL;
      *length=0;
      return -1;
    }

  temp=(double *) obstack_alloc(&glb_rate_stack,(l+1) * sizeof(double));
  if(smearing==GLB_PRE)
    for(i=0;i<l;i++) 
      temp[i]= ((struct glb_experiment *) glb_experiment_list[exp])->chrb[channel][i];
  else
    for(i=0;i<l;i++) 
      temp[i]= ((struct glb_experiment *) glb_experiment_list[exp])->chra[channel][i];

  temp[l]=-1.0;

  *data=temp;
  *length=l;
  return 0;
}



int
glbGetUserData(double **data, size_t *length, 
		   int exp, int channel, int smearing, int bgeff)
{
  size_t l;
  int i;
  double *temp;
  if(exp>=glb_num_of_exps) 
    {
      glb_error("Error in glbGetUserData: 5th argument is larger than"
		" glb_num_of_exps");
      *data=NULL;
      *length=0;
      return -1;
    }
  if(channel>= ((struct glb_experiment *) glb_experiment_list[exp])->numofchannels) 
    {
      glb_error("Error in glbGetUserData: 6th argument is larger than"
		" numofchannels");
      *data=NULL;
      *length=0;
      return -1;
    }

  if(exp<0) 
    {
      glb_error("Error in glbGetUserData: 5th argument must be larger than"
		" 0");
      *data=NULL;
      *length=0;
      return -1;
    }
 if(channel<0) 
    {
      glb_error("Error in glbGetUserData: 6th argument must be larger than"
		" 0");
      *data=NULL;
      *length=0;
      return -1;
    }
 if((smearing!=GLB_PRE) && (smearing!=GLB_POST))
   {
     glb_error("Error in glbGetUserData: 3rd argument must be either"
	       " GLB_PRE or GLB_POST");
     *data=NULL;
     *length=0;
     return -1;
   }
 if((bgeff!=GLB_EFF) && (bgeff!=GLB_BG))
   {
     glb_error("Error in glbGetUserData: 4th argument must be either"
	       " GLB_EFF or GLB_BG");
     *data=NULL;
     *length=0;
     return -1;
   }


  if(smearing==GLB_PRE) l=(size_t)
			  ((struct glb_experiment *) glb_experiment_list[exp])->simbins;
  else l=(size_t) 
	 ((struct glb_experiment *) glb_experiment_list[exp])->numofbins;
  if(l<=0)
    {
      glb_error("Error in glbGetUserData: Could not determine rate vector"
		" length");
      *data=NULL;
      *length=0;
      return -1;
    }

  temp=(double *) obstack_alloc(&glb_rate_stack,(l+1) * sizeof(double));
  if(smearing==GLB_PRE)
    {
      if(bgeff==GLB_EFF)
	for(i=0;i<l;i++) 
	  temp[i]= 
	    ((struct glb_experiment *) 
	     glb_experiment_list[exp])->user_pre_smearing_channel[channel][i];
      else
	for(i=0;i<l;i++) 
	  temp[i]= 
	    ((struct glb_experiment *) 
	     glb_experiment_list[exp])->user_pre_smearing_background[channel][i];
    }  
  else
    { 
      if(bgeff==GLB_EFF)
	for(i=0;i<l;i++) 
	  temp[i]= ((struct glb_experiment *) 
		    glb_experiment_list[exp])->user_post_smearing_channel[channel][i];
      else
	for(i=0;i<l;i++) 
	  temp[i]= ((struct glb_experiment *) 
		    glb_experiment_list[exp])->user_post_smearing_background[channel][i];
	
    }
  temp[l]=-1.0;

  *data=temp;
  *length=l;
  return 0;
}


void glbResetRateStack()
{
  obstack_free(&glb_rate_stack,NULL);
  obstack_init(&glb_rate_stack);
}

/* Some bound checking */

static int exp_range(int exp)
{
  if((exp>=0)&&(exp<glb_num_of_exps)) return 0;
  glb_error("Experiment index out of range");
  return -1;
}


static int channel_range(int exp, int channel)
{
  if(exp_range(exp)!=0) return -1;
  if((channel>=0)&&
     (channel<((struct glb_experiment *) glb_experiment_list[exp])->numofchannels)) return 0;
  glb_error("Channel index out of range");
  return -1;
}

static int signal_range(int signal)
{
  if((signal==GLB_SIG)||(signal==GLB_BG)) return 0;
  glb_error("Signal flag is neither GLB_SIG nor GLB_BG");
  return -1;
}



static int rule_range(int exp, int rule)
{
  if(exp_range(exp)!=0) return -1;
  if((rule>=0)&&
     (rule<((struct glb_experiment *) glb_experiment_list[exp])->numofrules)) return 0;
  glb_error("Rule index out of range");
  return -1;
}


static int pos_range(int exp, int rule, int pos, int signal)
{
  if(rule_range(exp,rule)!=0) return -1;
  if(signal_range(signal)!=0) return -1;
  if(signal==GLB_SIG) 
    if((pos>=0)&&
       (pos<((struct glb_experiment *) 
	    glb_experiment_list[exp])->lengthofrules[rule])) return 0;
  
  if(signal==GLB_BG) 
    if((pos>=0)&&
       (pos<((struct glb_experiment *) 
	     glb_experiment_list[exp])->lengthofbgrules[rule])) return 0;
  
  glb_error("Position index out of range");
  return -1;
}



/* Accessing channel and rule internals */

int glbGetNumberOfChannels(int exp)
{
  int i;
  if(exp_range(exp)!=0) return -1;
  i=((struct glb_experiment *) glb_experiment_list[exp])->numofchannels;
  return i;
}

int glbGetNumberOfRules(int exp)
{
  int i;
  if(exp_range(exp)!=0) return -1;
  i=((struct glb_experiment *) glb_experiment_list[exp])->numofrules; 
  return i;
}

int glbGetLengthOfRule(int exp, int rule, int signal)
{
  int i=-1;
  if(signal_range(signal)!=0) return -1;
  if(rule_range(exp,rule)!=0) return -1;

  if(signal==GLB_SIG)
    i=((struct glb_experiment *) glb_experiment_list[exp])->lengthofrules[rule];

  if(signal==GLB_BG)
    i=((struct glb_experiment *) glb_experiment_list[exp])->lengthofbgrules[rule];
  
  return i;
}

int glbGetChannelInRule(int exp, int rule, int pos, int signal)
{
  int i=-1;
  if(pos_range( exp,  rule,  pos,  signal)!=0) return -1;
 
  if(signal==GLB_SIG)
    i=((struct glb_experiment *) glb_experiment_list[exp])->rulechannellist[rule][pos];

  if(signal==GLB_BG)
    i=((struct glb_experiment *) glb_experiment_list[exp])->bgrulechannellist[rule][pos];

  return i;
} 

double glbGetCoefficientInRule(int exp, int rule, int pos, int signal)
{
  double res=-1.0;
  if(pos_range(exp,rule,pos,signal)!=0) return -1.0;

  if(signal==GLB_SIG)
    res=((struct glb_experiment *) glb_experiment_list[exp])->rulescoeff[rule][pos];

  if(signal==GLB_BG)
    res=((struct glb_experiment *) glb_experiment_list[exp])->bgrulescoeff[rule][pos];

  return res;
} 

double glbGetNormalizationInRule(int exp, int rule, int signal)
{
  double res=-1.0;
  if(signal_range(signal)!=0) return -1.0;
  if(rule_range(exp,rule)!=0) return -1.0;  

  if(signal==GLB_SIG)
    res=1.0;

  if(signal==GLB_BG)
    res=((struct glb_experiment *) glb_experiment_list[exp])->bgcenter[0][rule];

  return res;
} 



static 
int channel_bg(const double *ch, const double *bg, double **res, int flag)
{
  size_t lch,lbg;
  int i;
  double *temp;
  if((ch==NULL)||(bg==NULL)) 
    {
      glb_error("Rate or background vector may not be NULL");
      return -1;
    }
  /* FIXME That`s not very safe */
  lch=0;
  while(ch[lch]!=-1.0) lch++;
  lbg=0;
  while(bg[lbg]!=-1.0) lbg++;
  if(lbg!=lch)
    {
      glb_error("Length of rate or background vector may not be different");
      return -1;
    }
  if((flag!=GLB_W_BG)&&(flag!=GLB_WO_BG))
    {
      glb_error("Background flag is neither glb_W_BG nor GLB_WO_BG");
      return -1;
    }
  temp=(double *) obstack_alloc(&glb_rate_stack,sizeof(double) * (lch+1));
  
  if(flag==GLB_W_BG) 
    for(i=0;i<lch;i++) temp[i]=ch[i];
  
  if(flag==GLB_WO_BG)
    for(i=0;i<lch;i++) temp[i]=ch[i]-bg[i];
  
  temp[lch]=-1.0;
  *res=temp;
  return 0;
}

static int 
channel_eff(const double *ch, const double *bg, double **res, int flag)
{
  size_t lch,lbg;
  int i;
  double *temp;
  if((ch==NULL)||(bg==NULL)) 
    {
      glb_error("Rate or background vector may not be NULL");
      return -1;
    }
  /* FIXME That`s not very safe */
  lch=0;
  while(ch[lch]!=-1.0) lch++;
  lbg=0;
  while(bg[lbg]!=-1.0) lbg++;
  if(lbg!=lch)
    {
      glb_error("Length of rate or background vector may not be different");
      return -1;
    }
  if((flag!=GLB_W_EFF)&&(flag!=GLB_WO_EFF))
    {
      glb_error("Background flag is neither GLB_W_EFF nor GLB_WO_EFF");
      return -1;
    }
  temp=(double *) obstack_alloc(&glb_rate_stack,sizeof(double) * (lch+1));
  
  if(flag==GLB_W_EFF) 
    for(i=0;i<lch;i++) temp[i]=ch[i];
  
  if(flag==GLB_WO_EFF)
    for(i=0;i<lch;i++) temp[i]=ch[i]/bg[i];
  temp[lch]=-1.0;
  *res=temp;
  return 0;
}


static int 
rule_coeff(const double *ch, double coeff, double **res, int flag)
{
  size_t lch;
  int i;
  double *temp;
  if((ch==NULL)) 
    {
      glb_error("Rate or background vector may not be NULL");
      return -1;
    }
  /* FIXME That`s not very safe */
  lch=0;
  while(ch[lch]!=-1.0) lch++;

  if((flag!=GLB_W_COEFF)&&(flag!=GLB_WO_COEFF))
    {
      glb_error("Coefficient flag is neither glb_W_COEFF nor glb_WO_COEFF");
      return -1;
    }
  temp=(double *) obstack_alloc(&glb_rate_stack,sizeof(double) * (lch+1));
  
  if(flag==GLB_W_COEFF) 
    for(i=0;i<lch;i++) temp[i]=ch[i]*coeff;
  
  if(flag==GLB_WO_COEFF)
    for(i=0;i<lch;i++) temp[i]=ch[i];
  temp[lch]=-1.0;
  *res=temp;
  return 0;
}



static glb_smear 
*get_channel_smear_data(int exp, int channel)
{
  glb_smear *p;
  size_t w;
  if(channel_range(exp,channel)!=0) return NULL;
  w=(size_t) ((struct glb_experiment *) glb_experiment_list[exp])->listofchannels[5][channel];
  p=((struct glb_experiment *) glb_experiment_list[exp])->smear_data[w];
  return p;
}

void (*glb_channel_print_function)(FILE *stream,const double *energy, 
				   double **res, size_t l, size_t c);

int 
glbShowChannelRates(FILE *stream,
			int exp, int channel, int smearing, int effi, int bgi)
{
  int i,s,cc;
  size_t k,l,m,c;
  double *ch,*bg,*eff,*temp,**res,*energy,sum;
  glb_smear *smt;
  s=0;
  sum=0.0;
  cc=(size_t) glbGetNumberOfChannels(exp);
  if(channel!=GLB_ALL)
    {
      res=(double **) obstack_alloc(&glb_rate_stack,sizeof(double *)*2);
      s+=glbGetChannelRates(&ch,&l,exp,channel,smearing);
      s+=glbGetUserData(&bg,&k,exp,channel,smearing,GLB_BG);
      s+=glbGetUserData(&eff,&m,exp,channel,smearing,GLB_EFF);
      if(s!=0) return -1;
      if(!((l==k)&&(k==m))) return -1;
 
      s+=channel_bg(ch,bg,&temp,bgi);
      s+=channel_eff(temp,eff,&res[0],effi);
      res[1]=NULL;
      c=1;
    }
  else
    {
      res=(double **) obstack_alloc(&glb_rate_stack,sizeof(double *)*(cc+1));
      for(i=0;i<cc;i++)
	{
	  s+=glbGetChannelRates(&ch,&l,exp,i,smearing);
	  s+=glbGetUserData(&bg,&k,exp,i,smearing,GLB_BG);
	  s+=glbGetUserData(&eff,&m,exp,i,smearing,GLB_EFF);
	  if(s!=0) return -1;
	  if(!((l==k)&&(k==m))) return -1;
	  
	  s+=channel_bg(ch,bg,&temp,bgi);
	  s+=channel_eff(temp,eff,&res[i],effi);
	  
	}
      res[cc]=NULL;
      c=cc;
    }

  if(s!=0) return -1;
  if(channel!=GLB_ALL) smt=get_channel_smear_data(exp,channel);
  if(channel==GLB_ALL) smt=get_channel_smear_data(exp,0);
  energy=(double *) obstack_alloc(&glb_rate_stack,sizeof(double)*(l+1));
  if(smearing==GLB_PRE)
      for(i=0;i<l;i++) energy[i]=glb_sbin_center(i,smt);
  if(smearing==GLB_POST)
    for(i=0;i<l;i++) energy[i]=glb_bin_center(i,smt);
  energy[l]=-1.0;

  glb_channel_print_function(stream,energy,res,l,c);

  obstack_free(&glb_rate_stack,ch);
  return 0;
}



int 
glbShowRuleRates(FILE *stream,
		 int exp, int rule, int pos,
		 int effi, int bgi, int coeffi, 
		 int signal)
{
  int i,s,cc,channel;
  size_t k,l,m,c;
  double *ch,*bg,*eff,*temp,*ceff,**res,*energy,sum,coeff;
  glb_smear *smt;
  s=0;
  sum=0.0;

  cc=(size_t) glbGetLengthOfRule(exp,rule,signal);
  if(pos!=GLB_ALL)
    {
      res=(double **) obstack_alloc(&glb_rate_stack,sizeof(double *)*2);
      channel=glbGetChannelInRule(exp,rule,pos,signal);
      coeff=glbGetCoefficientInRule(exp,rule,pos,signal)*
	glbGetNormalizationInRule(exp,rule,signal);
      s+=glbGetChannelRates(&ch,&l,exp,channel,GLB_POST);
      s+=glbGetUserData(&bg,&k,exp,channel,GLB_POST,GLB_BG);
      s+=glbGetUserData(&eff,&m,exp,channel,GLB_POST,GLB_EFF);
      if(s!=0) return -1;
      if(!((l==k)&&(k==m))) return -1;
      s+=channel_bg(ch,bg,&temp,bgi);
      s+=channel_eff(temp,eff,&ceff,effi);
      s+=rule_coeff(ceff,coeff,&res[0],coeffi);
      res[1]=NULL;
      c=1;
    }
  else
    { 
      res=(double **) obstack_alloc(&glb_rate_stack,sizeof(double *)*(cc+1));
      for(i=0;i<cc;i++)
	{
	  channel=glbGetChannelInRule(exp,rule,i,signal);
	  coeff=glbGetCoefficientInRule(exp,rule,i,signal)*
	    glbGetNormalizationInRule(exp,rule,signal);
	 
	  s+=glbGetChannelRates(&ch,&l,exp,channel,GLB_POST);
	  s+=glbGetUserData(&bg,&k,exp,channel,GLB_POST,GLB_BG);
	  s+=glbGetUserData(&eff,&m,exp,channel,GLB_POST,GLB_EFF);
	  if(s!=0) return -1;
	  if(!((l==k)&&(k==m))) return -1;
	  s+=channel_bg(ch,bg,&temp,bgi);
	  s+=channel_eff(temp,eff,&ceff,effi);
	  s+=rule_coeff(ceff,coeff,&res[i],coeffi);
	}
      res[cc]=NULL;
      c=cc;
    }


  if(s!=0) return -1;
  smt=get_channel_smear_data(exp,0);
  
  energy=(double *) obstack_alloc(&glb_rate_stack,sizeof(double)*(l+1));
 
  for(i=0;i<l;i++) energy[i]=glb_bin_center(i,smt);
  energy[l]=-1.0;

  glb_channel_print_function(stream,energy,res,l,c);

  obstack_free(&glb_rate_stack,ch);
  return 0;
}




double 
glbTotalRuleRate(
		 int exp, int rule, int pos,
		 int effi, int bgi, int coeffi, 
		 int signal)
{
  double out;
  int i,s,cc,channel;
  size_t k,l,m,c,bins;
  double *ch,*bg,*eff,*temp,*ceff,**res,*energy,sum,coeff;
  glb_smear *smt;
  s=0;
  sum=0.0;
  bins=glb_experiment_list[exp]->numofbins;
  cc=(size_t) glbGetLengthOfRule(exp,rule,signal);
  if(pos!=GLB_ALL)
    {
      res=(double **) obstack_alloc(&glb_rate_stack,sizeof(double *)*2);
      channel=glbGetChannelInRule(exp,rule,pos,signal);
      coeff=glbGetCoefficientInRule(exp,rule,pos,signal)*
	glbGetNormalizationInRule(exp,rule,signal);
      s+=glbGetChannelRates(&ch,&l,exp,channel,GLB_POST);
      s+=glbGetUserData(&bg,&k,exp,channel,GLB_POST,GLB_BG);
      s+=glbGetUserData(&eff,&m,exp,channel,GLB_POST,GLB_EFF);
      if(s!=0) return -1;
      if(!((l==k)&&(k==m))) return -1;
      s+=channel_bg(ch,bg,&temp,bgi);
      s+=channel_eff(temp,eff,&ceff,effi);
      s+=rule_coeff(ceff,coeff,&res[0],coeffi);
      res[1]=NULL;
      c=1;
    }
  else
    { 
      res=(double **) obstack_alloc(&glb_rate_stack,sizeof(double *)*(cc+1));
      for(i=0;i<cc;i++)
	{
	  channel=glbGetChannelInRule(exp,rule,i,signal);
	  coeff=glbGetCoefficientInRule(exp,rule,i,signal)*
	    glbGetNormalizationInRule(exp,rule,signal);
	 
	  s+=glbGetChannelRates(&ch,&l,exp,channel,GLB_POST);
	  s+=glbGetUserData(&bg,&k,exp,channel,GLB_POST,GLB_BG);
	  s+=glbGetUserData(&eff,&m,exp,channel,GLB_POST,GLB_EFF);
	  if(s!=0) return -1;
	  if(!((l==k)&&(k==m))) return -1;
	  s+=channel_bg(ch,bg,&temp,bgi);
	  s+=channel_eff(temp,eff,&ceff,effi);
	  s+=rule_coeff(ceff,coeff,&res[i],coeffi);
	}
      res[cc]=NULL;
      c=cc;
    }


  if(s!=0) return -1;
  out=0;
  for(i=0;i<c;i++)
    for(k=0;res[i][k]!=-1;k++) out += res[i][k];

  obstack_free(&glb_rate_stack,ch);
  return out;
}


static char *printf_left=NULL;
static char *printf_middle=NULL;
static char *printf_right=NULL;

void glbSetPrintDelimiters(const char *left,const char *middle,
			   const char *right)
{


  glb_free((char *) printf_left);
  glb_free((char *) printf_right);
  glb_free((char *) printf_middle);

  printf_left=(char *) strdup(left);
  printf_middle=(char *) strdup(middle);
  printf_right=(char *) strdup(right);

}

void  glbPrintDelimiter(FILE *stream, int character)
{

  if(character=='l') fprintf(stream,"%s",printf_left);
  if(character=='m') fprintf(stream,"%s",printf_middle);
  if(character=='r') fprintf(stream,"%s",printf_right);
  return;

}


static void glb_builtin_channel_printf(FILE *stream,
				       const double *energy,
				       const double **res, size_t l,size_t c)
{
  int i,k;
  double *sum;
  sum=(double *) glb_malloc(sizeof(double)*c);
  for(k=0;k<c;k++) sum[k]=0.0;
  fprintf(stream,"\n");
  fprintf(stream,"%s",printf_left);
  for(i=0;i<l;i++)
    { 
      fprintf(stream,"%s",printf_left);
      fprintf(stream,"%6.4g%s",energy[i],printf_middle);
	for(k=0;k<c;k++)
	  {
	    fprintf(stream,"%12.6g%s",res[k][i],printf_middle);
	    sum[k]+=res[k][i];
	  }
	fprintf(stream,"%s",printf_right);

    }
   
  fprintf(stream,"----------------------");
  for(k=1;k<c;k++)
    {
      fprintf(stream,"----------------");
    }
  fprintf(stream,"\nTotal:\t");
  
  for(k=0;k<c;k++)
    {
      fprintf(stream,"%12.6g%s",sum[k],printf_middle);
    }

  fprintf(stream,"%s",printf_right);
  glb_free(sum);
}

void 
*glbSetChannelPrintFunction(void *fp)
{
  if(fp==NULL) return glb_channel_print_function;
  glb_channel_print_function=fp;
  return glb_channel_print_function;
}

void glb_clean_up()
{
  int i;
 for(i=0;i<32;i++) glbFreeExp(glb_experiment_list[i]);
 glb_clean_parser();
 glb_lexer_cleanup();
 obstack_free(&glb_rate_stack,NULL);
 glb_free((char *) printf_left);
 glb_free((char *) printf_right);
 glb_free((char *) printf_middle);
 glb_free((char *) glb_prog_name);
 for(i=0;i<glb_path_vector_length;i++)
   glb_free((char *) glb_path_vector[i]);
 glb_free((char **) glb_path_vector);
 
}



/* The all important init function */

void 
glb_init(char *name)
{
  int i;
  
 
  glb_prog_name_init(name);
  glb_setup_path();
 
  for(i=0;i<32;i++) glb_experiment_list[i]=glbAllocExp();
  glb_num_of_exps=0;
  glb_rule_number=0;
  obstack_init(&glb_rate_stack);
  glbSetPrintDelimiters("","\t","\n");
  glbSetChannelPrintFunction(glb_builtin_channel_printf);
  /* this has gone to glbSuperInit 
   * glbLoadPrior("prior-template");
   */
}

/* Toggle Systemtatics */

static int
glb_switch_systematics(int experiment, int rule, int on_off)
{
  int i;
  struct glb_experiment *in;
  if((on_off!=GLB_ON)&&(on_off!=GLB_OFF)) {
    glb_error("Invalid value for on_off");
    return -1;
  }

  in=(struct glb_experiment *) glb_experiment_list[experiment];
  if(rule==GLB_ALL)
    {
      for(i=0;i<in->numofrules;i++) {
	if(on_off==GLB_ON) in->errordim[i]=in->errordim_sys_on[i];
	if(on_off==GLB_OFF) in->errordim[i]=in->errordim_sys_off[i];
      }
    }
  else if((rule >= 0)&&(rule < in->numofrules ))
    {
      i=rule;
      if(on_off==GLB_ON) in->errordim[i]=in->errordim_sys_on[i];
      if(on_off==GLB_OFF) in->errordim[i]=in->errordim_sys_off[i];
    }
  else 
    {
      glb_error("Invalid value for rule number");
      return -1;
    }
  return 0;
}

int
glbSwitchSystematics(int experiment, int rule, int on_off)
{
  int i,s=0;
  if(experiment==GLB_ALL)
    {
      for(i=0;i<glb_num_of_exps;i++) 
	s+=glb_switch_systematics(i, rule, on_off);
    }
  else if((experiment >= 0)&&(experiment < glb_num_of_exps))
    {
	s+=glb_switch_systematics(experiment, rule, on_off);
    }
  else 
    {
      glb_error("Invalid value for experiment number");
      return -1;
    }
  return s;
}

static int
glb_set_systematics(int experiment, int rule, int on_off, int value)
{
  int i;
  struct glb_experiment *in;
  if((on_off!=GLB_ON)&&(on_off!=GLB_OFF)) {
    glb_error("Invalid value for on_off");
    return -1;
  }

  in=(struct glb_experiment *) glb_experiment_list[experiment];
  if(rule==GLB_ALL)
    {
      for(i=0;i<in->numofrules;i++) {
	if(on_off==GLB_ON) in->errordim_sys_on[i]=value;
	if(on_off==GLB_OFF) in->errordim_sys_off[i]=value;
      }
    }
  else if((rule >= 0)&&(rule < in->numofrules ))
    {
      i=rule;
      if(on_off==GLB_ON) in->errordim_sys_on[i]=value;
      if(on_off==GLB_OFF) in->errordim_sys_off[i]=value;
    }
  else 
    {
      glb_error("Invalid value for rule number");
      return -1;
    }
  return 0;
}


int 
glbSetErrorDim(int experiment, int rule, int on_off, int value)
{
  int i,s=0;
  if(experiment==GLB_ALL)
    {
      for(i=0;i<glb_num_of_exps;i++) 
	s+=glb_set_systematics(i, rule, on_off, value);
    }
  else if((experiment >= 0)&&(experiment < glb_num_of_exps))
    {
	s+=glb_set_systematics(experiment, rule, on_off,value);
    }
  else 
    {
      glb_error("Invalid value for experiment number");
      return -1;
    }
  return s;
}

int 
glbGetErrorDim(int experiment, int rule, int on_off)
{
  int s=0;
  struct glb_experiment *in;
  if((on_off!=GLB_ON)&&(on_off!=GLB_OFF)) {
    glb_error("Invalid value for on_off");
    return -1;
  }
  if((experiment >= 0)&&(experiment < glb_num_of_exps))
    {
      in=(struct glb_experiment *) glb_experiment_list[experiment];
      if((rule >= 0)&&(rule < in->numofrules )) 
	{
	  if(on_off==GLB_ON) s=in->errordim_sys_on[rule];
	  if(on_off==GLB_OFF) s=in->errordim_sys_off[rule];
	}
    }
  else 
    {
      glb_error("Invalid value for experiment number");
      return -1;
    }
  return s;
}


/* Here comes a bunch of set/get functions.
 * They all have the following ordering of their arguments:
 *    glbSetXXX(experiment, rule, flags, XXX)
 */


int 
glbSetTargetMass(int experiment,double mass)
{
  struct glb_experiment *in;
  int i;
  if(experiment==GLB_ALL)
    {
      if(mass > 0)
	{
	  for(i=0;i<glb_num_of_exps;i++)
	    {
	      in=(struct glb_experiment *) glb_experiment_list[i];
	      in->targetmass = mass;
	    }  
	}
      else
	{
	  glb_error("Target mass has to be positive");
	  return -1;
	}
      return 0;
    }

  if((experiment >= 0)&&(experiment < glb_num_of_exps))
    {
      in=(struct glb_experiment *) glb_experiment_list[experiment];
      if(mass > 0)
	{
	  in->targetmass = mass;  
	}
      else
	{
	  glb_error("Target mass has to be positive");
	  return -1;
	}
    }
  else 
    {
      glb_error("Invalid value for experiment number");
      return -1;
    }
  return 0;

  
}

double 
glbGetTargetMass(int experiment)
{
  struct glb_experiment *in;
  double out;
  if((experiment >= 0)&&(experiment < glb_num_of_exps))
    {
      in=(struct glb_experiment *) glb_experiment_list[experiment];
      out=in->targetmass;
      return out;
    }
  else 
    {
      glb_error("Invalid value for experiment number");
      return -1.0;
    }
  return -1.0;
}

int glbSetSignalErrors(int experiment, int rule, double norm, double tilt)
{
  struct glb_experiment *in;
  int i,k;
  /* Testing the arguments */
  if((norm <= 0) || (tilt <= 0)) { glb_error("Errors have to be positive");
  return -1;}

  /* Testing the experiment number */
  if(!(((experiment >= 0)&&(experiment < glb_num_of_exps))
       ||(experiment==GLB_ALL))) { 
    glb_error("Invalid value for experiment number");
    return -1;}

  for(i=0;i<glb_num_of_exps;i++)
    {
      if(experiment!=GLB_ALL) i=experiment;
      in=(struct glb_experiment *) glb_experiment_list[i];
      /* Testing the rule number */
      if(!(((rule >= 0)&&(rule < in->numofrules))
	   ||(rule==GLB_ALL))) { 
	glb_error("Invalid value for rule number");
	return -1;}     
      for(k=0;k<in->numofrules;k++)
	{
	   if(rule!=GLB_ALL) k=rule;
	   /* Here should come the assignment */
	  
	   in->signalruleerror[0][k]=norm;
	   in->signalruleerror[1][k]=tilt;
	   
	   if(rule!=GLB_ALL) break;
	}

      if(experiment!=GLB_ALL) break;
    }
  return 0;
}


int 
glbGetSignalErrors(int experiment, int rule, double *norm, double *tilt)
{
  struct glb_experiment *in;
  int i,k;
  /* Testing the arguments */
  if((norm == NULL) || (tilt == NULL)) { 
    glb_error("Input pointers may not be NULL");
  return -1;}

  /* Testing the experiment number */
  if(!(((experiment >= 0)&&(experiment < glb_num_of_exps)))) { 
    glb_error("Invalid value for experiment number");
    return -1;}

  for(i=0;i<glb_num_of_exps;i++)
    {
      if(experiment!=GLB_ALL) i=experiment;
      in=(struct glb_experiment *) glb_experiment_list[i];
      /* Testing the rule number */
      if(!(((rule >= 0)&&(rule < in->numofrules)))) { 
	glb_error("Invalid value for rule number");
	return -1;}     
      for(k=0;k<in->numofrules;k++)
	{
	   if(rule!=GLB_ALL) k=rule;
	   /* Here should come the assignment */
	  
	   *norm=in->signalruleerror[0][k];
	   *tilt=in->signalruleerror[1][k];
	   
	   if(rule!=GLB_ALL) break;
	}

      if(experiment!=GLB_ALL) break;
    }
  return 0;
}

/* same for BG errors */

int glbSetBGErrors(int experiment, int rule, double norm, double tilt)
{
  struct glb_experiment *in;
  int i,k;
  /* Testing the arguments */
  if((norm <= 0) || (tilt <= 0)) { glb_error("Errors have to be positive");
  return -1;}

  /* Testing the experiment number */
  if(!(((experiment >= 0)&&(experiment < glb_num_of_exps))
       ||(experiment==GLB_ALL))) { 
    glb_error("Invalid value for experiment number");
    return -1;}

  for(i=0;i<glb_num_of_exps;i++)
    {
      if(experiment!=GLB_ALL) i=experiment;
      in=(struct glb_experiment *) glb_experiment_list[i];
      /* Testing the rule number */
      if(!(((rule >= 0)&&(rule < in->numofrules))
	   ||(rule==GLB_ALL))) { 
	glb_error("Invalid value for rule number");
	return -1;}     
      for(k=0;k<in->numofrules;k++)
	{
	   if(rule!=GLB_ALL) k=rule;
	   /* Here should come the assignment */
	   in->bgerror[0][k]=norm;
	   in->bgerror[1][k]=tilt;
	   
	   
	   if(rule!=GLB_ALL) break;
	}

      if(experiment!=GLB_ALL) break;
    }
  return 0;
}


int 
glbGetBGErrors(int experiment, int rule, double *norm, double *tilt)
{
  struct glb_experiment *in;
  int i,k;
  /* Testing the arguments */
  if((norm == NULL) || (tilt == NULL)) { 
    glb_error("Input pointers may not be NULL");
  return -1;}

  /* Testing the experiment number */
  if(!(((experiment >= 0)&&(experiment < glb_num_of_exps)))) { 
    glb_error("Invalid value for experiment number");
    return -1;}

  for(i=0;i<glb_num_of_exps;i++)
    {
      if(experiment!=GLB_ALL) i=experiment;
      in=(struct glb_experiment *) glb_experiment_list[i];
      /* Testing the rule number */
      if(!(((rule >= 0)&&(rule < in->numofrules)))) { 
	glb_error("Invalid value for rule number");
	return -1;}     
      for(k=0;k<in->numofrules;k++)
	{
	   if(rule!=GLB_ALL) k=rule;
	   /* Here should come the assignment */
	  
	   *norm=in->bgerror[0][k];
	   *tilt=in->bgerror[1][k];
	   
	   if(rule!=GLB_ALL) break;
	}

      if(experiment!=GLB_ALL) break;
    }
  return 0;
}

/* same for BG errors */

int glbSetBGCenters(int experiment, int rule, double norm, double tilt)
{
  struct glb_experiment *in;
  int i,k;
  /* Testing the arguments */
  if((norm <= 0) || (tilt <= 0)) { glb_error("Errors have to be positive");
  return -1;}

  /* Testing the experiment number */
  if(!(((experiment >= 0)&&(experiment < glb_num_of_exps))
       ||(experiment==GLB_ALL))) { 
    glb_error("Invalid value for experiment number");
    return -1;}

  for(i=0;i<glb_num_of_exps;i++)
    {
      if(experiment!=GLB_ALL) i=experiment;
      in=(struct glb_experiment *) glb_experiment_list[i];
      /* Testing the rule number */
      if(!(((rule >= 0)&&(rule < in->numofrules))
	   ||(rule==GLB_ALL))) { 
	glb_error("Invalid value for rule number");
	return -1;}     
      for(k=0;k<in->numofrules;k++)
	{
	   if(rule!=GLB_ALL) k=rule;
	   /* Here should come the assignment */
	   in->bgcenter[0][k]=norm;
	   in->bgcenter[1][k]=tilt;
	   
	   
	   if(rule!=GLB_ALL) break;
	}

      if(experiment!=GLB_ALL) break;
    }
  return 0;
}


int 
glbGetBGCenters(int experiment, int rule, double *norm, double *tilt)
{
  struct glb_experiment *in;
  int i,k;
  /* Testing the arguments */
  if((norm == NULL) || (tilt == NULL)) { 
    glb_error("Input pointers may not be NULL");
  return -1;}

  /* Testing the experiment number */
  if(!(((experiment >= 0)&&(experiment < glb_num_of_exps)))) { 
    glb_error("Invalid value for experiment number");
    return -1;}

  for(i=0;i<glb_num_of_exps;i++)
    {
      if(experiment!=GLB_ALL) i=experiment;
      in=(struct glb_experiment *) glb_experiment_list[i];
      /* Testing the rule number */
      if(!(((rule >= 0)&&(rule < in->numofrules)))) { 
	glb_error("Invalid value for rule number");
	return -1;}     
      for(k=0;k<in->numofrules;k++)
	{
	   if(rule!=GLB_ALL) k=rule;
	   /* Here should come the assignment */
	  
	   *norm=in->bgcenter[0][k];
	   *tilt=in->bgcenter[1][k];
	   
	   if(rule!=GLB_ALL) break;
	}

      if(experiment!=GLB_ALL) break;
    }
  return 0;
}

const char *glbVersionOfExperiment(int experiment)
{
  struct glb_experiment *in;
  if(!(((experiment >= 0)&&(experiment < glb_num_of_exps)))) { 
    glb_error("Invalid value for experiment number");
    return NULL;}
  in=(struct glb_experiment *) glb_experiment_list[experiment];
  return in->version;
}


/* Flux access routines */

double
glbFlux(int experiment, int flux_ident, 
	double energy, double distance, int flavour, int anti)
{
  double out;
  struct glb_experiment *in;
  if(!(((experiment >= 0)&&(experiment < glb_num_of_exps)))) { 
    glb_error("Invalid value for experiment number");
    return -1;}
  in=(struct glb_experiment *) glb_experiment_list[experiment];
  if(!(((flux_ident >= 0)&&(flux_ident < in->num_of_fluxes)))) { 
    glb_error("Invalid value for flux number");
    return -1;}  

  


  out=glb_flux_calc(energy, distance,0,flavour, anti,in->fluxes[flux_ident]); 

  return out;
}

double glbXSection(int experiment, int xsec_ident, 
	double energy, int flavour, int anti)
{
  double out;
  struct glb_experiment *in;
  if(!(((experiment >= 0)&&(experiment < glb_num_of_exps)))) { 
    glb_error("Invalid value for experiment number");
    return -1;}
  in=(struct glb_experiment *) glb_experiment_list[experiment];
  if(!(((xsec_ident >= 0)&&(xsec_ident < in->num_of_xsecs)))) { 
    glb_error("Invalid value for X-section number");
    return -1;}  

  out=glb_xsec_calc(energy,flavour, anti,in->xsecs[xsec_ident]); 

  return out;
}

int
glbSetSourcePower(int experiment, int flux_ident, double power)
{
  int i,k;
  struct glb_experiment *in;

  /* Testing the experiment number */
  if(!(((experiment >= 0)&&(experiment < glb_num_of_exps))
       ||(experiment==GLB_ALL))) { 
    glb_error("Invalid value for experiment number");
    return -1;}
  if(power<=0)  {glb_error("Power must be positive");return -1;}
  
  for(i=0;i<glb_num_of_exps;i++)
    {
      if(experiment!=GLB_ALL) i=experiment;
      in=(struct glb_experiment *) glb_experiment_list[i];
      /* Testing theflux number */
      if(!(((flux_ident >= 0)&&(flux_ident < in->num_of_fluxes))
	   ||(flux_ident==GLB_ALL))) { 
	glb_error("Invalid value for flux number");
	return -1;}     
      for(k=0;k<in->num_of_fluxes;k++)
	{
	   if(flux_ident!=GLB_ALL) k=flux_ident;
	   /* Here should come the assignment */
	   (in->fluxes[k])->target_power=power;
	   
	   if(flux_ident!=GLB_ALL) break;
	}

      if(experiment!=GLB_ALL) break;
    }

 
  return 0;
}

double
glbGetSourcePower(int experiment, int flux_ident)
{
  double out;
  struct glb_experiment *in;
  if(!(((experiment >= 0)&&(experiment < glb_num_of_exps)))) { 
    glb_error("Invalid value for experiment number");
    return -1;}
  in=(struct glb_experiment *) glb_experiment_list[experiment];
  if(!(((flux_ident >= 0)&&(flux_ident < in->num_of_fluxes)))) { 
    glb_error("Invalid value for flux number");
    return -1;}
  
  out=(in->fluxes[flux_ident])->target_power;
  
  return out;
}

int
glbSetRunningTime(int experiment, int flux_ident, double time)
{
  
 int i,k;
  struct glb_experiment *in;

  /* Testing the experiment number */
  if(!(((experiment >= 0)&&(experiment < glb_num_of_exps))
       ||(experiment==GLB_ALL))) { 
    glb_error("Invalid value for experiment number");
    return -1;}
  if(time<=0)  {glb_error("Time must be positive");return -1;}
  
  for(i=0;i<glb_num_of_exps;i++)
    {
      if(experiment!=GLB_ALL) i=experiment;
      in=(struct glb_experiment *) glb_experiment_list[i];
      /* Testing theflux number */
      if(!(((flux_ident >= 0)&&(flux_ident < in->num_of_fluxes))
	   ||(flux_ident==GLB_ALL))) { 
	glb_error("Invalid value for flux number");
	return -1;}     
      for(k=0;k<in->num_of_fluxes;k++)
	{
	   if(flux_ident!=GLB_ALL) k=flux_ident;
	   /* Here should come the assignment */
	   (in->fluxes[k])->time=time;
	   
	   if(flux_ident!=GLB_ALL) break;
	}

      if(experiment!=GLB_ALL) break;
    }

 
  return 0;
}

double
glbGetRunningTime(int experiment, int flux_ident)
{
  double out;
  struct glb_experiment *in;
  if(!(((experiment >= 0)&&(experiment < glb_num_of_exps)))) { 
    glb_error("Invalid value for experiment number");
    return -1;}
  in=(struct glb_experiment *) glb_experiment_list[experiment];
  if(!(((flux_ident >= 0)&&(flux_ident < in->num_of_fluxes)))) { 
    glb_error("Invalid value for flux number");
    return -1;}
  
  out=(in->fluxes[flux_ident])->time;
 
  return out;
}

/* Access to the filter feature */

int
glbSetFilterState(int on_off)
{
  if (on_off==GLB_ON) {glb_switch_filter(GLB_ON);}
  if (on_off==GLB_OFF) {glb_switch_filter(GLB_OFF);}
  if((on_off!=GLB_ON)&&(on_off!=GLB_OFF)) 
    {glb_error("On/Off was not GLB_ON/GLB_OFF");return -1;}
  return 0;
}

int
glbGetFilterState()
{
  int i;
  i=glb_check_filter();
  return i;
}


int
glbSetFilterStateInExperiment(int experiment,int on_off)
{
  int i;
  struct glb_experiment *in;
  /* Testing the experiment number */
  if(!(((experiment >= 0)&&(experiment < glb_num_of_exps))
       ||(experiment==GLB_ALL))) { 
    glb_error("Invalid value for experiment number");
    return -1;}
  if((on_off!=GLB_ON)&&(on_off!=GLB_OFF))  
    {glb_error("on_off must be either GLB_ON or GLB_OFF");return -1;}
  
  for(i=0;i<glb_num_of_exps;i++)
    {
      if(experiment!=GLB_ALL) i=experiment;
      in=(struct glb_experiment *) glb_experiment_list[i];
      
      in->filter_state=on_off;

      if(experiment!=GLB_ALL) break;
    }
  return 0;
}

int
glbGetFilterStateInExperiment(int experiment)
{
  int i,out=-1;
  struct glb_experiment *in;
  /* Testing the experiment number */
  if(!(((experiment >= 0)&&(experiment < glb_num_of_exps)))) { 
    glb_error("Invalid value for experiment number");
    return -1;}
  
  for(i=0;i<glb_num_of_exps;i++)
    {
      if(experiment!=GLB_ALL) i=experiment;
      in=(struct glb_experiment *) glb_experiment_list[i];
      
      out=in->filter_state;

      if(experiment!=GLB_ALL) break;
    }
  return out;
}


int
glbSetFilter(double filter)
{
  if(filter<0) {glb_error("Filter must be positive");return -1;}
  glb_set_filter(filter);
  return 0;
}

double
glbGetFilter()
{
  double out;
  out=glb_get_filter();
  return out;
}

int
glbSetFilterInExperiment(int experiment,double filter)
{
  int i;
  struct glb_experiment *in;
  /* Testing the experiment number */
  if(!(((experiment >= 0)&&(experiment < glb_num_of_exps))
       ||(experiment==GLB_ALL))) { 
    glb_error("Invalid value for experiment number");
    return -1;}
  if(filter<0) {glb_error("Filter must be positive");return -1;}
  
  for(i=0;i<glb_num_of_exps;i++)
    {
      if(experiment!=GLB_ALL) i=experiment;
      in=(struct glb_experiment *) glb_experiment_list[i];
      
      in->filter_value=filter;

      if(experiment!=GLB_ALL) break;
    }
  return 0;
}


double
glbGetFilterInExperiment(int experiment)
{
  int i,out=-1;
  struct glb_experiment *in;
  /* Testing the experiment number */
  if(!(((experiment >= 0)&&(experiment < glb_num_of_exps)))) { 
    glb_error("Invalid value for experiment number");
    return -1;}

  for(i=0;i<glb_num_of_exps;i++)
    {
      if(experiment!=GLB_ALL) i=experiment;
      in=(struct glb_experiment *) glb_experiment_list[i];
      
      out=in->filter_value;

      if(experiment!=GLB_ALL) break;
    }
  return out;
}

/* Accessing parser meta-information */

int 
glbNameToValue(int exp,const char* context, const char *name)
{
  glb_naming *ptr;
  struct glb_experiment *in;
  /* Testing the experiment number */
  if(!(((exp >= 0)&&(exp < glb_num_of_exps)))) { 
    glb_error("Invalid value for experiment number");
    return -1;}
  in=(struct glb_experiment *) glb_experiment_list[exp];
  for (ptr = in->names; ptr != (glb_naming *) NULL;
       ptr = (glb_naming *)ptr->next)
    if (strcmp (ptr->name,name) == 0 && strcmp(ptr->context,context)==0)
	return (ptr->value)-1;
  return -1;
}


const char
*glbValueToName(int exp,const char* context, int value)
{
  glb_naming *ptr;
  struct glb_experiment *in;
  /* Testing the experiment number */
  if(!(((exp >= 0)&&(exp < glb_num_of_exps)))) { 
    glb_error("Invalid value for experiment number");
    return NULL;}
  in=(struct glb_experiment *) glb_experiment_list[exp];
  for (ptr = in->names; ptr != (glb_naming *) NULL;
       ptr = (glb_naming *)ptr->next)
    {
        if (ptr->value == value + 1  && strcmp(ptr->context,context)==0)
	  return ptr->name;
    }
  return NULL;
}


/* Acces to the more complicated probabilities */

double glbProfileProbability(int exp,int initial_flavour, int final_flavour,
			    int panti, double energy)
{
  double res;
  if(exp<0||exp>=glb_num_of_exps) {glb_error("Experiment index out of range");
  return -1;}

  glbSetExperiment(glb_experiment_list[exp]);

  res=glb_profile_probability(initial_flavour,final_flavour,panti,energy);
  return res;
}

/* Number of fluxes in an experiment */

int glbGetNumberOfFluxes(int exp)
{
  int s;
  /* Testing the experiment number */
  if(!(((exp >= 0)&&(exp < glb_num_of_exps)))) { 
    glb_error("Invalid value for experiment number");
    return -1;}

  s=glb_experiment_list[exp]->num_of_fluxes;
  return s;

}

