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
#if HAVE_CONFIG_H   /* config.h should come before any other includes */
#  include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <globes/globes.h>

#include "glb_probability.h"
#include "glb_fluxes.h"
#include "glb_rate_engine.h"

#include "glb_error.h"

#include "glb_minimize.h"
#include "glb_smear.h"
#include "glb_types.h"
#include "glb_multiex.h"
#include "glb_version.h"
#include "glb_sys.h"
#include "glb_wrapper.h"

#define FLOAT double

#define PI 3.1415

/* global variables */
int glb_current_exp = -1;  // Indicates which experiment is currently
                           // being processed by the rate engine
int glb_ignore_invalid_chi2 = 0; /* If != 0, use chiZero instead of throwing an error
                                    if unknown chi^2 function is used. This is necessary
                                    when an AEDL file using user-defined systematics is
                                    run through the standalone globes binary */

static struct glb_experiment MInitMemory0(struct glb_experiment in);

/* Functions for dealing with glb_naming */

static void
glb_free_names(glb_naming *stale)
{
    glb_naming *ptr;
    glb_naming *dummy;
    ptr=stale;
    while(ptr != (glb_naming *) NULL)
      {
        glb_free(ptr->name);
        glb_free(ptr->context);
        dummy=ptr->next;
        glb_free(ptr);
        ptr=dummy;
      }
    stale=NULL;
}


/***************************************************************************
 * Function glb_alloc_nuisance                                             *
 ***************************************************************************
 * Allocate a nuisance parameter data structure and initialize it to dummy *
 * values                                                                  *
 ***************************************************************************/
glb_nuisance *glb_alloc_nuisance()
{
  glb_nuisance *n = (glb_nuisance *) glb_malloc(sizeof(glb_nuisance));
  if (n)
  {
    n->name        = NULL;
    n->systype     = 0;
    n->error       = GLB_NAN;
    n->a           = GLB_NAN;
    n->n_energies  = 0;
    n->energy_list = NULL;
    n->error_list  = NULL;
    n->a_list      = NULL;
    n->ref_count   = 1;
  }

  return n;
}


/***************************************************************************
 * Function glb_copy_nuisance                                              *
 ***************************************************************************
 * Duplicates a nuisance parameter data structure. dest must already       *
 * exist. Note that the reference counter is *not* copied.                 *
 ***************************************************************************/
int glb_copy_nuisance(glb_nuisance *dest, glb_nuisance *src)
{
  if (src && dest)
  {
    dest->name        = strdup(src->name);
    dest->error       = src->error;
    dest->a           = src->a; 
    dest->systype     = src->systype;
    dest->n_energies  = src->n_energies;
    dest->energy_list = glb_duplicate_array(src->energy_list, src->n_energies*sizeof(double));
    dest->error_list  = glb_duplicate_array(src->error_list, src->n_energies*sizeof(double));
    dest->a_list      = glb_duplicate_array(src->a_list, src->n_energies*sizeof(double));

    return GLB_SUCCESS;
  }
  else
  {
    glb_error("glb_copy_nuisance: NULL is not a valid glb_nuisance structure");
    return GLBERR_UNINITIALIZED;
  }
}


/***************************************************************************
 * Function glb_free_nuisance                                              *
 ***************************************************************************
 * Release memory associated with a nuisance parameter data structure      *
 ***************************************************************************/
int glb_free_nuisance(glb_nuisance *n)
{
  if (n)
  {
    if (n->ref_count <= 1)
    {
      glb_free(n->name);          n->name        = NULL;
      glb_free(n->energy_list);   n->energy_list = NULL;
      glb_free(n->error_list);    n->error_list  = NULL;
      glb_free(n->a_list);        n->a_list      = NULL;
      n->error      = GLB_NAN;
      n->a          = GLB_NAN;
      n->n_energies = 0;
      n->ref_count  = 0;
      n->systype = 0;
      glb_free(n);
    }
    else
      n->ref_count--;
  }

  return GLB_SUCCESS;
}


/***************************************************************************
 * Function glbAllocExp                                                    *
 ***************************************************************************
 * Allocate a new glb_experiment data structure and initialize with dummy  *
 * values.                                                                 *
 ***************************************************************************/
glb_exp glbAllocExp()
{
  struct glb_experiment *in;
  in = (struct glb_experiment *) glb_malloc(sizeof(struct glb_experiment));
  glbInitExp((glb_exp) in);
  return (glb_exp) in;
}


/***************************************************************************
 * Function glbInitExp                                                     *
 ***************************************************************************
 * Initialize a glb_experiment data structure with safe dummy values.      *
 * This must be done before the experiment definition is read from an      *
 * AEDL file.                                                              *
 ***************************************************************************/
void glbInitExp(glb_exp ins)
{
  int i, j;
  struct glb_experiment *in = (struct glb_experiment *) ins;
  in->version=NULL;
  in->citation=NULL;
  in->parent=NULL;
  in->n_children = 0;
  for (i=0; i < GLB_MAX_EXP; i++)
    in->children[i] = NULL;
  in->filename=NULL;
  in->names=NULL;
  in->num_of_fluxes=-1;
  in->density_profile_type=-1;
  /* FIXME - a potential memory leak */
  for(i=0;i<GLB_MAX_FLUXES;i++) in->fluxes[i]=NULL;
  in->num_of_xsecs=-1;
  /* FIXME - a potential memory leak */
  for(i=0;i<GLB_MAX_XSECS;i++) in->xsecs[i]=NULL;

  in->binsize=NULL;
  in->simbinsize=NULL;
  in->baseline=-1;
  in->emin=-1;
  in->emax=-1;
  in->numofbins=-1;
  in->targetmass=-1;

  in->numofrules = -1;
  for (i=0; i < GLB_MAX_RULES; i++)
  {
    in->lengthofrules[i]         = -1;
    in->rulescoeff[i]            = NULL;
    in->rulechannellist[i]       = NULL;
    in->lengthofbgrules[i]       = -1;
    in->bgrulescoeff[i]          = NULL;
    in->bgrulechannellist[i]     = NULL;
                                 
    in->sys_on[i]                = NULL;
    in->sys_off[i]               = NULL;
    in->sys_on_strings[i]        = NULL;
    in->sys_off_strings[i]       = NULL;
    in->sys_on_off[i]            = GLB_ON;
    in->sys_on_errors[i]         = NULL;
    in->sys_on_startvals[i]      = NULL;
    in->sys_off_errors[i]        = NULL;
    in->sys_off_startvals[i]     = NULL;
    for (j=0; j < GLB_MAX_CHANNELS; j++)
    {
      in->sys_on_n_nuis_sig[i][j]          = 0;
      in->sys_on_n_nuis_bg[i][j]           = 0;
      in->sys_off_n_nuis_sig[i][j]         = 0;
      in->sys_off_n_nuis_bg[i][j]          = 0;
      in->sys_on_multiex_errors_sig[i][j]  = NULL;
      in->sys_on_multiex_errors_bg[i][j]   = NULL;
      in->sys_off_multiex_errors_sig[i][j] = NULL;
      in->sys_off_multiex_errors_bg[i][j]  = NULL;
    }
                                 
    in->bg_centers[0][i]         =  1.0;
    in->bg_centers[1][i]         =  0.0;
    in->signal_errors[0][i]      = -1.0;
    in->signal_errors[1][i]      = -1.0;
    in->signal_startvals[0][i]   =  0.0;
    in->signal_startvals[1][i]   =  0.0;
    in->bg_errors[0][i]          = -1.0;
    in->bg_errors[1][i]          = -1.0;
    in->bg_startvals[0][i]       =  0.0;
    in->bg_startvals[1][i]       =  0.0;

    in->energy_window[i][0]      = -1.0;
    in->energy_window[i][1]      = -1.0;
    in->energy_window_bins[i][0] = -1;
    in->energy_window_bins[i][1] = -1;

    in->SignalRates[i]           = NULL;
    in->BackgroundRates[i]       = NULL;
    in->rates0[i]                = NULL;
    in->rates1[i]                = NULL;
    in->rates1BG[i]              = NULL;
    /* PH 01/10/19 */
    in->data_on_off[i]=GLB_OFF;
    in->data[i]=NULL;
  }

  in->numofchannels=-1;
  for (i=0; i < GLB_MAX_CHANNELS; i++)
  {
    in->chrb_0[i]                        = NULL;
    in->chrb_1[i]                        = NULL;
    in->chra_0[i]                        = NULL;
    in->chra_1[i]                        = NULL;
    in->chr_template[i]                  = NULL; 
    in->user_pre_smearing_channel[i]     = NULL;
    in->user_post_smearing_channel[i]    = NULL;
    in->user_pre_smearing_background[i]  = NULL;
    in->user_post_smearing_background[i] = NULL;
  }
  for(i=0; i < 6; i++)
    in->listofchannels[i] = NULL;

  in->num_of_sm=-1;
  in->n_nuisance = -1;
  for(i=0; i < GLB_MAX_NUISANCE; i++)
    in->nuisance_params[i] = NULL;
  for(i=0;i<GLB_MAX_SMEAR;i++) {in->smear_data[i]=NULL;}
  for(i=0;i<GLB_MAX_SMEAR;i++) in->smear[i]=NULL;
  for(i=0;i<GLB_MAX_SMEAR;i++) in->lowrange[i]=NULL;
  for(i=0;i<GLB_MAX_SMEAR;i++) in->uprange[i]=NULL;
  in->simtresh=-1;
  in->simbeam=-1;
  in->simbins=-1;

  in->filter_state=-1;
  in->filter_value=-1;
  in->psteps=-1;
  in->lengthtab=NULL;
  in->densitytab=NULL;
  in->densitybuffer=NULL;
  in->energy_tab=NULL;

  /*
    in->probability_matrix         = NULL;
    in->set_oscillation_parameters = NULL;
    in->get_oscillation_parameters = NULL;
    in->probability_user_data      = NULL;
  */
  
  /* named oscillation engine setup */

  (in->osc_engine).num_of_params=-1;
  (in->osc_engine).matrix_function=NULL;
  (in->osc_engine).get_function=NULL;
  (in->osc_engine).set_function=NULL;
  (in->osc_engine).name=NULL; /* memory leak ? */
  (in->osc_engine).user_data=NULL;

}


static void my_free(void *ptr)
{
  if(ptr==NULL) return;
  glb_free(ptr);
  ptr=NULL;
}


/***************************************************************************
 * Function glbInitExpFromParent                                           *
 ***************************************************************************
 * Initialize a glb_experiment data structure by copying certain properties*
 * from another detector. This is used when the AEDL parser encounters a   *
 * #DETECTOR# directive. Note that rules inherited from a parent experiment*
 * will only be used if the child does not have its own rule defined in    *
 * the AEDL file.                                                          *
 ***************************************************************************/
void glbInitExpFromParent(struct glb_experiment *exp, struct glb_experiment *p)
{
  int i, j, k;

  if (!exp)
    glb_error("glbInitExpFromParent: NULL is not a valid experiment");
  if (!p)
    glb_error("glbInitExpFromParent: NULL is not a valid parent experiment");

  /* Start with a clean data structure */
  glbInitExp(exp);

  exp->parent     = p;
  glbExpAddChild(p, exp);
  exp->names      = NULL; /* Namespace will be copied by parser once all
                             AEDL declarations for this experiment are parsed */

  if (p->version)   exp->version  = strdup(p->version);
  if (p->citation)  exp->citation = strdup(p->citation);
  if (p->filename)  exp->filename = strdup(p->filename);
  /* I believe this is necessary so that named oscillation engines are
  inherited from the parent, and that clearly should be default
  behavior - PH 01/10/19 */
  if (p->osc_engine.name) exp->osc_engine.name = strdup(p->osc_engine.name);

  /* Binning */
  exp->emin       = p->emin;
  exp->emax       = p->emax;
  exp->numofbins  = p->numofbins;
  exp->binsize    = glb_duplicate_array(p->binsize, exp->numofbins * sizeof(double));
  exp->simtresh   = p->simtresh;
  exp->simbeam    = p->simbeam;
  exp->simbins    = p->simbins;
  exp->simbinsize = glb_duplicate_array(p->simbinsize, exp->simbins * sizeof(double));

  /* Detector parameteres */
  exp->targetmass           = p->targetmass;
  exp->baseline             = p->baseline;
  exp->density_profile_type = p->density_profile_type;
  exp->psteps               = p->psteps;
  /* FIXME What if lengthtab and densitytab have unequal lengths? */
  exp->lengthtab            = glb_duplicate_array(p->lengthtab, exp->psteps * sizeof(double));
  exp->densitytab           = glb_duplicate_array(p->lengthtab, exp->psteps * sizeof(double));
  exp->filter_state         = p->filter_state;
  exp->filter_value         = p->filter_value;

  /* Copy fluxes from parent */
  exp->num_of_fluxes = p->num_of_fluxes;
  for (i=0; i < exp->num_of_fluxes; i++)
    exp->fluxes[i] = p->fluxes[i];

  /* Copy cross sections from parent */
  exp->num_of_xsecs = p->num_of_xsecs;
  for (i=0; i < exp->num_of_xsecs; i++)
    exp->xsecs[i] = p->xsecs[i];

  /* Copy smearing definitions from parent */
  exp->num_of_sm = p->num_of_sm;
  for (i=0; i < exp->num_of_sm; i++)
  {
    if (p->smear_data[i])   /* Parametrized smearing */
    {
      exp->smear_data[i] = glb_smear_alloc();
      glb_copy_smear(exp->smear_data[i], p->smear_data[i]);
    }
    if (p->smear[i])        /* Explicit smearing matrix given in AEDL file */
    {
      exp->lowrange[i] = glb_malloc((p->numofbins+1) * sizeof(int));
      exp->uprange[i]  = glb_malloc((p->numofbins+1) * sizeof(int));
      exp->smear[i]    = glb_malloc((p->numofbins+1) * sizeof(double *));
      memset(exp->lowrange[i], 0, (p->numofbins+1) * sizeof(int));
      memset(exp->uprange[i],  0, (p->numofbins+1) * sizeof(int));
      memset(exp->smear[i],    0, (p->numofbins+1) * sizeof(double *));
      if (exp->lowrange[i]  &&  exp->uprange[i]  &&  exp->smear[i])
      {
        for (j=0; j < p->numofbins+1; j++)
        {
          exp->lowrange[i][j] = p->lowrange[i][j];
          exp->uprange[i][j]  = p->uprange[i][j];
          if (exp->lowrange[i][j] < 0 || exp->uprange[i][j] < 0)
            break; /* Just in case the smearing matrix has incorrect dimensions */
          if (p->smear[i][j])
          {
            exp->smear[i][j] = glb_malloc((exp->uprange[i][j] - exp->lowrange[i][j] + 2)
                                 * sizeof(double));
            if (exp->smear[i][j] != NULL)
              for (k=0; k < exp->uprange[i][j] - exp->lowrange[i][j] + 2; k++)
              {
                exp->smear[i][j][k] = p->smear[i][j][k];
                if (exp->smear[i][j][k] < 0)
                  break; /* Just in case the smearing matrix has incorrect dimensions */
              }
          }
          else
            exp->smear[i][j] = NULL;
        }
      }
      else
      {
        glb_free(exp->lowrange[i]);  exp->lowrange[i] = NULL;
        glb_free(exp->uprange[i]);   exp->uprange[i]  = NULL;
        glb_free(exp->smear[i]);     exp->smear[i]    = NULL;
      }
    }
  }

  /* Copy channels from parent */
  exp->numofchannels = p->numofchannels;
  for (i=0; i < 6; i++)
    exp->listofchannels[i] = glb_duplicate_array(p->listofchannels[i],
                                                 sizeof(int) * exp->numofchannels);
  for (i=0; i < exp->numofchannels; i++)
  {
    int n1 = sizeof(double) * (p->simbins + 1);
    int n2 = sizeof(double) * (p->numofbins + 1);
    exp->user_pre_smearing_channel[i]
      = glb_duplicate_array(p->user_pre_smearing_channel[i], n1);
    exp->user_post_smearing_channel[i]
      = glb_duplicate_array(p->user_post_smearing_channel[i], n2);
    exp->user_pre_smearing_background[i]
      = glb_duplicate_array(p->user_pre_smearing_background[i], n1);
    exp->user_post_smearing_background[i]
      = glb_duplicate_array(p->user_post_smearing_background[i], n2);
  }

  /* Copy rules from parent */
  /* JK - 2012-05-11 Removed because it was very confusing */
  exp->numofrules = 0;

//  exp->numofrules = p->numofrules;
//  for (i=0; i < exp->numofrules; i++)
//  {
//    exp->lengthofrules[i]         = p->lengthofrules[i];
//    exp->rulescoeff[i]            = glb_duplicate_array(p->rulescoeff[i],
//                                      exp->lengthofrules[i] * sizeof(double));
//    exp->rulechannellist[i]       = glb_duplicate_array(p->rulechannellist[i],
//                                      exp->lengthofrules[i] * sizeof(int));
//
//    exp->lengthofbgrules[i]       = p->lengthofbgrules[i];
//    exp->bgrulescoeff[i]          = glb_duplicate_array(p->bgrulescoeff[i],
//                                      exp->lengthofbgrules[i] * sizeof(double));
//    exp->bgrulechannellist[i]     = glb_duplicate_array(p->bgrulechannellist[i],
//                                      exp->lengthofbgrules[i] * sizeof(int));
//
//    exp->energy_window[i][0]      = p->energy_window[i][0];
//    exp->energy_window[i][1]      = p->energy_window[i][1];
//    exp->energy_window_bins[i][0] = p->energy_window_bins[i][0];
//    exp->energy_window_bins[i][1] = p->energy_window_bins[i][1];
//
//    exp->sys_on_off[i]            = p->sys_on_off[i];
//    exp->sys_on[i]                = p->sys_on[i];
//    exp->sys_off[i]               = p->sys_off[i];
//    exp->sys_on_strings[i]        = strdup(p->sys_on_strings[i]);
//    exp->sys_off_strings[i]       = strdup(p->sys_off_strings[i]);
//    if (p->sys_on_errors[i])
//    {
//      for (j=0; p->sys_on_errors[i][j] > 0.0; j++)
//        ;
//      exp->sys_on_errors[i] = glb_duplicate_array(p->sys_on_errors[i], (j+1)*sizeof(double));
//    }
//    if (p->sys_on_startvals[i])
//    {
//      for (j=0; p->sys_on_startvals[i][j] > 0.0; j++)
//        ;
//      exp->sys_on_startvals[i] = glb_duplicate_array(p->sys_on_startvals[i], (j+1)*sizeof(double));
//    }
//    if (p->sys_off_errors[i])
//    {
//      for (j=0; p->sys_off_errors[i][j] > 0.0; j++)
//        ;
//      exp->sys_off_errors[i] = glb_duplicate_array(p->sys_off_errors[i], (j+1)*sizeof(double));
//    }
//    if (p->sys_off_startvals[i])
//    {
//      for (j=0; p->sys_off_startvals[i][j] > 0.0; j++)
//        ;
//      exp->sys_off_startvals[i] = glb_duplicate_array(p->sys_off_startvals[i],(j+1)*sizeof(double));
//    }
//    for (j=0; j < 2; j++)
//    {
//      exp->signal_errors[j][i]    = p->signal_errors[j][i];
//      exp->signal_startvals[j][i] = p->signal_startvals[j][i];
//      exp->bg_errors[j][i]        = p->bg_errors[j][i];
//      exp->bg_startvals[j][i]     = p->bg_startvals[j][i];
//      exp->bg_centers[j][i]       = p->bg_centers[j][i];
//    }
//    for (j=0; j < GLB_MAX_CHANNELS; j++)
//    {
//      exp->sys_on_n_nuis_sig[i][j]  = p->sys_on_n_nuis_sig[i][j];
//      exp->sys_on_n_nuis_bg[i][j]   = p->sys_on_n_nuis_bg[i][j];
//      exp->sys_off_n_nuis_sig[i][j] = p->sys_off_n_nuis_sig[i][j];
//      exp->sys_off_n_nuis_bg[i][j]  = p->sys_off_n_nuis_bg[i][j];
//
//      /* Don't rely on sys_on/off_n_nuis_sig/bg for the copying since they
//       * may not be initialized yet (that happens only in glbDefaultExp, not
//       * in the parser */
//      if (p->sys_on_multiex_errors_sig[i][j])
//      {
//        for (k=0; p->sys_on_multiex_errors_sig[i][j][k] >= 0; k++)
//          ;
//        exp->sys_on_multiex_errors_sig[i][j]
//          = glb_duplicate_array(p->sys_on_multiex_errors_sig[i][j], (k+1)*sizeof(int));
//      }
//      if (p->sys_on_multiex_errors_bg[i][j])
//      {
//        for (k=0; p->sys_on_multiex_errors_bg[i][j][k] >= 0; k++)
//          ;
//        exp->sys_on_multiex_errors_bg[i][j]
//          = glb_duplicate_array(p->sys_on_multiex_errors_bg[i][j], (k+1)*sizeof(int));
//      }
//      if (p->sys_off_multiex_errors_sig[i][j])
//      {
//        for (k=0; p->sys_off_multiex_errors_sig[i][j][k] >= 0; k++)
//          ;
//        exp->sys_off_multiex_errors_sig[i][j]
//          = glb_duplicate_array(p->sys_off_multiex_errors_sig[i][j], (k+1)*sizeof(int));
//      }
//      if (p->sys_off_multiex_errors_bg[i][j])
//      {
//        for (k=0; p->sys_off_multiex_errors_bg[i][j][k] >= 0; k++)
//          ;
//        exp->sys_off_multiex_errors_bg[i][j]
//          = glb_duplicate_array(p->sys_off_multiex_errors_bg[i][j], (k+1)*sizeof(int));
//      }
//    }
//  }

  /* Copy nuisance parameter info */
  exp->n_nuisance = p->n_nuisance;
  for (i=0; i < exp->n_nuisance; i++)
  {
    exp->nuisance_params[i] = p->nuisance_params[i];
    exp->nuisance_params[i]->ref_count++;
  }
}


/***************************************************************************
 * Function glbCorrelateSys                                                *
 ***************************************************************************
 * Correlates the systematic uncertainties between two experiments. This   *
 * is done by looking for nuisance parameters with identical names, and    *
 * identifying them with each other                                        *
 ***************************************************************************/
int glbCorrelateSys(struct glb_experiment *e1, struct glb_experiment *e2)
{
  int i, j;
  for (i=0; i < e1->n_nuisance; i++)
    for (j=0; j < e2->n_nuisance; j++)
    {
      if (e1->nuisance_params[i]  &&  e2->nuisance_params[j])
        if (strcmp(e1->nuisance_params[i]->name, e2->nuisance_params[j]->name) == 0)
        {
//          printf("Correlating %s\n", e1->nuisance_params[i]->name);
          glb_free_nuisance(e2->nuisance_params[j]);
          e2->nuisance_params[j] = e1->nuisance_params[i];
          e1->nuisance_params[i]->ref_count++;
        }
    }

  return GLB_SUCCESS;
}


/***************************************************************************
 * Function glbExpAddChild                                                 *
 ***************************************************************************
 * Adds a child experiment to an experiment.                               *
 ***************************************************************************/
void glbExpAddChild(struct glb_experiment *parent, struct glb_experiment *child)
{
  if (parent && child)
  {
    child->parent = parent;
    parent->children[parent->n_children++] = child;
  }
  else
    glb_error("glbExpAddChild: Invalid pointer");
}


/***************************************************************************
 * Function glbExpRemoveChild                                              *
 ***************************************************************************
 * Removes a child experiment from an experiment. The child experiment     *
 * is not automatically destroyed, but this should be done manually.       *
 * (The child and the parent share common memory, which will be relased    *
 * when the parent is destroyed.)                                          *
 ***************************************************************************/
void glbExpRemoveChild(struct glb_experiment *parent, struct glb_experiment *child)
{
  if (parent && child)
  {
    int i;
    for (i=0; i < parent->n_children && parent->children[i] != child; i++)
      ;
    if (i == parent->n_children)
      glb_error("glbExpRemoveChild: Child not found");
    else
    {
      parent->n_children--;
      child->parent = NULL;
      while (i < parent->n_children)
        parent->children[i] = parent->children[i+1];
      parent->children[i] = NULL;
    }
  }
  else
    glb_error("glbExpRemoveChild: Invalid pointer");
}


/***************************************************************************
 * Function glbResetExp                                                    *
 ***************************************************************************
 * Release dynamically allocated memory for a given experiment and its     *
 * children and re-initialize the glb_experiment structure using           *
 * glbInitExp.                                                             *
 ***************************************************************************/
void glbResetExp(struct glb_experiment *in)
{
  int i, j;
  struct glb_experiment *p = in->parent;

  if (in == NULL)
    return;

  /* Destroy child experiments, starting from the back of the list (otherwise,
   * pointers to child experiments would change positions in the list while
   * being consecutively removed. We do not destroy the glb_experiment data
   * structures of the children themselves since there may still be pointers
   * to them in glb_experiment_list */
  for (i=in->n_children-1; i >= 0; i--)
    glbResetExp(in->children[i]);
  if (in->n_children != 0)
    glb_error("Something went wrong in the experiment management.\n");
  in->n_children = 0;

  glb_free(in->version);
  glb_free(in->citation);
  if (p)
    glbExpRemoveChild(in->parent, in);
  in->n_children = 0;
  glb_free(in->filename);

  /* Fluxes, cross sections and the namespace are shared between parents and
   * children -> remove only those belonging to this experiment. Nuisance parameters
   * have their own reference counter, so no special treatment is required for them. */
  if (!p)
  {
    glb_free_names(in->names);
    for(i=0; i < GLB_MAX_FLUXES; i++)
    {
      glb_free_flux(in->fluxes[i]);
      in->fluxes[i] = NULL;
    }
    for(i=0; i < GLB_MAX_XSECS; i++)
    {
      glb_free_xsec(in->xsecs[i]);
      in->xsecs[i] = NULL;
    }
  }
  else
  {
    for(i=MAX(0, p->num_of_fluxes); i < GLB_MAX_FLUXES; i++)
    {
      glb_free_flux(in->fluxes[i]);
      in->fluxes[i] = NULL;
    }
    for(i=MAX(0, p->num_of_xsecs); i < GLB_MAX_XSECS; i++)
    {
      glb_free_xsec(in->xsecs[i]);
      in->xsecs[i] = NULL;
    }
  }

  for(i=0; i < GLB_MAX_NUISANCE; i++)
  {
    glb_free_nuisance(in->nuisance_params[i]);
    in->nuisance_params[i] = NULL;
  }
  in->n_nuisance = -1;

  for(i=0;i<6;i++) my_free(in->listofchannels[i]);

  my_free(in->binsize);
  in->binsize=NULL;
  my_free(in->simbinsize);
  in->simbinsize=NULL;

  for(i=0;i<in->numofrules;i++)
    {
      my_free(in->rulescoeff[i]);
      my_free(in->rulechannellist[i]);
      my_free(in->bgrulescoeff[i]);
      my_free(in->bgrulechannellist[i]);
      my_free(in->data[i]);
    }

  for(i=0;i<in->numofrules;i++)
  {
    my_free(in->sys_on_strings[i]);
    my_free(in->sys_off_strings[i]);
  }


  for(i=0;i<GLB_MAX_SMEAR;i++) glb_smear_free(in->smear_data[i]);
  for(i=0;i<in->num_of_sm;i++)
    {
      my_free(in->lowrange[i]);
      my_free(in->uprange[i]);
      if (in->smear != NULL  &&  in->smear[i] != NULL)
      {
        for(j=0; in->smear[i][j] != NULL; j++) /* NULL signals end of list - length can be */
          my_free(in->smear[i][j]);            /* different from numofbins if smearing matrix */
                                               /* provided by user had too few lines */
        my_free(in->smear[i]);
      }
    }
  for(i=0;i<in->numofchannels;i++)
    {
      my_free(in->chr_template[i]);
      my_free(in->chrb_0[i]);
      my_free(in->chrb_1[i]);
      my_free(in->chra_0[i]);
      my_free(in->chra_1[i]);
      my_free(in->user_pre_smearing_channel[i]);
      my_free(in->user_post_smearing_channel[i]);
      my_free(in->user_pre_smearing_background[i]);
      my_free(in->user_post_smearing_background[i]);
    }
  my_free(in->lengthtab);
  my_free(in->densitytab);
  my_free(in->densitybuffer);

  for(i=0;i<in->numofrules;i++)
    {
      my_free(in->SignalRates[i]);
      my_free(in->BackgroundRates[i]);
      my_free(in->rates0[i]);
      my_free(in->rates1[i]);
      my_free(in->rates1BG[i]);
      my_free(in->sys_on_errors[i]);
      for(j=0;j<in->numofchannels;j++)
      {
        in->sys_on_n_nuis_sig[i][j]          = 0;
        in->sys_on_n_nuis_bg[i][j]           = 0;
        in->sys_off_n_nuis_sig[i][j]         = 0;
        in->sys_off_n_nuis_bg[i][j]          = 0;
        my_free(in->sys_on_multiex_errors_sig[i][j]);
        my_free(in->sys_on_multiex_errors_bg[i][j]);
        my_free(in->sys_off_multiex_errors_sig[i][j]);
        my_free(in->sys_off_multiex_errors_bg[i][j]);
      }
      my_free(in->sys_on_startvals[i]);
      my_free(in->sys_off_errors[i]);
      my_free(in->sys_off_startvals[i]);
    }

  my_free(in->energy_tab);
  my_free((in->osc_engine).name);
  
  /* Re-initialize */
  glbInitExp(in);
}


/***************************************************************************
 * Function glbFreeExp                                                     *
 ***************************************************************************
 * Release dynamically allocated memory for a given experiment and its     *
 * children, and destroy the glb_experiment data structure.                *
 ***************************************************************************/
void glbFreeExp(struct glb_experiment *in)
{
  glbResetExp(in);
  my_free(in);
}



static int setup_density_profile(glb_exp ins)
{
  int s=0;
  double *lb,*db;
  size_t l;
  double ttl=0;
  int i;
  struct glb_experiment *in;
  in=(struct glb_experiment *) ins;
  if(ins->density_profile_type==-1)
    {glb_exp_error(in, "No profile type specified");s=-1;}
  if(ins->density_profile_type==1)
    {
      glb_free(ins->lengthtab);
      glb_free(ins->densitytab);
      ins->psteps=1;
      if(ins->baseline > 0)
        {
          s+=glbAverageDensityProfile(ins->baseline,&lb,&db);
          ins->densitytab=db;
          ins->lengthtab=lb;
        }
      else
        {
          glb_exp_error(in, "Baseline must be a positive number");
          s=-1;
        }
    }

  if(ins->density_profile_type==2)
    {
      glb_free(ins->lengthtab);
      glb_free(ins->densitytab);

      if(ins->baseline > 0)
        {
          if(ins->psteps > 0)
            {
              l=(size_t) ins->psteps;
              s+=glbStaceyProfile(ins->baseline,l,&lb,&db);
              ins->densitytab=db;
              ins->lengthtab=lb;
            }
          else
            {
              glb_exp_error(in, "Densitysteps must be a positive number");
              s=-1;
            }

        }
      else
        {
          glb_exp_error(in, "Baseline must be a positive number");
          s=-1;
        }
    }

  if(ins->density_profile_type==3)
    {
      if(ins->lengthtab!=NULL && ins->psteps > 0)
      {
        for(i=0;i<ins->psteps;i++) ttl+=ins->lengthtab[i];
      }
      if(ttl!=ins->baseline&&ins->baseline!=-1)
        glb_warning("Override baseline with sum of lengthtab");
      ins->baseline=ttl;
    }

  return s;
}


/***************************************************************************
 * Function glbDefaultExp                                                  *
 ***************************************************************************
 * This function is called after the basic properties of the experiment    *
 * have been read from an AEDL file. It prepares the glb_experiment        *
 * structure for productive use by allocating memory (as specified in      *
 * the AEDL file), performing numerous consistency checks, and assigning   *
 * default values where appropriate.                                       *
 ***************************************************************************/
int glbDefaultExp(glb_exp ins)
{
  int i,j,k,status,def,ct;
  int sys_dim;

  struct glb_experiment *in;
  in=(struct glb_experiment *) ins;
  double *tmp_errorlist;
  double tmp;

  status=0;
  def=0;

  status+=setup_density_profile(ins);
  if(in->version==NULL)
    {
      glb_fatal("Missing version in AEDL file");
      def=-1;
      in->version=(char *) strdup(glb_release_version);
    }

  if(glbTestReleaseVersion(in->version)<0) {
    glb_error("AEDL file has a more recent version number than the"
                " installed globes package");}

  if(in->num_of_xsecs<1)  {glb_exp_error(in, "No X-section selected!");status=-1;}
  if(in->num_of_xsecs>31)  {glb_exp_error(in, "To many X-sections!");status=-1;}
  if(in->num_of_fluxes<1)  {glb_exp_error(in, "No flux selected!");status=-1;}
  if(in->num_of_fluxes>31)  {glb_exp_error(in, "To many fluxes!");status=-1;}

  
  /* Initialize flux tables */
  if(in->num_of_fluxes>0&&in->num_of_fluxes<GLB_MAX_FLUXES)
  {
    for(i=0;i<in->num_of_fluxes;i++)
    {
      if(in->fluxes[i] == NULL)
      {
        glb_exp_error(in, "Missing flux specification");
        status = -1;
      }
      else
      {
        if (glb_init_flux(in->fluxes[i]) != GLB_SUCCESS)
          status=-1;
      }
    }
  }

  /* Load cross sections */
  if(in->num_of_xsecs > 0 && in->num_of_xsecs < GLB_MAX_XSECS)
  {
    for(i=0;i < in->num_of_xsecs; i++)
    {
      if(in->xsecs[i] == NULL)
      {
        glb_exp_error(in, "Missing cross section specification");
        status = -1;
      }
      else
      {
        if (in->xsecs[i]->builtin >= 0)
        {
          glb_exp_error(in, "No builtin cross sections available. "
                            "Please specify cross section file");
          status = -1;
        }
        else
          glb_init_xsec(in->xsecs[i]);
      }
    }
  }

  if(in->baseline==-1){glb_exp_error(in, "No baseline specified!");status=-1;}
  if(in->emin==-1){glb_exp_error(in, "No emin specified!");status=-1;}
  if(in->emax==-1){glb_exp_error(in, "No emax specified!");status=-1;}
  if(in->emin>in->emax){glb_exp_error(in, "emin must be less than emax!");status=-1;}
  if(in->numofbins<=0){glb_exp_error(in, "numofbins is not set!");status=-1;}
  if(in->numofbins<=0) { glb_exp_error(in, "Too few bins defined!");status=-1;}
  if(in->targetmass==-1){glb_exp_error(in, "No targetmass specified!");status=-1;}

  if(in->numofchannels==-1){
    glb_exp_error(in, "numofchannels not specified!");status=-1;}

  for(i=0;i<6;i++){
    if(in->listofchannels[i]==NULL)
      {  glb_exp_error(in, "listofchannels not specified!"); status=-1;}
  }
  if(in->numofrules==-1){glb_exp_error(in, "numofrules not specified!");status=-1;}
  for(i=0;i<in->numofrules;i++)
    {
      if(in->lengthofrules[i]==-1)
        { glb_rule_error(in, i, "Invalid rule specification!"); status=-1; }
      if(in->rulescoeff[i]==NULL)
        { glb_rule_error(in, i, "No rulescoeff specified!"); status=-1; }
      if(in->rulechannellist[i]==NULL)
        { glb_rule_error(in, i, "No rulechannellist specified!"); status=-1; }
      if(in->lengthofbgrules[i]==-1)
        { glb_rule_error(in, i, "No lengthofbgrules specified!"); status=-1; }
      if(in->bgrulescoeff[i]==NULL)
        {  glb_rule_error(in, i, "No bgrulescoeff specified!"); status=-1; }
      if(in->bgrulechannellist[i]==NULL)
        { glb_rule_error(in, i, "No bgruloechannellist specified!"); status=-1; }
    }


  /* Check definitions of nuisance parameters for global/multi-experiment systematics */
  /* -------------------------------------------------------------------------------- */
  if (in->n_nuisance < 0)  in->n_nuisance = 0; /* The dummy value -1 is needed only in the parser */
  for (i=0; i < in->n_nuisance; i++)
  {
    if (in->nuisance_params[i] == NULL)
      { glb_exp_error(in, "Nuisance parameter #%d undefined", i);  status=-1; }
    else if (in->nuisance_params[i]->name == NULL)
      { glb_exp_error(in, "Nuisance parameter #%d has no name", i);  status=-1; }
    else if (in->nuisance_params[i]->n_energies <= 0)   /* Energy-independent nuisance parameter */
    {
      if (GLB_ISNAN(in->nuisance_params[i]->error)) {
        glb_exp_error(in, "No error specified for nuisance parameter %s",
                      in->nuisance_params[i]->name);
        status=-1;
      }
    }
    else if (in->nuisance_params[i]->n_energies == 1)  /* Only one energy given */
    {
      glb_exp_error(in, "@energy_list need >= 2 entries in nuisance parameter %s.\n"
                        "Use @error for energy-independent nuisance parameters",
                    in->nuisance_params[i]->name);
      status=-1;
    }
    else                                               /* Energy-dependent nuisance parameter */
    {
      if (!GLB_ISNAN(in->nuisance_params[i]->error)) {
        glb_exp_error(in, "Nuisance parameter %s cannot contain energy-dependent and "
                          "energy-independent parts", in->nuisance_params[i]->name);
        status = -1;
      }
      else if (in->nuisance_params[i]->energy_list == NULL ||
               in->nuisance_params[i]->error_list == NULL) {
        glb_exp_error(in, "Definition of nuisance parameter %s is incomplete",
                          in->nuisance_params[i]->name);
        status = -1;
      }
      else {
        in->nuisance_params[i]->a_list = glb_malloc(sizeof(double)
                                           * in->nuisance_params[i]->n_energies);
      }
    }
  }

  int old_verbosity = glbGetVerbosityLevel();
  int old_status    = status;
  if (glb_ignore_invalid_chi2)
    glbSetVerbosityLevel(0);


  /* Systematics definitions for systematics ON */
  /* ------------------------------------------ */
  for (i=0; i < in->numofrules; i++)
  {
    /* chi^2 function defined? */
    if (in->sys_on_strings[i] == NULL)
      { glb_rule_error(in, i, "No chi^2 function specified"); status=-1; }
    else
    {
      /* Special case: chiMultiExp */
      if (strcmp(in->sys_on_strings[i], "chiMultiExp") == 0)
      {
        /* Reject normal systematics definitions / require @sys_on_multiex_errors_sig/bg */
        if (in->sys_on_errors[i] != NULL)
        {
          glb_rule_error(in, i, "chiMultiExp is incompatible with @sys_on_errors\n"
                                "Use @sys_on_multiex_errors_sig and @sys_on_multiex_errors_bg");
          status = -1;
        }
        if (in->lengthofrules[i] >= GLB_MAX_CHANNELS  || /* Limitation due to fixed size arrays */
            in->lengthofbgrules[i] >= GLB_MAX_CHANNELS)
        {
          glb_rule_error(in, i, "chiMultiExp works only with rules of at most %d entries",
                         GLB_MAX_CHANNELS-1);
          status = -1;
        }

        /* Check length of @sys_on_multiex_errors_sig and check for invalid entries */
        for (j=0; j < in->lengthofrules[i]; j++)
          if (in->sys_on_multiex_errors_sig[i][j] == NULL) /* Definition too short? */
          {
            glb_rule_error(in, i, "@sys_on_multiex_errors_sig definition missing or too short");
            status = -1;
            j = in->lengthofrules[i] - 1;
          }
          else
          {
            for (k=0; in->sys_on_multiex_errors_sig[i][j][k] >= 0; k++)
              if (glbValueToNameByPointer(in,"sys",in->sys_on_multiex_errors_sig[i][j][k]) == NULL)
              {
                glb_rule_error(in, i, "Invalid entry in @sys_on_multiex_errors_sig");
                status = -1;
              }
            in->sys_on_n_nuis_sig[i][j] = k; /* Remember # of nuisance params for this ch */
          }
        if (in->sys_on_multiex_errors_sig[i][j] != NULL)   /* Definition too long? */
        {
          glb_rule_error(in, i, "@sys_on_multiex_errors_sig definition too long");
          status = -1;
          break;
        }

        /* Check length of @sys_on_multiex_errors_bg and check for invalid entries */
        for (j=0; j < in->lengthofbgrules[i]; j++)
          if (in->sys_on_multiex_errors_bg[i][j] == NULL) /* Definition too short? */
          {
            glb_rule_error(in, i, "@sys_on_multiex_errors_bg definition missing or too short");
            status = -1;
            j = in->lengthofbgrules[i] - 1;
          }
          else
          {
            for (k=0; in->sys_on_multiex_errors_bg[i][j][k] >= 0; k++)
              if (glbValueToNameByPointer(in,"sys",in->sys_on_multiex_errors_bg[i][j][k]) == NULL)
              {
                glb_rule_error(in, i, "Invalid entry in @sys_on_multiex_errors_bg");
                status = -1;
              }
            in->sys_on_n_nuis_bg[i][j] = k;  /* Remember # of nuisance params for this ch */
          }
        if (in->sys_on_multiex_errors_bg[i][j] != NULL)   /* Definition too long? */
        {
          glb_rule_error(in, i, "@sys_on_multiex_errors_bg definition too long");
          status = -1;
          break;
        }
      }
      /* All other chi^2 functions: Check that @sys_on_multiex_errors_sig/bg are absent */
      else
      {
        for (j=0; j < in->lengthofrules[i]; j++)
          if (in->sys_on_multiex_errors_sig[i][j] != NULL)
          {
            glb_rule_error(in, i, "@sys_on_multiex_errors_sig is allowed only in "
                                  "conjunction with chiMultiExp");
            status = -1;
            break;
          }
        for (j=0; j < in->lengthofbgrules[i]; j++)
          if (in->sys_on_multiex_errors_bg[i][j] != NULL)
          {
            glb_rule_error(in, i, "@sys_on_multiex_errors_bg is allowed only in "
                                  "conjunction with chiMultiExp");
            status = -1;
            break;
          }
      }

      /* ALL chi^2 functions: Check length of error list and make sure that all sigma are > 0.0 */
      tmp_errorlist = in->sys_on_errors[i];
      in->sys_on_errors[i] = NULL;
      sys_dim = glbGetSysDim(in->sys_on_strings[i]);
      if (sys_dim >= 0)
      {
        if (tmp_errorlist != NULL)
        {
          for (k=0; tmp_errorlist[k] > 0.0; k++)
            ;
          if (k != sys_dim)
          {
            glb_rule_error(in, i, "Invalid systematical error list @sys_on_errors");
            status=-1;
          }
        }
      }
      else
        { glb_rule_error(in, i, "Invalid @sys_on_function"); status=-1; }

      /* Get pointer to appropriate chi^2 function */
      if (glbSetChiFunctionInExperiment(in, i, GLB_ON, in->sys_on_strings[i],
                                        tmp_errorlist) != 0)
        { glb_rule_error(in, i, "Invalid @sys_on_function/@sys_on_errors"); status=-1; }

      glb_free(tmp_errorlist);
    }

    if (glb_ignore_invalid_chi2 && status != 0)
    {
      if (glbSetChiFunctionInExperiment(in, i, GLB_ON, "chiZero", NULL) == 0)
        status = old_status;
    }
  }

  /* Treatment of parameters for systematics OFF is equivalent to systematics ON */
  /* --------------------------------------------------------------------------- */
  for (i=0; i < in->numofrules; i++)
  {
    if (in->sys_off_strings[i] == NULL)
      { glb_rule_error(in, i, "No chi^2 function specified"); status=-1; }
    else
    {
      if (strcmp(in->sys_off_strings[i], "chiMultiExp") == 0)
      {
        if (in->sys_off_errors[i] != NULL)
        {
          glb_rule_error(in, i, "chiMultiExp is incompatible with @sys_off_errors.\n"
                                "Use @sys_off_multiex_errors_sig and @sys_off_multiex_errors_bg");
          status = -1;
        }
        if (in->lengthofrules[i] >= GLB_MAX_CHANNELS  ||
            in->lengthofbgrules[i] >= GLB_MAX_CHANNELS)
        {
          glb_rule_error(in, i, "chiMultiExp works offly with rules of at most %d entries",
                         GLB_MAX_CHANNELS-1);
          status = -1;
        }
        for (j=0; j < in->lengthofrules[i]; j++)
          if (in->sys_off_multiex_errors_sig[i][j] == NULL)
          {
            glb_rule_error(in, i, "@sys_off_multiex_errors_sig definition missing or too short");
            status = -1;
            j = in->lengthofrules[i] - 1;
          }
          else
          {
            for (k=0; in->sys_off_multiex_errors_sig[i][j][k] >= 0; k++)
              if (glbValueToNameByPointer(in,"sys",in->sys_off_multiex_errors_sig[i][j][k]) == NULL)
              {
                glb_rule_error(in, i, "Invalid entry in @sys_off_multiex_errors_sig");
                status = -1;
              }
            in->sys_off_n_nuis_sig[i][j] = k;/* Remember # of nuisance params for this ch */
          }
        if (in->sys_off_multiex_errors_sig[i][j] != NULL)
        {
          glb_rule_error(in, i, "@sys_off_multiex_errors_sig definition too long");
          status = -1;
          break;
        }
        for (j=0; j < in->lengthofbgrules[i]; j++)
          if (in->sys_off_multiex_errors_bg[i][j] == NULL)
          {
            glb_rule_error(in, i, "@sys_off_multiex_errors_bg definition missing or too short");
            status = -1;
            j = in->lengthofbgrules[i] - 1;
          }
          else
          {
            for (k=0; in->sys_off_multiex_errors_bg[i][j][k] >= 0; k++)
              if (glbValueToNameByPointer(in,"sys",in->sys_off_multiex_errors_bg[i][j][k]) == NULL)
              {
                glb_rule_error(in, i, "Invalid entry in @sys_off_multiex_errors_bg");
                status = -1;
              }
            in->sys_off_n_nuis_bg[i][j] = k; /* Remember # of nuisance params for this ch */
          }
        if (in->sys_off_multiex_errors_bg[i][j] != NULL)
        {
          glb_rule_error(in, i, "@sys_off_multiex_errors_bg definition too long");
          status = -1;
          break;
        }
      }
      else
      {
        for (j=0; j < in->lengthofrules[i]; j++)
          if (in->sys_off_multiex_errors_sig[i][j] != NULL)
          {
            glb_rule_error(in, i, "@sys_off_multiex_errors_sig is allowed only in "
                                  "conjunction with chiMultiExp");
            status = -1;
            break;
          }
        for (j=0; j < in->lengthofbgrules[i]; j++)
          if (in->sys_off_multiex_errors_bg[i][j] != NULL)
          {
            glb_rule_error(in, i, "@sys_off_multiex_errors_bg is allowed only in "
                                  "conjunction with chiMultiExp");
            status = -1;
            break;
          }
      }

      tmp_errorlist = in->sys_off_errors[i];
      in->sys_off_errors[i] = NULL;

      /* Check length of error list and make sure that all sigma are > 0.0 */
      sys_dim = glbGetSysDim(in->sys_off_strings[i]);
      if (sys_dim >= 0)
      {
        if (tmp_errorlist != NULL)
        {
          for (k=0; tmp_errorlist[k] > 0.0; k++)
            ;
          if (k != sys_dim)
          {
            glb_rule_error(in, i, "Invalid systematical error list @sys_off_errors");
            status=-1;
          }
        }
      }
      else
        { glb_rule_error(in, i, "Invalid @sys_off_function"); status=-1; }

      if (glbSetChiFunctionInExperiment(in, i, GLB_OFF, in->sys_off_strings[i],
                                        tmp_errorlist) != 0)
        { glb_rule_error(in, i, "Invalid @sys_off_function/@sys_off_errors"); status=-1; }

      glb_free(tmp_errorlist);
    }

    if (glb_ignore_invalid_chi2 && status != 0)
    {
      if (glbSetChiFunctionInExperiment(in, i, GLB_OFF, "chiZero", NULL) == 0)
        status = old_status;
    }
  }
  if (glb_ignore_invalid_chi2)
    glbSetVerbosityLevel(old_verbosity);


  /* If this experiment is a subdetector in a multi-detector setup, and if it has
   * its own rules, the parent's rules should not be used */
  /* JK - 2012-05-11 Removed because it was very confusing */
//  if (in->parent)
//  {
//    if (in->numofrules > in->parent->numofrules)
//    {
//      for (i=0; i < in->parent->numofrules; i++)
//      {
//        glbSetChiFunctionInRule(in, i, GLB_ON, "chiZero", NULL);
//        glbSetChiFunctionInRule(in, i, GLB_OFF, "chiZero", NULL);
//      }
//    }
//    else
//    {
//      for (i=0; i < in->numofrules; i++)
//      {
//        if (strcmp(in->sys_on[i]->name, "chiSpectrumTilt") != 0 &&
//            strcmp(in->sys_on[i]->name, "chiNoSysSpectrum") != 0 &&
//            strcmp(in->sys_on[i]->name, "chiSpectrumOnly") != 0 &&
//            strcmp(in->sys_on[i]->name, "chiTotalRatesTilt") != 0 &&
//            strcmp(in->sys_on[i]->name, "chiNoSysTotalRates") != 0 &&
//            strcmp(in->sys_on[i]->name, "chiSpectrumCalib") != 0 &&
//            strcmp(in->sys_on[i]->name, "chiZero") != 0)
//        {
//          glb_warning("in rule %d: Rule with non-trivial systematics treatment %s inherited\n"
//                      "  from parent experiment. Make sure correlations among systematical "
//                      "errors are treated correctly", i, in->sys_on[i]->name);
//        }
//
//        if (strcmp(in->sys_off[i]->name, "chiSpectrumTilt") != 0 &&
//            strcmp(in->sys_off[i]->name, "chiNoSysSpectrum") != 0 &&
//            strcmp(in->sys_off[i]->name, "chiSpectrumOnly") != 0 &&
//            strcmp(in->sys_off[i]->name, "chiTotalRatesTilt") != 0 &&
//            strcmp(in->sys_off[i]->name, "chiNoSysTotalRates") != 0 &&
//            strcmp(in->sys_off[i]->name, "chiSpectrumCalib") != 0 &&
//            strcmp(in->sys_off[i]->name, "chiZero") != 0)
//        {
//          glb_warning("in rule %d: Rule with non-trivial systematics treatment %s inherited\n"
//                      "  from parent experiment. Make sure correlations among systematical "
//                      "errors are treated correctly", i, in->sys_off[i]->name);
//        }
//      }
//    }
//  }


  if(in->filter_state==-1){in->filter_state=1;def=-1;}
  if(in->filter_value==-1){in->filter_value=0;def=-1;}


  /* Energy resolution/detector response function (``smearing'') */
  if (in->num_of_sm == -1) { glb_exp_error(in, "No smearing data specified!"); status=-1; }
  for(i=0; i < in->num_of_sm; i++)
  {
    /* Smearing matrix explicitly given in AEDL file? */
    if (in->smear[i] != NULL)
    {
      glb_set_up_smear_data(in->smear_data[i],in);
      glb_default_smear(in->smear_data[i],in);

      /* Cross-checks, relevant if the user redefines the binning in the AEDL file */
      if (in->lowrange[i] == NULL) {
        glb_exp_error(in, "No lowrange defined!"); status=-1;
      }
      else
      {
        int j;
        for (j=0; in->lowrange[i][j] >= 0; j++)
          if (in->lowrange[i][j] < 0 || in->lowrange[i][j] > in->simbins-1) {
            glb_exp_error(in, "in smearing matrix %d, line %d: Invalid lower "
                          "sampling point index: %d", i, j, in->lowrange[i][j]); status=-1;
          }
        if (j != in->numofbins) {
          glb_exp_error(in, "number of lines in smearing matrix %d (%d) does not match number "
                        "of bins (%d) in experiment", i, j, in->numofbins); status=-1;
          }
      }

      if (in->uprange[i] == NULL) {
        glb_exp_error(in, "No uprange defined!"); status=-1;
      }
      else
      {
        int j;
        for (j=0; in->uprange[i][j] >= 0; j++)
          if (in->uprange[i][j] < 0 || in->uprange[i][j] > in->simbins-1) {
            glb_exp_error(in, "in smearing matrix %d, line %d: Invalid upper "
                          "sampling point index: %d", i, j, in->uprange[i][j]); status=-1;
          }
        if (j != in->numofbins) {
          glb_exp_error(in, "number of lines in smearing matrix %d (%d) does not match number "
                        "of bins (%d) in experiment", i, j, in->numofbins); status=-1;
          }
      }
    } /* if (smear[i] != NULL) */

    /* If not, generate it according to the parameters from the AEDL file */
    else
    {
      /* Right now we do not support different binnings
       * in the reconstructed energy */
      in->smear_data[i]->numofbins=in->numofbins;
      in->smear_data[i]->e_min=in->emin;
      in->smear_data[i]->e_max=in->emax;

      if(glb_default_smear(in->smear_data[i],in)==1) def=-1;

      glb_compute_smearing_matrix(&in->smear[i],
                            &in->lowrange[i],&in->uprange[i],
                            in->smear_data[i],in);

      in->simtresh=in->smear_data[i]->e_sim_min;
      in->simbeam=in->smear_data[i]->e_sim_max;
      in->simbins=in->smear_data[i]->simbins;

      if (in->smear[i] == NULL) {
        glb_exp_error(in, "No smear matrix defined!"); status=-1;
      }
      if (in->lowrange[i] == NULL) {
        glb_exp_error(in, "No lowrange defined!"); status=-1;
      }
      if (in->uprange[i] == NULL) {
        glb_exp_error(in, "No uprange defined!"); status=-1;
      }
    } /* else */
  } /* for(i) */

  if(in->simtresh==-1){glb_exp_error(in, "No simtresh defined!");status=-1;}
  if(in->simbeam==-1){glb_exp_error(in, "No simbeam defined!");status=-1;}
  if(in->simbins==-1){glb_exp_error(in, "No simbins defined!");status=-1;}
  if(in->simbins<in->numofbins){glb_exp_error(in, "Less sampling points than bins");
                                status=-1;}

//---------------------------------------------------------

  *in=MInitMemory0(*in);

//-------------------------------------------------------


  if(in->psteps==-1){glb_exp_error(in, "psteps not defined!");status=-1;}
  if(in->lengthtab==NULL){glb_exp_error(in, "lengthtab not allocated!");status=-1;}
  if(in->densitytab==NULL){glb_exp_error(in, "densitytab not allocated!");status=-1;}
  if(in->densitybuffer==NULL){glb_exp_error(in, "densitybuffer not allocated!");
  status=-1;}
  // here we use lot of defaults -- and its not that clear how many
  // there are numofrules ? numofchannels ?
  for(i=0;i<in->numofchannels;i++)
    {
      if(in->simbins<=0) glb_exp_error(in, "Too few simbins defined!");
      else
        {

          if(in->user_pre_smearing_channel[i]==NULL)
            {
              in->user_pre_smearing_channel[i]=(double*)
                glb_malloc(sizeof(double)* in->simbins);
              for(k=0;k<in->simbins;k++)
                in->user_pre_smearing_channel[i][k]=1;
              def=-1;
            }
          else
            {
              for(ct=0;in->user_pre_smearing_channel[i][ct]!=-1;ct++) ct=ct;
              if(ct!=in->simbins)
              {
                glb_exp_error(in, "user_pre_smearing_channel has not simbins elements");
                status=-1;
              }
            }

          if(in->user_pre_smearing_background[i]==NULL)
            {
              in->user_pre_smearing_background[i]=(double*)
                glb_malloc(sizeof(double)*in->simbins);
              for(k=0;k<in->simbins;k++)
                {in->user_pre_smearing_background[i][k]=0;}
              def=-1;
            }
          else
            {
              for(ct=0;in->user_pre_smearing_background[i][ct]!=-1;ct++) ct=ct;
              if(ct!=in->simbins)
              {
                glb_exp_error(in, "user_pre_smearing_background has not simbins elements");
                status=-1;
              }
            }

          if(in->user_post_smearing_channel[i]==NULL)
            {
              in->user_post_smearing_channel[i]=(double*)
                glb_malloc(sizeof(double)*in->numofbins);
              for(k=0;k<in->numofbins;k++)
                in->user_post_smearing_channel[i][k]=1;
              def=-1;
            }
          else
            {
              for(ct=0;in->user_post_smearing_channel[i][ct]!=-1;ct++) ct=ct;
              if(ct!=in->numofbins)
              {
                glb_exp_error(in, "user_post_smearing_channel has not numofbins elements");
                status=-1;
              }
            }
          if(in->user_post_smearing_background[i]==NULL)
            {
              in->user_post_smearing_background[i]=(double*)
                glb_malloc(sizeof(double)*in->numofbins);
              for(k=0;k<in->numofbins;k++)
                {in->user_post_smearing_background[i][k]=0;}
              def=-1;
            }
          else
            {
              for(ct=0;in->user_post_smearing_background[i][ct]!=-1;ct++)
                ct=ct;
              if(ct!=in->numofbins)
              {
                glb_exp_error(in, "user_post_smearing_background has not numofbins elements");
                status=-1;
              }
            }
        }
    }

    for(i=0;i<in->numofrules;i++)
    {
      if (in->energy_window[i][0] == -1  ||  in->energy_window[i][0] < in->emin)
      {
        in->energy_window[i][0] = in->emin;
        def=-1;
      }
      if (in->energy_window[i][1] == -1  ||  in->energy_window[i][1] > in->emax)
      {
        in->energy_window[i][1]=in->emax;
        def=-1;
      }

      if (in->energy_window[i][0] > in->energy_window[i][1])
      {
        tmp                     = in->energy_window[i][1];
        in->energy_window[i][1] = in->energy_window[i][0];
        in->energy_window[i][0] = tmp;
      }

      /* Calculate bin ranges corresponding to the energy window,
       * trying to reproduce the behaviour of the old glb_window_function */
      k = 0;
      while (in->smear_data[0]->bincenter[k] <= in->energy_window[i][0])
        k++;
      in->energy_window_bins[i][0] = k;
      while (k < in->numofbins && in->smear_data[0]->bincenter[k] < in->energy_window[i][1])
        k++;
      in->energy_window_bins[i][1] = k-1;
    }
  for(i=0;i<in->numofrules;i++)
    {
      if(in->SignalRates[i]==NULL||
         in->BackgroundRates[i]==NULL||
         in->rates0[i]==NULL||
         in->rates1[i]==NULL||
         in->rates1BG[i]==NULL)
      {
        glb_rule_error(in, i, "No memory for ratevectors allocated!");
        status=-1;
      }
    }

  if(in->data[i]!=NULL )
    {
      for(ct=0;in->data[i][ct]!=-1;ct++) ct=ct;
      if(ct!=in->numofbins)
	{
	  glb_exp_error(in, "data has not numofbins elements");
	  status=-1;
	}
      
    }

  for(i=0; i < in->numofchannels; i++)
  {
    if(in->chrb_0[i]==NULL || in->chrb_1[i]==NULL ||
       in->chra_0[i]==NULL || in->chra_1[i]==NULL || in->chr_template[i]==NULL)
    {
      glb_exp_error(in, "No memory for convolution allocated!");
      status=-1;
    }
  }

  if(in->energy_tab==NULL){glb_exp_error(in, "energy_tab not allocated!");
  status=-1;}

  /* This is where the name of an oscillation engine is looked up in
     the list osc_engine list and the pointers are reset appropriately
     PH 01/10/19
  */
  
  glb_osc_engine *tmp_osce;

  if(in->osc_engine.name!=NULL) {

    tmp_osce =   glbFindOscEngineByName(in->osc_engine.name);

    if(tmp_osce==NULL) {status=-1;}
    else {
      in->osc_engine.matrix_function=tmp_osce->matrix_function; 
      in->osc_engine.get_function=tmp_osce->get_function;
      in->osc_engine.set_function=tmp_osce->set_function;
      in->osc_engine.num_of_params=tmp_osce->num_of_params;
      in->osc_engine.user_data= tmp_osce->user_data;

    }
  }

  
  /* Final report */
  if(status != 0)
    glb_fatal("Incompletely or incorrectly defined experiment!");
  if(def != 0)
    glb_warning("Incompletely defined experiment! Defaults will be used!");

  return def+status;
};


/***************************************************************************
 * Function glbPrintDblArray                                               *
 ***************************************************************************
 * Helper function for glbPrintExp - prints an array of doubles            *
 ***************************************************************************/
int glbPrintDblArray(double *x, int n)
{
  int j;

  printf("{ ");
  if (x)
  {
    for (j=0; j < n; j++)
    {
      printf("%10.7g", x[j]);
      if (j < n - 1)  printf(", ");
      if (j % 6 == 5)          printf("\n    ");
    }
  }
  printf(" }\n");

  return GLB_SUCCESS;
}


/***************************************************************************
 * Function glbPrintExp                                                    *
 ***************************************************************************
 * Prints selected experiment parameters to stdout in AEDL format          *
 ***************************************************************************/
int glbPrintExpByPointer(struct glb_experiment *e)
{
  int i, j, k;
  int status = GLB_SUCCESS;

  printf("%%!GLoBES\n");
  printf("\n");
  printf("// Experiment loaded from file %s\n", e->filename ? e->filename : "NULL");
  if (e->parent)
    printf("// Based on parent experiment from file %s\n",
           e->parent->filename ? e->parent->filename : "NULL");
  printf("// -------------------------------------------------------------------\n");
  printf("\n");
  if (e->version)
    printf("$version = ""%s""\n", e->version);
  printf("\n");

  printf("// Fluxes\n");
  for (i=0; i < e->num_of_fluxes; i++)
  {
    glb_flux *f = e->fluxes[i];
    printf("nuflux(%s)<\n", glbValueToNameByPointer(e, "nuflux", i));
    if (f->builtin > 0)        printf("  @builtin       = %d\n", f->builtin);
    if (f->time > 0.)          printf("  @time          = %g\n", f->time);
    if (f->target_power > 0.)  printf("  @power         = %g\n", f->target_power);
    if (f->parent_energy > 0.) printf("  @parent_energy = %g\n", f->parent_energy);
    if (f->stored_muons > 0.)  printf("  @stored_muons  = %g\n", f->stored_muons);
    if (f->gamma > 0.)         printf("  @gamma         = %g\n", f->gamma);
    if (f->end_point > 0.)     printf("  @end_point     = %g\n", f->end_point);
    if (f->file_name)          printf("  @file_name     = %s\n", f->file_name);
    if (f->norm > 0.)          printf("  @norm          = %g\n", f->norm);
    printf(">\n");
    printf("\n");
  }
  printf("\n");

  printf("// Cross sections\n");
  for (i=0; i < e->num_of_xsecs; i++)
  {
    glb_xsec *x = e->xsecs[i];
    printf("cross(%s)<\n", glbValueToNameByPointer(e, "cross", i));
    if (x->file_name)          printf("  @cross_file    = %s\n", x->file_name);
    printf(">\n");
    printf("\n");
  }
  printf("\n");

  printf("// Energy resolution\n");
  for (i=0; i < e->num_of_sm; i++)
  {
    printf("energy(%s)<\n", glbValueToNameByPointer(e, "energy", i));
    printf("  @energy = \n");
    for (j=0; j < e->numofbins; j++)
    {
      printf("  { %4d, %4d, ", e->lowrange[i][j], e->uprange[i][j]);
      for (k=0; k < e->uprange[i][j] - e->lowrange[i][j] + 1; k++)
      {
        printf("%10.7g", e->smear[i][j][k]);
        if (k < e->uprange[i][j] - e->lowrange[i][j])  printf(", ");
        if (k % 6 == 5) printf("\n                ");
      }
      if (j < e->numofbins-1)  printf(" }:\n");
      else                     printf(" };\n");
    }
    printf(">\n");
    printf("\n");
  }
  printf("\n");

  printf("// Detector properties\n");
  printf("$target_mass     = %g\n", e->targetmass);
  printf("\n");
  printf("$profiletype     = 3\n");
  printf("$lengthtab       = ");
  glbPrintDblArray(e->lengthtab, e->psteps);
  printf("$densitytab      = ");
  glbPrintDblArray(e->densitytab, e->psteps);
//  printf("$profile_type    = %d\n", e->density_profile_type);
//  printf("$baseline        = %g\n", e->baseline);
  printf("\n");

  printf("// Energy ranges\n");
  printf("$emin            = %g\n", e->emin);
  printf("$emax            = %g\n", e->emax);
  printf("$bins            = %d\n", e->numofbins);
  printf("\n");
  printf("$sampling_min    = %g\n", e->simtresh);
  printf("$sampling_max    = %g\n", e->simbeam);
  printf("$sampling_points = %d\n", e->simbins);
  printf("\n");
  printf("\n");

  printf("// Nuisance parameters\n");
  for (i=0; i < e->n_nuisance; i++)
  {
    printf("sys(%s)<\n", glbValueToNameByPointer(e, "sys", i));
    if (e->nuisance_params[i]->energy_list  &&  e->nuisance_params[i]->error_list)
    {
      printf("  @energy_list = ");
      glbPrintDblArray(e->nuisance_params[i]->energy_list, e->nuisance_params[i]->n_energies);
      printf("  @error_list = ");
      glbPrintDblArray(e->nuisance_params[i]->error_list, e->nuisance_params[i]->n_energies);
    }
    else
      printf("  @error = %g\n", e->nuisance_params[i]->error);
    printf(">\n");
    printf("\n");
  }
  printf("\n");
  printf("\n");


  printf("// Channels\n");
  printf("$filter_state    = %d\n", e->filter_state);
  printf("$filter_value    = %g\n", e->filter_value);
  for (i=0; i < e->numofchannels; i++)
  {
    char fi[10]=" ", ff[10]=" ";
    printf("channel(%s)<\n", glbValueToNameByPointer(e, "channel", i));
    switch(e->listofchannels[2][i])
    {
      case 1:   strcpy(fi, "e"); break;
      case 2:   strcpy(fi, "m"); break;
      case 3:   strcpy(fi, "t"); break;
      case 11:  strcpy(fi, "NOSC_e"); break;
      case 12:  strcpy(fi, "NOSC_m"); break;
      case 13:  strcpy(fi, "NOSC_t"); break;
    }
    switch(e->listofchannels[3][i])
    {
      case 1:   strcpy(ff, "e"); break;
      case 2:   strcpy(ff, "m"); break;
      case 3:   strcpy(ff, "t"); break;
      case 11:  strcpy(ff, "NOSC_e"); break;
      case 12:  strcpy(ff, "NOSC_m"); break;
      case 13:  strcpy(ff, "NOSC_t"); break;
    }
    printf("  @channel = %10s : %c : %6s : %6s : %10s : %10s\n",
           glbValueToNameByPointer(e, "nuflux", e->listofchannels[0][i]),
           e->listofchannels[1][i] > 0 ? '+' : '-', fi, ff,
           glbValueToNameByPointer(e, "cross", e->listofchannels[4][i]),
           glbValueToNameByPointer(e, "energy", e->listofchannels[5][i]));
    printf("  @pre_smearing_efficiencies = ");
    glbPrintDblArray(e->user_pre_smearing_channel[i], e->simbins);
    printf("  @post_smearing_efficiencies = ");
    glbPrintDblArray(e->user_post_smearing_channel[i], e->numofbins);
    printf("  @pre_smearing_background = ");
    glbPrintDblArray(e->user_pre_smearing_background[i], e->simbins);
    printf("  @post_smearing_background = ");
    glbPrintDblArray(e->user_post_smearing_background[i], e->numofbins);
    printf(">\n");
    printf("\n");
  }
  printf("\n");

  printf("// Rules\n");
  for (i=0; i < e->numofrules; i++)
  {
    printf("rule(%s)<\n", glbValueToNameByPointer(e, "rule", i));
    printf("  @signal            = ");
    for (j=0; j < e->lengthofrules[i]; j++)
    {
      printf("%10.7g @ %s", e->rulescoeff[i][j],
             glbValueToNameByPointer(e, "channel", e->rulechannellist[i][j]));
      if (j < e->lengthofrules[i]-1)  printf(" : ");
    }
    printf("\n");
    printf("  @background        = ");
    for (j=0; j < e->lengthofbgrules[i]; j++)
    {
      printf("%10.7g @ %s", e->bgrulescoeff[i][j],
             glbValueToNameByPointer(e, "channel", e->bgrulechannellist[i][j]));
      if (j < e->lengthofbgrules[i]-1)  printf(" : ");
    }
    printf("\n");
    printf("  @sys_on_function   = %s\n", e->sys_on_strings[i]);
    printf("  @sys_off_function  = %s\n", e->sys_off_strings[i]);
    printf("  @sys_on_errors     = ");
    glbPrintDblArray(e->sys_on_errors[i], glbGetSysDim(e->sys_on_strings[i]));
    printf("  @sys_off_errors    = ");
    glbPrintDblArray(e->sys_off_errors[i], glbGetSysDim(e->sys_off_strings[i]));
    if (e->sys_on_multiex_errors_sig[i][0])
    {
      printf("  @sys_on_multiex_errors_sig  = ");
      for (j=0; e->sys_on_multiex_errors_sig[i][j] != NULL; j++)
      {
        if (j > 0) printf(" : ");
        printf("{ ");
        for (k=0; k < e->sys_on_n_nuis_sig[i][j]; k++)
        {
          if (k > 0)  printf(", ");
          printf("%s", glbValueToNameByPointer(e, "sys", e->sys_on_multiex_errors_sig[i][j][k]));
        }
        printf("}");
      }
      printf("\n");
    }
    if (e->sys_on_multiex_errors_bg[i][0])
    {
      printf("  @sys_on_multiex_errors_bg   = ");
      for (j=0; e->sys_on_multiex_errors_bg[i][j] != NULL; j++)
      {
        if (j > 0) printf(" : ");
        printf("{ ");
        for (k=0; k < e->sys_on_n_nuis_bg[i][j]; k++)
        {
          if (k > 0)  printf(", ");
          printf("%s", glbValueToNameByPointer(e, "sys", e->sys_on_multiex_errors_bg[i][j][k]));
        }
        printf("}");
      }
      printf("\n");
    }
    if (e->sys_off_multiex_errors_sig[i][0])
    {
      printf("  @sys_off_multiex_errors_sig = ");
      for (j=0; e->sys_off_multiex_errors_sig[i][j] != NULL; j++)
      {
        if (j > 0) printf(" : ");
        printf("{ ");
        for (k=0; k < e->sys_off_n_nuis_sig[i][j]; k++)
        {
          if (k > 0)  printf(", ");
          printf("%s", glbValueToNameByPointer(e, "sys", e->sys_off_multiex_errors_sig[i][j][k]));
        }
        printf("}");
      }
      printf("\n");
    }
    if (e->sys_off_multiex_errors_bg[i][0])
    {
      printf("  @sys_off_multiex_errors_bg  = ");
      for (j=0; e->sys_off_multiex_errors_bg[i][j] != NULL; j++)
      {
        if (j > 0) printf(" : ");
        printf("{ ");
        for (k=0; k < e->sys_off_n_nuis_bg[i][j]; k++)
        {
          if (k > 0)  printf(", ");
          printf("%s", glbValueToNameByPointer(e, "sys", e->sys_off_multiex_errors_bg[i][j][k]));
        }
        printf("}");
      }
      printf("\n");
    }
    printf(">\n");
    printf("\n");
  }
  printf("\n");

  return status;
}




// this function moves the pointers to the ratevectors to the correct address
// allocated with MInitFile().
// this is crucial for multiple experiment support
// without reprogramming everything !
// here we pass a pointer for efficiency since this function is used in every step

static void MMovePointers(struct glb_experiment *in)
{
  int k;
  
  for (k=0;k<in->numofrules;k++)
    {
      glb_calc_data_on_off[k]=in->data_on_off[k];
      glb_calc_data[k] = in->data[k];
      glb_calc_rates_0[k] = in->rates0[k];
      glb_calc_rates_1[k] = in->rates1[k];
      glb_calc_rates_1BG[k] = in->rates1BG[k];
      glb_calc_signal_rates[k] = in->SignalRates[k];
      glb_calc_bg_rates[k] = in->BackgroundRates[k];
      /* FIXME -- wrong number */
    }
  for(k=0;k<in->num_of_sm;k++) glb_calc_smear_data[k]=in->smear_data[k];
  for(k=0;k<in->numofchannels;k++)
    {
      glb_calc_chrb_0[k] = in->chrb_0[k];
      glb_calc_chrb_1[k] = in->chrb_1[k];
      glb_calc_chra_0[k] = in->chra_0[k];
      glb_calc_chra_1[k] = in->chra_1[k];
      glb_calc_chr_template[k] = in->chr_template[k];
      glb_calc_user_pre_sm_channel[k]=in->user_pre_smearing_channel[k];
      glb_calc_user_post_sm_channel[k]= in->user_post_smearing_channel[k];
      glb_calc_user_pre_sm_background[k]=in->user_pre_smearing_background[k];
      glb_calc_user_post_sm_background[k]=in-> user_post_smearing_background[k];
    }

  for(k=0;k<in->num_of_fluxes;k++) glb_calc_fluxes[k]=in->fluxes[k];
  for(k=0;k<in->num_of_xsecs;k++) glb_calc_xsecs[k]=in->xsecs[k];
  glb_calc_energy_tab = in->energy_tab;
  return;
}

// memory allocation for all the ratevecvtors

static struct glb_experiment MInitMemory0(struct glb_experiment in)
{
 struct glb_experiment out;
  int k;
  int len,l2;
  out=in;
  len=out.numofbins;
  if(out.numofbins>=out.simbins) l2=out.numofbins;
  else l2=out.simbins;
  if(len<=0)
    {glb_exp_error(&in, "Too few bins defined!");}
  else
    {
      for (k=0;k<out.numofrules;k++)
        {
          out.rates0[k] =  (double*) glb_malloc( len*sizeof(double));
          out.rates1[k] = (double*) glb_malloc( len*sizeof(double));
          out.rates1BG[k] = (double*) glb_malloc( len*sizeof(double));
          if(out.simbins<=0) glb_exp_error(&in, "Too few simbins defined!");
          else
            {
              out.SignalRates[k] = (double*) glb_malloc(len*sizeof(double));
              out.BackgroundRates[k] = (double*) glb_malloc(len*sizeof(double));
            }
        }
      out.energy_tab= (double*) glb_malloc(len*sizeof(double));
    }
  for(k=0;k<out.numofchannels;k++)
    {
      if (out.simbins <= 0)
        glb_exp_error(&in, "Too few simbins defined!");
      else
      {
        out.chrb_0[k] = (double*) glb_malloc(out.simbins*sizeof(double));
        out.chrb_1[k] = (double*) glb_malloc(out.simbins*sizeof(double));
        out.chr_template[k] = (double*) glb_malloc(out.simbins*sizeof(double));
      }
      if (out.numofbins <= 0)
        glb_exp_error(&in, "Too few bins defined!");
      else
      {
        out.chra_0[k] = (double*) glb_malloc(out.numofbins*sizeof(double));
        out.chra_1[k] = (double*) glb_malloc(out.numofbins*sizeof(double));
      }
    }
  if(out.psteps<=0) glb_exp_error(&in, "Too few density steps defined!");
  else if (out.densitybuffer==NULL) out.densitybuffer=(double*)
                                      glb_malloc(out.psteps*sizeof(double));
  return out;
}






static void glb_set_profile_scaling_sec(struct glb_experiment *in)
{
  int k;
//  glb_set_length_ptr(in->lengthtab);
//  glb_set_psteps(in->psteps);

  for(k=0;k<in->psteps;k++)
    in->densitybuffer[k]=(double) 1.0*in->densitytab[k];
//  glb_set_density_ptr(in->densitybuffer);

}


// This function actually sets all switches and parameters to the values specified
// in an struct glb_experiment object and calls MMovePointers().
// here we pass a pointer for efficiency since this function is used in every step

void glbSetExperiment(glb_exp in)
{
  int i;
  struct glb_experiment *s;
  s=(struct glb_experiment *) in;

  glb_rule_number=0;
  glb_current_exp = -1;              // Find experiment number
  for (i=0; i < glb_num_of_exps; i++)
    if (glb_experiment_list[i] == s)
    {
      glb_current_exp = i;
      break;
    }


  // pointer moving for density profile
  glb_set_baseline(s->baseline);
  //density=in->densitytab;
  glb_set_profile_scaling_sec(s);
  //analog of SetEnergy

  glb_set_energy_treshold(s->emin);
  glb_set_max_energy(s->emax);
  glb_set_number_of_bins(s->numofbins);


  //analog of SetLuminosity

  glb_set_target_mass(s->targetmass);


  glb_set_num_of_channels(s->numofchannels);
  for (i=0;i<s->numofchannels;i++) glb_set_channel(i,s->listofchannels[0][i],
                                               s->listofchannels[1][i],
                                               s->listofchannels[2][i],
                                               s->listofchannels[3][i],
                                               s->listofchannels[4][i],
                                               s->listofchannels[5][i]
                                               );

  // analog of SetRules


  glb_set_number_of_rules(s->numofrules);
  for (i=0;i<s->numofrules;i++)
    {
      glb_set_rule(i,s->lengthofrules[i],
                   (int*) s->rulechannellist[i],s->rulescoeff[i]);
      glb_set_bg_rule(i,s->lengthofbgrules[i],
                              s->bgrulechannellist[i],s->bgrulescoeff[i]);
      glb_set_bg_center(i,(double) s->bg_centers[0][i],(double) s->bg_centers[1][i]);

      // this is preliminary (0 should be replaced by i)
      glb_calc_simbins=s->simbins;
      glb_calc_simtresh=s->simtresh;
      glb_calc_simbeam=s->simbeam;
    }

  for(i=0;i<s->num_of_sm;i++)
    {
      glb_calc_smearing[i]=s->smear[i];
      glb_calc_uprange[i]=s->uprange[i];
      glb_calc_lowrange[i]=s->lowrange[i];
    }

  /* set the current oscillation engine - PH 01/10/19 */
  glb_switch_osc_engine(&(in->osc_engine));


//  glb_switch_filter(s->filter_state);
//  glb_set_filter(s->filter_value);



  MMovePointers(s);

  return;
}

// redefinition of SetRates()

void glbSetRates()
{
  int i;
  for (i=0;i<glb_num_of_exps;i++)
    {
      glbSetExperiment(glb_experiment_list[i]);
      glb_set_rates();
    }
}

// redefinition of SetNewRates()

void glbSetNewRates()
{
  int i;
  for (i=0;i<glb_num_of_exps;i++)
    {
      glbSetExperiment(glb_experiment_list[i]);
      glb_set_new_rates();
    }
}

// a bunch of comands which allow to modify
// certain elements of a struct glb_experiment object
// in glb_experiment_list


void glb_set_profile_scaling(double scale,int i)
{
  int k;
//  glb_set_length_ptr(glb_experiment_list[i]->lengthtab);

//  glb_set_psteps(glb_experiment_list[i]->psteps);
  for(k=0;k<glb_experiment_list[i]->psteps;k++)
    {
      glb_experiment_list[i]->densitybuffer[k]=(double)
        scale*glb_experiment_list[i]->densitytab[k];
    }
//  glb_set_density_ptr(glb_experiment_list[i]->densitybuffer);

}


