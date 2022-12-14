/* ---------------------------------------------------------------------------- */
/* Finding degeneracies in neutrino oscillations                                */
/* ---------------------------------------------------------------------------- */
/* Author: Joachim Kopp (partly based on ideas by P. Huber and W. Winter)       */
/* ---------------------------------------------------------------------------- */
/* Note:                                                                        */
/*  - This file is written in C99. To compile it in gcc, use the command line   */
/*    option -std=gnu99                                                         */
/*  - When you use this code to produce a scientific publication, please        */
/*    cite the following reference:                                             */
/*                                                                              */
/*      @Article{Kopp:2008ds,                                                   */
/*        author    = "Kopp, Joachim and Ota, Toshihiko and Winter, Walter",    */
/*        title     = "{Neutrino factory optimization for non-standard          */
/*                     interactions}",                                          */
/*        journal   = "Phys. Rev.",                                             */
/*        volume    = "D78",                                                    */
/*        year      = "2008",                                                   */
/*        pages     = "053007",                                                 */
/*        eprint    = "0804.2261",                                              */
/*        archivePrefix = "arXiv",                                              */
/*        primaryClass  =  "hep-ph",                                            */
/*        doi       = "10.1103/PhysRevD.78.053007",                             */
/*        memo      = "JK's degfinder developed and first used in this paper",  */
/* ---------------------------------------------------------------------------- */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <globes/globes.h>
#include "degfinder.h"
/*#include "nsi_probability.h"*/ /* Uncomment if you want to use degfinder with JK's NSI engine */

#ifndef FIRST_NSI_INDEX
#define FIRST_NSI_INDEX -100
#define LAST_NSI_INDEX  -101
#endif

#define SQR(x)      ((x)*(x))                              // x^2 
#define POW10(x)    (exp(M_LN10*(x)))                      // 10^x
#define MIN(X,Y)    ( ((X) < (Y)) ? (X) : (Y) )
#define MAX(X,Y)    ( ((X) > (Y)) ? (X) : (Y) )
#define SIGN(a,b)   ( ((b) > 0.0) ? (fabs(a)) : (-fabs(a)) )
#define SGN(a)      ( ((a) >= 0.0) ? (1) : (-1) )
#define ROUND(x)    ( (int)((x) + 0.5) )

/* Global variables */
extern const char *param_strings[];
const int debug_level = 0;


/* ---------------------------------------------------------------------------- */
double ChiNPWrapper(glb_params base_values, double th12, double th13, double th23,
                    double delta, double dm21, double dm31, int hierarchy,
                    glb_params fit_values)
/* ---------------------------------------------------------------------------- */
/* Perform minimization for a certain set of test values                        */
/* ---------------------------------------------------------------------------- */
{
  double result;
  glb_params tv     = glbAllocParams();      /* Test values             */
  glb_params cv     = glbAllocParams();      /* Central values          */
  glb_params old_cv = glbAllocParams();      /* Previous central values */
  glbCopyParams(base_values, tv);
  glbGetCentralValues(cv);
  glbGetCentralValues(old_cv);

  if (!isnan(th12))
    glbSetOscParams(tv, th12, GLB_THETA_12);
  if (!isnan(th13))
    glbSetOscParams(tv, th13, GLB_THETA_13);
  if (!isnan(th23))
    glbSetOscParams(tv, th23, GLB_THETA_23);
  if (!isnan(delta))
    glbSetOscParams(tv, delta, GLB_DELTA_CP);
  if (!isnan(dm21))
    glbSetOscParams(tv, dm21, GLB_DM_21);
  if (!isnan(dm31))
    glbSetOscParams(tv, dm31, GLB_DM_31);

  // Note: Converting the dm31^2 for the normal hierarchy into a dm31^2 for the
  // inverted hierarchy that would give similar oscillation probabilities has bee
  // discussed in hep-ph/0503079 (see also hep-ph/0509359 and 0812.1879).
  // While a suitable transformation can be found approximately, it is different
  // for different oscillation channels. In all cases, it is similar to 2*s12^2*dm21^2
  // or 2*c12^2*dm21^2, which, very roughly, is of O(dm21^2), so that's what
  // we use here.
  double dm31_tv = glbGetOscParams(tv,GLB_DM_31);
  double dm31_cv = glbGetOscParams(cv,GLB_DM_31);
  if (hierarchy == HIERARCHY_NORMAL)
  {
    if (dm31_tv < 0)  glbSetOscParams(tv, -dm31_tv + glbGetOscParams(tv,GLB_DM_21), GLB_DM_31);
    if (dm31_cv < 0)  glbSetOscParams(cv, -dm31_cv + glbGetOscParams(cv,GLB_DM_21), GLB_DM_31);
    glbSetCentralValues(cv);
    result = glbChiNP(tv, fit_values, GLB_ALL);
  }
  else if (hierarchy == HIERARCHY_INVERTED)
  {
    if (dm31_tv > 0)  glbSetOscParams(tv, -dm31_tv + glbGetOscParams(tv,GLB_DM_21), GLB_DM_31);
    if (dm31_cv > 0)  glbSetOscParams(cv, -dm31_cv + glbGetOscParams(cv,GLB_DM_21), GLB_DM_31);
    glbSetCentralValues(cv);
    result = glbChiNP(tv, fit_values, GLB_ALL);
  }
  else
  {
    fprintf(stderr, "ChiNPWrapper: Please specify HIERARCHY_NORMAL or HIERARCHY_INVERTED!\n");
    result = -1.0;
  }

  glbSetCentralValues(old_cv);
  glbFreeParams(old_cv);
  glbFreeParams(cv);
  glbFreeParams(tv);

  return result;
}


/* ---------------------------------------------------------------------------- */
int degfinder(const glb_params base_values, const int n_prescan_params,
      const int *prescan_params, const double *prescan_min,
      const double *prescan_max, const int *prescan_steps,
      const glb_projection prescan_proj, const glb_projection fit_proj,
      int *n_deg, glb_params *deg_pos, double *deg_chi2, const long flags)
/* ---------------------------------------------------------------------------- */
/* Input parameters:                                                            */
/*   base_values: The oscillation parameters                                    */
/*   n_prescan_params: Number of parameters to perform prescan on               */
/*   prescan_params: Indices of the parameters for the prescan                  */
/*   prescan_min, prescan_max, prescan_steps: Minimum/Maximum values and        */
/*     numbers of steps for the prescan                                         */
/*   prescan_proj, fit_proj: Projections for the prescan and for the final fit  */
/*   flags: A combination of the DEG_XXX flags                                  */
/* Output parameters:                                                           */
/*   n_deg: Input:  Max. number of degenerate solutions to return               */
/*          Output: Number of degenerate solutions found                        */
/*   deg_pos: Positions of degeneracies in parameter space                      */
/*   deg_chi2: chi^2 values of the degenerate solutions                         */
/* ---------------------------------------------------------------------------- */
{
  /* Copy input parameters to private data structures */
  long private_flags = flags;
  int n_p_params = n_prescan_params;
  int p_params[n_p_params];
  double p_max[n_p_params], p_min[n_p_params];
  int p_steps[n_p_params];
  int p_flags[n_p_params];
  for (int i=0; i < n_p_params; i++)
  {
    p_params[i] = prescan_params[i];
    p_min[i]    = prescan_min[i];
    p_max[i]    = prescan_max[i];
    p_steps[i]  = prescan_steps[i];
    p_flags[i]  = 0;

  }
  glb_projection private_prescan_proj = glbAllocProjection();
  glb_projection private_fit_proj     = glbAllocProjection();
  glbCopyProjection(prescan_proj, private_prescan_proj);
  glbCopyProjection(fit_proj, private_fit_proj);
  
  /* Determine which parameters should be scanned on a log scale */
  for (int i=0; i < n_p_params; i++)
  {
    if (p_params[i] == GLB_THETA_13)
      p_flags[i] |= DEG_LOGSCALE;
    else if (p_params[i] >= FIRST_NSI_INDEX && p_params[i] <= LAST_NSI_INDEX &&
             strstr(param_strings[p_params[i]], "ARG") == NULL)
    {
      p_flags[i] |= (DEG_LOGSCALE | DEG_PLUS_MINUS);

      /* If phase of this parameter is explicitly included -> OK. Otherwise, use
       * DEG_PLUS_MINUS flag to at least include the sign ambiguity */
      for (int j=0; j < n_p_params; j++)
      {
        if (p_params[j] == p_params[i]+1  &&  strstr(param_strings[p_params[j]], "ARG") != NULL)
          p_flags[i] &= (~DEG_PLUS_MINUS);
      }
    }
  }

  /* For runs without systematics, switch them off now */
  int old_sys_state[glb_num_of_exps][128];
  for (int j=0; j < glb_num_of_exps; j++)
    for (int k=0; k < glbGetNumberOfRules(j); k++)
      old_sys_state[j][k] = glbGetSysOnOffState(j, k);
  if (private_flags & DEG_NO_SYS)
    glbSwitchSystematics(GLB_ALL, GLB_ALL, GLB_OFF);

  /* For runs without degeneracies, omit prescan and wrong hierarchy solutions */
  if (private_flags & DEG_NO_DEG)
  {
    n_p_params = 0;
    if (glbGetOscParams(base_values, GLB_DM_31) > 0)
      private_flags |= DEG_NO_IH;
    else
      private_flags |= DEG_NO_NH;
  }
  
  /* For runs without NSI degeneracies, remove all NSI parameters from
   * parameter list */
  if (private_flags & DEG_NO_NSI_DEG)
  {
    for (int i=0; i < n_p_params; i++)
    {
      if (p_params[i] >= FIRST_NSI_INDEX && p_params[i] <= LAST_NSI_INDEX)
      {
        for (int j=i; j < n_p_params-1; j++)
        {
          p_params[j] = p_params[j+1];
          p_min[j]    = p_min[j+1];
          p_max[j]    = p_max[j+1];
        }
        n_p_params--;
      }
    }
  }

  /* For runs without correlations, switch them off */
  if (private_flags & DEG_NO_CORR)
  {
    for (int i=0; i < glbGetNumOfOscParams(); i++)
    {
      glbSetProjectionFlag(private_prescan_proj, GLB_FIXED, i);
      glbSetProjectionFlag(private_fit_proj, GLB_FIXED, i);
    }
    glbSetDensityProjectionFlag(private_prescan_proj, GLB_FIXED, GLB_ALL);
    glbSetDensityProjectionFlag(private_fit_proj, GLB_FIXED, GLB_ALL);
  }


  /* Preparations for prescan */
  /* ------------------------ */
  
  /* For cyclic parameters (complex phases), make sure the upper boundary
   * of the scan region is omitted if it would imply double counting */
  for (int i=0; i < n_p_params; i++)
  {
    p_steps[i] = prescan_steps[i];
    p_min[i]   = prescan_min[i];
    if ( p_params[i] == GLB_DELTA_CP ||
         (p_params[i] >= FIRST_NSI_INDEX && p_params[i] <= LAST_NSI_INDEX &&
          strstr(param_strings[p_params[i]], "ARG") != NULL) )
    {
      if (fabs(prescan_min[i]) < 1e-10 && fabs(prescan_max[i] - 2*M_PI) < 1e-10)
      {
        p_max[i]    = prescan_max[i] - (prescan_max[i] - prescan_min[i])/prescan_steps[i];
        p_steps[i] -= 1;
      }
      else
        p_max[i] = prescan_max[i];
    }
    else
      p_max[i] = prescan_max[i];
  }

  /* If we are asked to use log scale, but consider positive and negative parameter
   * values, we have to use twice as many prescan points */
  for (int i=0; i < n_p_params; i++)
    if (p_flags[i] & DEG_LOGSCALE  &&  p_flags[i] & DEG_PLUS_MINUS)
      p_steps[i] = 2*p_steps[i] + 1;
 
  /* Compute number of points in prescan grid */
  unsigned long n_prescan_points = 1;
  for (int i=0; i < n_p_params; i++)
    n_prescan_points *= p_steps[i] + 1;

  /* Function for converting one-dimensional indices to a multi-dimensional
   * index for (d+1)-th dimension */
  int convert_index(int j, int d)
  {
    int k = n_prescan_points;
    for (int i=0; i <= d; i++)
    {
      j %= k;
      k /= p_steps[i] + 1;
    }
    j /= k;
    return j;
  }

  if (debug_level > 0)
    printf("#   Using %lu prescan points.\n", n_prescan_points);
  
  /* Create data structures */
  double chi2_table_NH[n_prescan_points];
  double chi2_table_IH[n_prescan_points];
  glb_params Fit_NH = glbAllocParams();
  glb_params Fit_IH = glbAllocParams();
  glb_params param_table_NH[n_prescan_points];
  glb_params param_table_IH[n_prescan_points];
  for (int j=0; j < n_prescan_points; j++)
  {
    param_table_NH[j] = glbAllocParams();
    param_table_IH[j] = glbAllocParams();    
  }
  memset(chi2_table_NH, 0.0, sizeof(chi2_table_NH[0]) * n_prescan_points);
  memset(chi2_table_IH, 0.0, sizeof(chi2_table_NH[0]) * n_prescan_points);

  if (debug_level > 0)
  {
    printf("#   Parameters used in prescan: ");
    for (int i=0; i < n_p_params; i++)
    {
      printf("%s ", param_strings[p_params[i]]);
      if (p_flags[i] & DEG_LOGSCALE)
        printf("(log) ");
      if (p_flags[i] & DEG_PLUS_MINUS)
        printf("(+-) ");
    }
    printf("\n");
  }

  
  /* Prescan */
  /* ------- */
 
  glbSetProjection(private_prescan_proj); 
  glbCopyParams(base_values, Fit_NH);
  glbCopyParams(base_values, Fit_IH);
  for (int j=0; j < n_prescan_points; j++)
  {
    /* Compute parameter values at current grid point */
    double prescan_test_values[n_p_params];
    for (int i=0; i < n_p_params; i++)
    {
      if (p_flags[i] & DEG_LOGSCALE)        /* Use log scale for this param? */
      {
        double x;
       
        if (p_flags[i] & DEG_PLUS_MINUS)
        {
          int k = convert_index(j,i);
          if (k <= (p_steps[i]-1) / 2)
            x =  -POW10( p_max[i] - k * (p_max[i]-p_min[i])/(p_steps[i]/2) );
          else
            x =   POW10( p_min[i] + (k - (p_steps[i]+1)/2) * (p_max[i]-p_min[i])/(p_steps[i]/2) );
        }
        else
          x = POW10( p_min[i] + convert_index(j,i) * (p_max[i]-p_min[i])/p_steps[i] );

        if (p_params[i] == GLB_THETA_13)
          prescan_test_values[i] = asin(sqrt(x))/2.0;
        else
          prescan_test_values[i] = x;
      }
      else                                  /* Otherwise: linear distribution */
      {
        prescan_test_values[i] = p_min[i]
          + convert_index(j,i) * (p_max[i]-p_min[i])/p_steps[i];
      }

      glbSetOscParams(Fit_NH, prescan_test_values[i], p_params[i]);
      glbSetOscParams(Fit_IH, prescan_test_values[i], p_params[i]);
    }

    /* Compute chi^2 _without_ systematics for both hierarchies */
    int old_sys_state2[glb_num_of_exps][128];
    for (int j=0; j < glb_num_of_exps; j++)
      for (int k=0; k < glbGetNumberOfRules(j); k++)
        old_sys_state2[j][k] = glbGetSysOnOffState(j, k);
    glbSwitchSystematics(GLB_ALL, GLB_ALL, GLB_OFF);
    if ( !(private_flags & DEG_NO_NH) )
    {
      chi2_table_NH[j] = ChiNPWrapper(Fit_NH, NAN, NAN, NAN, NAN, NAN, NAN, HIERARCHY_NORMAL, Fit_NH);
      glbCopyParams(Fit_NH, param_table_NH[j]);
    }
    if ( !(private_flags & DEG_NO_IH) )
    {
      chi2_table_IH[j] = ChiNPWrapper(Fit_IH, NAN, NAN, NAN, NAN, NAN, NAN, HIERARCHY_INVERTED, Fit_IH);
      glbCopyParams(Fit_IH, param_table_IH[j]);
    }
    for (int j=0; j < glb_num_of_exps; j++)
      for (int k=0; k < glbGetNumberOfRules(j); k++)
        glbSwitchSystematics(j, k, old_sys_state2[j][k]);

    if (debug_level > 1)
    {
      printf("#   Test params:");
      for (int i=0; i < n_p_params; i++)
        printf(" %8.5g", prescan_test_values[i]);
      printf(", chi2_NH = %12.7g, chi2_IH = %12.7g\n", chi2_table_NH[j], chi2_table_IH[j]);
    }
  }

  
  /* Determine locations of the degeneracies from chi^2 tables and perform fit */
  /* ------------------------------------------------------------------------- */

  glbSetProjection(private_fit_proj);
  int max_deg = *n_deg;
  *n_deg = 0;
  for (int j=0; j < n_prescan_points; j++)
  {
    /* Compute differences to the chi^2 values from neighbouring points to
     * find minima */
    int k = n_prescan_points;
    int min_NH = 1, min_IH = 1;
    for (int i=0; i < n_p_params; i++)
    {
      int k2 = k;
      k /= (p_steps[i] + 1);
      int s = (j+k) / k2 - j / k2;
      int t = (j+k2) / k2 - (j+k2-k) / k2;
                          /* This is 1 if wraparound happens, and 0 otherwise */
      
      /* Treat phases as cyclic prameters */
      if ( p_params[i] == GLB_DELTA_CP ||
           (p_params[i] >= FIRST_NSI_INDEX && p_params[i] <= LAST_NSI_INDEX &&
            strstr(param_strings[p_params[i]], "ARG") != NULL) )
      {
        if (chi2_table_NH[j] - chi2_table_NH[j + k - s*k2] > 0 ||
            chi2_table_NH[j] - chi2_table_NH[j - k + t*k2] > 0)
          min_NH = 0;
        if (chi2_table_IH[j] - chi2_table_IH[j + k - s*k2] > 0 ||
            chi2_table_IH[j] - chi2_table_IH[j - k + t*k2] > 0)
          min_IH = 0;
      }
      else
      {
        if (chi2_table_NH[j] - chi2_table_NH[MIN(n_prescan_points - (k-j%k), j+k)] > 0 ||
            chi2_table_NH[j] - chi2_table_NH[MAX(0, j-k)] > 0)
          min_NH = 0;
        if (chi2_table_IH[j] - chi2_table_IH[MIN(n_prescan_points - (k-j%k), j+k)] > 0 ||
            chi2_table_IH[j] - chi2_table_IH[MAX(0, j-k)] > 0)
          min_IH = 0;
      }
    }
    
    if ( !(private_flags & DEG_NO_NH) && min_NH )  /* We have found a minimum for the normal hierarchy */
    {
      if (debug_level > 0)
      {
        printf("#   Degeneracy NH at: ");
        for (int k=0; k < n_p_params; k++)
          printf("%d (%s = %g)    ", convert_index(j, k), param_strings[p_params[k]],
                 glbGetOscParams(param_table_NH[j], p_params[k]));
        printf("\n");
        printf("#     Prescan: ");
        for (int k=0; k < 6; k++)
          printf(" %7.4g", glbGetOscParams(param_table_NH[j], k));
        for (int k=0; k < n_p_params; k++)
          if (p_params[k] >= 6)
            printf(" %7.4g", glbGetOscParams(param_table_NH[j], p_params[k]));
        printf(", chi2 = %10.5g\n", chi2_table_NH[j]);
      }

      if (*n_deg >= max_deg)
      {
        fprintf(stderr, "degfinder: Too many degeneracies fond (max. is %d).\n", max_deg);
        break;
      }
      deg_chi2[*n_deg] = ChiNPWrapper(param_table_NH[j], NAN, NAN, NAN, NAN, NAN, NAN,
                                      HIERARCHY_NORMAL, deg_pos[*n_deg]);

      if (debug_level > 0)
      {
        printf("#     Fit:     ");
        for (int k=0; k < 6; k++)
          printf(" %7.4g", glbGetOscParams(deg_pos[*n_deg], k));
        for (int k=0; k < n_p_params; k++)
          if (p_params[k] >= 6)
            printf(" %7.4g", glbGetOscParams(deg_pos[*n_deg], p_params[k]));
        printf(", chi2 = %10.5g\n", deg_chi2[*n_deg]);
      }
      (*n_deg)++;
    }
    
    if ( !(private_flags & DEG_NO_IH) && min_IH )  /* We have found a minimum for the inverted hierarchy */
    {
      if (debug_level > 0)
      {
        printf("#   Degeneracy IH at: ");
        for (int k=0; k < n_p_params; k++)
          printf("%d (%s = %g)    ", convert_index(j, k), param_strings[p_params[k]],
                 glbGetOscParams(param_table_IH[j], p_params[k]));
        printf("\n");
        printf("#     Prescan: ");
        for (int k=0; k < 6; k++)
          printf(" %7.4g", glbGetOscParams(param_table_IH[j], k));
        for (int k=0; k < n_p_params; k++)
          if (p_params[k] >= 6)
            printf(" %7.4g", glbGetOscParams(param_table_IH[j], p_params[k]));
        printf(", chi2 = %10.5g\n", chi2_table_IH[j]);
      }

      if (*n_deg >= max_deg)
      {
        fprintf(stderr, "degfinder: Too many degeneracies fond (max. is %d).\n", max_deg);
        break;
      }
      deg_chi2[*n_deg] = ChiNPWrapper(param_table_IH[j], NAN, NAN, NAN, NAN, NAN, NAN,
                                      HIERARCHY_INVERTED, deg_pos[*n_deg]);

      if (debug_level > 0)
      {
        printf("#     Fit:     ");
        for (int k=0; k < 6; k++)
          printf(" %7.4g", glbGetOscParams(deg_pos[*n_deg], k));
        for (int k=0; k < n_p_params; k++)
          if (p_params[k] >= 6)
            printf(" %7.4g", glbGetOscParams(deg_pos[*n_deg], p_params[k]));
        printf(", chi2 = %10.5g\n", deg_chi2[*n_deg]);
      }
      (*n_deg)++;
    }
  } // for(j)
  
  /* Clean up */
  if (private_flags & DEG_NO_SYS)
  {
    for (int j=0; j < glb_num_of_exps; j++)
      for (int k=0; k < glbGetNumberOfRules(j); k++)
        glbSwitchSystematics(j, k, old_sys_state[j][k]);
  }
  for (int j=0; j < n_prescan_points; j++)
  {
    if (param_table_NH[j] != NULL)
      glbFreeParams(param_table_NH[j]);
   if (param_table_IH[j] != NULL)
      glbFreeParams(param_table_IH[j]);
  }
  glbFreeParams(Fit_NH);
  glbFreeParams(Fit_IH);
  glbFreeProjection(private_prescan_proj);
  glbFreeProjection(private_fit_proj);
  
  return 0;
}


