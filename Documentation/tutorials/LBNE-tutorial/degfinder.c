/* ---------------------------------------------------------------------------- */
/* Finding degeneracies in neutrino oscillations                                */
/* ---------------------------------------------------------------------------- */
/* Author: Joachim Kopp (partly based on ideas by P. Huber and W. Winter)       */
/* ---------------------------------------------------------------------------- */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <globes/globes.h>
#include "degfinder.h"

/* Global variables */
const int debug_level = 0;

#define MIN(X,Y)    ( ((X) < (Y)) ? (X) : (Y) )
#define MAX(X,Y)    ( ((X) > (Y)) ? (X) : (Y) )


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

  if (!GLB_ISNAN(th12))
    glbSetOscParams(tv, th12, GLB_THETA_12);
  if (!GLB_ISNAN(th13))
    glbSetOscParams(tv, th13, GLB_THETA_13);
  if (!GLB_ISNAN(th23))
    glbSetOscParams(tv, th23, GLB_THETA_23);
  if (!GLB_ISNAN(delta))
    glbSetOscParams(tv, delta, GLB_DELTA_CP);
  if (!GLB_ISNAN(dm21))
    glbSetOscParams(tv, dm21, GLB_DM_21);
  if (!GLB_ISNAN(dm31))
    glbSetOscParams(tv, dm31, GLB_DM_31);

  if (hierarchy == HIERARCHY_NORMAL)
  {
    glbSetOscParams(tv, fabs(glbGetOscParams(tv,GLB_DM_31)), GLB_DM_31);
    glbSetOscParams(cv, fabs(glbGetOscParams(cv,GLB_DM_31)), GLB_DM_31);
    glbSetCentralValues(cv);
    result = glbChiNP(tv, fit_values, GLB_ALL);
  }
  else if (hierarchy == HIERARCHY_INVERTED)
  {
    glbSetOscParams(tv, -fabs(glbGetOscParams(tv,GLB_DM_31)), GLB_DM_31);
    glbSetOscParams(cv, -fabs(glbGetOscParams(cv,GLB_DM_31)), GLB_DM_31);
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
int convert_index(int j, int d, int n, double *p_steps)
/* ---------------------------------------------------------------------------- */
/* Function for converting one-dimensional indices to a multi-dimensional       */
/* index for (d+1)-th dimension                                                 */
/* ---------------------------------------------------------------------------- */
/* Parameters:                                                                  */
/*   j: Index of prescan point for parameter d                                  */
/*   d: Index of parameter                                                      */
/*   n: Total number of points scanned in prescan over multi-dim. space         */
/*   p_steps: Number of prescan steps for each parameter in prescan             */
/* ---------------------------------------------------------------------------- */
{
  int i;
  for (i=0; i <= d; i++)
  {
    j %= n;
    n /= p_steps[i] + 1;
  }
  j /= n;
  return j;
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
/*   n_deg: Number of degenerate solutions                                      */
/*   deg_pos: Positions of degeneracies in parameter space                      */
/*   deg_chi2: chi^2 values of the degenerate solutions                         */
/* ---------------------------------------------------------------------------- */
{
  int i, j, k, m;

  /* Copy input parameters to private data structures */
  long private_flags = flags;
  int n_p_params = n_prescan_params;
  int p_params[n_p_params];
  double p_max[n_p_params], p_min[n_p_params], p_steps[n_p_params];
  for (i=0; i < n_p_params; i++)
  {
    p_params[i] = prescan_params[i];
    p_min[i]    = prescan_min[i];
    p_max[i]    = prescan_max[i];
    p_steps[i]  = prescan_steps[i];
  }
  glb_projection private_prescan_proj = glbAllocProjection();
  glb_projection private_fit_proj     = glbAllocProjection();
  glbCopyProjection(prescan_proj, private_prescan_proj);
  glbCopyProjection(fit_proj, private_fit_proj);
  
  /* For runs without systematics, switch them off now */
  int old_sys_state[glb_num_of_exps][128];
  for (j=0; j < glb_num_of_exps; j++)
    for (k=0; k < glbGetNumberOfRules(j); k++)
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
  
  /* For runs without correlations, switch them off */
  if (private_flags & DEG_NO_CORR)
  {
    for (i=0; i < glbGetNumOfOscParams(); i++)
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
  for (i=0; i < n_p_params; i++)
  {
    p_steps[i] = prescan_steps[i];
    p_min[i]   = prescan_min[i];
    if (p_params[i] == GLB_DELTA_CP)
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
 
  /* Compute number of points in prescan grid */
  unsigned long n_prescan_points = 1;
  for (i=0; i < n_p_params; i++)
    n_prescan_points *= p_steps[i] + 1;

  if (debug_level > 0)
    printf("#   Using %lu prescan points.\n", n_prescan_points);
  
  /* Create data structures */
  double chi2_table_NH[n_prescan_points];
  double chi2_table_IH[n_prescan_points];
  glb_params Fit_NH = glbAllocParams();
  glb_params Fit_IH = glbAllocParams();
  glb_params param_table_NH[n_prescan_points];
  glb_params param_table_IH[n_prescan_points];
  for (j=0; j < n_prescan_points; j++)
  {
    param_table_NH[j] = glbAllocParams();
    param_table_IH[j] = glbAllocParams();    
  }
  memset(chi2_table_NH, 0.0, sizeof(chi2_table_NH[0]) * n_prescan_points);
  memset(chi2_table_IH, 0.0, sizeof(chi2_table_NH[0]) * n_prescan_points);

  if (debug_level > 0)
  {
    printf("#   Parameters used in prescan: ");
    for (i=0; i < n_p_params; i++)
    {
      switch (p_params[i])
      {
        case GLB_THETA_12:
          printf("th12 ");
          break;
        case GLB_THETA_13:
          printf("th13 ");
          break;
        case GLB_THETA_23:
          printf("th23 ");
          break;
        case GLB_DELTA_CP:
          printf("delta ");
          break;
        case GLB_DM_21:
          printf("dm21 ");
          break;
        case GLB_DM_31:
          printf("dm31 ");
          break;
        default:
          printf("unknown ");
          break;
      }
    }
    printf("\n");
  }

  
  /* Prescan */
  /* ------- */
 
  glbSetProjection(private_prescan_proj); 
  glbCopyParams(base_values, Fit_NH);
  glbCopyParams(base_values, Fit_IH);
  for (j=0; j < n_prescan_points; j++)
  {
    /* Compute parameter values at current grid point */
    double prescan_test_values[n_p_params];
    for (i=0; i < n_p_params; i++)
    {
      /* Use log-distribution for sin^2 2 th13 */
      if (p_params[i] == GLB_THETA_13)
      {
        double log_param = p_min[i]
          + convert_index(j,i,n_prescan_points,p_steps) * (p_max[i]-p_min[i])/p_steps[i];
        prescan_test_values[i] = asin(sqrt(pow(10.0, log_param)))/2.0;
      }
      /* Use linear distribution for all other parameters */
      else
      {
        prescan_test_values[i] = p_min[i]
          + convert_index(j,i,n_prescan_points,p_steps) * (p_max[i]-p_min[i])/p_steps[i];
      }

      glbSetOscParams(Fit_NH, prescan_test_values[i], p_params[i]);
      glbSetOscParams(Fit_IH, prescan_test_values[i], p_params[i]);
    }

    /* Compute chi^2 _without_ systematics for both hierarchies */
    int old_sys_state2[glb_num_of_exps][128];
    for (m=0; m < glb_num_of_exps; m++)
      for (k=0; k < glbGetNumberOfRules(m); k++)
        old_sys_state2[m][k] = glbGetSysOnOffState(m, k);
    glbSwitchSystematics(GLB_ALL, GLB_ALL, GLB_OFF);
    if ( !(private_flags & DEG_NO_NH) )
    {
      chi2_table_NH[j] = ChiNPWrapper(Fit_NH, GLB_NAN, GLB_NAN, GLB_NAN,
                                      GLB_NAN, GLB_NAN, GLB_NAN, HIERARCHY_NORMAL, Fit_NH);
      glbCopyParams(Fit_NH, param_table_NH[j]);
    }
    if ( !(private_flags & DEG_NO_IH) )
    {
      chi2_table_IH[j] = ChiNPWrapper(Fit_IH, GLB_NAN, GLB_NAN, GLB_NAN,
                                      GLB_NAN, GLB_NAN, GLB_NAN, HIERARCHY_INVERTED, Fit_IH);
      glbCopyParams(Fit_IH, param_table_IH[j]);
    }
    for (m=0; m < glb_num_of_exps; m++)
      for (k=0; k < glbGetNumberOfRules(m); k++)
        glbSwitchSystematics(m, k, old_sys_state2[m][k]);

    if (debug_level > 1)
    {
      printf("#   Test params:");
      for (i=0; i < n_p_params; i++)
        printf(" %8.5g", prescan_test_values[i]);
      printf(", chi2_NH = %12.7g, chi2_IH = %12.7g\n", chi2_table_NH[j], chi2_table_IH[j]);
    }
  }

  
  /* Determine locations of the degeneracies from chi^2 tables and perform fit */
  /* ------------------------------------------------------------------------- */

  glbSetProjection(private_fit_proj);
  *n_deg = 0;
  for (j=0; j < n_prescan_points; j++)
  {
    /* Compute differences to the chi^2 values from neighbouring points to
     * find minima */
    int k = n_prescan_points;
    int min_NH = 1, min_IH = 1;
    for (i=0; i < n_p_params; i++)
    {
      int k2 = k;
      k /= (p_steps[i] + 1);
      int s = (j+k) / k2 - j / k2;
      int t = (j+k2) / k2 - (j+k2-k) / k2;
                          /* This is 1 if wraparound happens, and 0 otherwise */
      
      /* Treat phases as cyclic prameters */
      if (p_params[i] == GLB_DELTA_CP)
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
        for (k=0; k < n_p_params; k++)
          printf("%d ", convert_index(j, k, n_prescan_points, p_steps));
        printf("\n");
        printf("#     Prescan: ");
        for (k=0; k < 6; k++)
          printf(" %10.5g", glbGetOscParams(param_table_NH[j], k));
        printf(", chi2 = %10.5g\n", chi2_table_NH[j]);
      }

      deg_chi2[*n_deg] = ChiNPWrapper(param_table_NH[j], GLB_NAN, GLB_NAN, GLB_NAN,
                                      GLB_NAN, GLB_NAN, GLB_NAN,
                                      HIERARCHY_NORMAL, deg_pos[*n_deg]);

      if (debug_level > 0)
      {
        printf("#     Fit:     ");
        for (k=0; k < 6; k++)
          printf(" %10.5g", glbGetOscParams(deg_pos[*n_deg], k));
        printf(", chi2 = %10.5g\n", deg_chi2[*n_deg]);
      }
      (*n_deg)++;
    }
    
    if ( !(private_flags & DEG_NO_IH) && min_IH )  /* We have found a minimum for the inverted hierarchy */
    {
      if (debug_level > 0)
      {
        printf("#   Degeneracy IH at: ");
        for (k=0; k < n_p_params; k++)
          printf("%d ", convert_index(j, k, n_prescan_points, p_steps));
        printf("\n");
        printf("#     Prescan: ");
        for (k=0; k < 6; k++)
          printf(" %10.5g", glbGetOscParams(param_table_IH[j], k));
        printf(", chi2 = %10.5g\n", chi2_table_IH[j]);
      }

      deg_chi2[*n_deg] = ChiNPWrapper(param_table_IH[j], GLB_NAN, GLB_NAN, GLB_NAN,
                                      GLB_NAN, GLB_NAN, GLB_NAN,
                                      HIERARCHY_INVERTED, deg_pos[*n_deg]);

      if (debug_level > 0)
      {
        printf("#     Fit:     ");
        for (k=0; k < 6; k++)
          printf(" %10.5g", glbGetOscParams(deg_pos[*n_deg], k));
        printf(", chi2 = %10.5g\n", deg_chi2[*n_deg]);
      }
      (*n_deg)++;
    }
  }
  

  /* Clean up */
  if (private_flags & DEG_NO_SYS)
  {
    for (j=0; j < glb_num_of_exps; j++)
      for (k=0; k < glbGetNumberOfRules(j); k++)
        glbSwitchSystematics(j, k, old_sys_state[j][k]);
  }
  for (j=0; j < n_prescan_points; j++)
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


