#ifndef __DEGFINDER_H
#define __DEGFINDER_H

#include <math.h>
#include <float.h>
#include <globes/globes.h>

/* Mass hierarchies */
#define HIERARCHY_NORMAL     1
#define HIERARCHY_INVERTED  -1

/* Options for degfinder */
#define DEG_NO_NH      0x01   /* Omit normal hierarchy fits                    */
#define DEG_NO_IH      0x02   /* Omit inverted hierarchy fits                  */
#define DEG_NO_SYS     0x04   /* Switch off systematics                        */
#define DEG_NO_CORR    0x10   /* Switch off correlations                       */
#define DEG_NO_DEG     0x20   /* Switch off all degeneracies                   */

/* Define NaN (not a number) */
#ifdef NAN
  #define GLB_NAN      NAN
  #define GLB_ISNAN(x) isnan(x)
#else
  #define GLB_NAN      DBL_MAX
  #define GLB_ISNAN(x) ((x) == DBL_MAX)
#endif


/* degfinder.c */
int degfinder(const glb_params base_values, const int n_prescan_params,
      const int *prescan_params, const double *prescan_min,
      const double *prescan_max, const int *prescan_steps,
      const glb_projection prescan_proj, const glb_projection fit_proj,
      int *n_deg, glb_params *deg_pos, double *deg_chi2, const long flags);

#endif

