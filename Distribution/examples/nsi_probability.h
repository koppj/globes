#ifndef __NSI_PROBABILITY_H
#define __NSI_PROBABILITY_H

#include <globes/globes.h>

/* Names of non-standard parameters */
#define FIRST_NSI_INDEX  6
enum epsilon { ABS_EPS_S_EE=FIRST_NSI_INDEX, ARG_EPS_S_EE,
               ABS_EPS_S_MUE,    ARG_EPS_S_MUE,
               ABS_EPS_S_TAUE,   ARG_EPS_S_TAUE,
               ABS_EPS_S_EMU,    ARG_EPS_S_EMU, 
               ABS_EPS_S_MUMU,   ARG_EPS_S_MUMU,
               ABS_EPS_S_TAUMU,  ARG_EPS_S_TAUMU,
               ABS_EPS_S_ETAU,   ARG_EPS_S_ETAU,
               ABS_EPS_S_MUTAU,  ARG_EPS_S_MUTAU,
               ABS_EPS_S_TAUTAU, ARG_EPS_S_TAUTAU,

               EPS_M_EE, 
               ABS_EPS_M_EMU,    ARG_EPS_M_EMU,
               ABS_EPS_M_ETAU,   ARG_EPS_M_ETAU,
               EPS_M_MUMU,
               ABS_EPS_M_MUTAU,  ARG_EPS_M_MUTAU,
               EPS_M_TAUTAU,

               ABS_EPS_D_EE,     ARG_EPS_D_EE,
               ABS_EPS_D_MUE,    ARG_EPS_D_MUE,
               ABS_EPS_D_TAUE,   ARG_EPS_D_TAUE,
               ABS_EPS_D_EMU,    ARG_EPS_D_EMU, 
               ABS_EPS_D_MUMU,   ARG_EPS_D_MUMU,
               ABS_EPS_D_TAUMU,  ARG_EPS_D_TAUMU,
               ABS_EPS_D_ETAU,   ARG_EPS_D_ETAU,
               ABS_EPS_D_MUTAU,  ARG_EPS_D_MUTAU,
               ABS_EPS_D_TAUTAU, ARG_EPS_D_TAUTAU };
#define LAST_NSI_INDEX ARG_EPS_D_TAUTAU

/* Effective electron numbers for calculation of matter profile */
#define Ne_MANTLE       0.497
#define Ne_CORE         0.468

#define REARTH           6371  /* km */
#define RCORE            3480  /* km */

/* Macros */
#define SQR(x)      ((x)*(x))                              // x^2 
#define POW10(x)    (exp(M_LN10*(x)))                      // 10^x
#define MIN(X,Y)    ( ((X) < (Y)) ? (X) : (Y) )
#define MAX(X,Y)    ( ((X) > (Y)) ? (X) : (Y) )
#define SIGN(a,b)   ( ((b) > 0.0) ? (fabs(a)) : (-fabs(a)) )
#define SGN(a)      ( ((a) >= 0.0) ? (1) : (-1) )
#define ROUND(x)    ( (int)((x) + 0.5) )

/* Function declarations */
/* --------------------- */

/* probability.c */               
int nsi_init_probability_engine();
int nsi_free_probability_engine();
int nsi_set_oscillation_parameters(glb_params p, void *user_data);
int nsi_get_oscillation_parameters(glb_params p, void *user_data);
int nsi_probability_matrix(double P[3][3], int cp_sign, double E,
        int psteps, const double *length, const double *density,
        double filter_sigma, void *user_data);
int nsi_filtered_probability_matrix_cd(double P[3][3], double E, double L, double V,
                                       double sigma, int cp_sign);

/* prem.c */
int LoadPREMProfile(const char *prem_file);
double GetPREMDensity(double t, double L);
double GetAvgPREMDensity(double L_tot, double L1, double L2);
int GetPREM3LayerApprox(double L, int *n_layers, double *lengths,
                        double *densities);

#endif

