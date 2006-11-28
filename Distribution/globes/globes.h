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
#ifndef __GLOBES_H
#define __GLOBES_H 1

#include <stdio.h>

#ifdef __cplusplus
#  define BEGIN_C_DECLS extern "C" {
#  define END_C_DECLS   }
#else /* !__cplusplus */
#  define BEGIN_C_DECLS
#  define END_C_DECLS
#endif /* __cplusplus */


/* Constants */
/* --------- */

/* Activate compatibility mode, i.e. allow deprecated symbol names */
#define GLB_COMPAT

/* Number of neutrino flavours */
#define GLB_NU_FLAVOURS             3

/* Number of oscillation parameters                                         */
/* ---> May need to be modified when implementing non-standard physics <--- */
#define GLB_OSCP 6


// JK - I prefer enums over #defines since they are parsed by the compiler,
// which is more intelligent than the pre-processor.
enum glb_enum_oscillation_parameters
       { GLB_THETA_12 = 0, GLB_THETA_13 = 1, GLB_THETA_23 = 2,
         GLB_DELTA_CP = 3, GLB_DM_21 = 4, GLB_DM_31 = 5 };
enum glb_enum_switches
       { GLB_OFF=0, GLB_ON=1 };
enum glb_enum_rate_specs
       { GLB_W_EFF = 1, GLB_WO_EFF = 2,
         GLB_W_BG = 3, GLB_WO_BG = 4,
         GLB_W_COEFF = 5, GLB_WO_COEFF = 6 };
enum glb_enum_param_fixed_free
       { GLB_FIXED = 0, GLB_FREE = 1 };
enum glb_enum_efficiency_types
       { GLB_PRE = 1, GLB_POST = 2 };

#define GLB_EFF    1
#define GLB_BG     2
#define GLB_SIG    3

#define GLB_ALL   -1

#define GLB_EARTH_RADIUS 6371.0 /* km */

/* Unit conversion */
//#define EV_TO_KM_FACTOR  1.973e-10
#define EV_TO_KM_FACTOR  1.9747235e-10
#define KM_TO_EV(x)      ((x) / EV_TO_KM_FACTOR)
#define EV_TO_KM(x)      ((x) * EV_TO_KM_FACTOR)





/* In the application software we only need a handle for the
 * any lt_dlopen`ed module ...
 */
typedef struct lt_dlhandle_struct *glb_dlhandle;

typedef struct glb_projection_type *glb_projection;
typedef struct glb_params_type  *glb_params;



/* handle for the struct experiment */

typedef struct glb_experiment *glb_exp;

extern char **glb_path_vector;
extern size_t glb_path_vector_length;
extern int glb_num_of_exps;

extern int glb_single_experiment_number;
#ifdef SWIG
/* I really have no clue why I have to 
 * declare this one to SWIG and the other global 
 * variables not
 */
glb_exp glb_experiment_list[32];
#else
extern glb_exp glb_experiment_list[32];
#endif
extern int glb_rule_number;

/* Function declarations */
/* --------------------- */
BEGIN_C_DECLS

/* Initialization */
void glbInit(char *name);


/* Loading of AEDL files */


/* Version control / verbosity */
int glbTestReleaseVersion(const char *version);
int glbTestLibraryVersion(const char *version);
const char *glbVersionOfExperiment(int experiment);
int glbSetVerbosityLevel(int level);


/* Handling of oscillation parameter data structure */
glb_params glbAllocParams();
void glbFreeParams(glb_params stale);
glb_params glbDefineParams(glb_params in, double theta12, double theta13,
                           double theta23, double delta, double dm21, double dm31);
glb_params glbCopyParams(const glb_params source, glb_params dest);
glb_params glbSetOscParams(glb_params in, double osc, int which);
glb_params glbSetDensityParams(glb_params in,double dens, int which);
glb_params glbSetIteration(glb_params in, int iter);
double glbGetOscParams(glb_params in, int which);
double glbGetDensityParams(glb_params in, int which);
int glbGetIteration(const glb_params in);
void glbPrintParams(FILE *stream, const glb_params in);

int glbSetOscillationParameters(const glb_params in);
int glbGetOscillationParameters(glb_params in);


/* Screen output */


/* Event rate calculation */


/* chi^2 projections */
glb_projection glbAllocProjection();
void glbFreeProjection(glb_projection stale);
glb_projection glbDefineProjection(glb_projection in, int theta12,
    int theta13,int theta23, int delta, int dm21, int dm31);
glb_projection glbCopyProjection(const glb_projection source, glb_projection dest);
glb_projection glbSetProjectionFlag(glb_projection in, int flag, int which);
glb_projection glbSetDensityProjectionFlag(glb_projection in, int flag, int which);
int glbGetProjectionFlag(const glb_projection in, int which);
int glbGetDensityProjectionFlag(const glb_projection in, int which);
void glbPrintProjection(FILE *stream, const glb_projection in);

int glbSetProjection(const glb_projection in);
int glbGetProjection(glb_projection in);

double glbChiSys(const glb_params in,int experiment, int rule);
double glbChiTheta13(const glb_params in, glb_params out, int exp);
double glbChiTheta12(const glb_params in, glb_params out, int exp);
double glbChiTheta23(const glb_params in, glb_params out, int exp);
double glbChiDelta(const glb_params in, glb_params out, int exp);
double glbChiDm21(const glb_params in, glb_params out, int exp);
double glbChiDm31(const glb_params in, glb_params out, int exp);
double glbChiTheta13Delta(const glb_params in, glb_params out, int exp);
double glbChiNP(const glb_params in, glb_params out, int exp);
double glbChiAll(const glb_params in, glb_params out, int exp);


/* Matter profile handling */


/* Access to miscellaneous experiment parameters */
int glbSetTargetMass(int experiment, double mass);
int glbSetSourcePower(int experiment, int flux_ident, double power);
int glbSetRunningTime(int experiment, int flux_ident, double time);
int glbSetFilterStateInExperiment(int experiment,int on_off);
int glbSetFilterInExperiment(int experiment,double filter);

double glbGetTargetMass(int experiment);
double glbGetSourcePower(int experiment, int flux_ident);
double glbGetRunningTime(int experiment, int flux_ident);
int glbGetFilterStateInExperiment(int experiment); 
double glbGetFilterInExperiment(int experiment);
int glbGetNumberOfSimBins(int exp);
int glbGetNumberOfBins(int exp);
int glbGetNumberOfRules(int exp);
int glbGetNumberOfChannels(int exp);
int glbGetLengthOfRule(int exp, int rule, int signal);
double glbGetNormalizationInRule(int exp, int rule, int signal);
int glbGetChannelInRule(int exp, int rule, int pos, int signal);
double glbGetCoefficientInRule(int exp, int rule, int pos, int signal);
int glbGetNumberOfFluxes(int exp);
double glbFlux(int experiment, int flux_ident, 
        double energy, double distance, int flavour, int anti);
double glbXSection(int experiment, int xsec_ident,double energy,int flavour,
                   int anti);


/* User-defined Systematics */


/* Modules and user-defined priors */
int glbProbeModule(const char *module_name, int verbosity);
glb_dlhandle glbOpenModule(const char *module_name);
int glbCloseModule(glb_dlhandle stale);
void *glbSymModule(glb_dlhandle module,const char *symbol_name);

#ifndef SWIG
int glbRegisterPriorFunction(double (*prior)(const glb_params),
                             int (*starting)(const glb_params),
                             int (*error)(const glb_params));
#endif /* SWIG */
int glbUsePrior(glb_dlhandle module);


/* Access to event rate vectors */
double glbTotalRuleRate(int exp, int rule, int pos,int effi, 
			int bgi, int coeffi, int signal);
double *glbGetChannelRatePtr(int exp, int channel, int pre_post);
double *glbGetRuleRatePtr(int exp, int rule);
double *glbGetSignalRatePtr(int exp, int rule);
double *glbGetBGRatePtr(int exp, int rule);
double *glbGetChannelFitRatePtr(int exp, int channel, int pre_post);
//double *glbRuleFitRatePtr(int exp, int rule);
    // JK - new; this will not be implemented since the BG normalization is
    // usually free and thus only known during chi^2 calculation
double *glbGetSignalFitRatePtr(int exp, int rule);
double *glbGetBGFitRatePtr(int exp, int rule);

double *glbGetEfficiencyPtr(int exp, int channel, int pre_post);
double *glbGetBackgroundPtr(int exp, int channel, int pre_post);


/* Access to the probablity engine */
double glbVacuumProbability(int initial_flavour, int final_flavour,
                            int cp_sign, double E, double L);
double glbConstantDensityProbability(int initial_flavour, int final_flavour,
                            int cp_sign, double E, double L, double rho);
double glbProfileProbability(int exp,int initial_flavour, int final_flavour,
                            int panti, double energy);
double glbFilteredConstantDensityProbability(int exp,int initial_flavour,
                            int final_flavour, int panti, double energy);


/**********************************************************************/


void glbPrintDelimiter(FILE *stream, int character);



/* Setting the parameters */

int glbSetOscillationParameters(const glb_params in);
int glbGetOscillationParameters(glb_params in);


/* setting all the priors needed for projections */

int glbSetStartingValues(const glb_params in);
int glbSetInputErrors(const glb_params in);
int glbGetStartingValues(glb_params in);
int glbGetInputErrors(glb_params in);

/* Functions for user-defined chi^2 calculations */

void glbSetUserChi(int exp, int rule, double (*chi_func)(double x[]), int dim,
                   double params[], double errors[], char *info);
int glbGetCurrentExp();
int glbGetCurrentRule();
int glbShiftSignalEnergyScale(int exp, int rule, double amount);
int glbShiftBackgroundEnergyScale(int exp, int rule, double amount);

/* Interface to the projected chi^2 functions */

int glbSetProjection(const glb_projection in);
int glbGetProjection(glb_projection in);
double glbChiNP(const glb_params in, glb_params out, int exp);



void glbSetExperiment(glb_exp in);
int glbDefaultExp(glb_exp ins);
void glbInitExp(glb_exp ins);
glb_exp glbAllocExp();
void glbFreeExp(glb_exp ins);


void glbSetNewRates();
void glbSetRates();

int glbSetErrorDim(int experiment, int rule, int on_off, int value);
int glbGetErrorDim(int experiment, int rule, int on_off);

int glbSwitchSystematics(int experiment, int rule, int on_off);



int glbSetSignalErrors(int experiment, int rule, double norm, double tilt);
int glbGetSignalErrors(int experiment, int rule, double *norm, double *tilt);

int glbSetBGErrors(int experiment, int rule, double norm, double tilt);
int glbGetBGErrors(int experiment, int rule, double *norm, double *tilt);

int glbSetBGCenters(int experiment, int rule, double norm, double tilt);
int glbGetBGCenters(int experiment, int rule, double *norm, double *tilt);


/* gls_parser.c */
int glbInitExperiment(char *inf, glb_exp *in, int *counter);
void glbClearExperimentList();

/* Channel Access */
void glbResetRateStack();

int glbGetChannelRates(double **data, size_t *length, 
		       int exp, int channel,int smearing);
//FIXME Remove
//double glbGetChannelBin(int exp, int channel, int bin);
int glbGetUserData(double **data, size_t *length, 
		   int exp, int channel,int smearing, int bgeff);

/* Rule Access */
int glbShowRuleRates(FILE *stream,
		 int exp, int rule, int pos,
		 int effi, int bgi, int coeffi, int signal);
int glbShowChannelRates(FILE *stream,
			int exp, int channel, int smearing, int effi, int bgi);
int glbShowChannelProbs(FILE *stream,
			int exp, int channel, int smearing, int effi, int bgi);
void *glbSetChannelPrintFunction(void *fp);
void glbSetPrintDelimiters(const char *left,const char *middle,
			   const char *right);
//FIXME Remove
//double glbGetSignalBinInRule(int exp, int rule, int bin);
//double glbGetBackgroundBinInRule(int exp, int rule, int bin);
//double glbGetObservedBinInRule(int exp, int rule, int bin);


/* Matter profile access */
int glbGetProfileTypeInExperiment(int exp);
int glbLoadProfileData(const char* filename, size_t *layers, double **length, 
		       double **density );
int glbStaceyProfile(double baseline, size_t layers, double **length, 
		     double **density);
int glbAverageDensityProfile(double baseline, double **length, 
			     double **density);
int glbGetProfileData(size_t *layers, double** length, double** density);
int glbGetProfileDataInExperiment(int exp,size_t *layers, double** length, 
				  double** density);
int glbSetProfileDataInExperiment(int exp, size_t layers,const double* length, 
				  const double* density);
int glbSetBaselineInExperiment(int exp, double baseline);
double glbGetBaselineInExperiment(int exp);


void glbDefineAEDLVariable(const char *name, double value);
void glbClearAEDLVariables();

/* Filter feature access */


int glbSetFilter(double filter);
double glbGetFilter();
int glbSetFilterInExperiment(int experiment,double filter);
double glbGetFilterInExperiment(int experiment);

/* Translation of named variables to numbers and vice versa */
int glbNameToValue(int exp, const char* context, const char *name);
const char *glbValueToName(int exp,const char* context, int value);

/* Trying to get a tree-like linkage */
/* void glbInit(char *name); */


#ifdef GLB_EXPERIMENTAL
/* These functions are part of the user-defined chi^2 interface and will
 * not be documented in this release. I you want to use this feature
 * compile libglobes with '-DGLB_EXPERIMENTAL' and define GLB_EXPERIMENTAL
 * before you include globes.h into your program.
 */
  char* glbSetSys(int select, int rule, int experiment);
  void glbSetSysErrors(double *val, int rule, int experiment);
  void glbSetSysCenter(double *val, int rule, int experiment);
#endif /* GLB_EXPERIMENTAL */


/* Some symbols which are provided for compatibility with */
/* older versiond of GLoBES                               */
#ifdef GLB_COMPAT
  #define GLB_DM_SOL  GLB_DM_21   /* Redefined to avoid confusion of dm31 and dm32 */
  #define GLB_DM_ATM  GLB_DM_31 
    
  double glbChiTheta(const glb_params in, glb_params out, int exp);     /* Replaced by glbChiTheta13 */
  double glbChiDms(const glb_params in, glb_params out, int exp);       /* Replaced by glbChiDm21 */
  double glbChiDm(const glb_params in, glb_params out, int exp);        /* Replaced by glbChiDm31 */
  double glbChiThetaDelta(const glb_params in, glb_params out, int exp);/* Replaced by glbChiTheta13Delta */
#endif  /* #ifdef GLB_COMPAT */


END_C_DECLS

#endif /* GLS_GLOBES_H 1 */
