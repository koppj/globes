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



#ifndef GLS_GLOBES_H
#define GLS_GLOBES_H 1

#include <stdio.h>

#ifdef __cplusplus
#  define BEGIN_C_DECLS extern "C" {
#  define END_C_DECLS   }
#else /* !__cplusplus */
#  define BEGIN_C_DECLS
#  define END_C_DECLS
#endif /* __cplusplus */


// Number of neutrino flavours
#define GLB_NU_FLAVOURS             3


#define GLB_EFF 1
#define GLB_BG 2
#define GLB_SIG 3
#define GLB_PRE 1
#define GLB_POST 2
#define GLB_OSCP 6
#define GLB_ALL -1
#define GLB_W_EFF 1
#define GLB_WO_EFF 2
#define GLB_W_BG 3
#define GLB_WO_BG 4
#define GLB_W_COEFF 5
#define GLB_WO_COEFF 6
#define GLB_FIXED 0
#define GLB_FREE 1
#define GLB_THETA_12 0
#define GLB_THETA_13 1
#define GLB_THETA_23 2
#define GLB_DELTA_CP 3
#define GLB_DM_SOL 4 
#define GLB_DM_ATM 5
#define GLB_ON 1
#define GLB_OFF 0

#define GLB_EARTH_RADIUS 6371.0

// Unit conversion
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

BEGIN_C_DECLS

int glbSetVerbosityLevel(int level);
void glbPrintDelimiter(FILE *stream, int character);

glb_params glbSetIteration(glb_params in, int iter);
int glbGetIteration(const glb_params in);
glb_params glbCopyParams(const glb_params source, glb_params dest);
glb_params glbDefineParams(glb_params in,
			     double theta12, double theta13,double theta23,
			     double delta, 
			     double dms, double dma);
glb_params glbAllocParams();
void glbFreeParams(glb_params stale);
glb_params glbSetDensityParams(glb_params in,double dens, int which);
double glbGetDensityParams(glb_params in, int which);
glb_params glbSetOscParams(glb_params in, double osc, int which);
double glbGetOscParams(glb_params in, int which);
void glbPrintParams(FILE *stream, const glb_params in);



glb_projection glbSetProjectionFlag(glb_projection in, int flag, int which);
int glbGetProjectionFlag(const glb_projection in, int which);
glb_projection glbSetDensityProjectionFlag(glb_projection in, int flag, int which);
int glbGetDensityProjectionFlag(const glb_projection in, int which);
glb_projection glbCopyProjection(const glb_projection source, glb_projection dest);
glb_projection glbDefineProjection(glb_projection in,
			     int theta12, int theta13,int theta23,
			     int delta, 
			     int dms, int dma);
glb_projection glbAllocProjection();
void glbFreeProjection(glb_projection stale);
void glbPrintProjection(FILE *stream, const glb_projection in);



/* Setting the parameters */

int glbSetOscillationParameters(const glb_params in);
int glbGetOscillationParameters(glb_params in);


/* The probablities */
double glbVacuumProbability(int pl, int pm, int panti,double pen, 
			    double plength);
double glbProfileProbability(int exp,int initial_flavour, int final_flavour,
			     int panti, double energy);
double glbFilteredConstantDensityProbability(int exp,int initial_flavour, int final_flavour,
					     int panti, double energy);
/* Fluxes and X-sections */

double glbXSection(int experiment, int xsec_ident,double energy,int flavour, 
		   int anti);

/* Chi^2 including systematics */

double glbChiSys(const glb_params in,int experiment, int rule);

/* setting all the priors needed for projections */

int glbSetStartingValues(const glb_params in);
int glbSetInputErrors(const glb_params in);
int glbGetStartingValues(glb_params in);
int glbGetInputErrors(glb_params in);

double glbChiTheta(const glb_params in, glb_params out, int exp);
double glbChiTheta23(const glb_params in, glb_params out, int exp);
double glbChiDelta(const glb_params in, glb_params out, int exp);
double glbChiDms(const glb_params in, glb_params out, int exp);
double glbChiDm(const glb_params in, glb_params out, int exp);
double glbChiThetaDelta(const glb_params in, glb_params out, int exp);
double glbChiAll(const glb_params in, glb_params out, int exp);

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


int glbSetTargetMass(int experiment, double mass);
double glbGetTargetMass(int experiment);

int glbSetSignalErrors(int experiment, int rule, double norm, double tilt);
int glbGetSignalErrors(int experiment, int rule, double *norm, double *tilt);

int glbSetBGErrors(int experiment, int rule, double norm, double tilt);
int glbGetBGErrors(int experiment, int rule, double *norm, double *tilt);

int glbSetBGCenters(int experiment, int rule, double norm, double tilt);
int glbGetBGCenters(int experiment, int rule, double *norm, double *tilt);

int glbGetNumberOfFluxes(int exp);
double glbFlux(int experiment, int flux_ident, 
	double energy, double distance, int flavour, int anti);
int glbSetSourcePower(int experiment, int flux_ident, double power);
double glbGetSourcePower(int experiment, int flux_ident);
int glbSetRunningTime(int experiment, int flux_ident, double time);
double glbGetRunningTime(int experiment, int flux_ident);


/* gls_parser.c */
int glbInitExperiment(char *inf, glb_exp *in, int *counter);
void glbClearExperimentList();
void glbInit(char *name);

/* Channel Access */
void glbResetRateStack();

int glbGetNumberOfChannels(int exp);
int glbGetChannelRates(double **data, size_t *length, 
		       int exp, int channel,int smearing);
double glbGetChannelBin(int exp, int channel, int bin);
int glbGetUserData(double **data, size_t *length, 
		   int exp, int channel,int smearing, int bgeff);

/* Rule Access */
double glbTotalRuleRate(int exp, int rule, int pos,int effi, 
			int bgi, int coeffi, int signal);
int glbGetNumberOfRules(int exp);
int glbGetLengthOfRule(int exp, int rule, int signal);
double glbGetNormalizationInRule(int exp, int rule, int signal);
int glbGetChannelInRule(int exp, int rule, int pos, int signal);
double glbGetCoefficientInRule(int exp, int rule, int pos, int signal);
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
double glbGetSignalBinInRule(int exp, int rule, int bin);
double glbGetBackgroundBinInRule(int exp, int rule, int bin);
double glbGetObservedBinInRule(int exp, int rule, int bin);


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

int glbTestReleaseVersion(const char *version);
int glbTestLibraryVersion(const char *version);
const char *glbVersionOfExperiment(int experiment);

/* Filter feature access */


int glbSetFilterState(int on_off);
int glbGetFilterState();
int glbSetFilterStateInExperiment(int experiment,int on_off);
int glbGetFilterStateInExperiment(int experiment);
int glbSetFilter(double filter);
double glbGetFilter();
int glbSetFilterInExperiment(int experiment,double filter);
double glbGetFilterInExperiment(int experiment);

/* Translation of named variables to numbers and vice versa */
int glbNameToValue(int exp, const char* context, const char *name);
const char *glbValueToName(int exp,const char* context, int value);

/* Trying to get a tree-like linkage */
/* void glbInit(char *name); */

/* General Module support */
int glbProbeModule(const char *module_name, int verbosity);
glb_dlhandle glbOpenModule(const char *module_name);
int glbCloseModule(glb_dlhandle stale);
void *glbSymModule(glb_dlhandle module,const char *symbol_name);

/* Adding user defined priors */
#ifndef SWIG
int glbRegisterPriorFunction(double (*prior)(const glb_params),
			     int (*starting)(const glb_params),
			     int (*error)(const glb_params));
#endif /* SWIG */
int glbUsePrior(glb_dlhandle module);


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

END_C_DECLS

#endif /* GLS_GLOBES_H 1 */
