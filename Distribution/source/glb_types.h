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
#ifndef GLB_TYPES_H
#define GLB_TYPES_H 1

#if HAVE_CONFIG_H   /* config.h should come before any other includes */
#  include "config.h"
#endif

#include "globes/globes.h"

/* This is a close relative of glb_namerec, which is needed in order
 * to facilitate the cooperation between AEDL and C.
 */

/* this serves to handle the 5.2 fudge factor in versions prior to 3.0
 */
#define GLB_OLD_NORM -99

struct glb_naming
{
  char *name; /* name of symbol      */   
  char *context;  /* context rule, channel etc.           */
  int value;  /* e.g. rule number */
  struct glb_naming *next;    /* link field              */
};

typedef struct glb_naming glb_naming;

/* Data structure for handling fluxes and built-in fluxes */
#define GLB_FLUX_COLUMNS   7           /* Number of columns in flux file         */
typedef struct {
  int builtin;
  double time;
  double parent_energy;
  double stored_muons;
  double target_power;
  double norm;
  double gamma;
  double end_point;

  int n_lines;                         /* Number of columns in flux file         */
  char *file_name;                     /* Name of flux file                      */
  double *flux_data[GLB_FLUX_COLUMNS]; /* Flux data (7 x n_lines array)          */
} glb_flux;

/* Data structure for handling X-sections and built-in X-sections,
 * which will be added later.
 */
#define GLB_XSEC_COLUMNS   7           /* Number of columns in xsec file         */
typedef struct {
  int builtin;
  int n_lines;                         /* Number of lines in cross section file  */
  char *file_name;                     /* Name of cross section file             */
  double *xsec_data[GLB_XSEC_COLUMNS]; /* Cross section data (7 x n_lines array) */
} glb_xsec;



/* This are the additional data needed for type B matrices */

typedef struct {
  double corr_fac;
  double confidence_level;
  int offset;
  double low_bound;
  double up_bound;
} glb_option_type;


/* This is a pointer to the function which computes sigma(E,params) */
typedef double (*sigfun)(double, double* );

/* That is a first attempt for a data structure containing all
   the necessary information for computing a smear matrix. Right
   now it is sufficient for Type A matrices but lacks support
   for type B matrices */

#define GLB_TYPE_A 1
#define GLB_TYPE_B 2
#define GLB_TYPE_C 2 /* That is correct since there are now only two types */

typedef struct {
  int type;
  int numofbins;
  int simbins;
  double e_min; 
  double e_max;
  double e_sim_min;
  double e_sim_max;

    

  /* Pointer to the paramters for sig_f */
  double *sigma;
  int num_of_params;
  /* Here is the pointer to sigma(E,params) */
  sigfun sig_f;
  
  /* This is the binsize array -- unused so far.
   * In any case it is probably more useful to store the bin sizes
   * in a simple array.
   */
  
  double *binsize;

  double *simbinsize;

  /* Lookup tables */
  double *bincenter;
  double *simbincenter;

  /* Options for type B matrices */
  glb_option_type *options;
} glb_smear;


/* Data structure containing information about systematics functions */
typedef struct glb_systematic
{
  glb_chi_function chi_func;   /* Pointer to the chi^2 function                      */
  int dim;                     /* Number of systematics parameters                   */
  char *name;                  /* Unique name of this chi^2 routine                  */
  void *user_data;             /* Arbitrary user-defined parameter for chi^2 routine */
  struct glb_systematic *next; /* Pointer to next entry in glb_sys_list              */
} glb_systematic;




/* PH 01/10/19 -- in support of in-experiment oscillation engines  *
 * Data structure containing information about oscillation engines */

typedef struct glb_osc_engine
{
  glb_probability_matrix_function matrix_function;
  glb_set_oscillation_parameters_function set_function;
  glb_get_oscillation_parameters_function get_function;
  int num_of_params;            /* Number of oscillation parameters                   */
  char *name;                  /* Unique name of this oscillation engine                  */
  struct glb_osc_engine *next; /* Pointer to next entry in
				  glb_osc_engine_list             */
  void *user_data;             /* Arbitrary user-defined parameter for chi^2 routine */
} glb_osc_engine;


/* Data structure containing information about one nuisance parameter
 * (as defined in an AEDL "sys< >" block */
typedef struct glb_nuisance
{
  char *name;                  /* Name of this nuisance parameter                    */
  int systype; /* this allows to have different types of errors like calibration */

  double error;                /* The uncertainty on this nuisance parameter         */
  double a;                    /* The current value of the nuisance parameter        */
  double *a_list;              /* List of a values for energy-dependent parameters   */
  double *energy_list, *error_list; /* Energies and associated uncertainties if this */
                               /* nuisance parameter is energy-dependent             */
  int n_energies;              /* Number of entries in energy_list and error_list    */
  int ref_count;               /* Number of experiments using this nuisance param    */
} glb_nuisance;



/** This structure contains a large part of the information for defining
 * an experiment.
 * It has fixed size since it contains only pointers to dynamic objects.
 * You have to  make sure that the memory for dynamic objects is allocated
 * (and freed !) elsewhere !
 * Since this structure is large and will grow, care should be taken to pass
 * only pointers to it (esp. in functions which are called in each step or 
 * cycle).
 * The plan is, that it should at some point contain really all information
 * for describing an experiment like fluxes etc. ( or at least pointers to
 * the place where those things are)
 */

struct glb_experiment {

  /* Version string */
  char *version;

  /* A string containing citation information for this experiment */
  char *citation;

  /* Name of AEDL file on which this experiment is based */
  char *filename;

  /* A pointer to the primary detector in the case of multi-detector experiment
   * (defined using the #DETECTOR# diretive in AEDL) */
  struct glb_experiment *parent;
  struct glb_experiment *children[GLB_MAX_EXP];
  int n_children;
    /* The number of pointers that exist to this experiment, either from
     * glb_experiment_list, or from other experiments' parent pointers.
     * Shall be incremented/decremented by using glbIncrExpRefCount and
     * glbDecrExpRefCount, when ref_count reaches zero, the experimental
     * data structure will be destroyed. */

  /* This contains the parsing meta-information like names of rules etc. */
  glb_naming *names;
  
  /** The beam spectrum is loaded with glb_flux_loader(file_name, flux_ident)
   */

  /** That is a new way to define and load fluxes */
  glb_flux *fluxes[GLB_MAX_FLUXES];
  int num_of_fluxes;

  /** That is a new way to define and load x-sections */
  glb_xsec *xsecs[GLB_MAX_XSECS];
  int num_of_xsecs;

 
  /* Bin widths */
  double *binsize;
  double *simbinsize; 

  /** set the way a profile is computed or the baseline is changed */
  int density_profile_type;

  /** baseline is the baseline. It is a positive number. The unit
  * is km. For most cases it should not exceed 2*EARTHRADIUS.
  */
  double baseline;
  
  /** emin is the lower bound on the energy for the events used in the 
  * analysis.
  * It is a positive number and the unit is GeV.
  */ 
  double emin;

  /** emax is the upper bound on the energy for the  events in the analysis.
  * It is a positive number and the unit is GeV. It has to be larger than
  * emin.
  */
  double emax;

  /** Number of bins which divide the range from emin to emax for the analysis.
  * It is larger than zero.
  */
  int numofbins;


  /** Target mass in units of kt (=1000 tons). It is positive.
   */
  double targetmass;

  
  /** Number of channels. It is a number between 1 and GLB_MAX_CHANNELS.
   */
  int numofchannels;

  /** Strange construct. Holds the channel definition ( 5 int`s)
  * describing: beam polarity, initial state, final state, CP sign
  * cross section. The length was 5*numofchannnels. Not very nice!
  * will be changed for two reasons:
  * - it is ugly and the channel definition will be enlarged at least
  * by one number 
  * - the energy resolution function identifier.
  * This new defintion has now get numofchannels*sizeof(int) of memory.
  *
  * malloc!
  */
  int* listofchannels[6];

  /** Number of rules. Ranges from 1 to GLB_MAX_RULES.
   */
  int numofrules;

  /** Tells for each rule how many channels are added to form the signal.
  * Ranges from 1 to GLB_MAX_RULES.
  */
  int lengthofrules[GLB_MAX_RULES];

  /** Contains the pre-factor for each channel and rule which is applied to
  * each channel in computing the signal. Has length 
  * lengthofrules*sizeof(double). malloc!
  */
  double* rulescoeff[GLB_MAX_RULES];

  /** Contains the identifying number for each channel used in this rule.
  * Has got length lengthofrules*sizeof(int). Range is 0 to lengthofrule-1.
  * malloc!
  */
  int* rulechannellist[GLB_MAX_RULES];

  /** Tells for each rule how many channels are added to form the background.
  * Ranges from 1 to GLB_MAX_RULES.
  */
  int lengthofbgrules[GLB_MAX_RULES];

  /** Contains the pre-factor for each channel and rule which is applied to
  * each channel in computing the background. Has length 
  * lengthofrules*sizeof(double). malloc!
  */
  double* bgrulescoeff[GLB_MAX_RULES];

  /** Contains the identifying number for each channel used in this rule.
  * Has got length lengthofrules*sizeof(int). Range is 0 to lengthofrule-1.
  * malloc!
  */
  int* bgrulechannellist[GLB_MAX_RULES];

  /** Number of smearing data types stored */
  int num_of_sm;

  /** Meta information for computing the smear matrix for each rule */
  glb_smear *smear_data[GLB_MAX_SMEAR];

  /** Thus holds the pointer to the place where the energy resolution
  * function is stored. The energy resolution function is stored as
  * a matrix where only the non-zero elements are used. This and the
  * fact that there is a numer of bins in the analysis which can be
  * different from the bins used in the calculation makes the game
  * a little tricky. Be careful! smear contains for each different
  * energy resolution numofbins pointers which point each to vector
  * of length simbins and at this point in the memory you will find
  * a double.
  *
  * \todo This is now only the user-defined mode. Need implemenation
  * for computing this matrix and the other quantities. Also
  * the meta information for this process has to go in here.
  */
  double** smear[GLB_MAX_SMEAR];

  /** For each analysis bin (numofbins) exists a lower index in the
  * range given by simbins for which smear contains non-zero entries
  * (remember only those are stored). Thus lowrange is positive (including
  * zero) and smaller than simbins and smaller than uprange (see below)
  * \todo This is now only the user-defined mode. Need implemenation
  * for computing this matrix and the other quantities. Also
  * the meta information for this process has to go in here.
  */
  int* lowrange[GLB_MAX_SMEAR];

  /** For each analysis bin (numofbins) exists a upper index in the
  * range given by simbins for which smear contains non-zero entries
  * (remember only those are stored). Thus uprange is positive and 
  * smaller than simbins and bigger than uprange.
  * \todo This is now only the user-defined mode. Need implemenation
  * for computing this matrix and the other quantities. Also
  * the meta information for this process has to go in here.
  */
  int* uprange[GLB_MAX_SMEAR];
  
  /** simtresh has the same function for the computation of events than
  * emin for the analysis. simtresh has to be positive and smaller than 
  * emin. The unit is GeV. 
  */ 
  double simtresh;

  /** simbeam has the same function for the computation of events than
  * emax for the analysis. simtresh has to be positive and bigger than 
  * emax and bigger than simtresh. The unit is GeV.
  */ 
  double simbeam;

  /** simbins is the number of bins dividing (equi-distant) 
  * the range from simtresh to
  * simbeam. It is a positive number.
  * \todo This is now only the user-defined mode. Need implemenation
  * for computing this matrix and the other quantities. Also
  * the meta information for this process has to go in here.
  */
  int simbins; 

  /** Changed.
  * filter_state is either "GLB_ON" or "GLB_OFF".
  */
  int filter_state;

  /** Positive number.
   */
  double filter_value;

  /** The number of layers for the matter profile. It is a number greater
  * than 0. Default is 1.
  */
  int psteps;

  /** glb_List of the thickness of each matter layer. Unit is km and positive. 
  * Has length
  * psteps. The sum of all length in lengthtab should equal baseline.
  * Here a mechanism for ensuring this relation has to be implemented.
  * 
  */
  double* lengthtab;

  /** glb_List of the matter densities for each layer. It is a positive number
  * and the unit is g cm^-3. Has length psteps.
  * 
  */
  double* densitytab;

  /** A buffer needed in the computation of the matter profile uncertainty.
  * Does not appear in any config-file, but has to be malloced on
  * initialization. Has length psteps. Determined internally.
  */
  double* densitybuffer;

  /** New.
  * Contains a number for each simbin which is multiplied with each
  * channel before doing the smearing. Has length simbins. Positive
  * including zero.
  */
  double* user_pre_smearing_channel[GLB_MAX_CHANNELS];

  /** New.
  * Contains a number for each bin which is multiplied with each
  * channel after doing the smearing. Has length bins. Positive
  * including zero.
  */
  double* user_post_smearing_channel[GLB_MAX_CHANNELS];

  /** New.
  * Contains a number for each bin which is added to each
  * bin before doing the smearing. Has length simbins. Positive
  * including zero.
  */
  double* user_pre_smearing_background[GLB_MAX_CHANNELS];

  /** New.
  * Contains a number for each bin which is added to each
  * bin after doing the smearing. Has length bins. Positive
  * including zero.
  */
  double* user_post_smearing_background[GLB_MAX_CHANNELS];

  /** Energy range for each rule in which the events are used to compute
  * the chi^2. Lower limit and upper limit.
  */
  double energy_window[GLB_MAX_RULES][2];    /* Energy range */
  int energy_window_bins[GLB_MAX_RULES][2];  /* Bin range    */

  /** Has length numofbins */
  double* energy_tab;


  /** Now comes a bunch of pointers which finally are vectors containing
  * the different parts of event vectors needed during computation.
  * All mallocing has to be done at intialization of a given experiment!
  */
  double *chrb_0[GLB_MAX_CHANNELS], *chrb_1[GLB_MAX_CHANNELS]; /* True/fitted pre-sm. rates by ch  */
  double *chra_0[GLB_MAX_CHANNELS], *chra_1[GLB_MAX_CHANNELS]; /* True/fitted post-sm. rates by ch */
  double *chr_template[GLB_MAX_CHANNELS];   /* Products of fluxes, cross sections, and prefactors by ch */
  double* SignalRates[GLB_MAX_RULES];    /* "True" signal event rates for all rules */
  double* BackgroundRates[GLB_MAX_RULES];/* "True" background event rates for all rules */
  double* rates0[GLB_MAX_RULES];         /* "True" event rates for all rules */
  double* rates1[GLB_MAX_RULES];         /* Fitted signal rates for all rules */
  double* rates1BG[GLB_MAX_RULES];       /* Fitted background rates for all rules */

  /** Systematics functions and on/off states */
  int sys_on_off[GLB_MAX_RULES];         /* Systematics switch (GLB_ON/GLB_OFF)     */
  glb_systematic *sys_on[GLB_MAX_RULES]; /* The chi^2 functions of the rules        */
  glb_systematic *sys_off[GLB_MAX_RULES]; 
  char *sys_on_strings[GLB_MAX_RULES];   /* The errordim strings from the AEDL file */
  char *sys_off_strings[GLB_MAX_RULES];

  /** Systematical error and starting values (GLoBES 3.0 scheme) */
  double *sys_on_errors[GLB_MAX_RULES];     /* Systematical errors */
  double *sys_on_startvals[GLB_MAX_RULES];  /* Starting values for systematics minimization */
  double *sys_off_errors[GLB_MAX_RULES];    /* Systematical errors */
  double *sys_off_startvals[GLB_MAX_RULES]; /* Starting values for systematics minimization */

  /** Systematical error and starting values (GLoBES 2.0.11 scheme) */
  double signal_errors[2][GLB_MAX_RULES];   /* Signal norm/energy errors for old chi^2 functions */
  double signal_startvals[2][GLB_MAX_RULES];/* Signal starting values for old chi^2 functions */
  double bg_errors[2][GLB_MAX_RULES];       /* Background norm/energy errors for old chi^2 functions */
  double bg_startvals[2][GLB_MAX_RULES];    /* BG starting values for old chi^2 functions */
  double bg_centers[2][GLB_MAX_RULES];      /* BG central values for old chi^2 functions */

  /** Nuisance parameters for global/multi-experiment systematics */
  int n_nuisance;                           /* Number of nuisance parameters */
  glb_nuisance *nuisance_params[GLB_MAX_NUISANCE]; /* Pointers to nuisance params */

  int sys_on_n_nuis_sig[GLB_MAX_RULES][GLB_MAX_CHANNELS]; /* # of nuisance params applied to */
  int sys_on_n_nuis_bg[GLB_MAX_RULES][GLB_MAX_CHANNELS];  /* the channels making up a rule   */ 
  int sys_off_n_nuis_sig[GLB_MAX_RULES][GLB_MAX_CHANNELS];
  int sys_off_n_nuis_bg[GLB_MAX_RULES][GLB_MAX_CHANNELS];
  int *sys_on_multiex_errors_sig[GLB_MAX_RULES][GLB_MAX_CHANNELS]; /* Pointers to the nuisance */
  int *sys_on_multiex_errors_bg[GLB_MAX_RULES][GLB_MAX_CHANNELS];  /* parameters relevant to   */
  int *sys_off_multiex_errors_sig[GLB_MAX_RULES][GLB_MAX_CHANNELS];/* the channels making up   */
  int *sys_off_multiex_errors_bg[GLB_MAX_RULES][GLB_MAX_CHANNELS]; /* the rules                */


  /* added in version 3.3.18 */
  double * data[GLB_MAX_RULES]; /* contains the data if we try to fit a real
 			experiment */

  int data_on_off[GLB_MAX_RULES]; /* use data in data[] or use the simulated data
 			  in rates0[] */


  /* PH 01/10/19 this will be now all collected in one element
   * Probability engine for this experiment (if different from default) 
   *
   * glb_probability_matrix_function probability_matrix;
   * glb_set_oscillation_parameters_function set_oscillation_parameters;
   *  glb_get_oscillation_parameters_function get_oscillation_parameters;
   * void *probability_user_data;
   */

  glb_osc_engine osc_engine;
  
};

#endif /* GLB_TYPES_H 1 */
