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




/* --------------------------------------------------------
 * --- Comment on numerical accuracy ----------------------
 * --------------------------------------------------------
 * The whole code uses double precision, whatever this is
 * on Your machine. The following information was derived
 * on 32-bit Pentium III processor.
 * Even though double precision should give relative errors of
 * order <1E-12, functions involving minimization like ChiTheta
 * may be off in the order of 1E-8 (absolute). This may cause
 * ChiTheta and SingleChiTheta to yield results differing by 1E-8.
 * The same applies for all other ChiglbXSection and SingleChiglbXSection functions.
 * The crucial lines are:
 * erg2 = erg2
 *    + glb_prior(x[1],start[2],inp_errs[2])
 *    + glb_prior(x[3],start[4],inp_errs[4])
 *    + glb_prior(x[4],start[5],inp_errs[5])
 *    + glb_prior(x[5],(glb_experiment_list[glb_single_experiment_number]).density_center,
 *      (glb_experiment_list[glb_single_experiment_number]).density_error);
 * in the chi_xxx functions. This looks of course different
 * in the MDchi_xxx functions. Basically its a matter of how
 * and in which order the truncations are performed. This may
 * also change under compiler optimization. In any case the errors
 * should be very small (i.e. <<1E-6) and should not have any impact
 * on the physics results. 
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <setjmp.h>
#include <globes/globes.h>

#include "glb_probability.h"
#include "glb_fluxes.h"
#include "glb_rate_engine.h"
#include "glb_min_sup.h"
#include "glb_types.h"
#include "glb_multiex.h"
#include "glb_error.h"
#include "glb_wrapper.h"

/* global variabels */
int glb_single_experiment_number=0;



/* Warning-- fiddling around with these in general may
 * influence the stability of the minimization process !
 *
 * The latest Recator (D-CHOOZ as in hep-ph/0403068) setups 
 * have TOLSYS 1.0E-12, whereas all the other setups have TOLSYS 1.0E-5
 * and TOLSYS 1.0E-8 may be faster by 20% for Nufact (improves convergence
 * of the second layer of minimization)
 */

#define TOLSYS 1.0E-8
#define TOLOSC 1.0E-5


/* defines wether the program will run 
 * with (1) or w/o mathematica (0)
 */

#define GLB_MATHLINK 1

/* Defines wether the thing runs on UNIX (0) or not (1) */

#define UNIX 0

#define FLOAT double

//----------------------

/* this is for making the Chi... functions abortable */
static int okay_flag = 0;
/* this needed by the wrappers in order not to return garbage */
static jmp_buf env; 
/* this is need by setjmp() longjmp() */
//--------------------------------------------------------


static double *inp_errs = NULL;
static double *start = NULL;
//static double inp_errs[GLB_OSCP+1];
//static double start[GLB_OSCP+1];
static int count=0;
//int glb_single_experiment_number=0;

//This serves for ChiNP
static double *fix_params = NULL;
static int *para_tab = NULL;
static int *index_tab = NULL;
static int n_free;
static int n_fix;

//This serves for SingleChiNP
static double *s_fix_params = NULL;
static int *s_para_tab = NULL;
static int *s_index_tab = NULL;
static int s_n_free;
static int s_n_fix;

//-----------------------


//-------------------------------------------------------------
//------------------    MLAbort handling     ------------------
//-------------------------------------------------------------
// the whole stuff is being switched of if GLB_MATHLINK is set to 0.

// this function handels a caught abort message
// and jumps to the point where setjmp(env) has
// been called. usually from one of the ChiTheta(), 
// MDChiTheta() or Chi() or SingleChi(). this 
// function then will return a garbage value. this
// fact is comunicated to the wrappers with the 
// okay_flag. if they encounter an okay_flag=1 they
// will send the Abort[] mathematica comand and then
// return.

static void ml_abort_handler()
{
  fprintf(stderr,"Mathlink programe catched abort!\n");
  longjmp(env,1);
}

// this function checks wether a abort message is sent by the
// math-kernel, if so it calls the ml_abort_handler().
// otherwise it calls MLCallYieldFunction() in order to give
// the kernel a possibility to send its message (this should
// be needed only in Windows and on MacOS)

#ifdef MLabort
static void ml_abort_check(int flag)
{

  if(flag==0) return;
  if(!MLAbort)
    {
      if(UNIX==1) MLCallYieldFunction(MLYieldFunction(stdlink),
				      stdlink,(MLYieldParameters)0);
    }
  if(MLAbort)
    {
      ml_abort_handler();
    }
}
#else
static void  ml_abort_check(int flag)
{
  /* To silence gcc -Wall */
  int i;
  i=flag;
  return;
}

#endif


//--------------------------------------------------------------------------
//-------------------- Initialization of arrays ----------------------------
//--------------------------------------------------------------------------
int glb_init_minimizer()
{
  int i;
  
  glb_free_minimizer();
  
  /* General data structures */
  inp_errs = (double *) glb_malloc((glbGetNumOfOscParams()+1) * sizeof(inp_errs[0]));
  start = (double *) glb_malloc((glbGetNumOfOscParams()+1) * sizeof(start[0]));
  
  /* Data structures for internal_glbChiNP */
  fix_params = (double *) glb_malloc((glbGetNumOfOscParams()+32) * sizeof(fix_params[0]));
  para_tab = (int *) glb_malloc((glbGetNumOfOscParams()+32) * sizeof(para_tab[0]));
  index_tab = (int *) glb_malloc((glbGetNumOfOscParams()+32) * sizeof(index_tab[0]));

  /* Data structures for internal_glbSingleChiNP */
  s_fix_params = (double *) glb_malloc((glbGetNumOfOscParams()+1) * sizeof(s_fix_params[0]));
  s_para_tab = (int *) glb_malloc((glbGetNumOfOscParams()+1) * sizeof(s_para_tab[0]));
  s_index_tab = (int *) glb_malloc((glbGetNumOfOscParams()+1) * sizeof(s_index_tab[0]));

  /* Select default projection */
  for (i=0; i < glbGetNumOfOscParams()+32; i++)
    para_tab[i] = GLB_FREE;
  for (i=0; i < glbGetNumOfOscParams()+1; i++)
    s_para_tab[i] = GLB_FREE;
}


int glb_free_minimizer()
{
  if (inp_errs != NULL)     { glb_free(inp_errs);     inp_errs = NULL;     }
  if (start != NULL)        { glb_free(start);        start = NULL;        }
  if (fix_params != NULL)   { glb_free(fix_params);   fix_params = NULL;   }
  if (para_tab != NULL)     { glb_free(para_tab);     para_tab = NULL;     }
  if (index_tab != NULL)    { glb_free(index_tab);    index_tab = NULL;    }
  if (s_fix_params != NULL) { glb_free(s_fix_params); s_fix_params = NULL; }
  if (s_para_tab != NULL)   { glb_free(s_para_tab);   s_para_tab = NULL;   }
  if (s_index_tab != NULL)  { glb_free(s_index_tab);  s_index_tab = NULL;  }

  return 0;
}


//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//------------------------ Systematics -------------------------------------
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------


/* Callback function for the systematics minimizer;        */
/* Calls the actual (possibly user-defined) chi^2 function */
/* with the appropriate parameters                         */
// JK - FIXME - This should be solved more elegantly
static glb_chi_function glb_current_chi_function;
static int glb_current_n_sys;
double glb_chi_callback(double *params)
{
  return glb_current_chi_function(glb_current_exp, glb_rule_number,
                                  &(params[1]), glb_current_n_sys);
}


// this an init function -
// setting up the input matrix for minimization
// ie. the set of starting directions
static void init_mat(double **m, int dim)
{
  int i;
  int j;
  for (i=1 ;i<= dim;i++)
    {
      for (j=1;j<= dim;j++)
	{
	  if (i==j)
	    {
	      m[i][j]=1.0;
	    }
	  else
	    {
	      m[i][j]=0.0;
	    }	
	}
    }
}

// probably will be replaced by a more object oriented approach

// those two are needed in order to write unique
// chi^2 functions
double *glb_sys_errors;
double *glb_sys_centers;

//FIXME Remove
//static double chi_dispatch(int error_dim, int i)
//{     
//  // finally
//  double erg2=0;
//  
//  glb_rule_number = i;      // Tell the systematics functions which rule we are in
//  if(error_dim==0) //chi^2 with 4 sys. parameters
//    {
//      erg2=chi_sys_wrap(&glb_chi_sys_w_bg,4,i);
//      
//    }
//  else if(error_dim==1) //chi^2 with 6 sys. parameters
//    {
//      erg2=chi_sys_wrap(&glb_chi_sys_w_bg2,6,i); 
//    }
//  else if(error_dim==2) // chi^2 without sys. parameters
//    {
//      erg2=chi_sys_wrap(&glb_chi_sys_w_bg,0,i); 
//      // it seems that the minimizer correctly interprest zero dimensions
//      // as the function value at the starting point  
//    }
//  else if(error_dim==3) // very special for JHF-HK
//    {
//      if(i<4) erg2=chi_sys_wrap(&glb_chi_sys_w_bg,4,i); 
//      else erg2=chi_sys_wrap(&glb_chi_sys_w_bgtot,4,i); 
//    }
//  else if(error_dim==4) // total chi^2 with 4 sys parameters
//    {
//      erg2=chi_sys_wrap(&glb_chi_sys_w_bgtot2,4,i); 
//    }
//  else if(error_dim==5) // obsolete
//    {
//      fprintf(stderr,"Warning: obsolete errordim\n");
//      erg2=0;
//    }
//  else if(error_dim==6) // obsolete
//    {
//      fprintf(stderr,"Warning: obsolete errordim\n");
//      erg2=0;
//    }
//  else if(error_dim==7) // free normalization, i.e. spectrum only
//    {
//      erg2=chi_sys_wrap(&glb_chi_spec,1,i); 
//    }
//  else if(error_dim==8) //total chi^2 w/o systematics
//    {
//      erg2=chi_sys_wrap(&glb_chi_sys_w_bgtot2,0,i); 
//    }
//  else if(error_dim==9) // chi^2 with 4 syst params and 
//    // true energy calibration
//    {
//      erg2=chi_sys_wrap(&glb_chi_sys_w_bg_calib,4,i); 
//    }
//  else if(error_dim==10) // total chi^2 with 4 syst params 
//    {
//      erg2=chi_sys_wrap(&glb_chi_sys_w_bgtot,4,i); 
//    }
//  else if(error_dim==20) // User-defined systematics 
//    {
//      erg2=glb_evaluate_chi(&sys_calc[i]); 
//    }
//  else if(error_dim==21) // total chi^2 with 4 syst params 
//    {
//      erg2=sys_calc[i].evalf(&sys_calc[i]); 
//    }
//  else
//    {
//      erg2=0;
//    }
//return erg2;
//}

// this is the same as ChiSO()  but it allows access to a single rule !

static double ChiS0_Rule(int rule)
{
  double **mat; // contains the direction set
  double *sp;   // stores the coordinates of the minimum
  double res;   // stores the minimum value
  int it=0;     // counts the number of iterations
  glb_systematic *sys;
  struct glb_experiment *exp = glb_experiment_list[glb_current_exp];
  int i;
 
  res = 0.0;
  glb_rule_number = rule;
  if (exp->sys_on_off[rule] == GLB_ON)
    sys = exp->sys_on[rule];
  else
    sys = exp->sys_off[rule];
    
  mat  = glb_alloc_mat(1, sys->dim, 1, sys->dim);
  sp   = glb_alloc_vec(1, 6); // that is the maximal length //FIXME
  init_mat(mat, sys->dim);
  sp[1] = 1;
  sp[2] = 0;
  sp[3] = glb_bg_norm_center[rule];
  sp[4] = glb_bg_tilt_center[rule];
  sp[5] = glb_tre_null_center[rule];
  sp[6] = glb_tre_tilt_center[rule];
  glb_current_chi_function = sys->chi_func;
  glb_current_n_sys        = sys->dim;
  if (glb_powell(sp, mat, sys->dim, TOLSYS, &it, &res, &glb_chi_callback) != 0)
  {
    glb_warning("Systematics minimization failed.");
    return -res;
  }
  glb_free_vec(sp, 1, 6);
  glb_free_mat(mat, 1, sys->dim, 1, sys->dim);

  glb_rule_number = 0;
  return res;
}

static double ChiS0()
{
  int i;
  double res = 0.0;
  
//  ml_abort_check(GLB_MATHLINK);    /* RELICT??? WW */

  for (i=0; i < glb_num_of_rules; i++) 
    res += ChiS0_Rule(i);

  return res;
}

//---------------------------------------------------------------
//---------- MES begins here ------------------------------------
//---------------------------------------------------------------

// redfinition of ChiS

static double ChiS()
{
  double erg;
  int i;
  erg=0;
  for (i=0;i<glb_num_of_exps;i++)
    {
      glbSetExperiment(glb_experiment_list[i]);
      erg += ChiS0();
    }
  return erg;
}

//redefinition of Chi

static double Chi(double x[])
{
  int i;  
  double erg;

  glb_params p = glbAllocParams();
  for (i=0; i < glbGetNumOfOscParams(); i++)
    p->osc->osc_params[i] = x[i];
  glb_hook_set_oscillation_parameters(p);
  glbFreeParams(p);
    
//FIXME Remove
//  double nsp[glbGetNumOfOscParams()-6+1];
//  glb_set_c_vacuum_parameters(x[0],x[1],x[2],x[3]);
//  glb_set_c_squared_masses(0,x[4],x[5]);
//  if(glbGetNumOfOscParams()>6)
//  {
//	for(i=0;i<glbGetNumOfOscParams()-6;i++) nsp[i]=x[i+6];
//	glb_set_c_ns_params(nsp);
//  }
  
  for (i=0;i<glb_num_of_exps;i++)
    {
      glbSetExperiment(glb_experiment_list[i]);
      glb_set_profile_scaling(x[glbGetNumOfOscParams()+i],i);
      glb_set_new_rates();
    }
  if (setjmp(env)==1)
    {
      okay_flag=1;
      return erg;  
    }
  erg=ChiS();
  return erg;
}




// chi^2 with systematics for each Experiment

static double SingleChi(double x[glbGetNumOfOscParams()+1],int exp)
{
  int i;
  double erg;

  glb_params p = glbAllocParams();
  for (i=0; i < glbGetNumOfOscParams(); i++)
    p->osc->osc_params[i] = x[i];
  glb_hook_set_oscillation_parameters(p);
  glbFreeParams(p);

//FIXME remove
//  double nsp[glbGetNumOfOscParams()-6+1];
//  glb_set_c_vacuum_parameters(x[0], x[1],x[2],x[3]);
//  glb_set_c_squared_masses(0,x[4],x[5]);
//  if(glbGetNumOfOscParams()>6)
//  {
//	for(i=0;i<glbGetNumOfOscParams()-6;i++) nsp[i]=x[i+6];
//	glb_set_c_ns_params(nsp);
//  }

      glbSetExperiment(glb_experiment_list[exp]);
      glb_set_profile_scaling(x[glbGetNumOfOscParams()],exp);
      glb_set_new_rates();
     


  glbSetExperiment(glb_experiment_list[exp]);
  if (setjmp(env)==1) 
    {
      okay_flag=1;
      return erg;  
    }
  erg=ChiS0();
  return erg;
}

// chi^2 with systematics for each Experiment

static double SingleRuleChi(double x[glbGetNumOfOscParams()+1],int exp, int rule)
{
  int i;
  double erg;

  glb_params p = glbAllocParams();
  for (i=0; i < glbGetNumOfOscParams(); i++)
    p->osc->osc_params[i] = x[i];
  glb_hook_set_oscillation_parameters(p);
  glbFreeParams(p);

// FIXME Remove
//  double nsp[glbGetNumOfOscParams()-6+1];
//  glb_set_c_vacuum_parameters(x[0], x[1],x[2],x[3]);
//  glb_set_c_squared_masses(0,x[4],x[5]);
//  if(glbGetNumOfOscParams()>6)
//  {
//	for(i=0;i<glbGetNumOfOscParams()-6;i++) nsp[i]=x[i+6];
//	glb_set_c_ns_params(nsp);
//  }

      glbSetExperiment(glb_experiment_list[exp]);
      glb_set_profile_scaling(x[glbGetNumOfOscParams()],exp);
      glb_set_new_rates();
     

  glbSetExperiment(glb_experiment_list[exp]);
  erg=ChiS0_Rule(rule);
  return erg;
}

/* Wrappers for the API */

double glbChiSys(const glb_params in,int experiment, int rule)
{
  int i;
  double res,x[32+glbGetNumOfOscParams()];

  if(in==NULL)
    {
      glb_error("Failure in glbChiSys: Input pointer must be non-NULL");
      return -1;
    }
  for(i=0;i<glbGetNumOfOscParams();i++) x[i]=glbGetOscParams(in,i);
  if(experiment==GLB_ALL)
    {
      if(rule==GLB_ALL)
	{
	  for(i=0;i<glb_num_of_exps;i++) 
	    x[i+glbGetNumOfOscParams()]=glbGetDensityParams(in,i);
	  res=Chi(x);
	}
      else
	{
	  res=0;
	   for(i=0;i<glb_num_of_exps;i++)
	     {   
	       if(rule < glb_experiment_list[i]->numofrules)
		 {
		   x[glbGetNumOfOscParams()]=glbGetDensityParams(in,i);
		   res+=SingleRuleChi(x,i,rule);
		 }
	       else
		 {glb_error("Invalid rule number");return -1;}
	     }
	}
    }
  else
    {
      if(experiment >= glb_num_of_exps) 
	{
	  glb_error("Failure in glbChiSys: 2nd argument must be smaller than"
		    "glb_num_of_exps");
	  return -1;
	}
      x[glbGetNumOfOscParams()]=glbGetDensityParams(in,experiment);
      if(rule==GLB_ALL) res=SingleChi(x,experiment);
	else
	  {
	    if(rule >= glb_experiment_list[experiment]->numofrules)
	      {
		glb_error("Failure in glbChiSys: 3rd argument must be"
			  " smaller than numofrules");
		return -1;
	      }
	    res=SingleRuleChi(x,experiment,rule);
	  }
    }
  return res;
}


//-----------------------------------------------------------------
//------------------- END -----------------------------------------
//--------------- Systematics -------------------------------------
//-----------------------------------------------------------------


//--------------------------------------------------------------
//--------- Setting the glb_priors ---------------------------------
//--------------------------------------------------------------

int glb_set_input_errors(double a, double b, double c, double d, double e, double f)
{
 
  inp_errs[1]=a;
  inp_errs[2]=b;
  inp_errs[3]=c;
  inp_errs[4]=d;
  inp_errs[5]=e;
  inp_errs[6]=f;  /* Relict!!! Solar goes in Zero! */
  return 0;
}

int glb_set_starting_values(double a, double b, double c, double d, double e, double f)
{
  start[1]=a;
  start[2]=b;
  start[3]=c;
  start[4]=d;
  start[5]=e;
  start[6]=f;  /* Relict!!! Solar goes in Zero! */
  return 0;
}

int glb_set_ns_input_errors(double v[])
{ 
  int i;
  if(glbGetNumOfOscParams()>6)
	for(i=6;i<glbGetNumOfOscParams();i++)
		inp_errs[i]=v[i-6];
  return 0;
}

int glb_set_ns_starting_values(double v[])
{ 
  int i;
  if(glbGetNumOfOscParams()>6)
	for(i=6;i<glbGetNumOfOscParams();i++)
		start[i]=v[i-6];
  return 0;
}




// setting starting values for densities and errors on densitites
// for different experiments

void glbSetDensityPrior(double start, double error, int typ)
{
  glb_experiment_list[typ]->density_center=start;
  glb_experiment_list[typ]->density_error=error;
}

void glbSetDensityStartingValue(double start, int typ)
{
  glb_experiment_list[typ]->density_center=start;
}


void glbSetDensityInputError(double error, int typ)
{
  glb_experiment_list[typ]->density_error=error;
}


double* glb_return_input_errors()
{
  double* out;
  int i;
  i=glbGetNumOfOscParams()+glb_num_of_exps;
  out=(double*) glb_malloc(i*sizeof(double));
  for(i=0;i<glbGetNumOfOscParams();i++) out[i]=inp_errs[i];
  for(i=0;i<glb_num_of_exps;i++) out[i+glbGetNumOfOscParams()]=glb_experiment_list[i]->density_error;
  return out;
}

double* glb_return_input_values()
{
  double* out;
  int i;
  i=glbGetNumOfOscParams()+glb_num_of_exps;
  out=(double*) glb_malloc(i*sizeof(double));
  for(i=0;i<glbGetNumOfOscParams();i++) out[i]=start[i];
  for(i=0;i<glb_num_of_exps;i++) out[i+glbGetNumOfOscParams()]=glb_experiment_list[i]->density_center;
  return out;
}
  
 

//--------------------------------------------------------------
//--------- Setting the glb_priors for th12 --------------------
//--------------------------------------------------------------

int glb_set_solar_input_errors(double a)
{
  inp_errs[0]=a;
  return 0;
}

int glb_set_solar_starting_values(double a)
{
  start[0]=a;
  return 0;
}



//------------------------------------------------------------------
//---- Chi^2 with arbitrary number of free parameters --------------
//----------------------- 23.01.2004 -------------------------------
//------------------------------------------------------------------

static void SelectProjection(int *vec)
{
  int i,c,c1,c2;

  for(i=0;i<glbGetNumOfOscParams()+glb_num_of_exps;i++) para_tab[i]=vec[i];
  c=0;
  c1=0;
  c2=0;
  for(i=0;i<glbGetNumOfOscParams()+glb_num_of_exps;i++)
    {
      if(para_tab[i]==1) c++;
    }
  for(i=0;i<glbGetNumOfOscParams()+glb_num_of_exps;i++)
    {
      if(para_tab[i]==1) 
	{
	  index_tab[c1]=i;
	  c1++;
	}
      else if(para_tab[i]==0)
	{
	  index_tab[c+c2]=i;
	  c2++;
	}
      else
	{
	  glb_fatal("SelectProjection input error\n");
	} 
    }
  n_free=c;
  n_fix=c2;
  return;
}


static int CheckFree()
{
  int k;
  k=n_free;
  return k;
}

static int* CheckProjection()
{
  return &para_tab[0];
}


static double sglb_prior(double x, double center, double sigma)
{
  if(fabs(sigma-0)<1E-12) return 0;
  return (x-center)*(x-center)/sigma/sigma;
}  

// the pointer to the userdefined prior function
 double (*glb_user_defined_prior)(const glb_params);
 int (*glb_user_defined_starting_values)(const glb_params);
 int (*glb_user_defined_input_errors)(const glb_params);

static int my_default_sv(const glb_params in)
{
  return 0;
}

static int my_default_er(const glb_params in)
{
  return 0;
}

static double my_default_prior(const glb_params in)
{
  return 0;
}

// the user interface to register such a function ...
int glbRegisterPriorFunction(double (*prior)(const glb_params),
			     int (*starting)(const glb_params),
			     int (*error)(const glb_params))
{
  if(prior==NULL) 
    glb_user_defined_prior=my_default_prior;
  else
    glb_user_defined_prior=prior;

  if(starting==NULL) 
    glb_user_defined_starting_values=my_default_sv;
  else
     glb_user_defined_starting_values=starting;

  if(error==NULL) 
    glb_user_defined_input_errors=my_default_er;
  else
    glb_user_defined_input_errors=error;



  return 0;
}

// multi-experiment functions MDglbXSection
static double MD_chi_NP(double x[])
{
  glb_params prior_input;
  double erg2;
  double y[glbGetNumOfOscParams()+32];
  int i; 
//  double nsp[glbGetNumOfOscParams()-6+1];
  prior_input=glbAllocParams();
  count = count +1;
  for(i=0;i<n_free;i++) y[index_tab[i]]=x[i+1];
  // This basically is superflous, however it appears to be safer not
  // change a global (i.e to this file) parameter (fix_params) at this place
  for(i=n_free;i<n_free+n_fix;i++) y[index_tab[i]]=fix_params[index_tab[i]]; 

  glb_params p = glbAllocParams();
  for (i=0; i < glbGetNumOfOscParams(); i++)
    p->osc->osc_params[i] = y[i];
  glb_hook_set_oscillation_parameters(p);
  glbFreeParams(p);

// FIXME Remove
//  glb_set_c_vacuum_parameters(y[0],y[1],y[2],y[3]);
//  glb_set_c_squared_masses(0,y[4],y[5]);
//  if(glbGetNumOfOscParams()>6)
//  {
//	for(i=0;i<glbGetNumOfOscParams()-6;i++) nsp[i]=y[i+6];
//	glb_set_c_ns_params(nsp);
//  }
  
  for (i=0;i<glb_num_of_exps;i++)
    {
      glbSetExperiment(glb_experiment_list[i]);
      glb_set_profile_scaling(y[glbGetNumOfOscParams()+i],i);
      glb_set_new_rates();
    }
  
  erg2=ChiS();
  // adding  the user defined prior
  // shoufling the parameter vector y into an glb_params structure
  for (i=0;i<glbGetNumOfOscParams();i++) glbSetOscParams(prior_input,y[i],i);
  for (i=0;i<glb_num_of_exps;i++) glbSetDensityParams(prior_input,
						      y[i+glbGetNumOfOscParams()],i);
  glbSetIteration(prior_input,count);
  
  erg2 = erg2 + glb_user_defined_prior(prior_input); 

  // FIXME Remove
  // adding the glb_priors
  /*for(i=0;i<n_free;i++)
    {
      if(index_tab[i]<glbGetNumOfOscParams())
	{
	  erg2+=sglb_prior(y[index_tab[i]],
		      start[index_tab[i]],inp_errs[index_tab[i]]);  
	 
	}
      else
	{
	  erg2+=sglb_prior(y[index_tab[i]],
	  	     (glb_experiment_list[index_tab[i]-glbGetNumOfOscParams()])->density_center,
	       (glb_experiment_list[index_tab[i]-glbGetNumOfOscParams()])->density_error);
	}
	}
  */
  glbFreeParams(prior_input);
  return erg2;
}
 

// single-experiment functions ChiXXXSection

static void single_SelectProjection(int set)
{
  int i,c,c1,c2;

  for(i=0;i<glbGetNumOfOscParams();i++) s_para_tab[i]=para_tab[i];
  s_para_tab[glbGetNumOfOscParams()]=para_tab[glbGetNumOfOscParams()+set];
  c=0;
  c1=0;
  c2=0;
  for(i=0;i<glbGetNumOfOscParams()+1;i++)
    {
      if(s_para_tab[i]==1) c++;
    }
  for(i=0;i<glbGetNumOfOscParams()+1;i++)
    {
      if(s_para_tab[i]==1) 
	{
	  s_index_tab[c1]=i;
	  c1++;
	}
      else if(s_para_tab[i]==0)
	{
	  s_index_tab[c+c2]=i;
	  c2++;
	}
      else
	{
	  fprintf(stderr,"SelectProjection input error\n");
	} 
    }
 
  s_n_free=c;
  s_n_fix=c2;
  return;
}


static double chi_NP(double x[])
{
  glb_params prior_input;
  double erg2;
  double nsp[glbGetNumOfOscParams()-6+1];
  double y[glbGetNumOfOscParams()+1];
  int i;
  prior_input=glbAllocParams();
  count = count +1;
  for(i=0;i<s_n_free;i++) y[s_index_tab[i]]=x[i+1];
  // This basically is superflous, however it appears to be safer not
  // change a global (i.e to this file) parameter (fix_params) at this place
  for(i=s_n_free;i<s_n_free+s_n_fix;i++) y[s_index_tab[i]]
					 =s_fix_params[s_index_tab[i]];

  glb_params p = glbAllocParams();
  for (i=0; i < glbGetNumOfOscParams(); i++)
    p->osc->osc_params[i] = y[i];
  glb_hook_set_oscillation_parameters(p);
  glbFreeParams(p);

// FIXME Remove
//  glb_set_c_vacuum_parameters(y[0],y[1],y[2],y[3]);
//  glb_set_c_squared_masses(0,y[4],y[5]);
//  if(glbGetNumOfOscParams()>6)
//  {
//	for(i=0;i<glbGetNumOfOscParams()-6;i++) nsp[i]=y[i+6];
//	glb_set_c_ns_params(nsp);
//  }
  
  glbSetExperiment(glb_experiment_list[glb_single_experiment_number]);
  glb_set_profile_scaling(y[glbGetNumOfOscParams()],glb_single_experiment_number);
  glb_set_new_rates();
    
  
  erg2=ChiS0();
  
  // adding  the user defined prior
  // shoufling the parameter vector y into an glb_params structure
  for (i=0;i<glbGetNumOfOscParams();i++) glbSetOscParams(prior_input,y[i],i);
  glbSetDensityParams(prior_input,y[glbGetNumOfOscParams()],glb_single_experiment_number);
  glbSetIteration(prior_input,count);

  erg2 = erg2 + glb_user_defined_prior(prior_input); 

  //FIXME Remove  
  // adding the glb_priors
  /*
   for(i=0;i<s_n_free;i++)
    {
      if(s_index_tab[i]<glbGetNumOfOscParams())
	{
	  erg2+=sglb_prior(y[s_index_tab[i]],
		      start[s_index_tab[i]],inp_errs[s_index_tab[i]]);  
	 
	}
      else
	{
	  erg2+=sglb_prior(y[glbGetNumOfOscParams()],
	  	     (glb_experiment_list[glb_single_experiment_number])->density_center,
	       (glb_experiment_list[glb_single_experiment_number])->density_error);
	}
	}
  */
  glbFreeParams(prior_input);
  return erg2;
}
 




/* This implemenst the API for the ChiXXX functions */


 
static double internal_glbSingleChiNP(const glb_params in, glb_params out, 
			       int exp)
{
  double *sp2;
  double **mat2;
  double er1;
  glb_projection fbuf,fnew; 

  double x[32+glbGetNumOfOscParams()];
  int it;
  int i;
  int dim;
  fbuf=glbAllocProjection();
  fnew=glbAllocProjection();
  if(exp >=  glb_num_of_exps)
    {
      glb_error("Failure in internal_glbSingleChiNP: exp must be smaller than"
		" glb_num_of_exps");
      return -1;
    }

  

  glb_single_experiment_number=exp;
  
  //creating memory 
  single_SelectProjection(exp);

  dim=s_n_free;
  /* declaring temporariliy all densities of all other experiments as fixed */
  glbGetProjection(fbuf);
  fnew=glbCopyProjection(fbuf,fnew);
  fnew=glbSetDensityProjectionFlag(fnew,GLB_FIXED,GLB_ALL);
  fnew=glbSetDensityProjectionFlag(fnew,glbGetDensityProjectionFlag(fbuf,exp)
			      ,exp);  
  glbSetProjection(fnew);  
  /* - finis - */
  mat2=glb_alloc_mat(1,dim,1,dim);
  sp2=glb_alloc_vec(1,dim);
  init_mat(mat2,dim);
  //initializing various things
  count=0;
 
  if(out!=NULL) {
    out=glbCopyParams(in,out);
    if(out==NULL) 
      {
	glb_error("Failure while copying input of glbChiNP");
	return -1;
      }
  }
  
  for(i=0;i<glbGetNumOfOscParams();i++) x[i]=glbGetOscParams(in,i);
  x[glbGetNumOfOscParams()]=glbGetDensityParams(in,exp);

  for(i=0;i<glbGetNumOfOscParams()+1;i++) s_fix_params[i]=x[i];
  for(i=0;i<s_n_free;i++) sp2[i+1]=x[s_index_tab[i]];
  if (setjmp(env)==1) 
    {
      okay_flag=1;
      return er1;  
    }
  if(glb_powell2(sp2,mat2,dim,TOLOSC,&it,&er1,chi_NP)!=0) count=-count;
  if(out!=NULL)
    {
      for(i=0;i<s_n_free;i++) 
	{
	  if(s_index_tab[i]<glbGetNumOfOscParams())
	    glbSetOscParams(out,sp2[i+1],s_index_tab[i]);
	  else
	    glbSetDensityParams(out,sp2[i+1],exp);
	}
      out=glbSetIteration(out,count);
    }
  
  glb_free_vec(sp2,1,dim);
  glb_free_mat(mat2,1,dim,1,dim);
  glbSetProjection(fbuf);
  glbFreeProjection(fnew);
  glbFreeProjection(fbuf);
  return er1;
}

 

static double internal_glbChiNP(const glb_params in, glb_params out)
{
  double *sp2;
  double **mat2;
  double er1;
 

  double x[32+glbGetNumOfOscParams()];
  int it;
  int i;
  int dim;
  //creating memory 
 
  
  dim=n_free;
  
  mat2=glb_alloc_mat(1,dim,1,dim);
  sp2=glb_alloc_vec(1,dim);
  init_mat(mat2,dim);
  //initializing various things
  count=0;
  

 if(out!=NULL) {
    out=glbCopyParams(in,out);
    if(out==NULL) 
      {
	glb_error("Failure while copying input of glbChiNP");
	return -1;
      }
  }
  
  for(i=0;i<glbGetNumOfOscParams();i++) x[i]=glbGetOscParams(in,i);
  for(i=0;i<glb_num_of_exps;i++) x[i+glbGetNumOfOscParams()]=glbGetDensityParams(in,i);

  for(i=0;i<glbGetNumOfOscParams()+glb_num_of_exps;i++) fix_params[i]=x[i];
  for(i=0;i<n_free;i++) sp2[i+1]=x[index_tab[i]];
  if (setjmp(env)==1) 
    {
      okay_flag=1;
      return er1;  
    }
  if(glb_powell2(sp2,mat2,dim,TOLOSC,&it,&er1,MD_chi_NP)!=0) count=-count;
  if(out!=NULL)
    {
      for(i=0;i<n_free;i++) 
	{
	  if(index_tab[i]<glbGetNumOfOscParams())
	    glbSetOscParams(out,sp2[i+1],index_tab[i]);
	  else
	    glbSetDensityParams(out,sp2[i+1],index_tab[i]-glbGetNumOfOscParams());
	}
      out=glbSetIteration(out,count);
    }
  glb_free_vec(sp2,1,dim);
  glb_free_mat(mat2,1,dim,1,dim);  
  return er1;
}

double glbChiNP(const glb_params in, glb_params out, int exp)
{
  double res;
 
  if(in==NULL)
    {
      glb_error("Failure in glbChiNP: Input pointer must be non-NULL");
      return -1;
    }
  
  if(exp==GLB_ALL) res=internal_glbChiNP(in,out);
  else res=internal_glbSingleChiNP(in,out,exp);
  return res; 
}


/* Convenience wrappers for glbChiNP 
 *
 * First the 1-d wrappers
 */

static double glbChi1P(const glb_params in, 
	       	      glb_params out, int exp, int n)
{
  double res;
  int i,*b;
  int swit[32+glbGetNumOfOscParams()]; 
  int buff[32+glbGetNumOfOscParams()];
  b=CheckProjection();
  for(i=0;i<glbGetNumOfOscParams()+glb_num_of_exps;i++) buff[i]=b[i];
  for(i=0;i<glbGetNumOfOscParams()+glb_num_of_exps;i++) swit[i]=GLB_FREE;
  swit[n]=GLB_FIXED;
  SelectProjection(swit);
  if(in==NULL)
    {
      glb_error("Failure in glbChi1P: Input pointer must be non-NULL");
      return -1;
    }
  
  if(exp==GLB_ALL) res=internal_glbChiNP(in,out);
  else res=internal_glbSingleChiNP(in,out,exp);
  SelectProjection(buff);
  return res; 
}

/* Projection of chi^2 onto theta-12 */
double glbChiTheta12(const glb_params in,glb_params out, int exp)
{
  double res;
  res=glbChi1P(in,out,exp,GLB_THETA_12);
  return res;
}

/* Projection of chi^2 onto theta-13 */
double glbChiTheta13(const glb_params in,glb_params out, int exp)
{
  double res;
  res=glbChi1P(in,out,exp,GLB_THETA_13);
  return res;
}

/* Projection of chi^2 onto theta-23 */
double glbChiTheta23(const glb_params in,glb_params out, int exp)
{
  double res;
  res=glbChi1P(in,out,exp,GLB_THETA_23);
  return res;
}

/* Projection of chi^2 onto delta_CP */
double glbChiDelta(const glb_params in,glb_params out, int exp)
{
  double res;
  res=glbChi1P(in,out,exp,GLB_DELTA_CP);
  return res;
}

/* Projection of chi^2 onto Delta m_{21}^2 */
double glbChiDm21(const glb_params in,glb_params out, int exp)
{
  double res;
  res=glbChi1P(in,out,exp,GLB_DM_21);
  return res;
}

/* Projection of chi^2 onto Delta m_{31}^2 */
double glbChiDm31(const glb_params in,glb_params out, int exp)
{
  double res;
  res=glbChi1P(in,out,exp,GLB_DM_31);
  return res;
}

/* 2-dim Projection of chi^2 onto the theta-13 - delta_CP plane */
double glbChiTheta13Delta(const glb_params in,glb_params out, int exp)
{
  double res;
  int i,*b;
  int swit[32+glbGetNumOfOscParams()]; 
  int buff[32+glbGetNumOfOscParams()];
  b=CheckProjection();
  for(i=0;i<glbGetNumOfOscParams()+glb_num_of_exps;i++) buff[i]=b[i];
  for(i=0;i<glbGetNumOfOscParams()+glb_num_of_exps;i++) swit[i]=GLB_FREE;
  swit[GLB_THETA_13]=GLB_FIXED;
  swit[GLB_DELTA_CP]=GLB_FIXED;
  SelectProjection(swit);
  if(in==NULL)
    {
      glb_error("Failure in glbChiNP: Input pointer must be non-NULL");
      return -1;
    }
  
  if(exp==GLB_ALL) res=internal_glbChiNP(in,out);
  else res=internal_glbSingleChiNP(in,out,exp);
  SelectProjection(b);
  return res;
}

/* Wrapper for ChiAll */
double glbChiAll(const glb_params in,glb_params out, int exp)
{
  double res;
  int i,*b;
  int swit[32+glbGetNumOfOscParams()]; 
  int buff[32+glbGetNumOfOscParams()];
  b=CheckProjection();
  for(i=0;i<glbGetNumOfOscParams()+glb_num_of_exps;i++) buff[i]=b[i];
  for(i=0;i<glbGetNumOfOscParams()+glb_num_of_exps;i++) swit[i]=GLB_FREE;
  SelectProjection(swit);
  if(in==NULL)
    {
      glb_error("Failure in glbChiNP: Input pointer must be non-NULL");
      return -1;
    }
  
  if(exp==GLB_ALL) res=internal_glbChiNP(in,out);
  else res=internal_glbSingleChiNP(in,out,exp);
  SelectProjection(b);
  return res;
}


/* New functions dealing with setting the projection flags */


int 
glbSetProjection(const glb_projection in)
{
  int i;
  int *buf;
  if(in==NULL) return -1;
  buf=glb_malloc(sizeof(int)*(glbGetNumOfOscParams()+glb_num_of_exps));

  for(i=0;i<glbGetNumOfOscParams();i++) buf[i]=glbGetProjectionFlag(in,i);
  for(i=0;i<glb_num_of_exps;i++) buf[i+glbGetNumOfOscParams()]=
				   glbGetDensityProjectionFlag(in,i);
  SelectProjection(buf);
  glb_free(buf);
  return 0;

}

int
glbGetProjection(glb_projection in)
{
  int i;
  if(in==NULL) return -1;
  for(i=0;i<glbGetNumOfOscParams();i++) in=glbSetProjectionFlag(in,para_tab[i],i);
  for(i=0;i<glb_num_of_exps;i++) 
    in=glbSetDensityProjectionFlag(in,para_tab[i+glbGetNumOfOscParams()],i);
  return 0;
}

