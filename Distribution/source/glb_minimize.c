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


static double inp_errs[7];
static double start[7];
static int count=0;
static int errordim[32];
//int glb_single_experiment_number=0;
static double sysglb_calc_buffer[6];

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


// functions for returning the bfp of systematic parameters

static void to_glb_calc_buffer(double sys[],int len)
{
  int i;
  for (i=0;i<len;i++) sysglb_calc_buffer[i]=sys[i+1];
  return;
}

static void clear_glb_calc_buffer()
{
  int i;
  for (i=0;i<6;i++) sysglb_calc_buffer[i]=0.0;
  return;
}

static double* read_glb_calc_buffer()
{
  return &sysglb_calc_buffer[0];
}

//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//------------------------ Systematics -------------------------------------
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------


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
	     
	      m[i][j]=0.1;
	    }
	  else
	    {
	      m[i][j]=0.0;
	    }	
	}
    }
}

//This is a wraper for accessing the minimizer

static double chi_sys_wrap(double (*chi_func)(), int dimension, int alpha_in)
{
  double **mat; // contains the direction set
  double *sp; // stores the coordinates of the minimum
  double result=0; // stores the minimum value
  int it=0; // counts the number of iterations
  mat=glb_alloc_mat(1,dimension,1,dimension);
  sp=glb_alloc_vec(1,6); // that is the maximal length
  glb_rule_number=alpha_in;
  init_mat(mat,dimension);

  sp[1]=1;
  sp[2]=0;
  sp[3]=glb_bg_norm_center[alpha_in];
  sp[4]=glb_bg_tilt_center[alpha_in];
  sp[5]=glb_tre_null_center[alpha_in];
  sp[6]=glb_tre_tilt_center[alpha_in];
  glb_powell(sp,mat,dimension,TOLSYS,&it,&result,chi_func);  
  glb_free_vec(sp,1,6);
  glb_free_mat(mat,1,dimension,1,dimension);
  return result;
}
//Chi^2 where the systematics are integrated out

void glb_set_errordim(int typ, int rule)
{
  errordim[rule]=typ;
}

int glb_check_errordim(int rule)
{
  return errordim[rule];
}



// probably will be replaced by a more object oriented approach

// those two are needed in order to write unique
// chi^2 functions
double *glb_sys_errors;
double *glb_sys_centers;

struct glb_systematic glb_init_systematic(double (*chi_func)(),int dimension,
			    double *sp, double *errors, double (*evalf)(),
			    char info[])
{
  struct glb_systematic out;
  double *si,*se;
  int i;
  si=(double *) glb_malloc(sizeof(double)*dimension);
  se=(double *) glb_malloc(sizeof(double)*dimension);
  for(i=0;i<dimension;i++) si[i]=sp[i];
  for(i=0;i<dimension;i++) se[i]=errors[i];
  out.chi_func=chi_func;
  out.dimension=dimension;
  out.sp=si;
  out.errors=se;
  out.evalf=evalf;
  out.info=info;
  return out;
}



double glb_evaluate_chi(struct glb_systematic *in)
{
  int i;
  double **mat; // contains the direction set
  double *spi; // stores the coordinates of the minimum
  double result=0; // stores the minimum value
  int it=0; // counts the number of iterations
  mat=glb_alloc_mat(1,in->dimension,1,in->dimension);
  spi=glb_alloc_vec(1,in->dimension); // that is the maximal length
  init_mat(mat,in->dimension);
  for(i=1;i<in->dimension+1;i++) spi[i]=in->sp[i-1];
  // this part of the unified interface to this quantities
  glb_sys_centers=in->sp;
  glb_sys_errors=in->errors;
  //------------------------
  glb_powell(spi,mat,in->dimension,TOLSYS,&it,&result,in->chi_func);  
  glb_free_vec(spi,1,in->dimension);
  glb_free_mat(mat,1,in->dimension,1,in->dimension);
  return result;
}

// This is the central switch board for choosing the
// chi^2 function and  systematics treatment according
// to each value of errordim

static double chi_dispatch(int error_dim, int i)
{     
  // finally
  double erg2=0;
  
  if(error_dim==0) //chi^2 with 4 sys. parameters
    {
      erg2=chi_sys_wrap(&glb_chi_sys_w_bg,4,i);
      
    }
  else if(error_dim==1) //chi^2 with 6 sys. parameters
    {
      erg2=chi_sys_wrap(&glb_chi_sys_w_bg2,6,i); 
    }
  else if(error_dim==2) // chi^2 without sys. parameters
    {
      erg2=chi_sys_wrap(&glb_chi_sys_w_bg,0,i); 
      // it seems that the minimizer correctly interprest zero dimensions
      // as the function value at the starting point  
    }
  else if(error_dim==3) // very special for JHF-HK
    {
      if(i<4) erg2=chi_sys_wrap(&glb_chi_sys_w_bg,4,i); 
      else erg2=chi_sys_wrap(&glb_chi_sys_w_bgtot,4,i); 
    }
  else if(error_dim==4) // total chi^2 with 4 sys parameters
    {
      erg2=chi_sys_wrap(&glb_chi_sys_w_bgtot2,4,i); 
    }
  else if(error_dim==5) // obsolete
    {
      fprintf(stderr,"Warning: obsolete errordim\n");
      erg2=0;
    }
  else if(error_dim==6) // obsolete
    {
      fprintf(stderr,"Warning: obsolete errordim\n");
      erg2=0;
    }
  else if(error_dim==7) // free normalization, i.e. spectrum only
    {
      erg2=chi_sys_wrap(&glb_chi_spec,1,i); 
    }
  else if(error_dim==8) //total chi^2 w/o systematics
    {
      erg2=chi_sys_wrap(&glb_chi_sys_w_bgtot2,0,i); 
    }
  else if(error_dim==9) // chi^2 with 4 syst params and 
    // true energy calibration
    {
      erg2=chi_sys_wrap(&glb_chi_sys_w_bg_calib,4,i); 
    }
  else if(error_dim==10) // total chi^2 with 4 syst params 
    {
      erg2=chi_sys_wrap(&glb_chi_sys_w_bgtot,4,i); 
    }
  else if(error_dim==20) // total chi^2 with 4 syst params 
    {
      erg2=glb_evaluate_chi(&sys_calc[i]); 
    }
  else if(error_dim==21) // total chi^2 with 4 syst params 
    {
      erg2=sys_calc[i].evalf(&sys_calc[i]); 
    }
  else
    {
      erg2=0;
    }
return erg2;
}

static double ChiS0(int typ[])
{
  double erg2;
  int i,rul;
 
  ml_abort_check(GLB_MATHLINK);
 
  erg2=0;
  for (i=0;i<glb_num_of_rules;i++) 
    {
      erg2+=chi_dispatch(typ[i],i);
    }
  
  glb_rule_number=0;
  return erg2;
}

// this is the same as ChiSO()  but it allows access to a single rule !

static double ChiS0_Rule(int typ[],int rule)
{
  double erg2;
  clear_glb_calc_buffer();
  erg2=chi_dispatch(typ[rule],rule);
  glb_rule_number=0;
  return erg2;
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
      erg += ChiS0(glb_experiment_list[i]->errordim);
    }
  return erg;
}

//redefinition of Chi

static double Chi(double x[])
{
  int i;  
  double erg;
  glb_set_c_vacuum_parameters(x[0],x[1],x[2],x[3]);
  glb_set_c_squared_masses(0,x[4],x[5]); 
  for (i=0;i<glb_num_of_exps;i++)
    {
      glbSetExperiment(glb_experiment_list[i]);
      glb_set_profile_scaling(x[6+i],i);
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

static double SingleChi(double x[7],int exp)
{
  double erg;
  glb_set_c_vacuum_parameters(x[0], x[1],x[2],x[3]);
  glb_set_c_squared_masses(0,x[4],x[5]);

      glbSetExperiment(glb_experiment_list[exp]);
      glb_set_profile_scaling(x[6],exp);
      glb_set_new_rates();
     


  glbSetExperiment(glb_experiment_list[exp]);
  if (setjmp(env)==1) 
    {
      okay_flag=1;
      return erg;  
    }
  erg=ChiS0(errordim);
  return erg;
}

// chi^2 with systematics for each Experiment

static double SingleRuleChi(double x[7],int exp, int rule)
{
  double erg;
  glb_set_c_vacuum_parameters(x[0], x[1],x[2],x[3]);
  glb_set_c_squared_masses(0,x[4],x[5]);

      glbSetExperiment(glb_experiment_list[exp]);
      glb_set_profile_scaling(x[6],exp);
      glb_set_new_rates();
     


  glbSetExperiment(glb_experiment_list[exp]);
  erg=ChiS0_Rule(errordim,rule);
  return erg;
}

/* Wrappers for the API */

double glbChiSys(const glb_params in,int experiment, int rule)
{
  int i;
  double res,x[38];
  
  if(in==NULL) 
    {
      glb_error("Failure in glbChiSys: Input pointer must be non-NULL");
      return -1;
    }
  for(i=0;i<GLB_OSCP;i++) x[i]=glbGetOscParams(in,i);
  if(experiment==GLB_ALL)
    {
      if(rule==GLB_ALL)
	{
	  for(i=0;i<glb_num_of_exps;i++) 
	    x[i+GLB_OSCP]=glbGetDensityParams(in,i);
	  res=Chi(x);
	}
      else
	{
	  res=0;
	   for(i=0;i<glb_num_of_exps;i++)
	     {   
	       if(rule < glb_experiment_list[i]->numofrules)
		 {
		   x[GLB_OSCP]=glbGetDensityParams(in,i);
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
      x[GLB_OSCP]=glbGetDensityParams(in,experiment);
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
  inp_errs[6]=f;
  return 0;
}

int glb_set_starting_values(double a, double b, double c, double d, double e, double f)
{
  start[1]=a;
  start[2]=b;
  start[3]=c;
  start[4]=d;
  start[5]=e;
  start[6]=f;
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
  i=6+glb_num_of_exps;
  out=(double*) glb_malloc(i*sizeof(double));
  for(i=0;i<6;i++) out[i]=inp_errs[i];
  for(i=0;i<glb_num_of_exps;i++) out[i+6]=glb_experiment_list[i]->density_error;
  return out;
}

double* glb_return_input_values()
{
  double* out;
  int i;
  i=6+glb_num_of_exps;
  out=(double*) glb_malloc(i*sizeof(double));
  for(i=0;i<6;i++) out[i]=start[i];
  for(i=0;i<glb_num_of_exps;i++) out[i+6]=glb_experiment_list[i]->density_center;
  return out;
}
  
 

//--------------------------------------------------------------
//--------- Setting the glb_priors for th12 ------------------------
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



//----------------------------------------------------------------
//----------- shifted rate access --------------------------------
//----------------------------------------------------------------

static void ReturnShiftedRates(int rulenumber, int exp, double* inrates)
{
  int i;
  int bins;
  glbSetExperiment(glb_experiment_list[exp]);
  //glb_set_new_rates();
  //glbSetExperiment(&glb_experiment_list[exp]);
  //ChiS0_Rule(errordim,rulenumber);
  bins=rulenumber;
  bins=glb_get_number_of_bins();
  for(i=0;i<bins;i++) inrates[i]=glb_chirate[i];

}

//------------------------------------------------------------------
//---- Chi^2 with arbitrary number of free parameters --------------
//----------------------- 23.01.2004 -------------------------------
//------------------------------------------------------------------

//This serves for ChiNP

static double fix_params[37];
static int para_tab[37];
static int index_tab[37];
static int n_free;
static int n_fix;

//-----------------------


static void SelectProjection(int *vec)
{
  int i,c,c1,c2;

  for(i=0;i<6+glb_num_of_exps;i++) para_tab[i]=vec[i];
  c=0;
  c1=0;
  c2=0;
  for(i=0;i<6+glb_num_of_exps;i++)
    {
      if(para_tab[i]==1) c++;
    }
  for(i=0;i<6+glb_num_of_exps;i++)
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
static double (*glb_user_defined_prior)(const glb_params);
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
  double y[37];
  int i; 
  prior_input=glbAllocParams();
  count = count +1;
  for(i=0;i<n_free;i++) y[index_tab[i]]=x[i+1];  
  // This basically is superflous, however it appears to be safer not
  // change a global (i.e to this file) parameter (fix_params) at this place
  for(i=n_free;i<n_free+n_fix;i++) y[index_tab[i]]=fix_params[index_tab[i]]; 

  //fprintf(stderr,"x2 %f\n",x2[1]); 
  glb_set_c_vacuum_parameters(y[0],y[1],y[2],y[3]);
  glb_set_c_squared_masses(0,y[4],y[5]);
  for (i=0;i<glb_num_of_exps;i++)
    {
      glbSetExperiment(glb_experiment_list[i]);
      glb_set_profile_scaling(y[6+i],i);
      glb_set_new_rates();
    }
  
  erg2=ChiS();
  // adding  the user defined prior
  // shoufling the parameter vector y into an glb_params structure
  for (i=0;i<GLB_OSCP;i++) glbSetOscParams(prior_input,y[i],i);
  for (i=0;i<glb_num_of_exps;i++) glbSetDensityParams(prior_input,
						      y[i+GLB_OSCP],i);
  glbSetIteration(prior_input,count);
  
  erg2 = erg2 + glb_user_defined_prior(prior_input); 

  // adding the glb_priors
  /*for(i=0;i<n_free;i++)
    {
      if(index_tab[i]<6)
	{
	  erg2+=sglb_prior(y[index_tab[i]],
		      start[index_tab[i]],inp_errs[index_tab[i]]);  
	 
	}
      else
	{
	  erg2+=sglb_prior(y[index_tab[i]],
	  	     (glb_experiment_list[index_tab[i]-6])->density_center,
	       (glb_experiment_list[index_tab[i]-6])->density_error);
	}
	}
  */
  glbFreeParams(prior_input);
  return erg2;
}
 

// single-experiment functions ChiXXXSection
//This serves for ChiNP

static double s_fix_params[7];
static int s_para_tab[7];
static int s_index_tab[7];
static int s_n_free;
static int s_n_fix;

//-----------------------


static void single_SelectProjection(int set)
{
  int i,c,c1,c2;

  for(i=0;i<6;i++) s_para_tab[i]=para_tab[i];
  s_para_tab[6]=para_tab[6+set];
  c=0;
  c1=0;
  c2=0;
  for(i=0;i<6+1;i++)
    {
      if(s_para_tab[i]==1) c++;
    }
  for(i=0;i<6+1;i++)
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
  double y[7];
  int i; 
  prior_input=glbAllocParams();
  count = count +1;
  for(i=0;i<s_n_free;i++) y[s_index_tab[i]]=x[i+1];  
  // This basically is superflous, however it appears to be safer not
  // change a global (i.e to this file) parameter (fix_params) at this place
  for(i=s_n_free;i<s_n_free+s_n_fix;i++) y[s_index_tab[i]]
					 =s_fix_params[s_index_tab[i]]; 

  //fprintf(stderr,"x2 %f\n",x2[1]); 
  glb_set_c_vacuum_parameters(y[0],y[1],y[2],y[3]);
  glb_set_c_squared_masses(0,y[4],y[5]);
  
  glbSetExperiment(glb_experiment_list[glb_single_experiment_number]);
  glb_set_profile_scaling(y[6],glb_single_experiment_number);
  glb_set_new_rates();
    
  
  erg2=ChiS0(errordim);
  
  // adding  the user defined prior
  // shoufling the parameter vector y into an glb_params structure
  for (i=0;i<GLB_OSCP;i++) glbSetOscParams(prior_input,y[i],i);
  glbSetDensityParams(prior_input,y[GLB_OSCP],glb_single_experiment_number);
  glbSetIteration(prior_input,count);

  erg2 = erg2 + glb_user_defined_prior(prior_input); 

  
  // adding the glb_priors
  /*
   for(i=0;i<s_n_free;i++)
    {
      if(s_index_tab[i]<6)
	{
	  erg2+=sglb_prior(y[s_index_tab[i]],
		      start[s_index_tab[i]],inp_errs[s_index_tab[i]]);  
	 
	}
      else
	{
	  erg2+=sglb_prior(y[6],
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

  double x[38];
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
  
  for(i=0;i<GLB_OSCP;i++) x[i]=glbGetOscParams(in,i);
  x[GLB_OSCP]=glbGetDensityParams(in,exp);

  for(i=0;i<GLB_OSCP+1;i++) s_fix_params[i]=x[i];
  for(i=0;i<s_n_free;i++) sp2[i+1]=x[s_index_tab[i]];
  if (setjmp(env)==1) 
    {
      okay_flag=1;
      return er1;  
    }
  glb_powell2(sp2,mat2,dim,TOLOSC,&it,&er1,chi_NP);
  if(out!=NULL)
    {
      for(i=0;i<s_n_free;i++) 
	{
	  if(s_index_tab[i]<GLB_OSCP)
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
 

  double x[38];
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
  
  for(i=0;i<GLB_OSCP;i++) x[i]=glbGetOscParams(in,i);
  for(i=0;i<glb_num_of_exps;i++) x[i+GLB_OSCP]=glbGetDensityParams(in,i);

  for(i=0;i<6+glb_num_of_exps;i++) fix_params[i]=x[i];
  for(i=0;i<n_free;i++) sp2[i+1]=x[index_tab[i]];
  if (setjmp(env)==1) 
    {
      okay_flag=1;
      return er1;  
    }
  glb_powell2(sp2,mat2,dim,TOLOSC,&it,&er1,MD_chi_NP);
  if(out!=NULL)
    {
      for(i=0;i<n_free;i++) 
	{
	  if(index_tab[i]<GLB_OSCP)
	    glbSetOscParams(out,sp2[i+1],index_tab[i]);
	  else
	    glbSetDensityParams(out,sp2[i+1],index_tab[i]-GLB_OSCP);
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
  int swit[37]; 
  int buff[37];
  /*FIXME this should read buff= , that is the origin for the bug Campagne found */
  b=CheckProjection();
  for(i=0;i<6+glb_num_of_exps;i++) buff[i]=b[i];
  for(i=0;i<6+glb_num_of_exps;i++) swit[i]=GLB_FREE;
  swit[n-1]=GLB_FIXED;
  SelectProjection(swit);
  if(in==NULL)
    {
      glb_error("Failure in glbChiNP: Input pointer must be non-NULL");
      return -1;
    }
  
  if(exp==GLB_ALL) res=internal_glbChiNP(in,out);
  else res=internal_glbSingleChiNP(in,out,exp);
  SelectProjection(buff);
  return res; 
}


double glbChiTheta(const glb_params in,glb_params out, int exp)
{
  double res;
  res=glbChi1P(in,out,exp,2);
  return res;
}

double glbChiTheta23(const glb_params in,glb_params out, int exp)
{
  double res;
  res=glbChi1P(in,out,exp,3);
  return res;
}

double glbChiDelta(const glb_params in,glb_params out, int exp)
{
  double res;
  res=glbChi1P(in,out,exp,4);
  return res;
}

double glbChiDms(const glb_params in,glb_params out, int exp)
{
  double res;
  res=glbChi1P(in,out,exp,5);
  return res;
}

double glbChiDm(const glb_params in,glb_params out, int exp)
{
  double res;
  res=glbChi1P(in,out,exp,6);
  return res;
}

/* 2-d wrapper -- only ChiThetaDelta */

double glbChiThetaDelta(const glb_params in,glb_params out, int exp)
{
  double res;
  int i,*b;
  int swit[37]; 
  int buff[37];
  b=CheckProjection();
  for(i=0;i<6+glb_num_of_exps;i++) buff[i]=b[i];
  for(i=0;i<6+glb_num_of_exps;i++) swit[i]=GLB_FREE;
  swit[1]=GLB_FIXED;
  swit[3]=GLB_FIXED;
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
  int swit[37]; 
  int buff[37];
  b=CheckProjection();
  for(i=0;i<6+glb_num_of_exps;i++) buff[i]=b[i];
  for(i=0;i<6+glb_num_of_exps;i++) swit[i]=GLB_FREE;
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
  buf=glb_malloc(sizeof(int)*(GLB_OSCP+glb_num_of_exps));

  for(i=0;i<GLB_OSCP;i++) buf[i]=glbGetProjectionFlag(in,i);
  for(i=0;i<glb_num_of_exps;i++) buf[i+GLB_OSCP]=
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
  for(i=0;i<GLB_OSCP;i++) in=glbSetProjectionFlag(in,para_tab[i],i);
  for(i=0;i<glb_num_of_exps;i++) 
    in=glbSetDensityProjectionFlag(in,para_tab[i+GLB_OSCP],i);
  return 0;
}


