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
#include <string.h>
#include <globes/globes.h>
#include "glb_multiex.h"
#include "glb_rate_engine.h"
#include "glb_minimize.h"


// ----------------------------------------------------------------------------
void glbSetUserChi(int exp, int rule, double (*chi_func)(double x[]), int dim,
                   double params[], double errors[], char *info)
// ----------------------------------------------------------------------------
// Assign a user-defined chi^2 function to an experiment.
// Parameters:
//  exp:      Number of experiment for which the chi^2 function is intended
//  rule:     Number of rule for which the chi^2 function is intended
//  chi_func: A pointer to the actual function
//  dim:      Dimension of parameter space expected by this function
//  params:   Initial values for the parameters
//  errors:   Errors of the parameters
//  info:     An optional string describing the chi^2 function
// ----------------------------------------------------------------------------
{
  int i,j;
  
  if (exp == GLB_ALL)
  {
    for (i=0; i < glb_num_of_exps; i++)
      if (rule == GLB_ALL)
      {
        for (j=0; j < glbGetNumberOfRules(i); j++)
          glb_experiment_list[i]->sys[j] =
               glb_init_systematic(chi_func, dim, params, errors, &glb_evaluate_chi, info);
      }
      else
        glb_experiment_list[i]->sys[rule] =
             glb_init_systematic(chi_func, dim, params, errors, &glb_evaluate_chi, info);
  }
  else
  {
    if (rule == GLB_ALL)
    {
      for (j=0; j < glbGetNumberOfRules(exp); j++)
        glb_experiment_list[exp]->sys[j] =
             glb_init_systematic(chi_func, dim, params, errors, &glb_evaluate_chi, info);
    }
    else
      glb_experiment_list[exp]->sys[rule] =
           glb_init_systematic(chi_func, dim, params, errors, &glb_evaluate_chi, info);
  }
}


// ----------------------------------------------------------------------------
int glbGetCurrentExp()
// ----------------------------------------------------------------------------
// Returns the number of the experiment that is currently being processed.
// User-defined chi^2 function will need this information.
// ----------------------------------------------------------------------------
{
  return glb_current_exp;
}

// ----------------------------------------------------------------------------
int glbGetCurrentRule()
// ----------------------------------------------------------------------------
// Returns the number of the rule that is currently being processed.
// User-defined chi^2 function will need this information.
// ----------------------------------------------------------------------------
{
  return glb_rule_number;
}

// ----------------------------------------------------------------------------
double glbGetSignalBinInRule(int exp, int rule, int bin)
// ----------------------------------------------------------------------------
// Returns the signal rate for the given experiment, rule and bin
// ----------------------------------------------------------------------------
{
  return glb_experiment_list[exp]->rates1T[rule][bin] *
            glb_window_function(glb_experiment_list[exp]->energy_window[glb_rule_number][0],
                                glb_experiment_list[exp]->energy_window[glb_rule_number][1],bin);
}

// ----------------------------------------------------------------------------
double glbGetBackgroundBinInRule(int exp, int rule, int bin)
// ----------------------------------------------------------------------------
// Returns the background rate for the given experiment, rule and bin
// ----------------------------------------------------------------------------
{
  return glb_experiment_list[exp]->rates1BGT[rule][bin] *
            glb_window_function(glb_experiment_list[exp]->energy_window[glb_rule_number][0],
                                glb_experiment_list[exp]->energy_window[glb_rule_number][1],bin);
}
  
// ----------------------------------------------------------------------------
double glbGetChannelBin(int exp, int channel, int bin)
// ----------------------------------------------------------------------------
// Returns the signal rate for the given experiment, channel and bin
// ----------------------------------------------------------------------------
{
  return glb_experiment_list[exp]->chra[channel][bin] *
            glb_window_function(glb_experiment_list[exp]->energy_window[glb_rule_number][0],
                                glb_experiment_list[exp]->energy_window[glb_rule_number][1],bin);
}

// ----------------------------------------------------------------------------
double glbGetObservedBinInRule(int exp, int rule, int bin)
// ----------------------------------------------------------------------------
// Returns the measured event rate for the given experiment, rule and bin
// ----------------------------------------------------------------------------
{
  return glb_experiment_list[exp]->rates0[rule][bin];
}

// ----------------------------------------------------------------------------
int glbShiftSignalEnergyScale(int exp, int rule, double amount)
// ----------------------------------------------------------------------------
// Applies an energy calibration error of magnitude amount to the signal event
// rates of the given experiment and rule
// The algorithm is the same as in the internal function glb_shift_energy_scale
// ----------------------------------------------------------------------------
{
  double tmp_rates[glb_experiment_list[exp]->numofbins];
  int i;

  glb_shift_energy_scale(amount, glb_experiment_list[exp]->SignalRates[rule], tmp_rates);
  memcpy(glb_experiment_list[exp]->rates1T[rule], tmp_rates,
         glb_experiment_list[exp]->numofbins * sizeof(tmp_rates[0]));
  return 0;
}

// ----------------------------------------------------------------------------
int glbShiftBackgroundEnergyScale(int exp, int rule, double amount)
// ----------------------------------------------------------------------------
// Applies an energy calibration error of magnitude amount to the background
// event rates of the given experiment and rule
// The algorithm is the same as in the internal function glb_shift_energy_scale
// ----------------------------------------------------------------------------
{
  double tmp_rates[glb_experiment_list[exp]->numofbins];
  int i;

  glb_shift_energy_scale(amount, glb_experiment_list[exp]->BackgroundRates[rule], tmp_rates);
  memcpy(glb_experiment_list[exp]->rates1BGT[rule], tmp_rates,
         glb_experiment_list[exp]->numofbins * sizeof(tmp_rates[0]));
  return 0;
}


// The following part of this file contains highly experimental functions for
// handling user-defined systematics which are not documented.
// If you want to implement your own systematics, user the routines above and
// follow their documentation
#ifdef GLB_EXPERIMENTAL

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <globes/globes.h>
#include "glb_rate_engine.h"
#include "glb_multiex.h"
#include "glb_minimize.h"



double my_chi(double x[])
{     
  double erg; 
  double y[5],z[5];
  y[1]=x[1]+x[5];
  y[2]=x[2]+x[6];
  y[3]=x[3];
  y[4]=x[4];
  z[1]=x[7]+x[5];
  z[2]=x[8]+x[6];
  z[3]=x[9];
  z[4]=x[10];

  glbSetExperiment(glb_experiment_list[0]);
  erg=glb_chi_sys_w_bg_calib(y);
  glbSetExperiment(glb_experiment_list[1]);
  erg+=glb_chi_sys_w_bg_calib(z);
  erg+=glb_prior(x[5],glb_sys_centers[4],glb_sys_errors[4])+
    glb_prior(x[6],glb_sys_centers[5],glb_sys_errors[5]);
  return erg;
}

double my_chi2(double x[])
{     
  double erg; 
  int i;

  glb_shift_energy_scale(x[2]+x[6],glb_experiment_list[0]->rates1[glb_rule_number],
      glb_experiment_list[0]->rates1T[glb_rule_number]);
  glb_shift_energy_scale(x[4],glb_experiment_list[0]->rates1BG[glb_rule_number],
      glb_experiment_list[0]->rates1BGT[glb_rule_number]);

  for (i=0;i<glb_experiment_list[0]->numofbins;i++)
    {

     
      glb_chirate[i]= 
	((x[1]+x[5])*glb_experiment_list[0]->rates1T[glb_rule_number][i]+
	 x[3]*glb_experiment_list[0]->rates1BGT[glb_rule_number][i]) 

	* glb_window_function(glb_experiment_list[0]->energy_window[glb_rule_number][0],glb_experiment_list[0]->energy_window[glb_rule_number][1],i);
      
    }  
  erg = glb_list_likelihood(glb_experiment_list[0]->rates0[glb_rule_number],glb_chirate)+
    glb_prior(x[1],glb_sys_centers[0],glb_sys_errors[0])+
    glb_prior(x[2],glb_sys_centers[1],glb_sys_errors[1])+
    glb_prior(x[3],glb_sys_centers[2],glb_sys_errors[2])+
    glb_prior(x[4],glb_sys_centers[3],glb_sys_errors[3]);

  glb_shift_energy_scale(x[8]+x[6],glb_experiment_list[1]->rates1[glb_rule_number],
	    glb_experiment_list[1]->rates1T[glb_rule_number]);
  glb_shift_energy_scale(x[10],glb_experiment_list[1]->rates1BG[glb_rule_number],
	    glb_experiment_list[1]->rates1BGT[glb_rule_number]);

  for (i=0;i<glb_experiment_list[1]->numofbins;i++)
    {
      glb_chirate[i]= 
	((x[7]+x[5])*glb_experiment_list[1]->rates1T[glb_rule_number][i]+
	 x[9]*glb_experiment_list[1]->rates1BGT[glb_rule_number][i]) 

		* glb_window_function(glb_experiment_list[1]->energy_window[glb_rule_number][0],glb_experiment_list[1]->energy_window[glb_rule_number][1],i);
    }
  
  erg += glb_list_likelihood(glb_experiment_list[1]->rates0[glb_rule_number],glb_chirate)+
    glb_prior(x[7],glb_sys_centers[6],glb_sys_errors[6])+
    glb_prior(x[8],glb_sys_centers[7],glb_sys_errors[7])+
    glb_prior(x[9],glb_sys_centers[8],glb_sys_errors[8])+
    glb_prior(x[10],glb_sys_centers[9],glb_sys_errors[9]);
  erg+=glb_prior(x[5],glb_sys_centers[4],glb_sys_errors[4])+
    glb_prior(x[6],glb_sys_centers[5],glb_sys_errors[5]);
 
  return erg;
}

static char mch3[]="This a full blown reactor setup with near and far detector.\n"
"Additionally it also includes what is called shape error\n"
"in Nucl. Phys. B 665 (2003) 487\n";

double my_chi3(double x[])
{
  
  double erg; 
  int i;

 
 

  for (i=0;i<glb_experiment_list[0]->numofbins;i++)
    {

     
      glb_chirate[i]= 
	((x[1]+x[5]+x[11+i])*glb_experiment_list[0]->rates1[glb_rule_number][i]+
	 x[3]*glb_experiment_list[0]->rates1BG[glb_rule_number][i]) 

	* glb_window_function(glb_experiment_list[0]->energy_window[glb_rule_number][0],glb_experiment_list[0]->energy_window[glb_rule_number][1],i);
      
    }

  glb_shift_energy_scale(x[2]+x[6],glb_chirate,
	   glb_experiment_list[0]->rates1T[glb_rule_number]);  
  
  erg = glb_list_likelihood(glb_experiment_list[0]->rates0[glb_rule_number],
	   glb_experiment_list[0]->rates1T[glb_rule_number])+
    glb_prior(x[1],glb_sys_centers[0],glb_sys_errors[0])+
    glb_prior(x[2],glb_sys_centers[1],glb_sys_errors[1])+
    glb_prior(x[3],glb_sys_centers[2],glb_sys_errors[2])+
    glb_prior(x[4],glb_sys_centers[3],glb_sys_errors[3]);

/*   glb_shift_energy_scale(x[8]+x[6],glb_experiment_list[1]->rates1[glb_rule_number], */
/* 	    glb_experiment_list[1]->rates1T[glb_rule_number]); */
/*   glb_shift_energy_scale(x[10],glb_experiment_list[1]->rates1BG[glb_rule_number], */
/* 	    glb_experiment_list[1]->rates1BGT[glb_rule_number]); */

  for (i=0;i<glb_experiment_list[1]->numofbins;i++)
    {
      glb_chirate[i]= 
	((x[7]+x[5]+x[11+i])*glb_experiment_list[1]->rates1[glb_rule_number][i]+
	 x[9]*glb_experiment_list[1]->rates1BG[glb_rule_number][i]) 

		* glb_window_function(glb_experiment_list[1]->energy_window[glb_rule_number][0],glb_experiment_list[1]->energy_window[glb_rule_number][1],i);
    }
  glb_shift_energy_scale(x[8]+x[6],glb_chirate,
	    glb_experiment_list[1]->rates1T[glb_rule_number]);  
  erg += glb_list_likelihood(glb_experiment_list[1]->rates0[glb_rule_number]
			 ,glb_experiment_list[1]->rates1T[glb_rule_number])+
    glb_prior(x[7],glb_sys_centers[6],glb_sys_errors[6])+
    glb_prior(x[8],glb_sys_centers[7],glb_sys_errors[7])+
    glb_prior(x[9],glb_sys_centers[8],glb_sys_errors[8])+
    glb_prior(x[10],glb_sys_centers[9],glb_sys_errors[9]);
  erg+=glb_prior(x[5],glb_sys_centers[4],glb_sys_errors[4])+
    glb_prior(x[6],glb_sys_centers[5],glb_sys_errors[5]);
  for (i=0;i<glb_experiment_list[1]->numofbins;i++)
    {
      erg+=glb_prior(x[11+i],glb_sys_centers[10+i],glb_sys_errors[10+i]);
    } 

  return erg;
}

#define NBG 200

static double bgvec[5][NBG][2];

void BGLoader(char* filename, int number)
{
  int i;
  FILE *fp = fopen(filename,"r");
  if (fp!=NULL)
    {
      for (i=0; i<=NBG; i++) fscanf(fp,"%lf %lf",
				    &bgvec[number][i][0],
				    &bgvec[number][i][1]);
      fclose(fp);
    }
  else
    {
      fprintf(stderr,"Error: Could not open %s\n",filename); 
    }
  
}

double back(double en, int ident)
{
  double incre;
  double ergebnis;  
  int lowerind;
  int higherind;
  double energy;
  double norm[4]={1.2453,3.28348,6.19681,6.25757};
  energy=1000.0*en-0.8;//converting to Thomas' units 

  incre = (bgvec[ident][NBG-1][0]-bgvec[ident][0][0])/(NBG-1);
  lowerind = floor((energy-bgvec[ident][0][0])/incre); 
  higherind = lowerind + 1;
  
  if (lowerind<0 || higherind>NBG-1) ergebnis=0.0; 
  else
  {
    ergebnis=((bgvec[ident][higherind][1]-bgvec[ident][lowerind][1])*
	      (((energy-bgvec[ident][0][0])/incre)-lowerind)+
	      bgvec[ident][lowerind][1]);
  }
  //fprintf(stderr,"%lf %lf %d %lf\n",ergebnis, incre,lowerind,en);
  
  return ergebnis*norm[ident];
}

static char mch4[]="This a full blown reactor setup with near and far detector.\n"
"Additionally it also includes several backgrounds from Thomas\n";




double my_chi4(double x[])
{
  
  double erg; 
  int i,k;
  double *chi0, *chi1;

  
  // getting the norms right for the background
  double NORM_ND=pow(1.1,2)/pow(glb_experiment_list[0]->baseline,2)
    *glb_experiment_list[0]->targetmass/100.0;
  double NORM_FD=pow(1.1,2)/pow(glb_experiment_list[1]->baseline,2)
    *glb_experiment_list[0]->targetmass/100.0;
  // creating a buffer
  chi0=(double *) malloc(sizeof(double)*glb_experiment_list[0]->numofbins);
  chi1=(double *) malloc(sizeof(double)*glb_experiment_list[1]->numofbins);
  //adding the backgrounds to the zero vectors 
  for (i=0;i<glb_experiment_list[0]->numofbins;i++)
    {
      chi0[i]= glb_experiment_list[0]->rates0[glb_rule_number][i];
      chi1[i]= glb_experiment_list[1]->rates0[glb_rule_number][i];
    }
  for (i=0;i<glb_experiment_list[0]->numofbins;i++)
    {
      for(k=0;k<4;k++)
	{
	  chi0[i] += 
	    glb_sys_centers[k+10]*back(glb_calc_energy_tab[i],k)*NORM_ND;
	  chi1[i] += 
	    glb_sys_centers[k+10]*back(glb_calc_energy_tab[i],k)*NORM_FD;
	}

      // here the error lies we need a buffer !!!


    }

  glb_shift_energy_scale(x[2]+x[6],glb_experiment_list[0]->rates1[glb_rule_number],
      glb_experiment_list[0]->rates1T[glb_rule_number]);
  glb_shift_energy_scale(x[4],glb_experiment_list[0]->rates1BG[glb_rule_number],
      glb_experiment_list[0]->rates1BGT[glb_rule_number]);

  for (i=0;i<glb_experiment_list[0]->numofbins;i++)
    {

     
      glb_chirate[i]= 
	((x[1]+x[5])*glb_experiment_list[0]->rates1T[glb_rule_number][i]+
	 x[3]*glb_experiment_list[0]->rates1BGT[glb_rule_number][i]
	 +x[11]*NORM_ND*back(glb_calc_energy_tab[i],0)
	 +x[12]*NORM_ND*back(glb_calc_energy_tab[i],1)
	 +x[13]*NORM_ND*back(glb_calc_energy_tab[i],2)
	 +x[14]*NORM_ND*back(glb_calc_energy_tab[i],3)


) 

	* glb_window_function(glb_experiment_list[0]->energy_window[glb_rule_number][0],glb_experiment_list[0]->energy_window[glb_rule_number][1],i);
      
    }  
  erg = glb_list_likelihood(chi0,glb_chirate)+
    glb_prior(x[1],glb_sys_centers[0],glb_sys_errors[0])+
    glb_prior(x[2],glb_sys_centers[1],glb_sys_errors[1])+
    glb_prior(x[3],glb_sys_centers[2],glb_sys_errors[2])+
    glb_prior(x[4],glb_sys_centers[3],glb_sys_errors[3]);

  glb_shift_energy_scale(x[8]+x[6],glb_experiment_list[1]->rates1[glb_rule_number],
	    glb_experiment_list[1]->rates1T[glb_rule_number]);
  glb_shift_energy_scale(x[10],glb_experiment_list[1]->rates1BG[glb_rule_number],
	    glb_experiment_list[1]->rates1BGT[glb_rule_number]);

  for (i=0;i<glb_experiment_list[1]->numofbins;i++)
    {
      glb_chirate[i]= 
	((x[7]+x[5])*glb_experiment_list[1]->rates1T[glb_rule_number][i]+
	 x[9]*glb_experiment_list[1]->rates1BGT[glb_rule_number][i]
	 +x[11]*NORM_FD*back(glb_calc_energy_tab[i],0)
	 +x[12]*NORM_FD*back(glb_calc_energy_tab[i],1)
	 +x[13]*NORM_FD*back(glb_calc_energy_tab[i],2)
	 +x[14]*NORM_FD*back(glb_calc_energy_tab[i],3)
	 ) 

		* glb_window_function(glb_experiment_list[1]->energy_window[glb_rule_number][0],glb_experiment_list[1]->energy_window[glb_rule_number][1],i);
    }
  
  erg += glb_list_likelihood(chi1,glb_chirate)+
    glb_prior(x[7],glb_sys_centers[6],glb_sys_errors[6])+
    glb_prior(x[8],glb_sys_centers[7],glb_sys_errors[7])+
    glb_prior(x[9],glb_sys_centers[8],glb_sys_errors[8])+
    glb_prior(x[10],glb_sys_centers[9],glb_sys_errors[9]);
  erg+=glb_prior(x[5],glb_sys_centers[4],glb_sys_errors[4])+
    glb_prior(x[6],glb_sys_centers[5],glb_sys_errors[5]);
  erg+=glb_prior(x[11],glb_sys_centers[10],glb_sys_errors[10])
    +glb_prior(x[12],glb_sys_centers[11],glb_sys_errors[11])
    +glb_prior(x[13],glb_sys_centers[12],glb_sys_errors[12])
    +glb_prior(x[14],glb_sys_centers[13],glb_sys_errors[13]);
  free(chi1);
  free(chi0);  
  return erg;
  
}


double my_chi5(double x[])
{
  
  double erg; 
  int i,k;
  double *chi0=NULL;
  double *chi1=NULL;
  // getting the norms right for the background
  double NORM_ND=pow(1.1,2)/pow(glb_experiment_list[0]->baseline,2)
    *glb_experiment_list[0]->targetmass/100.0;
  double NORM_FD=pow(1.1,2)/pow(glb_experiment_list[1]->baseline,2)
    *glb_experiment_list[1]->targetmass/100.0;

 // creating a buffer
  chi0=(double *) malloc(sizeof(double)*glb_experiment_list[0]->numofbins);
  chi1=(double *) malloc(sizeof(double)*glb_experiment_list[1]->numofbins);

  //adding the backgrounds to the zero vectors

  for (i=0;i<glb_experiment_list[0]->numofbins;i++)
    {
      chi0[i]= glb_experiment_list[0]->rates0[glb_rule_number][i];
      chi1[i]= glb_experiment_list[1]->rates0[glb_rule_number][i];
    }
  for (i=0;i<glb_experiment_list[0]->numofbins;i++)
    {
      for(k=0;k<4;k++)
	{
	  chi0[i] += 
	    glb_sys_centers[k+10]*back(glb_calc_energy_tab[i],k)*NORM_ND;
	  chi1[i] += 
	    glb_sys_centers[k+10]*back(glb_calc_energy_tab[i],k)*NORM_FD;
	}
    }



  for (i=0;i<glb_experiment_list[0]->numofbins;i++)
    {
      glb_chirate[i]= 
	((x[1]+x[5])*glb_experiment_list[0]->rates1[glb_rule_number][i]+
	 x[3]*glb_experiment_list[0]->rates1BG[glb_rule_number][i] 
	 +x[11]*NORM_ND*back(glb_calc_energy_tab[i],0)
	 +x[12]*NORM_ND*back(glb_calc_energy_tab[i],1)
	 +x[13]*NORM_ND*back(glb_calc_energy_tab[i],2)
	 +x[14]*NORM_ND*back(glb_calc_energy_tab[i],3)
	 ) 	
	* glb_window_function(glb_experiment_list[0]->energy_window[glb_rule_number][0],
			  glb_experiment_list[0]->energy_window[glb_rule_number][1],i); 
    }  

  
  glb_shift_energy_scale(x[2]+x[6],glb_chirate,
	    glb_experiment_list[0]->rates1T[glb_rule_number]);
  erg = glb_list_likelihood(chi0,glb_experiment_list[0]->rates1T[glb_rule_number])+
    glb_prior(x[1],glb_sys_centers[0],glb_sys_errors[0])+
    glb_prior(x[2],glb_sys_centers[1],glb_sys_errors[1])+
    glb_prior(x[3],glb_sys_centers[2],glb_sys_errors[2])+
    glb_prior(x[4],glb_sys_centers[3],glb_sys_errors[3]);

  glb_shift_energy_scale(x[8]+x[6],glb_experiment_list[1]->rates1[glb_rule_number],
	    glb_experiment_list[1]->rates1T[glb_rule_number]);
  glb_shift_energy_scale(x[10],glb_experiment_list[1]->rates1BG[glb_rule_number],
	    glb_experiment_list[1]->rates1BGT[glb_rule_number]);

  for (i=0;i<glb_experiment_list[1]->numofbins;i++)
    {
      glb_chirate[i]= 
	((x[7]+x[5])*glb_experiment_list[1]->rates1[glb_rule_number][i]+
	 x[9]*glb_experiment_list[0]->rates1BG[glb_rule_number][i] 
	 +x[11]*NORM_FD*back(glb_calc_energy_tab[i],0)
	 +x[12]*NORM_FD*back(glb_calc_energy_tab[i],1)
	 +x[13]*NORM_FD*back(glb_calc_energy_tab[i],2)
	 +x[14]*NORM_FD*back(glb_calc_energy_tab[i],3)
	 ) 
	* glb_window_function(glb_experiment_list[1]->energy_window[glb_rule_number][0],
			  glb_experiment_list[1]->energy_window[glb_rule_number][1],i);
    }
  
  glb_shift_energy_scale(x[8]+x[6],glb_chirate,
	    glb_experiment_list[1]->rates1T[glb_rule_number]);
  erg += glb_list_likelihood(chi1,glb_experiment_list[1]->rates1T[glb_rule_number])+
    glb_prior(x[7],glb_sys_centers[6],glb_sys_errors[6])+
    glb_prior(x[8],glb_sys_centers[7],glb_sys_errors[7])+
    glb_prior(x[9],glb_sys_centers[8],glb_sys_errors[8])+
    glb_prior(x[10],glb_sys_centers[9],glb_sys_errors[9]);
  erg+=glb_prior(x[5],glb_sys_centers[4],glb_sys_errors[4])+
    glb_prior(x[6],glb_sys_centers[5],glb_sys_errors[5]);
  erg+=glb_prior(x[11],glb_sys_centers[10],glb_sys_errors[10])
    +glb_prior(x[12],glb_sys_centers[11],glb_sys_errors[11])
    +glb_prior(x[13],glb_sys_centers[12],glb_sys_errors[12])
    +glb_prior(x[14],glb_sys_centers[13],glb_sys_errors[13]);
  
  // for having acces to the BG vector
  for (i=0;i<glb_experiment_list[0]->numofbins;i++)
    {
      for(k=0;k<4;k++)
	{
	glb_chirate[i]= 
	    glb_sys_centers[k+10]*back(glb_calc_energy_tab[i],k)*NORM_FD;
	  
	}
    }

  free(chi0);
  free(chi1);

  return erg;
}



double my_chi6(double x[])
{
  
  double erg; 
  int i,k;
  double *chi0=NULL;
  double *chi1=NULL;
  // getting the norms right for the background
  double NORM_ND=pow(1.1,2)/pow(glb_experiment_list[0]->baseline,2)
    *glb_experiment_list[0]->targetmass/100.0;
  double NORM_FD=pow(1.1,2)/pow(glb_experiment_list[1]->baseline,2)
    *glb_experiment_list[1]->targetmass/100.0;

  // creating a buffer
  chi0=(double *) malloc(sizeof(double)*glb_experiment_list[0]->numofbins);
  chi1=(double *) malloc(sizeof(double)*glb_experiment_list[1]->numofbins);

  //adding the backgrounds to the zero vectors

  for (i=0;i<glb_experiment_list[0]->numofbins;i++)
    {
      chi0[i]= glb_experiment_list[0]->rates0[glb_rule_number][i];
      chi1[i]= glb_experiment_list[1]->rates0[glb_rule_number][i];
    }
  for (i=0;i<glb_experiment_list[0]->numofbins;i++)
    {
      for(k=0;k<4;k++)
	{
	  chi0[i] += 
	    glb_sys_centers[k+10]*back(glb_calc_energy_tab[i],k)*NORM_ND;
	  chi1[i] += 
	    glb_sys_centers[k+14]*back(glb_calc_energy_tab[i],k)*NORM_FD;
	}
    }



  for (i=0;i<glb_experiment_list[0]->numofbins;i++)
    {
      glb_chirate[i]= 
	((x[1]+x[5])*glb_experiment_list[0]->rates1[glb_rule_number][i]+
	 x[3]*glb_experiment_list[0]->rates1BG[glb_rule_number][i] 
	 +x[11]*NORM_ND*back(glb_calc_energy_tab[i],0)
	 +x[12]*NORM_ND*back(glb_calc_energy_tab[i],1)
	 +x[13]*NORM_ND*back(glb_calc_energy_tab[i],2)
	 +x[14]*NORM_ND*back(glb_calc_energy_tab[i],3)
	 ) 	
	* glb_window_function(glb_experiment_list[0]->energy_window[glb_rule_number][0],
			  glb_experiment_list[0]->energy_window[glb_rule_number][1],i); 
    }  

  
  glb_shift_energy_scale(x[2]+x[6],glb_chirate,
	    glb_experiment_list[0]->rates1T[glb_rule_number]);
  erg = glb_list_likelihood(chi0,glb_experiment_list[0]->rates1T[glb_rule_number])+
    glb_prior(x[1],glb_sys_centers[0],glb_sys_errors[0])+
    glb_prior(x[2],glb_sys_centers[1],glb_sys_errors[1])+
    glb_prior(x[3],glb_sys_centers[2],glb_sys_errors[2])+
    glb_prior(x[4],glb_sys_centers[3],glb_sys_errors[3]);

  glb_shift_energy_scale(x[8]+x[6],glb_experiment_list[1]->rates1[glb_rule_number],
	    glb_experiment_list[1]->rates1T[glb_rule_number]);
  glb_shift_energy_scale(x[10],glb_experiment_list[1]->rates1BG[glb_rule_number],
	    glb_experiment_list[1]->rates1BGT[glb_rule_number]);

  for (i=0;i<glb_experiment_list[1]->numofbins;i++)
    {
      glb_chirate[i]= 
	((x[7]+x[5])*glb_experiment_list[1]->rates1[glb_rule_number][i]+
	 x[9]*glb_experiment_list[0]->rates1BG[glb_rule_number][i] 
	 +x[15]*NORM_FD*back(glb_calc_energy_tab[i],0)
	 +x[16]*NORM_FD*back(glb_calc_energy_tab[i],1)
	 +x[17]*NORM_FD*back(glb_calc_energy_tab[i],2)
	 +x[18]*NORM_FD*back(glb_calc_energy_tab[i],3)
	 ) 
	* glb_window_function(glb_experiment_list[1]->energy_window[glb_rule_number][0],
			  glb_experiment_list[1]->energy_window[glb_rule_number][1],i);
    }
  
  glb_shift_energy_scale(x[8]+x[6],glb_chirate,
	    glb_experiment_list[1]->rates1T[glb_rule_number]);
  erg += glb_list_likelihood(chi1,glb_experiment_list[1]->rates1T[glb_rule_number])+
    glb_prior(x[7],glb_sys_centers[6],glb_sys_errors[6])+
    glb_prior(x[8],glb_sys_centers[7],glb_sys_errors[7])+
    glb_prior(x[9],glb_sys_centers[8],glb_sys_errors[8])+
    glb_prior(x[10],glb_sys_centers[9],glb_sys_errors[9]);
  erg+=glb_prior(x[5],glb_sys_centers[4],glb_sys_errors[4])+
    glb_prior(x[6],glb_sys_centers[5],glb_sys_errors[5]);
  erg+=glb_prior(x[11],glb_sys_centers[10],glb_sys_errors[10])
    +glb_prior(x[12],glb_sys_centers[11],glb_sys_errors[11])
    +glb_prior(x[13],glb_sys_centers[12],glb_sys_errors[12])
    +glb_prior(x[14],glb_sys_centers[13],glb_sys_errors[13]);
  erg+=glb_prior(x[15],glb_sys_centers[14],glb_sys_errors[14])
    +glb_prior(x[16],glb_sys_centers[15],glb_sys_errors[15])
    +glb_prior(x[17],glb_sys_centers[16],glb_sys_errors[16])
    +glb_prior(x[18],glb_sys_centers[17],glb_sys_errors[17]);
  
  // for having acces to the BG vector
  for (i=0;i<glb_experiment_list[0]->numofbins;i++)
    {
      for(k=0;k<4;k++)
	{
	glb_chirate[i]= 
	    glb_sys_centers[k+10]*back(glb_calc_energy_tab[i],k)*NORM_FD;
	  
	}
    }

  free(chi0);
  free(chi1);

  return erg;
}

static char mch7[]="This setup includes most of the things we have:\n"
"   a near and a far detector\n"
"   2% bin-to-bin shape uncertainty (near and far correlated)\n"
"   2.5% flux uncertainty\n"
"   0.6% relative normalization\n"
"   0.5% relative energy calibration\n"
"   four types of background:\n"
"      flat (0.004), accidental (0.002), cosmogenic I (0.002),"
" cosmogenic I (0.002)\n"
"   including 50% error on each background"; 

double my_chi7(double x[])
{
  
  double erg; 
  int i,k;
  double *chi0=NULL;
  double *chi1=NULL;
  // getting the norms right for the background
  double NORM_ND=pow(1.1,2)/pow(glb_experiment_list[0]->baseline,2)
    *glb_experiment_list[0]->targetmass/100.0;
  double NORM_FD=pow(1.1,2)/pow(glb_experiment_list[1]->baseline,2)
    *glb_experiment_list[1]->targetmass/100.0;

  // creating a buffer
  chi0=(double *) malloc(sizeof(double)*glb_experiment_list[0]->numofbins);
  chi1=(double *) malloc(sizeof(double)*glb_experiment_list[1]->numofbins);

  //adding the backgrounds to the zero vectors

  for (i=0;i<glb_experiment_list[0]->numofbins;i++)
    {
      chi0[i]= glb_experiment_list[0]->rates0[glb_rule_number][i];
      chi1[i]= glb_experiment_list[1]->rates0[glb_rule_number][i];
    }
  for (i=0;i<glb_experiment_list[0]->numofbins;i++)
    {
      for(k=0;k<4;k++)
	{
	  chi0[i] += 
	    glb_sys_centers[k+10]*back(glb_calc_energy_tab[i],k)*NORM_ND;
	  chi1[i] += 
	    glb_sys_centers[k+14]*back(glb_calc_energy_tab[i],k)*NORM_FD;
	}
    }



  for (i=0;i<glb_experiment_list[0]->numofbins;i++)
    {
      glb_chirate[i]= 
	((x[1]+x[5]+x[19+i])*glb_experiment_list[0]->rates1[glb_rule_number][i]+
	 x[3]*glb_experiment_list[0]->rates1BG[glb_rule_number][i] 
	 +x[11]*NORM_ND*back(glb_calc_energy_tab[i],0)
	 +x[12]*NORM_ND*back(glb_calc_energy_tab[i],1)
	 +x[13]*NORM_ND*back(glb_calc_energy_tab[i],2)
	 +x[14]*NORM_ND*back(glb_calc_energy_tab[i],3)
	 ) 	
	* glb_window_function(glb_experiment_list[0]->energy_window[glb_rule_number][0],
			  glb_experiment_list[0]->energy_window[glb_rule_number][1],i); 
    }  

  
  glb_shift_energy_scale(x[2]+x[6],glb_chirate,
	    glb_experiment_list[0]->rates1T[glb_rule_number]);
  erg = glb_list_likelihood(chi0,glb_experiment_list[0]->rates1T[glb_rule_number])+
    glb_prior(x[1],glb_sys_centers[0],glb_sys_errors[0])+
    glb_prior(x[2],glb_sys_centers[1],glb_sys_errors[1])+
    glb_prior(x[3],glb_sys_centers[2],glb_sys_errors[2])+
    glb_prior(x[4],glb_sys_centers[3],glb_sys_errors[3]);

  glb_shift_energy_scale(x[8]+x[6],glb_experiment_list[1]->rates1[glb_rule_number],
	    glb_experiment_list[1]->rates1T[glb_rule_number]);
  glb_shift_energy_scale(x[10],glb_experiment_list[1]->rates1BG[glb_rule_number],
	    glb_experiment_list[1]->rates1BGT[glb_rule_number]);

  for (i=0;i<glb_experiment_list[1]->numofbins;i++)
    {
      glb_chirate[i]= 
	((x[7]+x[5]+x[19+i])*glb_experiment_list[1]->rates1[glb_rule_number][i]+
	 x[9]*glb_experiment_list[0]->rates1BG[glb_rule_number][i] 
	 +x[15]*NORM_FD*back(glb_calc_energy_tab[i],0)
	 +x[16]*NORM_FD*back(glb_calc_energy_tab[i],1)
	 +x[17]*NORM_FD*back(glb_calc_energy_tab[i],2)
	 +x[18]*NORM_FD*back(glb_calc_energy_tab[i],3)
	 ) 
	* glb_window_function(glb_experiment_list[1]->energy_window[glb_rule_number][0],
			  glb_experiment_list[1]->energy_window[glb_rule_number][1],i);
    }
  
  glb_shift_energy_scale(x[8]+x[6],glb_chirate,
	    glb_experiment_list[1]->rates1T[glb_rule_number]);

  erg += glb_list_likelihood(chi1,glb_experiment_list[1]->rates1T[glb_rule_number])+
    glb_prior(x[7],glb_sys_centers[6],glb_sys_errors[6])+
    glb_prior(x[8],glb_sys_centers[7],glb_sys_errors[7])+
    glb_prior(x[9],glb_sys_centers[8],glb_sys_errors[8])+
    glb_prior(x[10],glb_sys_centers[9],glb_sys_errors[9]);
  erg+=glb_prior(x[5],glb_sys_centers[4],glb_sys_errors[4])+
    glb_prior(x[6],glb_sys_centers[5],glb_sys_errors[5]);
  erg+=glb_prior(x[11],glb_sys_centers[10],glb_sys_errors[10])
    +glb_prior(x[12],glb_sys_centers[11],glb_sys_errors[11])
    +glb_prior(x[13],glb_sys_centers[12],glb_sys_errors[12])
    +glb_prior(x[14],glb_sys_centers[13],glb_sys_errors[13]);
  erg+=glb_prior(x[15],glb_sys_centers[14],glb_sys_errors[14])
    +glb_prior(x[16],glb_sys_centers[15],glb_sys_errors[15])
    +glb_prior(x[17],glb_sys_centers[16],glb_sys_errors[16])
    +glb_prior(x[18],glb_sys_centers[17],glb_sys_errors[17]);
  
  for (i=0;i<glb_experiment_list[0]->numofbins;i++)
    {
      erg+=glb_prior(x[19+i],glb_sys_centers[18+i],glb_sys_errors[18+i]);
    }

  // for having acces to the BG vector
  for (i=0;i<glb_experiment_list[0]->numofbins;i++)
    {
      for(k=0;k<4;k++)
	{
	glb_chirate[i]= 
	    glb_sys_centers[k+10]*back(glb_calc_energy_tab[i],k)*NORM_FD;
	  
	}
    }

  free(chi0);
  free(chi1);

  return erg;
}

static char mch8[]="This setup includes most of the things we have:\n"
"   a near and a far detector\n"
"   2% bin-to-bin shape uncertainty (near and far correlated)\n"
"   2.5% flux uncertainty\n"
"   0.6% relative normalization\n"
"   0.5% relative energy calibration\n"
"   four types of background:\n"
"      flat (0.004), accidental (0.002), cosmogenic I (0.002),"
" cosmogenic I (0.002)\n"
"   including 50% error on each background\n"
"since this is not enough ;-) there is also\n"
"   a 0.005 uncorrelated flat background (Gaußian chi^2)\n"
"   with 50% uncertainty";

static char mch9[]="This setup includes most of the things we have:\n"
"   a near and a far detector\n"
"   2% bin-to-bin shape uncertainty (near and far correlated)\n"
"   2.5% flux uncertainty\n"
"   0.6% relative normalization\n"
"   0.5% relative energy calibration\n"
"   four types of background:\n"
"      flat (0.004), accidental (0.002), cosmogenic I (0.002),"
" cosmogenic I (0.002)\n"
"   including 50% error on each background\n"
"   same as 8 but w/o uncorrelated background";

// this does not cost much, we just need a Gaußian chi^2



double gauss(double ex, double th, double fbg)
{
  double ti;
  ti=fabs(fbg)+fabs(ex);
  if(ti!=0) return pow((ex-th),2)/ti;
  else return 0;
}

double list_gauss(double *ex, double *th, double fb)
{
  int i;
  double erg=0;
  for(i=0;i<glb_experiment_list[0]->numofbins;i++) erg+= gauss(ex[i],th[i],fb);
  return erg;
}


// missing target mass in the calculation of the flat background !
double my_chi8(double x[])
{
  
  double erg; 
  int i,k;
  double *chi0=NULL;
  double *chi1=NULL;

  //Here we have a nice example of how it should not look like
  // By-passing the usual flux system is not a very clever idea
  // obviously ...

  // getting the norms right for the background
  double NORM_ND=pow(1.1,2)/pow(glb_experiment_list[0]->baseline,2)
    *glb_experiment_list[0]->targetmass/100.0;
  double NORM_FD=pow(1.1,2)/pow(glb_experiment_list[1]->baseline,2)
    *glb_experiment_list[1]->targetmass/100.0;
  double fsk_ND=pow(
		    (glb_experiment_list[0]->targetmass/315.809)*
		    (pow(1.1,2)/pow(glb_experiment_list[0]->baseline,2)
		     *60000.0/glb_experiment_list[0]->numofbins*0.5*0.005),2);
  double fsk_FD=pow(
		    (glb_experiment_list[1]->targetmass/315.809)*
		    (pow(1.1,2)/pow(glb_experiment_list[1]->baseline,2)
		     *60000.0/glb_experiment_list[1]->numofbins*0.5*0.005),2);
  double frg_ND=sqrt(fsk_ND)/0.5;
  double frg_FD=sqrt(fsk_FD)/0.5;
 
 /*  fsk_ND= pow(60000.0/glb_experiment_list[0]->numofbins*0.5*0.005,2); */
/*   fsk_FD=fsk_ND; */
/*   frg_ND=sqrt(fsk_ND)/0.5; */
/*   frg_FD=sqrt(fsk_FD)/0.5; */
/*   fsk_FD=0; */
/*   fsk_ND=0; */
 
  // creating a buffer
  chi0=(double *) malloc(sizeof(double)*glb_experiment_list[0]->numofbins);
  chi1=(double *) malloc(sizeof(double)*glb_experiment_list[1]->numofbins);

  //adding the backgrounds to the zero vectors

  for (i=0;i<glb_experiment_list[0]->numofbins;i++)
    {
      chi0[i]= glb_experiment_list[0]->rates0[glb_rule_number][i]+frg_ND;
      chi1[i]= glb_experiment_list[1]->rates0[glb_rule_number][i]+frg_FD;
    }
  for (i=0;i<glb_experiment_list[0]->numofbins;i++)
    {
      for(k=0;k<4;k++)
	{
	  chi0[i] += 
	    glb_sys_centers[k+10]*back(glb_calc_energy_tab[i],k)*NORM_ND;
	  chi1[i] += 
	    glb_sys_centers[k+14]*back(glb_calc_energy_tab[i],k)*NORM_FD;
	}
    
    }



  for (i=0;i<glb_experiment_list[0]->numofbins;i++)
    {
      glb_chirate[i]= 
	((x[1]+x[5]+x[19+i])*glb_experiment_list[0]->rates1[glb_rule_number][i]+
	 x[3]*glb_experiment_list[0]->rates1BG[glb_rule_number][i] 
	 +x[11]*NORM_ND*back(glb_calc_energy_tab[i],0)
	 +x[12]*NORM_ND*back(glb_calc_energy_tab[i],1)
	 +x[13]*NORM_ND*back(glb_calc_energy_tab[i],2)
	 +x[14]*NORM_ND*back(glb_calc_energy_tab[i],3)
	+frg_ND) 	
	* glb_window_function(glb_experiment_list[0]->energy_window[glb_rule_number][0],
			  glb_experiment_list[0]->energy_window[glb_rule_number][1],i); 
    }  

  
  glb_shift_energy_scale(x[2]+x[6],glb_chirate,
	    glb_experiment_list[0]->rates1T[glb_rule_number]);
  erg = list_gauss(chi0,glb_experiment_list[0]->rates1T[glb_rule_number],fsk_ND)+
    glb_prior(x[1],glb_sys_centers[0],glb_sys_errors[0])+
    glb_prior(x[2],glb_sys_centers[1],glb_sys_errors[1])+
    glb_prior(x[3],glb_sys_centers[2],glb_sys_errors[2])+
    glb_prior(x[4],glb_sys_centers[3],glb_sys_errors[3]);

  glb_shift_energy_scale(x[8]+x[6],glb_experiment_list[1]->rates1[glb_rule_number],
	    glb_experiment_list[1]->rates1T[glb_rule_number]);
  glb_shift_energy_scale(x[10],glb_experiment_list[1]->rates1BG[glb_rule_number],
	    glb_experiment_list[1]->rates1BGT[glb_rule_number]);

  for (i=0;i<glb_experiment_list[1]->numofbins;i++)
    {
      glb_chirate[i]= 
	((x[7]+x[5]+x[19+i])*glb_experiment_list[1]->rates1[glb_rule_number][i]+
	 x[9]*glb_experiment_list[0]->rates1BG[glb_rule_number][i] 
	 +x[15]*NORM_FD*back(glb_calc_energy_tab[i],0)
	 +x[16]*NORM_FD*back(glb_calc_energy_tab[i],1)
	 +x[17]*NORM_FD*back(glb_calc_energy_tab[i],2)
	 +x[18]*NORM_FD*back(glb_calc_energy_tab[i],3)
	 +frg_FD
	 ) 
	* glb_window_function(glb_experiment_list[1]->energy_window[glb_rule_number][0],
			  glb_experiment_list[1]->energy_window[glb_rule_number][1],i);
    }
  
  glb_shift_energy_scale(x[8]+x[6],glb_chirate,
	    glb_experiment_list[1]->rates1T[glb_rule_number]);

  erg += list_gauss(chi1,glb_experiment_list[1]->rates1T[glb_rule_number],fsk_FD)+
    glb_prior(x[7],glb_sys_centers[6],glb_sys_errors[6])+
    glb_prior(x[8],glb_sys_centers[7],glb_sys_errors[7])+
    glb_prior(x[9],glb_sys_centers[8],glb_sys_errors[8])+
    glb_prior(x[10],glb_sys_centers[9],glb_sys_errors[9]);
  erg+=glb_prior(x[5],glb_sys_centers[4],glb_sys_errors[4])+
    glb_prior(x[6],glb_sys_centers[5],glb_sys_errors[5]);
  erg+=glb_prior(x[11],glb_sys_centers[10],glb_sys_errors[10])
    +glb_prior(x[12],glb_sys_centers[11],glb_sys_errors[11])
    +glb_prior(x[13],glb_sys_centers[12],glb_sys_errors[12])
    +glb_prior(x[14],glb_sys_centers[13],glb_sys_errors[13]);
  erg+=glb_prior(x[15],glb_sys_centers[14],glb_sys_errors[14])
    +glb_prior(x[16],glb_sys_centers[15],glb_sys_errors[15])
    +glb_prior(x[17],glb_sys_centers[16],glb_sys_errors[16])
    +glb_prior(x[18],glb_sys_centers[17],glb_sys_errors[17]);
  
  for (i=0;i<glb_experiment_list[0]->numofbins;i++)
    {
      erg+=glb_prior(x[19+i],glb_sys_centers[18+i],glb_sys_errors[18+i]);
    }

  // for having acces to the BG vector
  for (i=0;i<glb_experiment_list[0]->numofbins;i++)
    {
      for(k=0;k<4;k++)
	{
	glb_chirate[i]= 
	    glb_sys_centers[k+10]*back(glb_calc_energy_tab[i],k)*NORM_FD;
	  
	}
    }

  free(chi0);
  free(chi1);

  return erg;
}


double my_chi9(double x[])
{
  
  double erg; 
  int i,k;
  double *chi0=NULL;
  double *chi1=NULL;
  // getting the norms right for the background
  double NORM_ND=pow(1.1,2)/pow(glb_experiment_list[0]->baseline,2)
    *glb_experiment_list[0]->targetmass/100.0;
  double NORM_FD=pow(1.1,2)/pow(glb_experiment_list[1]->baseline,2)
    *glb_experiment_list[1]->targetmass/100.0;

  // creating a buffer
  chi0=(double *) malloc(sizeof(double)*glb_experiment_list[0]->numofbins);
  chi1=(double *) malloc(sizeof(double)*glb_experiment_list[1]->numofbins);

  //adding the backgrounds to the zero vectors

  for (i=0;i<glb_experiment_list[0]->numofbins;i++)
    {
      chi0[i]= glb_experiment_list[0]->rates0[glb_rule_number][i];
      chi1[i]= glb_experiment_list[1]->rates0[glb_rule_number][i];
    }
  for (i=0;i<glb_experiment_list[0]->numofbins;i++)
    {
      for(k=0;k<4;k++)
	{
	  chi0[i] += 
	    glb_sys_centers[k+10]*back(glb_calc_energy_tab[i],k)*NORM_ND;
	  chi1[i] += 
	    glb_sys_centers[k+14]*back(glb_calc_energy_tab[i],k)*NORM_FD;
	}
    }



  for (i=0;i<glb_experiment_list[0]->numofbins;i++)
    {
      glb_chirate[i]= 
	((x[1]+x[5]+x[19+i])*glb_experiment_list[0]->rates1[glb_rule_number][i]+
	 x[3]*glb_experiment_list[0]->rates1BG[glb_rule_number][i] 
	 +x[11]*NORM_ND*back(glb_calc_energy_tab[i],0)
	 +x[12]*NORM_ND*back(glb_calc_energy_tab[i],1)
	 +x[13]*NORM_ND*back(glb_calc_energy_tab[i],2)
	 +x[14]*NORM_ND*back(glb_calc_energy_tab[i],3)
	 ) 	
	* glb_window_function(glb_experiment_list[0]->energy_window[glb_rule_number][0],
			  glb_experiment_list[0]->energy_window[glb_rule_number][1],i); 
    }  

  
  glb_shift_energy_scale(x[2]+x[6],glb_chirate,
	    glb_experiment_list[0]->rates1T[glb_rule_number]);
  erg = list_gauss(chi0,glb_experiment_list[0]->rates1T[glb_rule_number],0)+
    glb_prior(x[1],glb_sys_centers[0],glb_sys_errors[0])+
    glb_prior(x[2],glb_sys_centers[1],glb_sys_errors[1])+
    glb_prior(x[3],glb_sys_centers[2],glb_sys_errors[2])+
    glb_prior(x[4],glb_sys_centers[3],glb_sys_errors[3]);

  glb_shift_energy_scale(x[8]+x[6],glb_experiment_list[1]->rates1[glb_rule_number],
	    glb_experiment_list[1]->rates1T[glb_rule_number]);
  glb_shift_energy_scale(x[10],glb_experiment_list[1]->rates1BG[glb_rule_number],
	    glb_experiment_list[1]->rates1BGT[glb_rule_number]);

  for (i=0;i<glb_experiment_list[1]->numofbins;i++)
    {
      glb_chirate[i]= 
	((x[7]+x[5]+x[19+i])*glb_experiment_list[1]->rates1[glb_rule_number][i]+
	 x[9]*glb_experiment_list[0]->rates1BG[glb_rule_number][i] 
	 +x[15]*NORM_FD*back(glb_calc_energy_tab[i],0)
	 +x[16]*NORM_FD*back(glb_calc_energy_tab[i],1)
	 +x[17]*NORM_FD*back(glb_calc_energy_tab[i],2)
	 +x[18]*NORM_FD*back(glb_calc_energy_tab[i],3)
	 ) 
	* glb_window_function(glb_experiment_list[1]->energy_window[glb_rule_number][0],
			  glb_experiment_list[1]->energy_window[glb_rule_number][1],i);
    }
  
  glb_shift_energy_scale(x[8]+x[6],glb_chirate,
	    glb_experiment_list[1]->rates1T[glb_rule_number]);

  erg += list_gauss(chi1,glb_experiment_list[1]->rates1T[glb_rule_number],0)+
    glb_prior(x[7],glb_sys_centers[6],glb_sys_errors[6])+
    glb_prior(x[8],glb_sys_centers[7],glb_sys_errors[7])+
    glb_prior(x[9],glb_sys_centers[8],glb_sys_errors[8])+
    glb_prior(x[10],glb_sys_centers[9],glb_sys_errors[9]);
  erg+=glb_prior(x[5],glb_sys_centers[4],glb_sys_errors[4])+
    glb_prior(x[6],glb_sys_centers[5],glb_sys_errors[5]);
  erg+=glb_prior(x[11],glb_sys_centers[10],glb_sys_errors[10])
    +glb_prior(x[12],glb_sys_centers[11],glb_sys_errors[11])
    +glb_prior(x[13],glb_sys_centers[12],glb_sys_errors[12])
    +glb_prior(x[14],glb_sys_centers[13],glb_sys_errors[13]);
  erg+=glb_prior(x[15],glb_sys_centers[14],glb_sys_errors[14])
    +glb_prior(x[16],glb_sys_centers[15],glb_sys_errors[15])
    +glb_prior(x[17],glb_sys_centers[16],glb_sys_errors[16])
    +glb_prior(x[18],glb_sys_centers[17],glb_sys_errors[17]);
  
  for (i=0;i<glb_experiment_list[0]->numofbins;i++)
    {
      erg+=glb_prior(x[19+i],glb_sys_centers[18+i],glb_sys_errors[18+i]);
    }

  // for having acces to the BG vector
  for (i=0;i<glb_experiment_list[0]->numofbins;i++)
    {
      for(k=0;k<4;k++)
	{
	glb_chirate[i]= 
	    glb_sys_centers[k+10]*back(glb_calc_energy_tab[i],k)*NORM_FD;
	  
	}
    }

  free(chi0);
  free(chi1);

  return erg;
}



double my_chi10(double x[])
{
  
  double erg; 
  int i,k;
  double b1,b2,b3,b4;
  double sl,th,thf;
  double *chi0=NULL;
  double *chi1=NULL;
  double *ri0=NULL;
  double *ri1=NULL;

  // getting the norms right for the background
  double NORM_ND=pow(1.1,2)/pow(glb_experiment_list[0]->baseline,2)
    *glb_experiment_list[0]->targetmass/100.0;
  double NORM_FD=pow(1.1,2)/pow(glb_experiment_list[1]->baseline,2)
    *glb_experiment_list[1]->targetmass/100.0;
  double fsk_ND=pow(
		    (glb_experiment_list[0]->targetmass/315.809)*
		    (pow(1.1,2)/pow(glb_experiment_list[0]->baseline,2)
		     *60000.0/glb_experiment_list[0]->numofbins*0.5*0.005),2);
  double fsk_FD=pow(
		    (glb_experiment_list[1]->targetmass/315.809)*
		    (pow(1.1,2)/pow(glb_experiment_list[1]->baseline,2)
		     *60000.0/glb_experiment_list[1]->numofbins*0.5*0.005),2);
  double frg_ND=sqrt(fsk_ND)/0.5;
  double frg_FD=sqrt(fsk_FD)/0.5;
 
  // creating a buffer
  chi0=(double *) malloc(sizeof(double)*glb_experiment_list[0]->numofbins);
  chi1=(double *) malloc(sizeof(double)*glb_experiment_list[1]->numofbins);
  // creating a buffer
  ri0=(double *) malloc(sizeof(double)*glb_experiment_list[0]->numofbins);
  ri1=(double *) malloc(sizeof(double)*glb_experiment_list[1]->numofbins);

  //adding the backgrounds to the zero vectors

  for (i=0;i<glb_experiment_list[0]->numofbins;i++)
    {
      chi0[i]= glb_experiment_list[0]->rates0[glb_rule_number][i]+frg_ND;
      chi1[i]= glb_experiment_list[1]->rates0[glb_rule_number][i]+frg_FD;
    }
  for (i=0;i<glb_experiment_list[0]->numofbins;i++)
    {
      for(k=0;k<4;k++)
	{
	  chi0[i] += 
	    glb_sys_centers[k+10]*back(glb_calc_energy_tab[i],k)*NORM_ND;
	  chi1[i] += 
	    glb_sys_centers[k+14]*back(glb_calc_energy_tab[i],k)*NORM_FD;
	}
    
    }


  for (i=0;i<glb_experiment_list[0]->numofbins;i++)
    {

      b2= glb_window_function(glb_experiment_list[0]->energy_window[glb_rule_number][0],
			  glb_experiment_list[0]->energy_window[glb_rule_number][1],i); 
      b1=x[3]*glb_experiment_list[0]->rates1BG[glb_rule_number][i] 
	+x[11]*NORM_ND*back(glb_calc_energy_tab[i],0)
	+x[12]*NORM_ND*back(glb_calc_energy_tab[i],1)
	+x[13]*NORM_ND*back(glb_calc_energy_tab[i],2)
	+x[14]*NORM_ND*back(glb_calc_energy_tab[i],3)
	+frg_ND;

      b3=glb_window_function(glb_experiment_list[1]->energy_window[glb_rule_number][0],
			  glb_experiment_list[1]->energy_window[glb_rule_number][1],i);
      
      b4= x[9]*glb_experiment_list[0]->rates1BG[glb_rule_number][i] 
	 +x[15]*NORM_FD*back(glb_calc_energy_tab[i],0)
	 +x[16]*NORM_FD*back(glb_calc_energy_tab[i],1)
	 +x[17]*NORM_FD*back(glb_calc_energy_tab[i],2)
	 +x[18]*NORM_FD*back(glb_calc_energy_tab[i],3)
	+frg_FD;
      th=glb_experiment_list[0]->rates1[glb_rule_number][i];
      thf=glb_experiment_list[1]->rates1[glb_rule_number][i];

   
      sl=
	(thf*(frg_ND + chi0[i])*(-b4 + chi1[i]) - 
	 th*(frg_FD + chi1[i])*(b1 - chi0[i] + th*(x[1] + x[5])) 
	 -(frg_ND + chi0[i])*(x[5] + x[7])*pow(thf,2))*
	pow(0.02,2)*
	pow((frg_ND + chi0[i])*pow(thf,2)*
	    pow(0.02,2) + 
	    (frg_FD + chi1[i])*(frg_ND + chi0[i] + 
			       pow(th,2)*pow(0.02,2)),-1);
      ri0[i]= 
	((x[1]+x[5]+sl)*th+b1)*b2;
      ri1[i]= 
	((x[7]+x[5]+sl)*thf+b4)*b3;	

    }  

  
  glb_shift_energy_scale(x[2]+x[6],ri0,
	    glb_experiment_list[0]->rates1T[glb_rule_number]);
  erg = list_gauss(chi0,glb_experiment_list[0]->rates1T[glb_rule_number],fsk_ND)+
    glb_prior(x[1],glb_sys_centers[0],glb_sys_errors[0])+
    glb_prior(x[2],glb_sys_centers[1],glb_sys_errors[1])+
    glb_prior(x[3],glb_sys_centers[2],glb_sys_errors[2])+
    glb_prior(x[4],glb_sys_centers[3],glb_sys_errors[3]);

  glb_shift_energy_scale(x[8]+x[6],glb_experiment_list[1]->rates1[glb_rule_number],
	    glb_experiment_list[1]->rates1T[glb_rule_number]);
  glb_shift_energy_scale(x[10],glb_experiment_list[1]->rates1BG[glb_rule_number],
	    glb_experiment_list[1]->rates1BGT[glb_rule_number]);


  glb_shift_energy_scale(x[8]+x[6],ri1,
	    glb_experiment_list[1]->rates1T[glb_rule_number]);

  erg += list_gauss(chi1,glb_experiment_list[1]->rates1T[glb_rule_number],fsk_FD)+
    glb_prior(x[7],glb_sys_centers[6],glb_sys_errors[6])+
    glb_prior(x[8],glb_sys_centers[7],glb_sys_errors[7])+
    glb_prior(x[9],glb_sys_centers[8],glb_sys_errors[8])+
    glb_prior(x[10],glb_sys_centers[9],glb_sys_errors[9]);
  erg+=glb_prior(x[5],glb_sys_centers[4],glb_sys_errors[4])+
    glb_prior(x[6],glb_sys_centers[5],glb_sys_errors[5]);
  erg+=glb_prior(x[11],glb_sys_centers[10],glb_sys_errors[10])
    +glb_prior(x[12],glb_sys_centers[11],glb_sys_errors[11])
    +glb_prior(x[13],glb_sys_centers[12],glb_sys_errors[12])
    +glb_prior(x[14],glb_sys_centers[13],glb_sys_errors[13]);
  erg+=glb_prior(x[15],glb_sys_centers[14],glb_sys_errors[14])
    +glb_prior(x[16],glb_sys_centers[15],glb_sys_errors[15])
    +glb_prior(x[17],glb_sys_centers[16],glb_sys_errors[16])
    +glb_prior(x[18],glb_sys_centers[17],glb_sys_errors[17]);
  
 

  // for having acces to the BG vector
  for (i=0;i<glb_experiment_list[0]->numofbins;i++)
    {
      for(k=0;k<4;k++)
	{
	glb_chirate[i]= 
	    glb_sys_centers[k+10]*back(glb_calc_energy_tab[i],k)*NORM_FD;
	  
	}
    }

  free(chi0);
  free(chi1);
  free(ri0);
  free(ri1);
  return erg;
}



double my_chi11(double x[])
{
  
  double erg; 
  int i,k;
  double *b1, *b2, *b3, *b4;
  double sl, *th, *thf;
  double *chi0=NULL;
  double *chi1=NULL;
  double *ri0=NULL;
  double *ri1=NULL;

  // getting the norms right for the background
  double NORM_ND=pow(1.1,2)/pow(glb_experiment_list[0]->baseline,2)
    *glb_experiment_list[0]->targetmass/100.0;
  double NORM_FD=pow(1.1,2)/pow(glb_experiment_list[1]->baseline,2)
    *glb_experiment_list[1]->targetmass/100.0;
  double fsk_ND=pow(
		    (glb_experiment_list[0]->targetmass/315.809)*
		    (pow(1.1,2)/pow(glb_experiment_list[0]->baseline,2)
		     *60000.0/glb_experiment_list[0]->numofbins*0.5*0.005),2);
  double fsk_FD=pow(
		    (glb_experiment_list[1]->targetmass/315.809)*
		    (pow(1.1,2)/pow(glb_experiment_list[1]->baseline,2)
		     *60000.0/glb_experiment_list[1]->numofbins*0.5*0.005),2);
  double frg_ND=sqrt(fsk_ND)/0.5;
  double frg_FD=sqrt(fsk_FD)/0.5;
 
  // creating a buffer
  chi0=(double *) malloc(sizeof(double)*glb_experiment_list[0]->numofbins);
  chi1=(double *) malloc(sizeof(double)*glb_experiment_list[1]->numofbins);
  // creating a buffer
  ri0=(double *) malloc(sizeof(double)*glb_experiment_list[0]->numofbins);
  ri1=(double *) malloc(sizeof(double)*glb_experiment_list[1]->numofbins);

  b1=(double *) malloc(sizeof(double)*glb_experiment_list[0]->numofbins);
  b2=(double *) malloc(sizeof(double)*glb_experiment_list[0]->numofbins);
  b3=(double *) malloc(sizeof(double)*glb_experiment_list[0]->numofbins);
  b4=(double *) malloc(sizeof(double)*glb_experiment_list[0]->numofbins);

  th=(double *) malloc(sizeof(double)*glb_experiment_list[0]->numofbins);
  thf=(double *) malloc(sizeof(double)*glb_experiment_list[0]->numofbins);
 
  //adding the backgrounds to the zero vectors

  for (i=0;i<glb_experiment_list[0]->numofbins;i++)
    {
      chi0[i]= glb_experiment_list[0]->rates0[glb_rule_number][i]+frg_ND;
      chi1[i]= glb_experiment_list[1]->rates0[glb_rule_number][i]+frg_FD;
    }
  for (i=0;i<glb_experiment_list[0]->numofbins;i++)
    {
      for(k=0;k<4;k++)
	{
	  chi0[i] += 
	    glb_sys_centers[k+10]*back(glb_calc_energy_tab[i],k)*NORM_ND;
	  chi1[i] += 
	    glb_sys_centers[k+14]*back(glb_calc_energy_tab[i],k)*NORM_FD;
	}
    
    }
  
  glb_shift_energy_scale(x[8]+x[6],glb_experiment_list[1]->rates1[glb_rule_number],
	    glb_experiment_list[1]->rates1T[glb_rule_number]);
  glb_shift_energy_scale(x[10],glb_experiment_list[1]->rates1BG[glb_rule_number],
	    glb_experiment_list[1]->rates1BGT[glb_rule_number]);


  for (i=0;i<glb_experiment_list[0]->numofbins;i++)
    {
      
      b2[i]= glb_window_function(glb_experiment_list[0]->energy_window[glb_rule_number][0],
			  glb_experiment_list[0]->energy_window[glb_rule_number][1],i); 
      b1[i]=(x[3]*glb_experiment_list[0]->rates1BG[glb_rule_number][i] 
	     +x[11]*NORM_ND*back(glb_calc_energy_tab[i],0)
	     +x[12]*NORM_ND*back(glb_calc_energy_tab[i],1)
	     +x[13]*NORM_ND*back(glb_calc_energy_tab[i],2)
	     +x[14]*NORM_ND*back(glb_calc_energy_tab[i],3)
	     +frg_ND)*b2[i];
      
      b3[i]=glb_window_function(glb_experiment_list[1]->energy_window[glb_rule_number][0],
			    glb_experiment_list[1]->energy_window[glb_rule_number][1],i);
      
      // this presumably not correct it should read
      // glb_experiment_list[1]->rates1BG...
      b4[i]= (x[9]*glb_experiment_list[0]->rates1BG[glb_rule_number][i] 
		    +x[15]*NORM_FD*back(glb_calc_energy_tab[i],0)
		    +x[16]*NORM_FD*back(glb_calc_energy_tab[i],1)
		    +x[17]*NORM_FD*back(glb_calc_energy_tab[i],2)
		    +x[18]*NORM_FD*back(glb_calc_energy_tab[i],3)
		    +frg_FD)*b3[i];
      th[i]=glb_experiment_list[0]->rates1[glb_rule_number][i]*b2[i];
      thf[i]=glb_experiment_list[1]->rates1[glb_rule_number][i]*b3[i];
      
    }
  
  glb_shift_energy_scale(x[2]+x[6],b1,glb_chirate);
  for (i=0;i<glb_experiment_list[0]->numofbins;i++) b1[i]=glb_chirate[i];
  glb_shift_energy_scale(x[2]+x[6],th,glb_chirate);
  for (i=0;i<glb_experiment_list[0]->numofbins;i++) th[i]=glb_chirate[i]; 


  glb_shift_energy_scale(x[8]+x[6],b4,glb_chirate);
  for (i=0;i<glb_experiment_list[0]->numofbins;i++) b4[i]=glb_chirate[i];
  glb_shift_energy_scale(x[8]+x[6],thf,glb_chirate);
  for (i=0;i<glb_experiment_list[0]->numofbins;i++) thf[i]=glb_chirate[i];
 


  glb_shift_energy_scale(x[2]+x[6],ri0,
	    glb_experiment_list[0]->rates1T[glb_rule_number]);
  glb_shift_energy_scale(x[8]+x[6],ri1,
	    glb_experiment_list[1]->rates1T[glb_rule_number]);

  
  for  (i=0;i<glb_experiment_list[0]->numofbins;i++)
    {
      sl=
	(thf[i]*(frg_ND + chi0[i])*(-b4[i] + chi1[i]) - 
	 th[i]*(frg_FD + chi1[i])*(b1[i] - chi0[i] + th[i]*(x[1] + x[5])) 
	 -(frg_ND + chi0[i])*(x[5] + x[7])*pow(thf[i],2))*
	pow(0.02,2)*
	pow((frg_ND + chi0[i])*pow(thf[i],2)*
	    pow(0.02,2) + 
	    (frg_FD + chi1[i])*(frg_ND + chi0[i] + 
				pow(th[i],2)*pow(0.02,2)),-1);
      ri0[i]= 
	((x[1]+x[5]+sl)* th[i]+b1[i]);
      ri1[i]= 
	((x[7]+x[5]+sl)* thf[i]+b4[i]);	
      
    }  

  
 
  erg = list_gauss(chi0,ri0,fsk_ND)+
    glb_prior(x[1],glb_sys_centers[0],glb_sys_errors[0])+
    glb_prior(x[2],glb_sys_centers[1],glb_sys_errors[1])+
    glb_prior(x[3],glb_sys_centers[2],glb_sys_errors[2])+
    glb_prior(x[4],glb_sys_centers[3],glb_sys_errors[3]);

  erg += list_gauss(chi1,ri1,fsk_FD)+
    glb_prior(x[7],glb_sys_centers[6],glb_sys_errors[6])+
    glb_prior(x[8],glb_sys_centers[7],glb_sys_errors[7])+
    glb_prior(x[9],glb_sys_centers[8],glb_sys_errors[8])+
    glb_prior(x[10],glb_sys_centers[9],glb_sys_errors[9]);
  erg+=glb_prior(x[5],glb_sys_centers[4],glb_sys_errors[4])+
    glb_prior(x[6],glb_sys_centers[5],glb_sys_errors[5]);
  erg+=glb_prior(x[11],glb_sys_centers[10],glb_sys_errors[10])
    +glb_prior(x[12],glb_sys_centers[11],glb_sys_errors[11])
    +glb_prior(x[13],glb_sys_centers[12],glb_sys_errors[12])
    +glb_prior(x[14],glb_sys_centers[13],glb_sys_errors[13]);
  erg+=glb_prior(x[15],glb_sys_centers[14],glb_sys_errors[14])
    +glb_prior(x[16],glb_sys_centers[15],glb_sys_errors[15])
    +glb_prior(x[17],glb_sys_centers[16],glb_sys_errors[16])
    +glb_prior(x[18],glb_sys_centers[17],glb_sys_errors[17]);
  
 

  // for having acces to the BG vector
  for (i=0;i<glb_experiment_list[0]->numofbins;i++)
    {
      for(k=0;k<4;k++)
	{
	glb_chirate[i]= 
	    glb_sys_centers[k+10]*back(glb_calc_energy_tab[i],k)*NORM_FD;
	  
	}
    }

  free(chi0);
  free(chi1);
  free(ri0);
  free(ri1);
  free(b1);
  free(b2);
  free(b3);
  free(b4);
  free(th);
  free(thf);
  return erg;
}

char mch10[]="Same as eight, but with analytic treatment of the shape error.";


void glb_add_sys()
{
  double sp[6];
  double e[6];
  int i;
  double s[100];
  double eq[100];
  double sp1[]={1,0,1,0,0,0,1,0,1,0};
  double e1[]={0.006,0.005,0.005,0.006,0.025,0.025,0.006,0.005,0.005,0.006};
  sp[0]=1;
  sp[1]=0;
  sp[2]=glb_bg_norm_center[0];
  sp[3]=glb_bg_tilt_center[0];
  
  sys_list[0]=glb_init_systematic(glb_chi_sys_w_bg,4,&sp[0],&e[0],glb_evaluate_chi,"");
 
  sys_list[1]=glb_init_systematic(my_chi,10,&sp1[0],&e1[0],glb_evaluate_chi,"");
  sys_list[2]=glb_init_systematic(my_chi2,10,&sp1[0],&e1[0],glb_evaluate_chi,"");
  for(i=10;i<100;i++) s[i]=0;
  for(i=10;i<100;i++) eq[i]=0.02;
  for(i=0;i<10;i++) s[i]=sp1[i];
  for(i=0;i<10;i++) eq[i]=e1[i];
  sys_list[3]=glb_init_systematic(my_chi3,10+glb_experiment_list[0]->numofbins,
			&s[0],&eq[0],glb_evaluate_chi,mch3); 
  
  // Loading the backgrounds 

  BGLoader("flat.dat",0);
  BGLoader("accidental.dat",1);
  BGLoader("cosmogenic1.dat",2);
  BGLoader("cosmogenic2.dat",3);
  for(i=10;i<14;i++) s[i]=0.002;
  s[10]=0.004;
  for(i=10;i<14;i++) eq[i]=s[i]*0.5;


  sys_list[4]=glb_init_systematic(my_chi4,14,&s[0],&eq[0],glb_evaluate_chi,mch4); 
  sys_list[5]=glb_init_systematic(my_chi5,14,&s[0],&eq[0],glb_evaluate_chi,mch4);
  for(i=14;i<18;i++) s[i]=0.002;
  s[14]=0.004;
  for(i=14;i<18;i++) eq[i]=s[i]*0.5;
  sys_list[6]=glb_init_systematic(my_chi6,18,&s[0],&eq[0],glb_evaluate_chi,mch4); 
  sys_list[7]=glb_init_systematic(my_chi7,18+glb_experiment_list[0]->numofbins
			,&s[0],&eq[0],glb_evaluate_chi,mch7); 
  sys_list[8]=glb_init_systematic(my_chi8,18+glb_experiment_list[0]->numofbins
		      ,&s[0],&eq[0],glb_evaluate_chi,mch8); 
  sys_list[9]=glb_init_systematic(my_chi9,18+glb_experiment_list[0]->numofbins
		      ,&s[0],&eq[0],glb_evaluate_chi,mch9); 
  sys_list[10]=glb_init_systematic(my_chi10,18
		      ,&s[0],&eq[0],glb_evaluate_chi,"obsolete"); 
  sys_list[11]=glb_init_systematic(my_chi11,18
		      ,&s[0],&eq[0],glb_evaluate_chi,mch10); 
}

#else /* GLB_EXPERIMENTAL */

void glb_add_sys()
{
  return;
}

#endif /* GLB_EXPERIMENTAL */
