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





#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <setjmp.h>
#include <globes/globes.h>


#include "glb_types.h"
#include "glb_smear.h"
#include "glb_probability.h"
#include "glb_multiex.h"
#include "glb_fluxes.h"
#include "glb_error.h"


#define PI 3.1415
#define EARTH_RADIUS 6371
#define NUCLEON_MASS 1.6726E-27
#define CROSS_SECTION 0.67E-45
#define CROSS_SECTION_ANTI 0.34E-45

#define FLOAT double




static double baseline = 1000 ;   // baseline in km
static double target_mass = 10;   // target mass in kt;
static double years = 1;          // running time in years
static int bins = 20;          // number of energy bins
static double treshold = 5.0;     // treshold energy

static double Probs[3][3];
static double ProbsAnti[3][3];
double** glb_calc_smearing[32];
int* glb_calc_uprange[32];
int* glb_calc_lowrange[32];
static int* rules[32];
static int rule_length[32];
static double* rule_coeff[32]; 
static int* BGrules[32];
static int BGrule_length[32];
static double* BGrule_coeff[32];
 double* glb_calc_signal_rates[32];
 double* glb_calc_bg_rates[32]; 
// speeding up the rate calculation
static double* signal_prefactors[32];
static double* background_prefactors[32];
// ---------------------------------
double* glb_calc_rates_0[32];
double* glb_calc_rates_1[32];
double* glb_calc_rates_1T[32];
double* glb_calc_rates_1BG[32];
double* glb_calc_rates_1BGT[32];
double* glb_calc_energy_tab;
double glb_glb_sig_norm_error[32];
double glb_sig_tilt_error[32];
double glb_bg_norm_error[32];
double glb_bg_tilt_error[32];
double glb_bg_norm_center[32];
double glb_bg_tilt_center[32];

double* glb_calc_ratevec[32];
double* glb_calc_glb_calc_ratevec_bg[32];

/* rate glb_calc_buffer for chi functions in order to allow
 * access to shifted rate vector
 */
double* glb_chirate;

static int errortype=1;
double glb_tre_null_center[32];
double glb_tre_tilt_center[32];
double glb_tre_null_error[32];
double glb_tre_tilt_error[32];

/* this ensures that we have some bins below and
 * above the analysis range for properly calculating
 * the effect of an energy calibration uncertainty
 */
double glb_calc_energy_window[32][2];

int glb_calc_simbins;
double glb_calc_simtresh;
double glb_calc_simbeam;

static double* smm;
static int exptype = 1;
static int channel_list[32][6];
static int num_of_ch = 1;
int glb_num_of_rules=1;


/* workspace of glb_calc_simbins length */
double *glb_calc_buffer;

/* storing the channel result before glb_calc_smearing (glb_calc_simbins) */
double *glb_calc_chrb[32];

/* storing the channel result after glb_calc_smearing (bins) */
double *glb_calc_chra[32];


//binwise glb_calc_efficiencies
double* glb_calc_efficiencies[32];
//constant background
double* glb_calc_const_background[32];


double *glb_calc_user_pre_sm_channel[32];
double *glb_calc_user_post_sm_channel[32];
double *glb_calc_user_pre_sm_background[32];
double *glb_calc_user_post_sm_background[32];

/* the new glb_calc_smearing structure */

glb_smear *glb_calc_smear_data[32];


glb_flux *glb_calc_fluxes[32];
glb_xsec *glb_calc_xsecs[32];



// ------------------------------------------
// -----      Initialize Parameters      ----
// ------------------------------------------
// some of these functions may be superflous

void glb_set_type(int type)
{
  exptype =type;
}

int glb_get_type()
{
  return exptype;
}


void glb_set_baseline(double l)
{
  struct glb_experiment *e = glb_experiment_list[glb_current_exp];
  double sum,*ll;
  size_t s;
  int i;
  s = e->psteps;
  ll = e->lengthtab;
  sum=0;
  if(s<=0)
    { 
      glb_error("No baseline data to change");
      return;
    }
  for(i=0;i<s;i++) sum += ll[i];
  if(fabs(sum-l)>1E-8)
    {
      glb_error("Sum of lengths of layers does not match baseline");
      return;
    }
  baseline = sum;  
}

double glb_check_baseline()
{
  return baseline;
}



/* in kt */
void glb_set_target_mass(double mass)
{
  target_mass = mass;
}

double glb_check_target_mass()
{
  return target_mass;
}

void glb_set_years(double year)
{
  years = year;
}

double glb_check_years()
{
  return years;
}



void glb_set_energy_treshold(double etresh)
{
  treshold = etresh;
}

double glb_check_energy_treshold()
{
  return treshold;
}


void glb_set_number_of_bins(int nr)
{
  bins = nr;
}

int glb_check_number_of_bins()
{
  return bins;
}

int glb_get_number_of_bins()
{
  return bins;
}

/* Initialize the controlstructure for signal and backgrounds 
 * Here I introduce the concept of channels. In general there are 
 * three flavours in a beam, with to CP-signs, thus there six initial 
 * states, which can oscillate to six final states. The reaction in 
 * the detector can be anything, e.g. CC. This are altogether 72 
 * possibilities, of which you usually need only a few. Therefore we 
 * define channels. 
 * A channel is defined by a flux, the CP-sign, an initial and final state, 
 * the cross section and the energy resolution function.
 */


void glb_set_channel(int i, int polarity, int anti, int l, int m, int cc, 
		int energysmear)
{
  channel_list[i][0]=polarity;
  channel_list[i][1]=anti;
  channel_list[i][2]=l;
  channel_list[i][3]=m;
  channel_list[i][4]=cc;
  channel_list[i][5]=energysmear;
}

double* glb_get_channel(int i)
{
  double* erg;
  int k;
  erg=(double*) glb_malloc(6 * sizeof(double));
  for (k=0;k<6;k++) erg[k]=channel_list[i][k];
  return erg;
}



void glb_set_num_of_channels(int cc)
{
 num_of_ch=cc;
}


int glb_calc_check_num_of_channels()
{
 return num_of_ch;
}

/* Here we define rules, how to convert the channels into the observable 
 * signal and background. A rule has variable length, a list of channel 
 * numbers and the coefficients for each channel, like the efficiency.
 */

void glb_set_number_of_rules(int cn)
{
 glb_num_of_rules=cn;
}

int glb_calc_check_num_of_rules()
{
  return glb_num_of_rules;
}

void glb_set_rule(int i, int cn, int *rule, double *coeff)
{
  rule_coeff[i]= coeff;
  rule_length[i] = cn;
  rules[i]=rule;
}



int glb_calc_glb_calc_check_rule_length(int i, int bg)
{
  if (bg==1)
    {
      return rule_length[i];
    }
  else
  {
    return BGrule_length[i];
  }
}

double* glb_calc_glb_calc_check_rule_coeff(int i, int bg)
{
  if (bg==1)
    {
      return rule_coeff[i];
    }
  else
    {
      return BGrule_coeff[i];
    }
}
int* glb_calc_check_rule(int i,int  bg)
{
  if (bg==1)
    {
      return rules[i];
    }
  else
    {
      return BGrules[i];
    }
}

void glb_set_bg_rule(int i, int cn, int *rule, double *coeff)
{
  BGrule_coeff[i]= coeff;
  BGrule_length[i] = cn;
  BGrules[i]=rule;
}



// ----------------------------------------
// -----      Binned Rate Vectors      ----
// ----------------------------------------




// adaption to Martins version of oszprob.cc

static void CalcAllProbs(double en, double baseline)
{
  struct glb_experiment *e = glb_experiment_list[glb_current_exp];
  int i, j;
  int status;
  
  if ((status=glb_probability_matrix(Probs, +1, en, e->psteps, e->lengthtab, e->densitybuffer,
          (e->filter_state == GLB_ON) ? e->filter_value : -1.0)) != GLB_SUCCESS)
    glb_error("Calculation of oscillation probabilities failed.");

  if ((status=glb_probability_matrix(ProbsAnti, -1, en, e->psteps, e->lengthtab, e->densitybuffer,
          (e->filter_state == GLB_ON) ? e->filter_value : -1.0)) != GLB_SUCCESS)
    glb_error("Calculation of oscillation probabilities failed.");
}

// Sum of all NC rates

static double RatesSXX(double en, double baseline, int polarity, int anti, int l, int m,int ident)
{
  double ergebnis;
  /*glbXSection(ident,en,m,anti)*/
  
  ergebnis=glb_xsec_calc(en,m,anti,glb_calc_xsecs[ident])*
    (glb_flux_calc(en,baseline,polarity,1,anti,glb_calc_fluxes[polarity])+
     glb_flux_calc(en,baseline,polarity,2,anti,glb_calc_fluxes[polarity])+
     glb_flux_calc(en,baseline,polarity,3,anti,glb_calc_fluxes[polarity]))
    *target_mass;
  
  return ergebnis;
}


static double RatesNOSC(double en, double baseline, 
		 int polarity, int anti, int l, int m,int ident)
{
  double ergebnis;
 
  ergebnis=glb_xsec_calc(en,m,anti,glb_calc_xsecs[ident])*
    glb_flux_calc(en,baseline,polarity,l,anti,glb_calc_fluxes[polarity])
    *target_mass;
 
  return ergebnis;
}

// Rates for arbitrary XSection

static double RatesXX(double en, double baseline, int polarity, int anti, int l, int m,int ident)
{
  double ergebnis;
  if (anti == 1) 
    {
      ergebnis=glb_xsec_calc(en,m,anti,glb_calc_xsecs[ident])*
	glb_flux_calc(en,baseline,polarity,l,anti,glb_calc_fluxes[polarity])
	*Probs[l-1][m-1]*target_mass;
    }
  else
    {
      ergebnis=glb_xsec_calc(en,m,anti,glb_calc_xsecs[ident])
	*glb_flux_calc(en,baseline,polarity,l,anti,glb_calc_fluxes[polarity])
	*ProbsAnti[l-1][m-1]*target_mass;
      
    }
  return ergebnis;
}


 
double glb_calc_channel(int i, double en, double baseline)
{
  double ergebnis;
  int final_flavour,initial_flavour;
  int nosc=0;
  if(channel_list[i][3]>9) {final_flavour=channel_list[i][3]-10;nosc=1;}
  else {final_flavour=channel_list[i][3];}
  if(channel_list[i][2]>9) {initial_flavour=channel_list[i][2]-10;nosc=1;}
  else {initial_flavour=channel_list[i][2];}


    if(nosc!=0)
      {
	if(final_flavour!=initial_flavour) {return 0;}
	ergebnis=RatesNOSC(en,baseline,channel_list[i][0],
			   channel_list[i][1],initial_flavour,
			   final_flavour,channel_list[i][4]);
      }
    else
      { 
	ergebnis=RatesXX(en,baseline,channel_list[i][0],
			 channel_list[i][1],initial_flavour,
			 final_flavour,channel_list[i][4]);
      }
  return ergebnis;
}

static void CalcAllChannels(int mode)
{
  int i,k,s;
     
  for (i=0; i<glb_calc_simbins; i++)
    {
      CalcAllProbs(glb_sbin_center(i,glb_calc_smear_data[0]),baseline);
      for(k=0;k<num_of_ch;k++)
	{
	  if(!((mode!=0)&&((channel_list[k][2]>9)||(channel_list[k][3]>9))))
	    {
	      s=channel_list[k][5];
	      glb_calc_chrb[k][i]=glb_calc_channel(k,glb_sbin_center(i,glb_calc_smear_data[s]),baseline)
		* glb_calc_user_pre_sm_channel[k][i] 
		* glb_calc_smear_data[s]->simbinsize[i]
		+ glb_calc_user_pre_sm_background[k][i];
	    } 
	}	
    }
}




static double Signal(int cn, double en, double baseline)
{
  int k;
  double ergebnis=0.0;
  for (k=0;k<rule_length[cn];k++)
    {
      ergebnis += rule_coeff[cn][k] * glb_calc_channel(rules[cn][k],en,baseline);
    
     
    }
  return ergebnis;
}

static double Background(int cn, double en, double baseline)
{
  int k;
  double ergebnis=0.0;
  for (k=0;k<BGrule_length[cn];k++)
    {
      ergebnis += BGrule_coeff[cn][k] * glb_calc_channel(BGrules[cn][k],en,baseline);
    }
  return ergebnis;
}


inline static double BinEnergy(int i)
{
  /* This fixed a bug with energy window stuff for variable bin width */
  return glb_bin_center(i,glb_calc_smear_data[0]);
}

void glb_set_error_function(int typ)
{
  errortype=typ;
}

int glb_check_error_function()
{
  return errortype;
}

inline double glb_window_function(double low,double up,int bin)
{
  double en = glb_calc_energy_tab[bin];

  if (en > low && en < up)
    return 1.0;
  else
    return 0.0;
}

inline double errorf(double x, double y, int bin)
{
  double xx;
  double erg;
  xx=glb_calc_energy_tab[bin];

  switch(errortype)
    {
    case 1:
      if(y<=0)
	{
	  if(xx<=x) 
	    {
	      erg=0;
	    }
	  else
	    {
	      erg=1;
	    }
	  
	}
      else
	{
	  erg= 0.5/y*xx + (0.5 - x*0.5/y);
	  if(erg<=0) erg=0;
	  if(erg>=1) erg=1;
	}
      break;
    case 2:
      erg=1/xx*y+1;
      break;
    default: erg=1; break;
    }
  return erg;
}

static void SmearAllChannels(int mode);

static void BinnedSignalB(int mode)
{
  int k,j,i;
  CalcAllChannels(mode);
  SmearAllChannels(mode);
  for(i=0;i<bins;i++)
    {
      for(k=0;k<glb_num_of_rules;k++)
	{
	  glb_calc_bg_rates[k][i]=0;
	  for (j=0;j<BGrule_length[k];j++)
	    {
	      glb_calc_bg_rates[k][i] += BGrule_coeff[k][j] * 
		glb_calc_chra[BGrules[k][j]][i];
	    }
	  
	  glb_calc_signal_rates[k][i]=0;
	  for (j=0;j<rule_length[k];j++)
	    {
	      glb_calc_signal_rates[k][i] += rule_coeff[k][j] * 
		glb_calc_chra[rules[k][j]][i];
	    } 
	}     
    }  
}


//-----------------------------------------
//-----   Smearing                 --------
//-----------------------------------------



//this does the actual folding

static double* SmearRates(double* rates, int pl)
{
  int i;
  int j;
  double* rateout;
  rateout = (double*) glb_malloc(bins*sizeof(double));
  
  for (i=0; i<bins; i++)
    {
      rateout[i]=0;
      for (j=glb_calc_lowrange[pl][i]; j<glb_calc_uprange[pl][i]+1; j++)
	{
	  rateout[i] += glb_calc_smearing[pl][i][j-glb_calc_lowrange[pl][i]]* rates[j];
	}
      rateout[i]=rateout[i]*glb_calc_smear_data[pl]->simbinsize[i];
    }
  
  for(i=0;i<bins;i++) rates[i]=rateout[i];  
  glb_free(rateout);
  return rates;
}


static void SmearAllChannels(int mode)
{
  int i,k,l,s;
  
  for(k=0;k<num_of_ch;k++)
    {
      if(!((mode!=0)&&((channel_list[k][2]>9)||(channel_list[k][3]>9))))
	{
	  s=channel_list[k][5];
	  for(i=0;i<bins;i++)
	    {
	      glb_calc_chra[k][i]=0;
	      for(l=glb_calc_lowrange[s][i]; l<glb_calc_uprange[s][i]+1; l++)
		{
		  glb_calc_chra[k][i] += glb_calc_smearing[s][i][l-glb_calc_lowrange[s][i]]
		    *glb_calc_chrb[k][l];
		}
	      glb_calc_chra[k][i]=glb_calc_chra[k][i] * glb_calc_user_post_sm_channel[k][i]
		+ glb_calc_user_post_sm_background[k][i];
	    }
	}    
    }
}

static void SmearAllChannelsTilt(int mode)
{
  int i,k,l,s;
  
  for(k=0;k<num_of_ch;k++)
    {
      if(!((mode!=0)&&((channel_list[k][2]>9)||(channel_list[k][3]>9))))
	{
	  s=channel_list[k][5];
	  for(i=0;i<bins;i++)
	    {
	      glb_calc_chra[k][i]=0;
	      for(l=glb_calc_lowrange[s][i]; l<glb_calc_uprange[s][i]+1; l++)
		{
		  glb_calc_chra[k][i] += glb_calc_smearing[s][i][l-glb_calc_lowrange[s][i]]
		    *glb_calc_chrb[k][l]
		    *glb_sbin_center(l,glb_calc_smear_data[s]);
		}
	      glb_calc_chra[k][i]=glb_calc_chra[k][i] * glb_calc_user_post_sm_channel[k][i]
		+ glb_calc_user_post_sm_background[k][i];
	    }
	}
    }
}


void glb_shift_energy_scale(double g,double* ratesin, double* ratesout)
{
  /* needed by user_chi.c */
  int i;
  int k;
  double ki;
  double beam_energy=glb_get_max_energy();
  double di=treshold/((beam_energy-treshold)/bins);
  for (i=0;i<bins; i++)
    {
      ki=g*(i+0.5+di)+i;
      k=(int) floor(ki);
      // fprintf(stderr,"k val: %d ki val: %lf\n",k,ki);
      if(k<0) 
	{
	  ratesout[i]=0.0;
	}
      else if (k>bins-1) 
	{
	  ratesout[i]=0.0;
	}
      else if (k==bins-1)   // This prevents reading beyond array boundaries below - JK
        {
          ratesout[i]=(1+g)*((0-ratesin[k])*(ki-k)+ratesin[k]);
        }
      else
	{
	  ratesout[i]=(1+g)*((ratesin[k+1]-ratesin[k])*(ki-k)+ratesin[k]);
	}
    }
}

/* Set rate vector(s) */


void glb_set_rates()
{
  int i;
  int j;

 
  BinnedSignalB(0);
  for(i=0;i<bins;i++) glb_calc_energy_tab[i]=BinEnergy(i);

  for (i=0;i<glb_num_of_rules;i++)
    {
      for (j=0; j<bins;j++)
	{
	  glb_calc_buffer[j]=(glb_calc_signal_rates[i][j]+ (glb_calc_bg_rates[i][j])*
		  (glb_bg_norm_center[i]));
	   
	}
    
      for(j=0;j<bins;j++) glb_calc_rates_0[i][j]=glb_calc_buffer[j] 
			    *errorf(glb_tre_null_center[i],glb_tre_tilt_center[i],j)
			    *glb_window_function(glb_calc_energy_window[i][0],
					     glb_calc_energy_window[i][1],j);

    }
 

}

void glb_set_new_rates()
{
  int i;
  int j;
 
  
  
  BinnedSignalB(1);
  for (i=0;i<glb_num_of_rules;i++)
    {
      for (j=0;j<bins;j++) glb_calc_rates_1[i][j]= glb_calc_signal_rates[i][j];

      for (j=0;j<bins;j++) glb_calc_rates_1BG[i][j]= glb_calc_bg_rates[i][j];
    }
 
  /* FIXME add. glb_calc_buffer for tiltet rates prob. safer */
  SmearAllChannelsTilt(1);
  for (i=0;i<glb_num_of_rules;i++)
    {
      for (j=0;j<bins;j++) glb_calc_rates_1T[i][j]= glb_calc_signal_rates[i][j];
      for (j=0;j<bins;j++) glb_calc_rates_1BGT[i][j]= glb_calc_bg_rates[i][j];
      
    }
 
}


/* // ---------------------------------------- */
/* // -----      Statistics Functions      ---- */
/* // ---------------------------------------- */

void glb_calc_set_energy_window(int i, double a, double b)
{
  glb_calc_energy_window[i][0]=a;
  glb_calc_energy_window[i][1]=b;
}

void check_glb_calc_energy_window(int i,double *in)
{
  
  in[0]=glb_calc_energy_window[i][0];
  in[1]=glb_calc_energy_window[i][1];
  
}

void glb_calc_set_tresh_center(int i, double a, double b)
{
  glb_tre_null_center[i]=a;
  glb_tre_tilt_center[i]=b;
}


void glb_calc_set_tresh_errors(int i, double a, double b)
{
  glb_tre_null_error[i]=a;
  glb_tre_tilt_error[i]=b;
}

void glb_calc_check_tresh_center(int i,double *in)
{
  
  in[0]=glb_tre_null_center[i];
  in[1]=glb_tre_tilt_center[i];
  
}

void glb_calc_check_tresh_errors(int i,double *in)
{
  
  in[0]=glb_tre_null_error[i];
  in[1]=glb_tre_tilt_error[i];
  
}


void glb_set_signal_errors(int i,double norm, double tilt)
{
  glb_glb_sig_norm_error[i]=norm;
  glb_sig_tilt_error[i]=tilt;
}

double* glb_calc_check_signal_errors(int i)
{
  double* erg;
  erg=(double*) glb_malloc(2*sizeof(double));
  erg[0]=glb_glb_sig_norm_error[i];
  erg[1]=glb_sig_tilt_error[i];
  return erg;
}

double* glb_calc_check_bg_errors(int i)
{
  double* erg;
  erg=(double*) glb_malloc(2*sizeof(double));
  erg[0]=glb_bg_norm_error[i];
  erg[1]=glb_bg_tilt_error[i];
  return erg;
}

double* glb_calc_check_bg_center(int i)
{
  double* erg;
  erg=(double*) glb_malloc(2*sizeof(double));
  erg[0]=glb_bg_norm_center[i];
  erg[1]=glb_bg_tilt_center[i];
  return erg;
}

void glb_set_bg_errors(int i,double norm, double tilt)
{
  glb_bg_norm_error[i]=norm;
  glb_bg_tilt_error[i]=tilt;
}

void glb_set_bg_center(int i,double norm, double tilt)
{
  glb_bg_norm_center[i]=norm;
  glb_bg_tilt_center[i]=tilt;
}



//likelihood funktion
inline double likelihood(double ntheo, double nobs)
{
  if(nobs>0)
      return 2*(ntheo-nobs) + 2*nobs * log(nobs/ntheo);
  else
    return 2*(ntheo-nobs);
}

/* double likelihood(double ntheo, double nobs) */
/* { */
/*   if(nobs>0) */
/*     return (ntheo-nobs)*(ntheo-nobs)/ntheo; */
/*   else */
/*     return 0; */
/* } */




inline double glb_list_likelihood(double* ratest, double* ratesm)
{
  /* needed by user_chi.c */
  register int i;

  double sum = 0;
  for (i=0; i< bins; i++)
    sum += likelihood(ratest[i],ratesm[i]);
  return sum;
}

inline double glb_prior(double x, double center, double sigma)
{
  double tmp = (x - center)/sigma;
  return tmp*tmp;
}  


/* now comes a whole bunch of functions which calculate different chi^2 for different situations */

// standard chi^2
double glb_chi_sys_w_bg(double x[5])
{
  
  double erg;
  int i;
  

  // The following code admittedly does not look very nice, but it's fast because
  // it decouples consecutive instructions, enabling modern processor to partially
  // parallelise them
  for (i=0; i < bins; i++)
    glb_chirate[i] = x[1]*glb_calc_rates_1[glb_rule_number][i];
  for (i=0; i < bins; i++)
    glb_chirate[i] += x[2]*glb_calc_rates_1T[glb_rule_number][i];
  for (i=0; i < bins; i++)
    glb_chirate[i] += x[3]*glb_calc_rates_1BG[glb_rule_number][i];
  for (i=0; i < bins; i++)
    glb_chirate[i] += x[4]*glb_calc_rates_1BGT[glb_rule_number][i];
  for (i=0; i < bins; i++)
    glb_chirate[i] *= errorf(glb_tre_null_center[glb_rule_number],
                             glb_tre_tilt_center[glb_rule_number], i);
  for (i=0; i < bins; i++)
    glb_chirate[i] *= glb_window_function(glb_calc_energy_window[glb_rule_number][0],
                                          glb_calc_energy_window[glb_rule_number][1],i);
  
  
/*  for (i=0;i<bins;i++)
    {
      glb_chirate[i]= 
	(x[1]*glb_calc_rates_1[glb_rule_number][i]+
	x[2]*glb_calc_rates_1T[glb_rule_number][i]+
	x[3]*glb_calc_rates_1BG[glb_rule_number][i]+
	x[4]*glb_calc_rates_1BGT[glb_rule_number][i]) * errorf(glb_tre_null_center[glb_rule_number],glb_tre_tilt_center[glb_rule_number],i)
	*glb_window_function(glb_calc_energy_window[glb_rule_number][0],glb_calc_energy_window[glb_rule_number][1],i);  
    }*/
 
  erg = glb_list_likelihood(glb_calc_rates_0[glb_rule_number],glb_chirate)+
    glb_prior(x[1],1,glb_glb_sig_norm_error[glb_rule_number])+
    glb_prior(x[2],0,glb_sig_tilt_error[glb_rule_number])+
    glb_prior(x[3],glb_bg_norm_center[glb_rule_number],glb_bg_norm_error[glb_rule_number])+
    glb_prior(x[4],glb_bg_tilt_center[glb_rule_number],glb_bg_tilt_error[glb_rule_number]);
  
  return erg;
  
}

// standard chi^2
double glb_chi_spec(double x[2])
{
  
  double erg;
  int i;
  
  
  
  for (i=0;i<bins;i++)
    {
      glb_chirate[i]= 
	(x[1]*glb_calc_rates_1[glb_rule_number][i]+
	 glb_bg_norm_center[glb_rule_number]*glb_calc_rates_1BG[glb_rule_number][i]) 
	* errorf(glb_tre_null_center[glb_rule_number],glb_tre_tilt_center[glb_rule_number],i)
	*glb_window_function(glb_calc_energy_window[glb_rule_number][0],glb_calc_energy_window[glb_rule_number][1],i);
    
    }
 
  erg = glb_list_likelihood(glb_calc_rates_0[glb_rule_number],glb_chirate);
  
  return erg;
  
}

// split chi^2 total rates + spectrum
double glb_chi_split(double x[6])
{
   
  double erg;
  double r0;
  double r1;
  int i;
  r1=0;
  r0=0;
  
 
  for (i=0;i<bins;i++)
    {
      r1 += 
	(1.0*glb_calc_rates_1[glb_rule_number][i]+
	x[2]*glb_calc_rates_1T[glb_rule_number][i]+
	x[3]*glb_calc_rates_1BG[glb_rule_number][i]+
	x[4]*glb_calc_rates_1BGT[glb_rule_number][i]) * 
	errorf(glb_tre_null_center[glb_rule_number],glb_tre_tilt_center[glb_rule_number],i)
	*glb_window_function(glb_calc_energy_window[glb_rule_number][0],glb_calc_energy_window[glb_rule_number][1],i);
     
      glb_chirate[i]= 
	(x[1]*glb_calc_rates_1[glb_rule_number][i]+
	x[2]*glb_calc_rates_1T[glb_rule_number][i]+
	x[3]*glb_calc_rates_1BG[glb_rule_number][i]+
	x[4]*glb_calc_rates_1BGT[glb_rule_number][i]) * 
	errorf(glb_tre_null_center[glb_rule_number],glb_tre_tilt_center[glb_rule_number],i)
	*glb_window_function(glb_calc_energy_window[glb_rule_number][0],glb_calc_energy_window[glb_rule_number][1],i); 
      
     
      
      r0 += glb_calc_rates_0[glb_rule_number][i];

    }


  
  
  erg =  
    likelihood(r0,r1)+
    glb_list_likelihood(glb_calc_rates_0[glb_rule_number],glb_chirate)+
    //glb_prior(x[1],1,glb_glb_sig_norm_error[glb_rule_number])+
    glb_prior(x[2],0,glb_sig_tilt_error[glb_rule_number])+
    glb_prior(x[3],glb_bg_norm_center[glb_rule_number],glb_bg_norm_error[glb_rule_number])+
    glb_prior(x[4],glb_bg_tilt_center[glb_rule_number],glb_bg_tilt_error[glb_rule_number]);
  
  return erg;
  
}

double glb_chi_sys_w_bg2(double x[7])
{
  
  double erg;
  int i;
 
  for (i=0;i<bins;i++)
    {
      glb_chirate[i]= 
	(x[1]*glb_calc_rates_1[glb_rule_number][i]+
	x[2]*glb_calc_rates_1T[glb_rule_number][i]+
	x[3]*glb_calc_rates_1BG[glb_rule_number][i]+
	x[4]*glb_calc_rates_1BGT[glb_rule_number][i]) * errorf(x[5],x[6] ,i)
	*glb_window_function(glb_calc_energy_window[glb_rule_number][0],glb_calc_energy_window[glb_rule_number][1],i);
    
    }
 
  erg = glb_list_likelihood(glb_calc_rates_0[glb_rule_number],glb_chirate)+
    glb_prior(x[1],1,glb_glb_sig_norm_error[glb_rule_number])+
    glb_prior(x[2],0,glb_sig_tilt_error[glb_rule_number])+
    glb_prior(x[3],glb_bg_norm_center[glb_rule_number],glb_bg_norm_error[glb_rule_number])+
    glb_prior(x[4],glb_bg_tilt_center[glb_rule_number],glb_bg_tilt_error[glb_rule_number])+
    glb_prior(x[5],glb_tre_null_center[glb_rule_number],glb_tre_null_error[glb_rule_number])+
    glb_prior(x[6],glb_tre_tilt_center[glb_rule_number],glb_tre_tilt_error[glb_rule_number]);

  
  return erg;
  
}

//total rates
double glb_chi_sys_w_bgtot(double x[5])
{
  
  double erg;
  double r0;
  double r1;
  int i;
 
  r1=0;
  r0=0;
  
  for (i=0;i<bins;i++)
    {
      glb_chirate[i]= 
	(x[1]*glb_calc_rates_1[glb_rule_number][i]+
	x[2]*glb_calc_rates_1T[glb_rule_number][i]+
	x[3]*glb_calc_rates_1BG[glb_rule_number][i]+
	x[4]*glb_calc_rates_1BGT[glb_rule_number][i])
	*glb_window_function(glb_calc_energy_window[glb_rule_number][0],glb_calc_energy_window[glb_rule_number][1],i);
    }
  for (i=0;i<bins;i++)
    {
      r1 += glb_chirate[i];
      r0 += glb_calc_rates_0[glb_rule_number][i];
    }
  
  erg = likelihood(r0,r1)+
    glb_prior(x[1],1,glb_glb_sig_norm_error[glb_rule_number])+
    glb_prior(x[2],0,glb_sig_tilt_error[glb_rule_number])+
    glb_prior(x[3],glb_bg_norm_center[glb_rule_number],glb_bg_norm_error[glb_rule_number])+
    glb_prior(x[4],glb_bg_tilt_center[glb_rule_number],glb_bg_tilt_error[glb_rule_number]);
  
  return erg;
  
}

double glb_chi_sys_w_bgtot2(double x[5])
{
  
  double erg;
  int i;
  double r0;
  double r1;
  r0=0;
  r1=0;

  
  
  for (i=0;i<bins;i++)
    {
      glb_chirate[i]= 
	(x[1]*glb_calc_rates_1[glb_rule_number][i]+
	x[2]*glb_calc_rates_1T[glb_rule_number][i]+
	x[3]*glb_calc_rates_1BG[glb_rule_number][i]+
	x[4]*glb_calc_rates_1BGT[glb_rule_number][i]) * 
	errorf(glb_tre_null_center[glb_rule_number],glb_tre_tilt_center[glb_rule_number],i)
	*glb_window_function(glb_calc_energy_window[glb_rule_number][0],glb_calc_energy_window[glb_rule_number][1],i);
     
    }
  
 for (i=0;i<bins;i++)
    {
      r1 += glb_chirate[i];
      r0 += glb_calc_rates_0[glb_rule_number][i];
    } 

 
 erg = likelihood(r0,r1)+
   glb_prior(x[1],1,glb_glb_sig_norm_error[glb_rule_number])+
    glb_prior(x[2],0,glb_sig_tilt_error[glb_rule_number])+
    glb_prior(x[3],glb_bg_norm_center[glb_rule_number],glb_bg_norm_error[glb_rule_number])+
    glb_prior(x[4],glb_bg_tilt_center[glb_rule_number],glb_bg_tilt_error[glb_rule_number]);
  
  return erg;
  
}

// new tilting function in chi^2
double glb_chi_sys_w_bg_calib(double x[5])
{
  
  double erg;
  int i;

  glb_shift_energy_scale(x[2],glb_calc_rates_1[glb_rule_number],glb_calc_rates_1T[glb_rule_number]);
  glb_shift_energy_scale(x[4],glb_calc_rates_1BG[glb_rule_number],glb_calc_rates_1BGT[glb_rule_number]);
  //fprintf(stderr,"Tilt x[2] %lf\n",x[2]);
  //fprintf(stderr,"first bin %lf %lf\n",glb_calc_rates_1[glb_rule_number][10],glb_calc_rates_1T[glb_rule_number][10]); 
  for (i=0;i<bins;i++)
    {
      glb_chirate[i]= 
	(x[1]*glb_calc_rates_1T[glb_rule_number][i]+
	 x[3]*glb_calc_rates_1BGT[glb_rule_number][i]) 
	* errorf(glb_tre_null_center[glb_rule_number],glb_tre_tilt_center[glb_rule_number],i)
	* glb_window_function(glb_calc_energy_window[glb_rule_number][0],glb_calc_energy_window[glb_rule_number][1],i);
      
    }
  
  erg = glb_list_likelihood(glb_calc_rates_0[glb_rule_number],glb_chirate)+
    glb_prior(x[1],1,glb_glb_sig_norm_error[glb_rule_number])+
    glb_prior(x[2],0,glb_sig_tilt_error[glb_rule_number])+
    glb_prior(x[3],glb_bg_norm_center[glb_rule_number],glb_bg_norm_error[glb_rule_number])+
    glb_prior(x[4],glb_bg_tilt_center[glb_rule_number],glb_bg_tilt_error[glb_rule_number]);
  
  return erg;
  
}

// reprogrammed on 12.09.02
//______________________________________________________
//Acces to glb_calc_ratevectors
//______________________________________________________

static double signal_parts(int cn, int k, double en, double baseline)
{
  double ergebnis;  
  ergebnis = rule_coeff[cn][k] * glb_calc_channel(rules[cn][k],en,baseline);
  return ergebnis;
}

static double background_parts(int cn, int k, double en, double baseline)
{
  double ergebnis;
  ergebnis = BGrule_coeff[cn][k] * glb_calc_channel(BGrules[cn][k],en,baseline);
  return ergebnis;
}

static void binned_signal_parts(int rulenumber,int partnumber,double* inrates)
{
  int i;
  double beam_energy = glb_get_max_energy(); 
  int k;
 
  double deltae2=(beam_energy-treshold)/bins;
  double deltae = (glb_calc_simbeam-glb_calc_simtresh)/glb_calc_simbins;
  
  k=rulenumber;
 
  for (i=0; i<glb_calc_simbins; i++)
    {
      CalcAllProbs(glb_calc_simtresh+(0.5+i)*deltae,baseline);
      inrates[i]=signal_parts(k,partnumber,glb_calc_simtresh + (0.5+i)*deltae,baseline)*deltae2; 
    }
}

static void binned_background_parts(int rulenumber,int partnumber,double* inrates)
{
  int i;
  double beam_energy = glb_get_max_energy(); 
  int k;
  
  double deltae2=(beam_energy-treshold)/bins;
  double deltae = (glb_calc_simbeam-glb_calc_simtresh)/glb_calc_simbins;
  
  
  k=rulenumber;
 
  for (i=0; i<glb_calc_simbins; i++)
    {
      CalcAllProbs(glb_calc_simtresh+(0.5+i)*deltae,baseline);
      inrates[i]=background_parts(k,partnumber,glb_calc_simtresh + (0.5+i)*deltae,baseline)*deltae2*glb_bg_norm_center[k]; 
    }
}

// bug fixed in Returnglb_calc_signal_rates and glb_calc_bg_rates on 25.01.03
// SmearRates was called with the second argument 1, but
// correctly its i (the rulenumber)

static void Returnglb_calc_signal_rates(int rulenumber, int partnumber, double* inrates)
{
  int i;
  int j;

  double* out;
  i=rulenumber;
  out=inrates;
 
  binned_signal_parts(i,partnumber,out);
  for(j=0;j<bins;j++) glb_calc_energy_tab[j]=BinEnergy(j);
  out=SmearRates(out,i);
  for(j=0;j<bins;j++) out[j]=out[j] * errorf(glb_tre_null_center[i],glb_tre_tilt_center[i],j)*glb_window_function(glb_calc_energy_window[i][0],glb_calc_energy_window[i][1],j);
}

static void Returnglb_calc_bg_rates(int rulenumber, int partnumber, double* inrates)
{
  int i;
  int j;

  double* out;
  i=rulenumber;
  out=inrates;
 
  binned_background_parts(i,partnumber,out);
  
  for(j=0;j<bins;j++) glb_calc_energy_tab[j]=BinEnergy(j);
  out=SmearRates(out,i);
  for(j=0;j<bins;j++) out[j]=out[j] * errorf(glb_tre_null_center[i],glb_tre_tilt_center[i],j)*glb_window_function(glb_calc_energy_window[i][0],glb_calc_energy_window[i][1],j);
 
}



// function for avoiding malloc problems

void glb_remove_calc_pointers()
{
  int k;
  for(k=0;k<32;k++)
    {
      glb_calc_smear_data[k]=NULL;
      glb_calc_smearing[k]=NULL;
      rules[k]=NULL;
      rule_coeff[k]=NULL; 
      BGrules[k]=NULL;
      BGrule_coeff[k]=NULL;
      glb_calc_signal_rates[k]=NULL;
      glb_calc_bg_rates[k]=NULL; 
      signal_prefactors[k]=NULL;
      background_prefactors[k]=NULL;
      glb_calc_rates_0[k]=NULL;
      glb_calc_rates_1[k]=NULL;
      glb_calc_rates_1T[k]=NULL;
      glb_calc_rates_1BG[k]=NULL;
      glb_calc_rates_1BGT[k]=NULL;
      glb_calc_energy_tab=NULL;
      glb_calc_ratevec[k]=NULL;
      glb_calc_glb_calc_ratevec_bg[k]=NULL;
      glb_calc_efficiencies[k]=NULL;
      glb_calc_const_background[k]=NULL;
      glb_calc_chra[k]=NULL;
      glb_calc_chrb[k]=NULL;
    }
  glb_calc_buffer=NULL;
  smm=NULL;
  glb_chirate=NULL;

}
