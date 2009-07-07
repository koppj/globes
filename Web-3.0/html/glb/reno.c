/* GLoBES -- General LOng Baseline Experiment Simulator
 * (C) 2002 - 2009,  The GLoBES Team
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
#include <math.h>
#include <globes/globes.h>   /* GLoBES library */

#define MAX_SYS 500

static double re_sys_errors[MAX_SYS];       /* Uncertainties of
					       systematics params */
static double re_sys_startval[MAX_SYS];     /* Starting values for
					       systematics
					       minimizer */
static double re_sigma_binbin = 0.0;   

#define RE_EXP_FAR1  0
#define RE_EXP_FAR2  1
#define RE_EXP_FAR3  2
#define RE_EXP_FAR4  3
#define RE_EXP_FAR5  4
#define RE_EXP_FAR6  5

#define RE_EXP_NEAR1 6
#define RE_EXP_NEAR2 7
#define RE_EXP_NEAR3 8
#define RE_EXP_NEAR4 9
#define RE_EXP_NEAR5 10
#define RE_EXP_NEAR6 11


/***************************************************************************
 *                        H E L P E R   F U N C T I O N S                  *
 ***************************************************************************/

/* Minimum of two numbers */
inline double min(double x, double y)
{
  if (x < y)
    return x;
  else
    return y;
}

/* Square of real number */
inline double square(double x)
{
  return x*x;
}

/* Gauss likelihood (this is sufficient here due to the large event
 * numbers in a reactor experiment; for other setups, one should use
 * Poisson statistics) */
inline double likelihood(double true_rate, double fit_rate, double sqr_sigma)
{
  if (sqr_sigma > 0)
    return square(true_rate - fit_rate) / sqr_sigma;
  else
    return 0.0;
}


double reno_distances[12]={1.52164, 1.43291, 1.38642, 1.38642, 1.43291, 
			   1.52164, 0.7033,0.482169, 0.318457, 0.318457, 
			   0.482169, 0.7033};


double re_errors[29]={0.02,0.02,0.02,0.02,0.02,0.02,
		      0.005,0.005,
		      0.001,0.001,
		      0.02,
		      0.02,0.02,0.02,
		      0.02,0.02,0.02,
		      0.02,0.02,0.02,
		      0.02,0.02,0.02,
		      0.02,0.02,0.02,
		      0.02,0.02,0.02};


/***************************************************************
 *                     RENO                                    *   
 ***************************************************************
 * x[0] - x[5] are the errors associated with each core        *
 * x[6] is the far detector norm error                         *
 * x[7] is the near decteor 1 norm error                       *
 * x[8] is the overall normalization error                     *
 * x[9] far energy scale                                       * 
 * x[10] near energy scale                                     *
 * x[11]-x[13] is the isotopic abundance error for core 1      *
 * x[14]-x[16] is the isotopic abundance error for core 2      *
 * x[17]-x[19] is the isotopic abundance error for core 3      * 
 * x[20]-x[22] is the isotopic abundance error for core 4      *
 * x[23]-x[25] is the isotopic abundance error for core 5      *
 * x[26]-x[28] is the isotopic abundance error for core 6      *
 ***************************************************************
 */

double chi_reno(int exp, int rule, int n_params, double *x, double *errors,
                     void *user_data)
{
  int n_bins = glbGetNumberOfBins(RE_EXP_FAR1);
  int j;
  double signal_norm_N1, signal_norm_F1;
  double signal_norm_N3, signal_norm_F3;
  double signal_norm_N2, signal_norm_F2;
  double signal_norm_N4, signal_norm_F4;
  double signal_norm_N5, signal_norm_F5;
  double signal_norm_N6, signal_norm_F6;
 
  double signal_N1[n_bins];
  double signal_F[n_bins];

  /* the new index 4 is to have room for the event rates from the four
   *  isotopes seperately 
   */
  double *true_rates_F1[4];
  double *true_rates_F2[4];
  double *true_rates_F3[4];
  double *true_rates_F4[4];
  double *true_rates_F5[4];
  double *true_rates_F6[4];

  double *true_rates_N1[4];
  double *true_rates_N2[4];
  double *true_rates_N3[4];
  double *true_rates_N4[4];
  double *true_rates_N5[4];
  double *true_rates_N6[4];

  double *signal_fit_rates_F1[4];
  double *signal_fit_rates_F2[4];
  double *signal_fit_rates_F3[4];
  double *signal_fit_rates_F4[4];
  double *signal_fit_rates_F5[4];
  double *signal_fit_rates_F6[4];

  double *signal_fit_rates_N1[4];
  double *signal_fit_rates_N2[4];
  double *signal_fit_rates_N3[4];
  double *signal_fit_rates_N4[4];
  double *signal_fit_rates_N5[4];
  double *signal_fit_rates_N6[4];


  for(j=0;j<4;j++) /* looping over the isotopes */
    {
      true_rates_F1[j]= glbGetRuleRatePtr(RE_EXP_FAR1, j);
      true_rates_F2[j]= glbGetRuleRatePtr(RE_EXP_FAR2, j);
      true_rates_F3[j]= glbGetRuleRatePtr(RE_EXP_FAR3, j);

      true_rates_F4[j]= glbGetRuleRatePtr(RE_EXP_FAR4, j);
      true_rates_F5[j]= glbGetRuleRatePtr(RE_EXP_FAR5, j);
      true_rates_F6[j]= glbGetRuleRatePtr(RE_EXP_FAR6, j);

      true_rates_N1[j]= glbGetRuleRatePtr(RE_EXP_NEAR1, j);
      true_rates_N2[j]= glbGetRuleRatePtr(RE_EXP_NEAR2, j);
      true_rates_N3[j]= glbGetRuleRatePtr(RE_EXP_NEAR3, j);

      true_rates_N4[j]= glbGetRuleRatePtr(RE_EXP_NEAR4, j);
      true_rates_N5[j]= glbGetRuleRatePtr(RE_EXP_NEAR5, j);
      true_rates_N6[j]= glbGetRuleRatePtr(RE_EXP_NEAR6, j);

      signal_fit_rates_F1[j] = glbGetSignalFitRatePtr(RE_EXP_FAR1,j);
      signal_fit_rates_F2[j] = glbGetSignalFitRatePtr(RE_EXP_FAR2,j);
      signal_fit_rates_F3[j] = glbGetSignalFitRatePtr(RE_EXP_FAR3,j);

      signal_fit_rates_F4[j] = glbGetSignalFitRatePtr(RE_EXP_FAR4,j);
      signal_fit_rates_F5[j] = glbGetSignalFitRatePtr(RE_EXP_FAR5,j);
      signal_fit_rates_F6[j] = glbGetSignalFitRatePtr(RE_EXP_FAR6,j);

      signal_fit_rates_N1[j] = glbGetSignalFitRatePtr(RE_EXP_NEAR1,j);
      signal_fit_rates_N2[j] = glbGetSignalFitRatePtr(RE_EXP_NEAR2,j);
      signal_fit_rates_N3[j] = glbGetSignalFitRatePtr(RE_EXP_NEAR3,j);

      signal_fit_rates_N4[j] = glbGetSignalFitRatePtr(RE_EXP_NEAR4,j);
      signal_fit_rates_N5[j] = glbGetSignalFitRatePtr(RE_EXP_NEAR5,j);
      signal_fit_rates_N6[j] = glbGetSignalFitRatePtr(RE_EXP_NEAR6,j);
    }

  int ew_low, ew_high;
  double emin, emax;
  double fit_rate,true_rate;
  double chi2 = 0.0;
  int i;

  /* Request simulated energy interval and analysis energy window */
  glbGetEminEmax(exp, &emin, &emax);
  glbGetEnergyWindowBins(exp, rule, &ew_low, &ew_high);

  signal_norm_F1 = 1.0 + x[0] + x[6];
  signal_norm_F2 = 1.0 + x[1] + x[6];
  signal_norm_F3 = 1.0 + x[2] + x[6];
  signal_norm_F4 = 1.0 + x[3] + x[6];
  signal_norm_F5 = 1.0 + x[4] + x[6];
  signal_norm_F6 = 1.0 + x[5] + x[6];

  signal_norm_N1 = 1.0 + x[0] + x[7]; 
  signal_norm_N2 = 1.0 + x[1] + x[7];
  signal_norm_N3 = 1.0 + x[2] + x[7];
  signal_norm_N4 = 1.0 + x[3] + x[7];
  signal_norm_N5 = 1.0 + x[4] + x[7];
  signal_norm_N6 = 1.0 + x[5] + x[7];

  /* Loop over all bins in energy window */
  for (i=0; i < n_bins; i++)
    {

      signal_F[i]=0.0;
      signal_N1[i]=0.0;
 
	  /* introducing the error on the isotope abundance, only
	     three since we assume that error on U235 is accounted for
	     by the normalization. This parametrization assumes that
	     the isotopic composition is uncorrelated between cores.
	  */

      j=0;
      signal_F[i]+=(1.0+x[8])*(signal_norm_F1*signal_fit_rates_F1[j][i]
			       +signal_norm_F2*signal_fit_rates_F2[j][i]
			       +signal_norm_F3*signal_fit_rates_F3[j][i]
			       +signal_norm_F4*signal_fit_rates_F4[j][i]
			       +signal_norm_F5*signal_fit_rates_F5[j][i]
			       +signal_norm_F6*signal_fit_rates_F6[j][i]);
      
      signal_N1[i]+=(1.0+x[8])*(signal_norm_N1*signal_fit_rates_N1[j][i]
				+signal_norm_N2*signal_fit_rates_N2[j][i]
				+signal_norm_N3*signal_fit_rates_N3[j][i]
				+signal_norm_N4*signal_fit_rates_N4[j][i]
				+signal_norm_N5*signal_fit_rates_N5[j][i]
				+signal_norm_N6*signal_fit_rates_N6[j][i]);
      
       for(j=1;j<4;j++) /* summing the events for the remaining 3 isotopes */ 
	{
	  signal_F[i]+=(1.0+x[8])
	    *(signal_norm_F1*signal_fit_rates_F1[j][i]*(1.0+x[10+j])
	      +signal_norm_F2*signal_fit_rates_F2[j][i]*(1.0+x[13+j])
	      +signal_norm_F3*signal_fit_rates_F3[j][i]*(1.0+x[16+j])
	      +signal_norm_F4*signal_fit_rates_F4[j][i]*(1.0+x[19+j])
	      +signal_norm_F5*signal_fit_rates_F5[j][i]*(1.0+x[22+j-1])
	      +signal_norm_F6*signal_fit_rates_F6[j][i]*(1.0+x[25+j-1]));
	  
	  signal_N1[i]+=(1.0+x[8])
	    *(signal_norm_N1*signal_fit_rates_N1[j][i]*(1.0+x[10+j])
	      +signal_norm_N2*signal_fit_rates_N2[j][i]*(1.0+x[13+j])
	      +signal_norm_N3*signal_fit_rates_N3[j][i]*(1.0+x[16+j])
	      +signal_norm_N4*signal_fit_rates_N4[j][i]*(1.0+x[19+j])
	      +signal_norm_N5*signal_fit_rates_N5[j][i]*(1.0+x[22+j])
	      +signal_norm_N6*signal_fit_rates_N6[j][i]*(1.0+x[25+j]));
	  
	}
    }

  /* energy scale error for each detector x[9]-x[10] */
  glbShiftEnergyScale(x[9],signal_F,signal_F, n_bins, emin, emax);
  glbShiftEnergyScale(x[10],signal_N1,signal_N1, n_bins, emin, emax);

  for (i=ew_low; i <= ew_high; i++)
  {
    /* chi^2 from far detector */
    fit_rate=signal_F[i]*(1.0+x[i+29]);
    true_rate=0.0;
    for(j=0;j<4;j++) /* summing the four isotopes */
      {
	true_rate+=true_rates_F1[j][i]+true_rates_F2[j][i]+true_rates_F3[j][i]
	  +true_rates_F4[j][i]+true_rates_F5[j][i]+true_rates_F6[j][i];
      }
    chi2 += likelihood( true_rate, fit_rate, 
			true_rate * (1.0 + true_rate
				     *square(re_sigma_binbin)) );

    /* chi^2 from near detector 1 */
    fit_rate=signal_N1[i]*(1.0+x[i+29]);
    true_rate=0.0;
    for(j=0;j<4;j++) /* summing the four isotopes */
      {
	true_rate+=true_rates_N1[j][i]+true_rates_N2[j][i]+true_rates_N3[j][i]
	  +true_rates_N4[j][i]+true_rates_N5[j][i]+true_rates_N6[j][i];
      }   
    chi2 += likelihood( true_rate, fit_rate, 
			true_rate * (1.0 + true_rate
				     *square(re_sigma_binbin)) );
  }

  /* Systematical part of chi^2 (= priors) */
  for (i=0; i < n_params; i++)
    chi2 += square(x[i] / errors[i]);

  /* Save the systematics parameters as starting values for the next step */
  for (i=0; i < n_params; i++)
    re_sys_startval[i] = x[i];

  return chi2;
}  


/* These functions need to be called in order to set up the reactor
 * experiments in the correct way.
 */


void init_reno(double *errors, double b2b, double shape)
{
  int i;
  
  glbDefineChiFunction(&chi_reno,29+62,"reno", NULL);
  for(i=0;i<12;i++)
    {
      glbInitExperiment("reno.glb", 
			&glb_experiment_list[0], 
			&glb_num_of_exps);
      
      glbSetBaselineInExperiment(i,reno_distances[i]);
      glbSetTargetMass(i,15);        
    }
  
  for (i=0; i < 29; i++) /* Normalization and energy calibration errors */
    { re_sys_errors[i] = errors[i];}
  for (i=29; i < 29+62; i++)
    re_sys_errors[i] = shape;  /* Spectral error */

  for(i=0;i<29+62;i++) re_sys_startval[i]=0;

  glbSetChiFunction(RE_EXP_FAR1, 0, GLB_ON, "reno", re_sys_errors);
  glbSetChiFunction(RE_EXP_FAR2, 0, GLB_ON, "chiZero", re_sys_errors);
  glbSetChiFunction(RE_EXP_FAR3, 0, GLB_ON, "chiZero", re_sys_errors);
  glbSetChiFunction(RE_EXP_FAR4, 0, GLB_ON, "chiZero", re_sys_errors);
  glbSetChiFunction(RE_EXP_FAR5, 0, GLB_ON, "chiZero", re_sys_errors);
  glbSetChiFunction(RE_EXP_FAR6, 0, GLB_ON, "chiZero", re_sys_errors);
  
  glbSetChiFunction(RE_EXP_NEAR1, 0, GLB_ON, "chiZero", re_sys_errors);
  glbSetChiFunction(RE_EXP_NEAR2, 0, GLB_ON, "chiZero", re_sys_errors);
  glbSetChiFunction(RE_EXP_NEAR3, 0, GLB_ON, "chiZero", re_sys_errors);
  glbSetChiFunction(RE_EXP_NEAR4, 0, GLB_ON, "chiZero", re_sys_errors);
  glbSetChiFunction(RE_EXP_NEAR5, 0, GLB_ON, "chiZero", re_sys_errors);
  glbSetChiFunction(RE_EXP_NEAR6, 0, GLB_ON, "chiZero", re_sys_errors);

  re_sigma_binbin=b2b;

  return;


}


#ifdef TEST





int main(int argc, char* argv[])
{
  

  glbInit(argv[0]);
 

  init_reno(&re_errors[0],0.003,0.02);  


  /* Set standard oscillation parameters arXiv: 0808.2016 */
  double theta12 = asin(sqrt(0.304));
  double theta13 = 0;
  double theta23 = M_PI/4;
  double deltacp = 0;
  double sdm = 7.65e-5;
  double ldm = 0.0024; 

  glb_params true_values  = glbAllocParams();
  glb_params test_values  = glbAllocParams();
  glb_params input_errors = glbAllocParams();

  /* Calculate "true" event rates */
  glbDefineParams(true_values,theta12,theta13,theta23,0,sdm,ldm);
  glbSetDensityParams(true_values,1.0,GLB_ALL);
  glbDefineParams(test_values,theta12,theta13,theta23,deltacp,sdm,ldm);  
  glbSetDensityParams(test_values,1.0,GLB_ALL);
  glbDefineParams(input_errors, 0.1*theta12, 0, 
		  0.15*theta23, 0, 0.05*sdm, 0.1*ldm);
  glbSetDensityParams(input_errors, 0.05, GLB_ALL);
  glbSetOscillationParameters(true_values);
  glbSetInputErrors(input_errors);
  glbSetRates();



  double x;
  for(x=-4;x<-1;x+=0.1)
    {
      double thetheta13, chi2;
      /* Set vector of test values */
      thetheta13 = asin(sqrt(pow(10,x)))/2;
      glbSetOscParams(test_values, thetheta13, GLB_THETA_13);
      
      /* Set starting values for systematics minimiyer to the coordinates of
       * minimum in the last iteration. This accelerates the minimization and
       * prevents convergence problems. */
      glbSetSysStartingValuesList(RE_EXP_FAR1, 0, GLB_ON, re_sys_startval);
     
      
       chi2 = glbChiSys(test_values, GLB_ALL, GLB_ALL);
       printf("%f\t%f\n",x,chi2);
    }
  return 0;
}

#endif /* TEST */
