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




#ifndef GLB_CALC_H
#define GLB_CALC_H 1

#include "glb_types.h"

int glb_num_of_rules;
double glb_glb_sig_norm_error[32];
double glb_sig_tilt_error[32];
double glb_bg_norm_error[32];
double glb_bg_tilt_error[32];
double glb_bg_norm_center[32];
double glb_bg_tilt_center[32];

double glb_tre_null_center[32];
double glb_tre_tilt_center[32];
double glb_tre_null_error[32];
double glb_tre_tilt_error[32];

int glb_calc_simbins;
double glb_calc_simtresh;
double glb_calc_simbeam;
double glb_calc_energy_window[32][2];


void glb_set_baseline(double l);
double glb_check_baseline();
void glb_set_target_mass(double mass);
double glb_check_target_mass();
void glb_set_years(double year);
double glb_check_years();
void glb_set_energy_treshold(double etresh);
double glb_check_energy_treshold();
void glb_set_number_of_bins(int nr);
int glb_check_number_of_bins();
int glb_get_number_of_bins();
void glb_set_type(int type);
int glb_get_type();

double glb_list_likelihood(double* ratest, double* ratesm);


void glb_set_channel(int i, int polarity, int anti, int l, int m, int cc,
		int energysmear);
void glb_set_num_of_channels(int cc);
double glb_calc_channel(int i, double en, double baseline);
void glb_set_number_of_rules(int cn);
void glb_set_rule(int i, int cn, int *rule, double *coeff);
void glb_set_bg_rule(int i, int cn, int *rule, double *coeff);


void glb_set_rates();
void glb_set_new_rates();
void glb_set_signal_errors(int i,double norm, double tilt);
void glb_set_bg_errors(int i,double norm, double tilt);
void glb_set_bg_center(int i,double norm, double tilt);
double glb_chi_sys_w_bg(double x[5]);




double glb_prior(double x, double center, double sigma);




double* glb_get_channel(int i);
int glb_calc_check_num_of_channels();
int glb_calc_glb_calc_check_rule_length(int i,int bg);
double* glb_calc_glb_calc_check_rule_coeff(int i,int bg);
int* glb_calc_check_rule(int i,int bg);
int glb_calc_check_num_of_rules();
double* glb_calc_check_signal_errors(int i);
double* glb_calc_check_bg_errors(int i);
double* glb_calc_check_bg_center(int i);
double glb_chi_sys_w_bgtot(double x[5]);
void glb_calc_set_tresh_center(int i, double a, double b);
void glb_calc_set_tresh_errors(int i, double a, double b);
void glb_calc_check_tresh_center(int i,double *in);
void glb_calc_check_tresh_errors(int i,double *in);
double glb_chi_sys_w_bg2(double x[7]);

void glb_set_error_function(int typ);
int glb_check_error_function();

void glb_remove_calc_pointers();
double glb_chi_sys_w_bgtot2(double x[5]);
double glb_chi_split(double x[6]);
double glb_chi_spec(double x[2]);

void glb_calc_set_energy_window(int i, double a, double b);
void check_glb_calc_energy_window(int i,double *in);
double glb_window_function(double low,double up,int bin);
double glb_chi_sys_w_bg_calib(double x[5]);
void glb_shift_energy_scale(double g,double* ratesin, double* ratesout);

double* glb_chirate;
double** glb_calc_smearing[32];
int* glb_calc_uprange[32];
int* glb_calc_lowrange[32];
double* glb_calc_signal_rates[32];
double* glb_calc_bg_rates[32]; 
double* glb_calc_rates_0[32];
double* glb_calc_rates_1[32];
double* glb_calc_rates_1T[32];
double* glb_calc_rates_1BG[32];
double* glb_calc_rates_1BGT[32];
double* glb_calc_energy_tab;
double* glb_calc_ratevec[32];
double* glb_calc_glb_calc_ratevec_bg[32];

glb_smear *glb_calc_smear_data[32];

glb_flux *glb_calc_fluxes[32];
glb_xsec *glb_calc_xsecs[32];

double *glb_calc_buffer;
double *glb_calc_chrb[32];
double *glb_calc_chra[32];

double *glb_calc_user_pre_sm_channel[32];
double *glb_calc_user_post_sm_channel[32];
double *glb_calc_user_pre_sm_background[32];
double *glb_calc_user_post_sm_background[32];

double* glb_calc_efficiencies[32];
double* glb_calc_const_background[32];


#endif /* GLB_CALC_H */