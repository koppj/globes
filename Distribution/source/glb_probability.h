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





#ifndef GLS_OSZPROB_H
#define GLS_OSZPROB_H 1


double glbVacuumProbability(int pl, int pm, int panti, double pen, double plength);
double glbProfileProbability(int pl, int pm, int panti, double pen);
void glb_probability_matrix(double prob[3][3], int panti, double pen);
void glb_set_parameters(double dm21, double dm31, double pt12, double pt13, double pt23, double pdelta1);
void glb_get_profile_data(double** lval, double** rval, int* n);
double* glb_get_vacuum_parameters();
double* glb_get_squared_masses();

void glb_set_c_squared_masses(double mq1, double mq2, double mq3);
void glb_set_c_vacuum_parameters(double pt12, double pt13, double pt23, double pdelta1);

 double *glb_get_density_ptr();
 double *glb_get_length_ptr();
 int glb_get_psteps();
 void glb_set_density_ptr(double *d);
 void glb_set_length_ptr(double *l);
 int glb_set_psteps(int p);


double glb_filtered_vac_probability(int pl, int pm, int panti, double pen, double plength);
double glb_filtered_const_density_probability(int pl, int pm, int panti,double pen);

void glb_switch_filter(int on);
int glb_check_filter();
void glb_set_filter(double x);
double glb_get_filter();


#endif /* GLS_OSZPROB_H */
