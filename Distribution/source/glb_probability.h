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

#include <gsl/gsl_complex.h>

//#define GLB_V_FACTOR        7.56e-14   // Conversion factor for matter potentials
#define GLB_V_FACTOR        7.5e-14   // Conversion factor for matter potentials
//#define GLB_Ne_MANTLE       0.497      // Effective electron numbers for calculation
#define GLB_Ne_MANTLE       0.5
#define GLB_Ne_CORE         0.468      //   of MSW potentials

BEGIN_C_DECLS

int glb_33_zheev(gsl_complex A[3][3], gsl_complex Q[3][3], double lambda[3]);
double glb_profile_probability(int pl, int pm, int panti, double pen);
int glb_probability_matrix(double P[3][3], int cp_sign, double E,
    int psteps, const double *length, const double *density,
    double filter_sigma);
int glb_filtered_probability_matrix_cd(double E, double L, double V, double sigma, int cp_sign);
int glb_set_oscillation_parameters(glb_params p);
int glb_get_oscillation_parameters(glb_params p);


// Macros
// ------
#define SQR(x)      ((x)*(x))                              // x^2 
#define SQR_ABS(x)  (SQR(GSL_REAL(x)) + SQR(GSL_IMAG(x)))  // |x|^2

END_C_DECLS

#endif /* GLS_OSZPROB_H */


