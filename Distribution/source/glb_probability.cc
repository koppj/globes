/* GLoBES -- General LOng Baseline Experiment Simulator
 * (C) 2002 - 2005,  The GLoBES Team
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
#include <complex>
#include <iostream>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <globes/globes.h>
#include "glb_wrapper.h"
#include "glb_probability.h"
using namespace std;



// Standard operations on complex numbers
// --------------------------------------

inline gsl_complex operator+(gsl_complex a, gsl_complex b)         // r = a + b
{
  gsl_complex r;
  GSL_REAL(r) = GSL_REAL(a) + GSL_REAL(b);
  GSL_IMAG(r) = GSL_IMAG(a) + GSL_IMAG(b);
  return r;
}

inline gsl_complex operator-(gsl_complex a, gsl_complex b)         // r = a - b
{
  gsl_complex r;
  GSL_REAL(r) = GSL_REAL(a) - GSL_REAL(b);
  GSL_IMAG(r) = GSL_IMAG(a) - GSL_IMAG(b);
  return r;
}

inline gsl_complex operator*(gsl_complex a, gsl_complex b)         // r = a * b
{
  gsl_complex r;
  GSL_REAL(r) = GSL_REAL(a)*GSL_REAL(b) - GSL_IMAG(a)*GSL_IMAG(b);
  GSL_IMAG(r) = GSL_REAL(a)*GSL_IMAG(b) + GSL_IMAG(a)*GSL_REAL(b);
  return r;
}

inline gsl_complex operator/(gsl_complex a, gsl_complex b)         // r = a / b
{
  double s = 1.0 / (SQR(GSL_REAL(b)) + SQR(GSL_IMAG(b)));
  gsl_complex r;
  GSL_REAL(r) = s * (GSL_REAL(a)*GSL_REAL(b) + GSL_IMAG(a)*GSL_IMAG(b));
  GSL_IMAG(r) = s * (GSL_IMAG(a)*GSL_REAL(b) - GSL_REAL(a)*GSL_IMAG(b));
  return r;
}

inline gsl_complex operator-(gsl_complex a)                        // r = -a
{
  gsl_complex r;
  GSL_REAL(r) = -GSL_REAL(a);
  GSL_IMAG(r) = -GSL_IMAG(a);
  return r;
}

inline gsl_complex conj(gsl_complex a)                             // r = a*
{
  gsl_complex r;
  GSL_REAL(r) = GSL_REAL(a);
  GSL_IMAG(r) = -GSL_IMAG(a);
  return r;
}


// Standard operations mixing complex and real numbers
inline gsl_complex operator+(gsl_complex a, double x)              // r = a + x
{
  gsl_complex r;
  GSL_REAL(r) = GSL_REAL(a) + x;
  GSL_IMAG(r) = GSL_IMAG(a);
  return r;
}

inline gsl_complex operator+(double x, gsl_complex a)              // r = x + a
{
  gsl_complex r;
  GSL_REAL(r) = GSL_REAL(a) + x;
  GSL_IMAG(r) = GSL_IMAG(a);
  return r;
}

inline gsl_complex operator-(gsl_complex a, double x)              // r = a - x
{
  gsl_complex r;
  GSL_REAL(r) = GSL_REAL(a) - x;
  GSL_IMAG(r) = GSL_IMAG(a);
  return r;
}

inline gsl_complex operator-(double x, gsl_complex a)              // r = x - a
{
  gsl_complex r;
  GSL_REAL(r) = x - GSL_REAL(a);
  GSL_IMAG(r) = -GSL_IMAG(a);
  return r;
}

inline gsl_complex operator*(gsl_complex a, double x)              // r = a * x
{
  gsl_complex r;
  GSL_REAL(r) = x * GSL_REAL(a);
  GSL_IMAG(r) = x * GSL_IMAG(a);
  return r;
}

inline gsl_complex operator*(double x, gsl_complex a)              // r = x * a
{
  gsl_complex r;
  GSL_REAL(r) = x * GSL_REAL(a);
  GSL_IMAG(r) = x * GSL_IMAG(a);
  return r;
}

inline gsl_complex operator/(gsl_complex a, double x)              // r = a / x
{
  double inv_x = 1.0/x;
  gsl_complex r;
  GSL_REAL(r) = inv_x * GSL_REAL(a);
  GSL_IMAG(r) = inv_x * GSL_IMAG(a);
  return r;
}

inline gsl_complex operator/(double x, gsl_complex a)              // r = x / a
{
  double s = x / (SQR(GSL_REAL(a)) + SQR(GSL_IMAG(a)));
  gsl_complex r;
  GSL_REAL(r) = s * GSL_REAL(a);
  GSL_IMAG(r) = -s * GSL_IMAG(a);
  return r; 
}


class glbNuOscWorkspace
{
  public:
    glbNuOscWorkspace()
    {
      U = gsl_matrix_complex_calloc(GLB_NU_FLAVOURS, GLB_NU_FLAVOURS);
      H = gsl_matrix_complex_calloc(GLB_NU_FLAVOURS, GLB_NU_FLAVOURS);
      Q = gsl_matrix_complex_calloc(GLB_NU_FLAVOURS, GLB_NU_FLAVOURS);
      lambda = gsl_vector_alloc(GLB_NU_FLAVOURS);
      S = gsl_matrix_complex_calloc(GLB_NU_FLAVOURS, GLB_NU_FLAVOURS);
    
      H0_template = gsl_matrix_complex_calloc(GLB_NU_FLAVOURS, GLB_NU_FLAVOURS);
      S1 = gsl_matrix_complex_calloc(GLB_NU_FLAVOURS, GLB_NU_FLAVOURS);
      T0 = gsl_matrix_complex_calloc(GLB_NU_FLAVOURS, GLB_NU_FLAVOURS);
    }
    virtual ~glbNuOscWorkspace()
    {
      gsl_matrix_complex_free(T0);
      gsl_matrix_complex_free(S1);
      gsl_matrix_complex_free(H0_template);
      
      gsl_matrix_complex_free(S);
      gsl_vector_free(lambda);
      gsl_matrix_complex_free(Q);
      gsl_matrix_complex_free(H);
      gsl_matrix_complex_free(U);
    }

    gsl_matrix_complex *U;      // The vacuum mixing matrix
    gsl_matrix_complex *H;      // Neutrino Hamiltonian
    gsl_matrix_complex *Q;      // Eigenvectors of Hamiltonian (= eff. mixing matrix)
    gsl_vector *lambda;         // Eigenvalues of Hamiltonian
    gsl_matrix_complex *S;      // The neutrino S-matrix

    gsl_matrix_complex *H0_template;  // Used in the construction of the (vacuum) Hamiltonian
    gsl_matrix_complex *S1, *T0;// Temporary matrix storage
};

static glbNuOscWorkspace w;

// Fundamental oscillation parameters
static double th12, th13, th23; // Mixing angles
static double delta;            // Dirac CP phase
static double mq[3];            // Squared masses
static double nsparams[GLB_OSCP-6+1];  // Non-standard parameters



extern "C"{ 
#undef GLS_ERROR_H
#include "glb_error.h"
}
//interestingly enough this does *not* work 


#define MAX_ITER_STEFFENSON    50   // Max. number of iterations in Steffenson root finder

#define V_THRESHOLD 0.001*GLB_V_FACTOR*GLB_Ne_MANTLE  // The minimum matter potential below
                                                      // which vacuum algorithms are used

// ----------------------------------------------------------------------------
int glb_33_zheev(gsl_complex A[3][3], gsl_complex Q[3][3], double lambda[3])
// ----------------------------------------------------------------------------
// Calculates the eigenvalues and eigenvectors for a hermitian 3x3 matrix A.
// A is destroyed in the calculation. Eigenvectors are normalized to unity.
// Only the upper triangle of A is accessed.
// ----------------------------------------------------------------------------
// Parameters:
//   A:      The hermitian input matrix
//   Q:      Storage buffer for eigenvectors
//   lambda: Storage buffer for eigenvalues
// ----------------------------------------------------------------------------
{
#if GLB_NU_FLAVOURS != 3
#  error Probability engine can only handle three-flavour systems at the moment
#endif
  
  gsl_complex u[3], p[3], q[3];   // Temporary vectors for Householder transformation
  gsl_complex omega, K;           // Temporary scalars for Householder transformation
  double householder_sign;        // The sign in the Householder transformation
  double norm, inv_norm, norm0, norm1;
  gsl_complex v2[3], v3[3];       // Second and third eigenvector
  double a[3];
  
  // Determine coefficients of characteristic poynomial. We write
  //       | a   d   f  |
  //  A =  | d*  b   e  |
  //       | f*  e*  c  |
  gsl_complex de = A[0][1] * A[1][2];                               // d * e
  double d_abs = SQR_ABS(A[0][1]);                                  // d * conj(d)
  double e_abs = SQR_ABS(A[1][2]);                                  // e * conj(e)
  double f_abs = SQR_ABS(A[0][2]);                                  // f * conj(f)
  a[2] = GSL_REAL(A[0][0]) + GSL_REAL(A[1][1]) + GSL_REAL(A[2][2]); // a + b + c
  a[1] = d_abs + e_abs + f_abs  // d*conj(d) + e*conj(e) + f*conj(f) - a*b - a*c -b*c
                      - GSL_REAL(A[0][0])*GSL_REAL(A[1][1])
                      - GSL_REAL(A[0][0])*GSL_REAL(A[2][2])
                      - GSL_REAL(A[1][1])*GSL_REAL(A[2][2]);
  a[0] = GSL_REAL(A[0][0])*GSL_REAL(A[1][1])*GSL_REAL(A[2][2])
                     - GSL_REAL(A[2][2])*d_abs - GSL_REAL(A[0][0])*e_abs - GSL_REAL(A[1][1])*f_abs
                     + 2.0 * (GSL_REAL(A[0][2])*GSL_REAL(de) + GSL_IMAG(A[0][2])*GSL_IMAG(de));
                             // a*b*c - c*d*conj(d) - a*e*conj(e) - b*f*conj(f) + 2*Re(conj(f)*d*e)

  // Find first eigenvalue using the Steffenson algorithm to locate a root
  // of the characteristic polynomial. The choice of starting value is inspired
  // by the Gerschgorin theorem.
  double lambda_max = max(fabs(GSL_REAL(A[0][0])) + fabs(GSL_REAL(A[0][1])) + fabs(GSL_IMAG(A[0][1]))
                                                  + fabs(GSL_REAL(A[0][2])) + fabs(GSL_IMAG(A[0][2])),
                      max(fabs(GSL_REAL(A[1][1])) + fabs(GSL_REAL(A[0][1])) + fabs(GSL_IMAG(A[0][1]))
                                                  + fabs(GSL_REAL(A[1][2])) + fabs(GSL_IMAG(A[1][2])),
                          fabs(GSL_REAL(A[2][2])) + fabs(GSL_REAL(A[0][2])) + fabs(GSL_IMAG(A[0][2]))
                                                  + fabs(GSL_REAL(A[1][2])) + fabs(GSL_IMAG(A[1][2])) ));
  double x=lambda_max, x1=0.0, x2=0.0;
  double f, df;
  int n = 1;
  double lambda_old = x;
  do
  {
    double z=x;
    f   = a[0] + a[1]*z;        // Evaluate function and derivative at the new position
    df  = a[1] + 2.0*a[2]*z;
    z  *= x;
    f  += a[2]*z;
    df -= 3.0*z;
    z  *= x;
    f  -= z;

    if (df == 0.0)
    {
      if (f == 0.0)
      {
        lambda[0] = x;
        break;
      }
      else
      {
        cout << "Nu33zheev: Steffenson root finder encountered zero derivative at x = "
             << x << "." << endl;
        return GLBERR_NO_CONVERGENCE;
      }
    }
    x2  = x1;
    x1  = x;
    x   = x - f/df;             // Perform Newton step

    if (n >= 3)
    {
      double u = x - x1;        // Use Aitken Delta-squared acceleration
      double v = x - 2*x1 + x2;
      if (v != 0)               // Avoid division by zero
        lambda[0] = x - u*u/v;
      else
        lambda[0] = x;
    }
    else
      lambda[0] = x;

    if (fabs(lambda[0] - lambda_old) <= 1e-9 * lambda_max)   // Test for convergence
      break;
    lambda_old = lambda[0];
      
    if (++n > MAX_ITER_STEFFENSON)
    {
      cout << "Nu33zheev: Limit of " << MAX_ITER_STEFFENSON << " iterations exceeded." << endl;
      return GLBERR_NO_CONVERGENCE;
    }
  } while (1);

  // Calculate eigenvector for first eigenvalue by the formula
  //   v = conj( (H - lambda[0]).e1 x (H - lambda[0]).e2 )
  A[0][0] = A[0][0] - lambda[0];
  A[1][1] = A[1][1] - lambda[0];
  Q[0][0] = A[0][1]*A[1][2] - A[0][2]*GSL_REAL(A[1][1]);
  Q[1][0] = A[0][2]*conj(A[0][1]) - A[1][2]*GSL_REAL(A[0][0]);
  Q[2][0] = gsl_complex_rect(GSL_REAL(A[0][0])*GSL_REAL(A[1][1]) - d_abs, 0.0);
  norm  = SQR_ABS(Q[0][0]) + SQR_ABS(Q[1][0]) + SQR_ABS(Q[2][0]);
  norm0 = SQR_ABS(A[0][0]) + SQR_ABS(A[0][1]) + SQR_ABS(A[0][2]);
  norm1 = SQR_ABS(A[0][1]) + SQR_ABS(A[1][1]) + SQR_ABS(A[1][2]);

  if (norm0 < 1e-15*SQR(lambda_max))
    norm = norm0 = 0.0;
  else if (norm1 < 1e-15*SQR(lambda_max))
    norm = norm1 = 0.0;

  if (norm != 0.0)
  {
    inv_norm = 1.0 / norm;
    if (norm0 * norm1 * inv_norm < 1e15) // Accept cross product only if the angle between
    {                                    // the two vectors was not too small
      inv_norm = sqrt(inv_norm);
      for (int i=0; i < 3; i++)
        Q[i][0] = Q[i][0] * inv_norm;
    }
    else                                 // ... but then (1, -A0/A1, 0) is an ev.
    {
      double max_A = -1.0;
      int i_max = 0;
      gsl_complex factor;
      for (int i=0; i < 3; i++)
        if (SQR_ABS(A[0][i]) > max_A)
        {
          max_A = SQR_ABS(A[0][i]);
          i_max = i;
        }
      if (i_max > 1)
        factor = -conj(A[0][i_max]/A[1][i_max]);
      else
        factor = -conj(A[0][i_max])/A[i_max][1];
      inv_norm = 1.0/sqrt(1 + SQR_ABS(factor));
      Q[0][0] = gsl_complex_rect(inv_norm, 0.0);
      Q[1][0] = factor * inv_norm;
      Q[2][0] = GSL_COMPLEX_ZERO;
    }
  }
  else if (norm0 == 0)                   // If the first column was zero, then (1,0,0) is an eigenvector
  {
    Q[0][0] = GSL_COMPLEX_ONE;
    Q[1][0] = GSL_COMPLEX_ZERO;
    Q[2][0] = GSL_COMPLEX_ZERO;
  }
  else if (norm1 == 0)                   // If the second column was zero, then (0,1,0) is an eigenvector
  {
    Q[0][0] = GSL_COMPLEX_ZERO;
    Q[1][0] = GSL_COMPLEX_ONE;
    Q[2][0] = GSL_COMPLEX_ZERO;
  }
  else                                   // If the first and second columns are linearly dependent,
  {                                      // then (1, -A0/A1, 0) is an eigenvector
    double max_A = -1.0;
    int i_max = 0;
    gsl_complex factor;
    for (int i=0; i < 3; i++)
      if (SQR_ABS(A[0][i]) > max_A)
      {
        max_A = SQR_ABS(A[0][i]);
        i_max = i;
      }
    if (i_max > 1)
      factor = -conj(A[0][i_max]/A[1][i_max]);
    else
      factor = -conj(A[0][i_max])/A[i_max][1];
    inv_norm = 1.0/sqrt(1 + SQR_ABS(factor));
    Q[0][0] = gsl_complex_rect(inv_norm, 0.0);
    Q[1][0] = factor * inv_norm;
    Q[2][0] = GSL_COMPLEX_ZERO;
  }

  // Prepare Householder transformation
  householder_sign = GSL_REAL(Q[0][0]) < 0 ? -1.0 : 1.0; // Choose sign such that roundoff error is minimal
  A[0][0] = A[0][0] + lambda[0];
  A[1][1] = A[1][1] + lambda[0];
  u[0] = Q[0][0] + householder_sign;
  u[1] = Q[1][0];
  u[2] = Q[2][0];
  omega = 1.0 + householder_sign*conj(Q[0][0]);          // = v1^\dag.u as |v1| = 1.0
  omega = (1.0 + omega/conj(omega)) / (SQR_ABS(u[0]) + SQR_ABS(u[1]) + SQR_ABS(u[2]));

  // Perform Householder transform on H. Note that we omit the calculation of the upper
  // triangle and of the complete first column.
  K = GSL_COMPLEX_ZERO;
  for (int i=0; i < 3; i++)
  {
    p[i] = conj(omega) * (lambda[0]*Q[i][0] + householder_sign*conj(A[0][i]));
    K    = K + conj(u[i]) * p[i];
  }
  K = K * 0.5*omega;
  for (int i=1; i < 3; i++)
    q[i] = p[i] - GSL_REAL(K) * u[i];
  for (int i=1; i < 3; i++)
    for (int j=i; j < 3; j++)
      A[i][j] = A[i][j] - q[i]*conj(u[j]) - u[i]*conj(q[j]);

  // Calculate eigenvalues of remaining 2x2 matrix
  double D = 0.5 * sqrt(fabs( SQR(GSL_REAL(A[1][1]) + GSL_REAL(A[2][2]))
                           - 4*( GSL_REAL(A[1][1])*GSL_REAL(A[2][2]) - SQR_ABS(A[1][2]) )));
  lambda[1]  = lambda[2] = 0.5 * (GSL_REAL(A[1][1]) + GSL_REAL(A[2][2]));
  lambda[1] += D;
  lambda[2] -= D;

  // Calculate eigenvectors of remaining 2x2 matrix
  A[1][1] = A[1][1] - lambda[1];
  A[2][2] = A[2][2] - lambda[1];
  norm1 = norm0 = SQR_ABS(A[1][2]);
  norm0 += SQR(GSL_REAL(A[1][1]));
  norm1 += SQR(GSL_REAL(A[2][2]));
  v2[0] = GSL_COMPLEX_ZERO;
  if (norm0 > norm1)
  {
    inv_norm = sqrt(1.0 / norm0);
    v2[1] = -A[1][2] * inv_norm;
    v2[2] = A[1][1] * inv_norm;
  }
  else
  {
    if (norm1 != 0.0)
    {
      inv_norm = sqrt(1.0 / norm1);
      v2[1] = -A[2][2] * inv_norm;
      v2[2] = conj(A[1][2]) * inv_norm;
    }
    else             // This branch will only be reached if norm0 == norm1 == 0.0
    {
      v2[1] = GSL_COMPLEX_ONE;
      v2[2] = GSL_COMPLEX_ZERO;
    }
  }
  v3[0] = GSL_COMPLEX_ZERO;
  v3[1] = conj(v2[2]);
  v3[2] = -conj(v2[1]);

  // Transform eigenvectors back by multiplying with hermitian conjugate of first Householder matrix
  // P^\dagger = 1 - conj(omega) u u^\dagger
  gsl_complex r = conj(omega) * (conj(u[1])*v2[1] + conj(u[2])*v2[2]);
  for (int i=0; i < 3; i++)
    Q[i][1] = v2[i] - r * u[i];
  r = conj(omega) * (conj(u[1])*v3[1] + conj(u[2])*v3[2]);
  for (int i=0; i < 3; i++)
    Q[i][2] = v3[i] - r * u[i];
  
  return GLB_SUCCESS;
}


// ---------------------------------------------------------------------------- 
extern "C" int glb_set_oscillation_parameters(glb_params p)
// ---------------------------------------------------------------------------- 
// Sets the fundamental oscillation parameters and precomputes the mixing
// matrix and part of the Hamiltonian
// ----------------------------------------------------------------------------
// Parameters:
//   p: Pointer to data structure containing the new oscillation parameters
// ----------------------------------------------------------------------------
{
  gsl_complex (*U)[GLB_NU_FLAVOURS]
    = (gsl_complex (*)[GLB_NU_FLAVOURS]) gsl_matrix_complex_ptr(w.U, 0, 0);

  // Copy parameters
  th12  = p->osc->osc_params[0];
  th13  = p->osc->osc_params[1];
  th23  = p->osc->osc_params[2];
  delta = p->osc->osc_params[3];
  mq[0] = abs(p->osc->osc_params[5]);
  mq[1] = abs(p->osc->osc_params[5]) + p->osc->osc_params[4];
  mq[2] = abs(p->osc->osc_params[5]) + p->osc->osc_params[5];
   
  // Compute vacuum mixing matrix
  U[0][0] = gsl_complex_rect(cos(th12)*cos(th13), 0.0);
  U[0][1] = gsl_complex_rect(sin(th12)*cos(th13), 0.0);
  U[0][2] = gsl_complex_polar(sin(th13), -delta);

  U[1][0] = -sin(th12)*cos(th23) - gsl_complex_polar(cos(th12)*sin(th23)*sin(th13), delta);
  U[1][1] =  cos(th12)*cos(th23) - gsl_complex_polar(sin(th12)*sin(th23)*sin(th13), delta);
  U[1][2] = gsl_complex_rect(sin(th23)*cos(th13), 0.0);

  U[2][0] =  sin(th12)*sin(th23) - gsl_complex_polar(cos(th12)*cos(th23)*sin(th13), delta);
  U[2][1] = -cos(th12)*sin(th23) - gsl_complex_polar(sin(th12)*cos(th23)*sin(th13), delta);
  U[2][2] = gsl_complex_rect(cos(th23)*cos(th13), 0.0);
 
  // Calculate energy independent matrix H0 * E
  gsl_matrix_complex_set_zero(w.H0_template);
  gsl_matrix_complex_set_zero(w.H);
  for (int i=0; i < GLB_NU_FLAVOURS; i++)
    gsl_matrix_complex_set(w.H0_template, i, i, gsl_complex_rect(0.5*mq[i], 0.0));

  gsl_matrix_complex *T = gsl_matrix_complex_alloc(GLB_NU_FLAVOURS, GLB_NU_FLAVOURS);
  gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, GSL_COMPLEX_ONE, w.H0_template, w.U, // T  = H0.U^\dagger
                 GSL_COMPLEX_ZERO, T);
  gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, w.U, T,               // H0 = U.T
                 GSL_COMPLEX_ZERO, w.H0_template);
  gsl_matrix_complex_free(T);
  return GLB_SUCCESS;
}


// ---------------------------------------------------------------------------- 
extern "C" int glb_get_oscillation_parameters(glb_params p)
// ---------------------------------------------------------------------------- 
// Retrieves the fundamental oscillation parameters
// ----------------------------------------------------------------------------
// Parameters:
//   p: Pointer to data structure where the oscillation parameters are stored
// ----------------------------------------------------------------------------
{
  glbDefineParams(p, th12, th13, th23, delta, mq[1] - mq[0], mq[2] - mq[0]);
  return GLB_SUCCESS;
}


// ----------------------------------------------------------------------------
static int glb_hamiltonian_cd(double E, double V, int cp_sign)
// ----------------------------------------------------------------------------
// Calculates the Hamiltonian for neutrinos or antineutrinos with energy E
// propagating in matter of density V and stores the result in w.H.
// ----------------------------------------------------------------------------
// Parameters:
//   E: Neutrino energy
//   V: Matter potential (must be > 0 even for antineutrinos, the minus sign
//      for antineutrinos is incorporated automatically by the function)
//   cp_sign: +1 for neutrinos, -1 for antineutrinos
// ----------------------------------------------------------------------------
{
  double inv_E = 1.0 / E;
  gsl_complex (*H)[GLB_NU_FLAVOURS]
    = (gsl_complex (*)[GLB_NU_FLAVOURS]) gsl_matrix_complex_ptr(w.H, 0, 0);
  gsl_complex (*H0_template)[GLB_NU_FLAVOURS]
    = (gsl_complex (*)[GLB_NU_FLAVOURS]) gsl_matrix_complex_ptr(w.H0_template, 0, 0);

  if (cp_sign > 0)
  {
    for (int i=0; i < GLB_NU_FLAVOURS; i++)
      for (int j=0; j < GLB_NU_FLAVOURS; j++)
        H[i][j] = H0_template[i][j] * inv_E;
  }
  else
  {
    for (int i=0; i < GLB_NU_FLAVOURS; i++)
      for (int j=0; j < GLB_NU_FLAVOURS; j++)
        H[i][j] = conj(H0_template[i][j] * inv_E);   // delta_CP -> -delta_CP
  }
 
  H[0][0] = H[0][0] + cp_sign*V;
  return GLB_SUCCESS;
}


// ----------------------------------------------------------------------------
int glb_S_matrix_cd(double E, double L, double V, int cp_sign)
// ----------------------------------------------------------------------------
// Calculates the S matrix for neutrino oscillations in matter of constant
// density using a fast eigenvalue solver optimized to 3x3 matrices
// ----------------------------------------------------------------------------
// Parameters:
//   E: Neutrino energy
//   L: Baseline
//   V: Matter potential (must be > 0 even for antineutrinos)
//   cp_sign: +1 for neutrinos, -1 for antineutrinos
// ----------------------------------------------------------------------------
{
  // Introduce some abbreviations
  gsl_complex (*S)[3]  = (gsl_complex (*)[3]) gsl_matrix_complex_ptr(w.S,0,0);
  gsl_complex (*Q)[3]  = (gsl_complex (*)[3]) gsl_matrix_complex_ptr(w.Q,0,0);
  gsl_complex (*T0)[3] = (gsl_complex (*)[3]) gsl_matrix_complex_ptr(w.T0,0,0);
  double *lambda = gsl_vector_ptr(w.lambda,0);
  int status;
  
  if (V < V_THRESHOLD)                             // Vacuum
  {
    // Use vacuum mixing angles and masses
    double inv_E = 0.5/E;
    for (int i=0; i < GLB_NU_FLAVOURS; i++)
      lambda[i] = mq[i] * inv_E;

    if (cp_sign > 0)
      gsl_matrix_complex_memcpy(w.Q, w.U);
    else
    {
      gsl_complex (*U)[3]  = (gsl_complex (*)[3]) gsl_matrix_complex_ptr(w.U,0,0);
      for (int i=0; i < GLB_NU_FLAVOURS; i++)
        for (int j=0; j < GLB_NU_FLAVOURS; j++)
          Q[i][j] = conj(U[i][j]);
    }
  }
  else                                             // Matter
  {
    gsl_complex (*H)[3] = (gsl_complex (*)[3]) gsl_matrix_complex_ptr(w.H,0,0);
    
    // Calculate neutrino Hamiltonian
    if ((status=glb_hamiltonian_cd(E, V, cp_sign)) != GLB_SUCCESS)
      return status;
    
    // Calculate eigenvalues of Hamiltonian
    if ((status=glb_33_zheev(H, Q, lambda)) != GLB_SUCCESS)
      return status;
  }

/*  cout << lambda[0] << "  " << lambda[1] << "  " << lambda[2] << endl; 
  for (int i=0; i < GLB_NU_FLAVOURS; i++)
  {
    for (int j=0; j < GLB_NU_FLAVOURS; j++)
      printf("%g + %g i     ", GSL_REAL(Q[i][j]), GSL_IMAG(Q[i][j]));
    printf("\n");
  }
  getchar();*/
  
  // Calculate S-Matrix in mass basis in matter ...
  double phase;
  gsl_matrix_complex_set_zero(w.S);
  for (int i=0; i < GLB_NU_FLAVOURS; i++)
  {
    phase = -L * lambda[i];
    GSL_REAL(S[i][i]) = cos(phase);
    GSL_IMAG(S[i][i]) = sin(phase);
  } 
  
  // ... and transform it to the flavour basis
  gsl_matrix_complex_set_zero(w.T0);
  gsl_complex *p = &T0[0][0];
  for (int i=0; i < GLB_NU_FLAVOURS; i++)          // T0 = S.Q^\dagger
    for (int j=0; j < GLB_NU_FLAVOURS; j++)
    {
      for (int k=0; k < GLB_NU_FLAVOURS; k++)
      {
        GSL_REAL(*p) += GSL_REAL(S[i][k])*GSL_REAL(Q[j][k]) + GSL_IMAG(S[i][k])*GSL_IMAG(Q[j][k]);
        GSL_IMAG(*p) += GSL_IMAG(S[i][k])*GSL_REAL(Q[j][k]) - GSL_REAL(S[i][k])*GSL_IMAG(Q[j][k]);
      }
      p++;
    }
  gsl_matrix_complex_set_zero(w.S);
  p = &S[0][0];
  for (int i=0; i < GLB_NU_FLAVOURS; i++)          // S = Q.T0
    for (int j=0; j < GLB_NU_FLAVOURS; j++)
    {
      for (int k=0; k < GLB_NU_FLAVOURS; k++)
      {
        GSL_REAL(*p) += GSL_REAL(Q[i][k])*GSL_REAL(T0[k][j]) - GSL_IMAG(Q[i][k])*GSL_IMAG(T0[k][j]);
        GSL_IMAG(*p) += GSL_IMAG(Q[i][k])*GSL_REAL(T0[k][j]) + GSL_REAL(Q[i][k])*GSL_IMAG(T0[k][j]);
      }
      p++;
    }

  return GLB_SUCCESS;
}


// ----------------------------------------------------------------------------
int glb_filtered_probability_matrix_cd(double P[3][3], double E, double L, double V,
    double sigma, int cp_sign)
// ----------------------------------------------------------------------------
// Calculates the probability matrix for neutrino oscillations in matter of
// constant density, including a low pass filter to suppress aliasing due to
// very fast oscillations.
// ----------------------------------------------------------------------------
// Parameters:
//   P: Storage buffer for the probability matrix
//   E: Neutrino energy
//   L: Baseline
//   V: Matter potential (must be > 0 even for antineutrinos)
//   sigma: Width of Gaussian filter
//   cp_sign: +1 for neutrinos, -1 for antineutrinos
// ----------------------------------------------------------------------------
{
  // Introduce some abbreviations
  gsl_complex (*Q)[3]  = (gsl_complex (*)[3]) gsl_matrix_complex_ptr(w.Q,0,0);
  gsl_complex (*T0)[3] = (gsl_complex (*)[3]) gsl_matrix_complex_ptr(w.T0,0,0);
  double *lambda = gsl_vector_ptr(w.lambda,0);
  int status;
 
  if (V < V_THRESHOLD)                             // Vacuum
  {
    // Use vacuum mixing angles and masses
    double inv_E = 0.5/E;
    for (int i=0; i < GLB_NU_FLAVOURS; i++)
      lambda[i] = mq[i] * inv_E;

    if (cp_sign > 0)
      gsl_matrix_complex_memcpy(w.Q, w.U);
    else
    {
      gsl_complex (*U)[3]  = (gsl_complex (*)[3]) gsl_matrix_complex_ptr(w.U,0,0);
      for (int i=0; i < GLB_NU_FLAVOURS; i++)
        for (int j=0; j < GLB_NU_FLAVOURS; j++)
          Q[i][j] = conj(U[i][j]);
    }
  }
  else                                             // Matter
  {
    gsl_complex (*H)[3] = (gsl_complex (*)[3]) gsl_matrix_complex_ptr(w.H,0,0);
    
    // Calculate neutrino Hamiltonian
    if ((status=glb_hamiltonian_cd(E, V, cp_sign)) != GLB_SUCCESS)
      return status;
    
    // Calculate eigenvalues of Hamiltonian
    if ((status=glb_33_zheev(H, Q, lambda)) != GLB_SUCCESS)
      return status;
  }

  // Calculate probability matrix (see GLoBES manual for a discussion of the algorithm)
  double phase, filter_factor;
  gsl_matrix_complex_set_zero(w.T0);
  for (int i=0; i < GLB_NU_FLAVOURS; i++)
    for (int j=i+1; j < GLB_NU_FLAVOURS; j++)
    {
      phase              = -L * (lambda[i] - lambda[j]);
      filter_factor      = exp(-0.5 * SQR(phase*sigma) / SQR(1.0e-9 * E));
      GSL_REAL(T0[i][j]) = cos(phase) * filter_factor;
      GSL_IMAG(T0[i][j]) = sin(phase) * filter_factor;
    }

  for (int k=0; k < GLB_NU_FLAVOURS; k++)
    for (int l=0; l < GLB_NU_FLAVOURS; l++)
    {
      P[k][l] = 0.0;
      for (int i=0; i < GLB_NU_FLAVOURS; i++)
      {
        for (int j=i+1; j < GLB_NU_FLAVOURS; j++)
          P[k][l] += 2 * GSL_REAL(Q[k][j] * conj(Q[l][j]) * conj(Q[k][i]) * Q[l][i] * T0[i][j]);
        P[k][l] += SQR_ABS(Q[k][i]) * SQR_ABS(Q[l][i]);
      }
    }
    
  return GLB_SUCCESS;
}


// ----------------------------------------------------------------------------
extern "C" int glb_probability_matrix(double P[3][3], int cp_sign, double E,
    int psteps, const double *length, const double *density,
    double filter_sigma)
// ----------------------------------------------------------------------------
// Calculates the neutrino oscillation probability matrix
// ----------------------------------------------------------------------------
// Parameters:
//   P: Buffer for the storage of the matrix
//   cp_sign: +1 for neutrinos, -1 for antineutrinos
//   E: Neutrino energy (in GeV)
//   psteps: Number of layers in the matter density profile
//   length: Lengths of the layers in the matter density profile
//   density: The matter densities
//   filter_sigma: Width of low-pass filter or <0 for no filter
// ----------------------------------------------------------------------------
{
  int status;

  // Convert energy to eV
  E *= 1.0e9;
//E *= 1.0e6;
  
//  cout << psteps << "  " << length[0] << "  " << density[0] << "  " << filter_sigma << endl;
  if (filter_sigma > 0.0)                     // With low-pass filter
  {
    if (psteps == 1)
      glb_filtered_probability_matrix_cd(P, E, KM_TO_EV(length[0]), density[0]*GLB_V_FACTOR*GLB_Ne_MANTLE,
          filter_sigma, cp_sign);
    else
      return GLBERR_NOT_IMPLEMENTED;
  }
  else                                        // Without low-pass filter
  {
    if (psteps > 1)
    {
      gsl_matrix_complex_set_identity(w.S1);                                   // S1 = 1
      for (int i=0; i < psteps; i++)
      {
        status = glb_S_matrix_cd(E, KM_TO_EV(length[i]), density[i]*GLB_V_FACTOR*GLB_Ne_MANTLE, cp_sign);
        if (status != GLB_SUCCESS)
          return status;
        gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, w.S, w.S1, // T0 = S.S1
                       GSL_COMPLEX_ZERO, w.T0);
        gsl_matrix_complex_memcpy(w.S1, w.T0);                                 // S1 = T0
      } 
      gsl_matrix_complex_memcpy(w.S, w.S1);                                    // S = S1
    }
    else
    {
      status = glb_S_matrix_cd(E, KM_TO_EV(length[0]), density[0]*GLB_V_FACTOR*GLB_Ne_MANTLE, cp_sign);
      if (status != GLB_SUCCESS)
        return status;
    }

    gsl_complex (*S)[3] = (gsl_complex (*)[3]) gsl_matrix_complex_ptr(w.S,0,0);
    for (int i=0; i < GLB_NU_FLAVOURS; i++)
      for (int j=0; j < GLB_NU_FLAVOURS; j++)
        P[j][i] = SQR_ABS(S[i][j]);
  }

/*  printf("Probabilities:\n");
  for (int i=0; i < GLB_NU_FLAVOURS; i++)
  {
    for (int j=0; j < GLB_NU_FLAVOURS; j++)
      printf("%g    ", P[i][j]);
    printf("\n");
  }
  getchar();*/

  return GLB_SUCCESS;
}





