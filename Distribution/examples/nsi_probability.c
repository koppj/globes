/***************************************************************************
 * Probability engine for simulating non-standard interactions             *
 ***************************************************************************
 * Author: Joachim Kopp, Toshihiko Ota                                     *
 ***************************************************************************
 * Note:                                                                   *
 *  - This file is written in C99, which has data types for comples        *
 *    number. To compile it in gcc, use the command line option -std=gnu99 *
 *  - When you use this code to produce a scientific publication, please   *
 *    cite the following references:                                       *
 *                                                                         *
 *      @Article{Kopp:2006wp,                                              *
 *        author    = "Kopp, Joachim",                                     *
 *        title     = "{Efficient numerical diagonalization of hermitian   *
 *                     $3 \times 3$ matrices}",                            *
 *        journal   = "Int. J. Mod. Phys.",                                *
 *        volume    = "C19",                                               *
 *        year      = "2008",                                              *
 *        pages     = "523-548",                                           *
 *        eprint    = "physics/0610206",                                   *
 *        archivePrefix = "arXiv",                                         *
 *        doi       = "10.1142/S0129183108012303",                         *
 *        SLACcitation  = "%%CITATION = PHYSICS/0610206;%%",               *
 *        note      = "Erratum ibid.\ {\bf C19} (2008) 845",               *
 *        memo      = "Algorithms for fast diagonalization of 3x3 matrices"*
 *      }                                                                  *
 *                                                                         *
 *      @Article{Kopp:2007ne,                                              *
 *        author    = "Kopp, Joachim and Lindner, Manfred and Ota,         *
 *                     Toshihiko and Sato, Joe",                           *
 *        title     = "{Non-standard neutrino interactions in reactor and  *
 *                     superbeam experiments}",                            *
 *        journal   = "Phys. Rev.",                                        *
 *        volume    = "D77",                                               *
 *        year      = "2008",                                              *
 *        pages     = "013007",                                            *
 *        eprint    = "0708.0152",                                         *
 *        archivePrefix = "arXiv",                                         *
 *        primaryClass  =  "hep-ph",                                       *
 *        doi       = "10.1103/PhysRevD.77.013007",                        *
 *        SLACcitation  = "%%CITATION = 0708.0152;%%",                     *
 *        memo      = "NSI oscillation engine developed in this work",     *
 *      }                                                                  *
 ***************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <complex.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include "globes/globes.h"
#include "nsi_probability.h"

/* Constants */
#define GLB_V_FACTOR        7.5e-14   /* Conversion factor for matter potentials */
#define GLB_Ne_MANTLE       0.5        /* Effective electron numbers for calculation */
#define GLB_Ne_CORE         0.468      /*   of MSW potentials                        */
#define V_THRESHOLD 0.001*GLB_V_FACTOR*GLB_Ne_MANTLE  /* The minimum matter potential below */
                                                      /* which vacuum algorithms are used   */
#define M_SQRT3  1.73205080756887729352744634151      /* sqrt(3) */

/* Macros */
#define SQR(x)      ((x)*(x))                        /* x^2   */
#define SQR_ABS(x)  (SQR(creal(x)) + SQR(cimag(x)))  /* |x|^2 */


/* Fundamental oscillation parameters */
static double th12, th13, th23; // Mixing angles
static double delta;            // Dirac CP phase
static double mq[3];            // Squared masses
static double complex epsilon_s_plus_1[3][3]; // NSI in the source
static double complex epsilon_m[3][3]; // NSI in the propagation
static double complex epsilon_d_plus_1[3][3]; // NSI in the detector

/* Names of NSI parameters */
const char *param_strings[] =
  { "TH12", "TH13", "TH23", "DELTA_CP", "DM21", "DM31",
    "ABS_EPS_S_EE",     "ARG_EPS_S_EE",
    "ABS_EPS_S_MUE",    "ARG_EPS_S_MUE",
    "ABS_EPS_S_TAUE",   "ARG_EPS_S_TAUE",
    "ABS_EPS_S_EMU",    "ARG_EPS_S_EMU",
    "ABS_EPS_S_MUMU",   "ARG_EPS_S_MUMU", 
    "ABS_EPS_S_TAUMU",  "ARG_EPS_S_TAUMU",
    "ABS_EPS_S_ETAU",   "ARG_EPS_S_ETAU",
    "ABS_EPS_S_MUTAU",  "ARG_EPS_S_MUTAU",
    "ABS_EPS_S_TAUTAU", "ARG_EPS_S_TAUTAU",

    "EPS_M_EE",
    "ABS_EPS_M_EMU",    "ARG_EPS_M_EMU",
    "ABS_EPS_M_ETAU",   "ARG_EPS_M_ETAU",
    "EPS_M_MUMU",
    "ABS_EPS_M_MUTAU",  "ARG_EPS_M_MUTAU",
    "EPS_M_TAUTAU",

    "ABS_EPS_D_EE",     "ARG_EPS_D_EE",
    "ABS_EPS_D_MUE",    "ARG_EPS_D_MUE",
    "ABS_EPS_D_TAUE",   "ARG_EPS_D_TAUE",
    "ABS_EPS_D_EMU",    "ARG_EPS_D_EMU",
    "ABS_EPS_D_MUMU",   "ARG_EPS_D_MUMU",
    "ABS_EPS_D_TAUMU",  "ARG_EPS_D_TAUMU",
    "ABS_EPS_D_ETAU",   "ARG_EPS_D_ETAU", 
    "ABS_EPS_D_MUTAU",  "ARG_EPS_D_MUTAU",
    "ABS_EPS_D_TAUTAU", "ARG_EPS_D_TAUTAU" };

/* 1\sigma bounds on NSI parameters, taken from 0907.0097 and 0704.1800               */
/* and converted to 1\sigma level according to \sigma_{68} = \sigma_{90} / \sqrt{2.7} */
double param_bounds[] = 
  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.025, 0.0,   0.016, 0.0,   0.073, 0.0,
    0.015, 0.0,   0.047, 0.0,   0.011, 0.0,
    0.025, 0.0,   0.008, 0.0,   0.079, 0.0,

    2.556,        0.201, 0.0,   1.826, 0.0,
                  0.041,        0.038, 0.0,
                                0.122,     

    0.025, 0.0,   0.015, 0.0,   0.025, 0.0,
    0.016, 0.0,   0.047, 0.0,   0.008, 0.0,
    0.073, 0.0,   0.011, 0.0,   0.079, 0.0 };


/* Internal temporary variables */
static gsl_matrix_complex *U=NULL; /* The vacuum mixing matrix                           */
static gsl_matrix_complex *H=NULL; /* Neutrino Hamiltonian                               */
static gsl_matrix_complex *Q=NULL; /* Eigenvectors of Hamiltonian (= eff. mixing matrix) */
static gsl_vector *lambda=NULL;    /* Eigenvalues of Hamiltonian                         */
static gsl_matrix_complex *S=NULL; /* The neutrino S-matrix                              */

static gsl_matrix_complex *H0_template=NULL;  /* Used in the construction of the vac. Hamiltonian */
static gsl_matrix_complex *S1=NULL, *T0=NULL; /* Temporary matrix storage                         */
static gsl_matrix_complex *Q1=NULL, *Q2=NULL; /* More temporary storage                           */

//extern int density_corr[];


/***************************************************************************
 *            3 x 3   E I G E N S Y S T E M   F U N C T I O N S            *
 ***************************************************************************
 * These functions are taken from physics/0610206.                         *
 ***************************************************************************/

// ----------------------------------------------------------------------------
static void zhetrd3(double complex A[3][3], double complex Q[3][3],
                    double d[3], double e[2])
// ----------------------------------------------------------------------------
// Reduces a hermitian 3x3 matrix to real tridiagonal form by applying
// (unitary) Householder transformations:
//            [ d[0]  e[0]       ]
//    A = Q . [ e[0]  d[1]  e[1] ] . Q^T
//            [       e[1]  d[2] ]
// The function accesses only the diagonal and upper triangular parts of
// A. The access is read-only.
// ---------------------------------------------------------------------------
{
  const int n = 3;
  double complex u[n], q[n];
  double complex omega, f;
  double K, h, g;
  
  // Initialize Q to the identitity matrix
#ifndef EVALS_ONLY
  for (int i=0; i < n; i++)
  {
    Q[i][i] = 1.0;
    for (int j=0; j < i; j++)
      Q[i][j] = Q[j][i] = 0.0;
  }
#endif

  // Bring first row and column to the desired form 
  h = SQR_ABS(A[0][1]) + SQR_ABS(A[0][2]);
  if (creal(A[0][1]) > 0)
    g = -sqrt(h);
  else
    g = sqrt(h);
  e[0] = g;
  f    = g * A[0][1];
  u[1] = conj(A[0][1]) - g;
  u[2] = conj(A[0][2]);
  
  omega = h - f;
  if (creal(omega) > 0.0)
  {
    omega = 0.5 * (1.0 + conj(omega)/omega) / creal(omega);
    K = 0.0;
    for (int i=1; i < n; i++)
    {
      f    = conj(A[1][i]) * u[1] + A[i][2] * u[2];
      q[i] = omega * f;                  // p
      K   += creal(conj(u[i]) * f);      // u* A u
    }
    K *= 0.5 * SQR_ABS(omega);

    for (int i=1; i < n; i++)
      q[i] = q[i] - K * u[i];
    
    d[0] = creal(A[0][0]);
    d[1] = creal(A[1][1]) - 2.0*creal(q[1]*conj(u[1]));
    d[2] = creal(A[2][2]) - 2.0*creal(q[2]*conj(u[2]));
    
    // Store inverse Householder transformation in Q
#ifndef EVALS_ONLY
    for (int j=1; j < n; j++)
    {
      f = omega * conj(u[j]);
      for (int i=1; i < n; i++)
        Q[i][j] = Q[i][j] - f*u[i];
    }
#endif

    // Calculate updated A[1][2] and store it in f
    f = A[1][2] - q[1]*conj(u[2]) - u[1]*conj(q[2]);
  }
  else
  {
    for (int i=0; i < n; i++)
      d[i] = creal(A[i][i]);
    f = A[1][2];
  }

  // Make (23) element real
  e[1] = cabs(f);
#ifndef EVALS_ONLY
  if (e[1] != 0.0)
  {
    f = conj(f) / e[1];
    for (int i=1; i < n; i++)
      Q[i][n-1] = Q[i][n-1] * f;
  }
#endif
}


// ----------------------------------------------------------------------------
static int zheevc3(double complex A[3][3], double w[3])
// ----------------------------------------------------------------------------
// Calculates the eigenvalues of a hermitian 3x3 matrix A using Cardano's
// analytical algorithm.
// Only the diagonal and upper triangular parts of A are accessed. The access
// is read-only.
// ----------------------------------------------------------------------------
// Parameters:
//   A: The hermitian input matrix
//   w: Storage buffer for eigenvalues
// ----------------------------------------------------------------------------
// Return value:
//   0: Success
//  -1: Error
// ----------------------------------------------------------------------------
{
  double m, c1, c0;
  
  // Determine coefficients of characteristic poynomial. We write
  //       | a   d   f  |
  //  A =  | d*  b   e  |
  //       | f*  e*  c  |
  double complex de = A[0][1] * A[1][2];                            // d * e
  double dd = SQR_ABS(A[0][1]);                                  // d * conj(d)
  double ee = SQR_ABS(A[1][2]);                                  // e * conj(e)
  double ff = SQR_ABS(A[0][2]);                                  // f * conj(f)
  m  = creal(A[0][0]) + creal(A[1][1]) + creal(A[2][2]);
  c1 = (creal(A[0][0])*creal(A[1][1])  // a*b + a*c + b*c - d*conj(d) - e*conj(e) - f*conj(f)
          + creal(A[0][0])*creal(A[2][2])
          + creal(A[1][1])*creal(A[2][2]))
          - (dd + ee + ff);
  c0 = creal(A[2][2])*dd + creal(A[0][0])*ee + creal(A[1][1])*ff
            - creal(A[0][0])*creal(A[1][1])*creal(A[2][2])
            - 2.0 * (creal(A[0][2])*creal(de) + cimag(A[0][2])*cimag(de));
                             // c*d*conj(d) + a*e*conj(e) + b*f*conj(f) - a*b*c - 2*Re(conj(f)*d*e)

  double p, sqrt_p, q, c, s, phi;
  p = SQR(m) - 3.0*c1;
  q = m*(p - (3.0/2.0)*c1) - (27.0/2.0)*c0;
  sqrt_p = sqrt(fabs(p));

  phi = 27.0 * ( 0.25*SQR(c1)*(p - c1) + c0*(q + 27.0/4.0*c0));
  phi = (1.0/3.0) * atan2(sqrt(fabs(phi)), q);
  
  c = sqrt_p*cos(phi);
  s = (1.0/M_SQRT3)*sqrt_p*sin(phi);

  w[1]  = (1.0/3.0)*(m - c);
  w[2]  = w[1] + s;
  w[0]  = w[1] + c;
  w[1] -= s;

  return 0;
}


// ----------------------------------------------------------------------------
static int zheevq3(double complex A[3][3], double complex Q[3][3], double w[3])
// ----------------------------------------------------------------------------
// Calculates the eigenvalues and normalized eigenvectors of a hermitian 3x3
// matrix A using the QL algorithm with implicit shifts, preceded by a
// Householder reduction to real tridiagonal form.
// The function accesses only the diagonal and upper triangular parts of A.
// The access is read-only.
// ----------------------------------------------------------------------------
// Parameters:
//   A: The hermitian input matrix
//   Q: Storage buffer for eigenvectors
//   w: Storage buffer for eigenvalues
// ----------------------------------------------------------------------------
// Return value:
//   0: Success
//  -1: Error (no convergence)
// ----------------------------------------------------------------------------
// Dependencies:
//   zhetrd3()
// ----------------------------------------------------------------------------
{
  const int n = 3;
  double e[3];                 // The third element is used only as temporary workspace
  double g, r, p, f, b, s, c;  // Intermediate storage
  double complex t;
  int nIter;
  int m;

  // Transform A to real tridiagonal form by the Householder method
  zhetrd3(A, Q, w, e);
  
  // Calculate eigensystem of the remaining real symmetric tridiagonal matrix
  // with the QL method
  //
  // Loop over all off-diagonal elements
  for (int l=0; l < n-1; l++)
  {
    nIter = 0;
    while (1)
    {
      // Check for convergence and exit iteration loop if off-diagonal
      // element e(l) is zero
      for (m=l; m <= n-2; m++)
      {
        g = fabs(w[m])+fabs(w[m+1]);
        if (fabs(e[m]) + g == g)
          break;
      }
      if (m == l)
        break;
      
      if (nIter++ >= 30)
        return -1;

      // Calculate g = d_m - k
      g = (w[l+1] - w[l]) / (e[l] + e[l]);
      r = sqrt(SQR(g) + 1.0);
      if (g > 0)
        g = w[m] - w[l] + e[l]/(g + r);
      else
        g = w[m] - w[l] + e[l]/(g - r);

      s = c = 1.0;
      p = 0.0;
      for (int i=m-1; i >= l; i--)
      {
        f = s * e[i];
        b = c * e[i];
        if (fabs(f) > fabs(g))
        {
          c      = g / f;
          r      = sqrt(SQR(c) + 1.0);
          e[i+1] = f * r;
          c     *= (s = 1.0/r);
        }
        else
        {
          s      = f / g;
          r      = sqrt(SQR(s) + 1.0);
          e[i+1] = g * r;
          s     *= (c = 1.0/r);
        }
        
        g = w[i+1] - p;
        r = (w[i] - g)*s + 2.0*c*b;
        p = s * r;
        w[i+1] = g + p;
        g = c*r - b;

        // Form eigenvectors
#ifndef EVALS_ONLY
        for (int k=0; k < n; k++)
        {
          t = Q[k][i+1];
          Q[k][i+1] = s*Q[k][i] + c*t;
          Q[k][i]   = c*Q[k][i] - s*t;
        }
#endif 
      }
      w[l] -= p;
      e[l]  = g;
      e[m]  = 0.0;
    }
  }

  return 0;
}


// ----------------------------------------------------------------------------
static int zheevh3(double complex A[3][3], double complex Q[3][3], double w[3])
// ----------------------------------------------------------------------------
// Calculates the eigenvalues and normalized eigenvectors of a hermitian 3x3
// matrix A using Cardano's method for the eigenvalues and an analytical
// method based on vector cross products for the eigenvectors. However,
// if conditions are such that a large error in the results is to be
// expected, the routine falls back to using the slower, but more
// accurate QL algorithm. Only the diagonal and upper triangular parts of A need
// to contain meaningful values. Access to A is read-only.
// ----------------------------------------------------------------------------
// Parameters:
//   A: The hermitian input matrix
//   Q: Storage buffer for eigenvectors
//   w: Storage buffer for eigenvalues
// ----------------------------------------------------------------------------
// Return value:
//   0: Success
//  -1: Error
// ----------------------------------------------------------------------------
// Dependencies:
//   zheevc3(), zhetrd3(), zheevq3()
// ----------------------------------------------------------------------------
// Version history:
//   v1.1: Simplified fallback condition --> speed-up
//   v1.0: First released version
// ----------------------------------------------------------------------------
{
#ifndef EVALS_ONLY
  double norm;          // Squared norm or inverse norm of current eigenvector
//  double n0, n1;        // Norm of first and second columns of A
  double error;         // Estimated maximum roundoff error
  double t, u;          // Intermediate storage
  int j;                // Loop counter
#endif

  // Calculate eigenvalues
  zheevc3(A, w);

#ifndef EVALS_ONLY
//  n0 = SQR(creal(A[0][0])) + SQR_ABS(A[0][1]) + SQR_ABS(A[0][2]);
//  n1 = SQR_ABS(A[0][1]) + SQR(creal(A[1][1])) + SQR_ABS(A[1][2]);
  
  t = fabs(w[0]);
  if ((u=fabs(w[1])) > t)
    t = u;
  if ((u=fabs(w[2])) > t)
    t = u;
  if (t < 1.0)
    u = t;
  else
    u = SQR(t);
  error = 256.0 * DBL_EPSILON * SQR(u);
//  error = 256.0 * DBL_EPSILON * (n0 + u) * (n1 + u);

  Q[0][1] = A[0][1]*A[1][2] - A[0][2]*creal(A[1][1]);
  Q[1][1] = A[0][2]*conj(A[0][1]) - A[1][2]*creal(A[0][0]);
  Q[2][1] = SQR_ABS(A[0][1]);

  // Calculate first eigenvector by the formula
  //   v[0] = conj( (A - w[0]).e1 x (A - w[0]).e2 )
  Q[0][0] = Q[0][1] + A[0][2]*w[0];
  Q[1][0] = Q[1][1] + A[1][2]*w[0];
  Q[2][0] = (creal(A[0][0]) - w[0]) * (creal(A[1][1]) - w[0]) - Q[2][1];
  norm    = SQR_ABS(Q[0][0]) + SQR_ABS(Q[1][0]) + SQR(creal(Q[2][0]));

  // If vectors are nearly linearly dependent, or if there might have
  // been large cancellations in the calculation of A(I,I) - W(1), fall
  // back to QL algorithm
  // Note that this simultaneously ensures that multiple eigenvalues do
  // not cause problems: If W(1) = W(2), then A - W(1) * I has rank 1,
  // i.e. all columns of A - W(1) * I are linearly dependent.
  if (norm <= error)
    return zheevq3(A, Q, w);
  else                      // This is the standard branch
  {
    norm = sqrt(1.0 / norm);
    for (j=0; j < 3; j++)
      Q[j][0] = Q[j][0] * norm;
  }
  
  // Calculate second eigenvector by the formula
  //   v[1] = conj( (A - w[1]).e1 x (A - w[1]).e2 )
  Q[0][1]  = Q[0][1] + A[0][2]*w[1];
  Q[1][1]  = Q[1][1] + A[1][2]*w[1];
  Q[2][1]  = (creal(A[0][0]) - w[1]) * (creal(A[1][1]) - w[1]) - creal(Q[2][1]);
  norm     = SQR_ABS(Q[0][1]) + SQR_ABS(Q[1][1]) + SQR(creal(Q[2][1]));
  if (norm <= error)
    return zheevq3(A, Q, w);
  else
  {
    norm = sqrt(1.0 / norm);
    for (j=0; j < 3; j++)
      Q[j][1] = Q[j][1] * norm;
  }
  
  // Calculate third eigenvector according to
  //   v[2] = conj(v[0] x v[1])
  Q[0][2] = conj(Q[1][0]*Q[2][1] - Q[2][0]*Q[1][1]);
  Q[1][2] = conj(Q[2][0]*Q[0][1] - Q[0][0]*Q[2][1]);
  Q[2][2] = conj(Q[0][0]*Q[1][1] - Q[1][0]*Q[0][1]);
#endif

  return 0;
}



/***************************************************************************
 *                  I N T E R N A L   F U N C T I O N S                    *
 ***************************************************************************/

/***************************************************************************
 * Function nsi_init_probability_engine                                    *
 ***************************************************************************
 * Allocates internal data structures for the probability engine.          *
 ***************************************************************************/
int nsi_init_probability_engine()
{
  nsi_free_probability_engine();
  
  U = gsl_matrix_complex_calloc(GLB_NU_FLAVOURS, GLB_NU_FLAVOURS);
  H = gsl_matrix_complex_calloc(GLB_NU_FLAVOURS, GLB_NU_FLAVOURS);
  Q = gsl_matrix_complex_calloc(GLB_NU_FLAVOURS, GLB_NU_FLAVOURS);
  lambda = gsl_vector_alloc(GLB_NU_FLAVOURS);
  S = gsl_matrix_complex_calloc(GLB_NU_FLAVOURS, GLB_NU_FLAVOURS);
    
  H0_template = gsl_matrix_complex_calloc(GLB_NU_FLAVOURS, GLB_NU_FLAVOURS);
  S1 = gsl_matrix_complex_calloc(GLB_NU_FLAVOURS, GLB_NU_FLAVOURS);
  T0 = gsl_matrix_complex_calloc(GLB_NU_FLAVOURS, GLB_NU_FLAVOURS);
  Q1 = gsl_matrix_complex_calloc(GLB_NU_FLAVOURS, GLB_NU_FLAVOURS);
  Q2 = gsl_matrix_complex_calloc(GLB_NU_FLAVOURS, GLB_NU_FLAVOURS);

  return 0;
}


/***************************************************************************
 * Function nsi_free_probability_engine                                    *
 ***************************************************************************
 * Destroys internal data structures of the probability engine.            *
 ***************************************************************************/
int nsi_free_probability_engine()
{
  if (Q2!=NULL)     { gsl_matrix_complex_free(Q2);  Q2 = NULL; }
  if (Q1!=NULL)     { gsl_matrix_complex_free(Q1);  Q1 = NULL; }
  if (T0!=NULL)     { gsl_matrix_complex_free(T0);  T0 = NULL; }
  if (S1!=NULL)     { gsl_matrix_complex_free(S1);  S1 = NULL; }
  if (H0_template!=NULL) { gsl_matrix_complex_free(H0_template);  H0_template = NULL; }
  
  if (S!=NULL)      { gsl_matrix_complex_free(S);   S = NULL; }
  if (lambda!=NULL) { gsl_vector_free(lambda);      lambda = NULL; }
  if (Q!=NULL)      { gsl_matrix_complex_free(Q);   Q = NULL; }
  if (H!=NULL)      { gsl_matrix_complex_free(H);   H = NULL; }
  if (U!=NULL)      { gsl_matrix_complex_free(U);   U = NULL; }

  return 0;
}


/***************************************************************************
 * Function nsi_set_oscillation_parameters                                 *
 ***************************************************************************
 * Sets the fundamental oscillation parameters and precomputes the mixing  *
 * matrix and part of the Hamiltonian.                                     *
 ***************************************************************************/
int nsi_set_oscillation_parameters(glb_params p, void *user_data)
{
  double complex (*_U)[GLB_NU_FLAVOURS]
    = (double complex (*)[GLB_NU_FLAVOURS]) gsl_matrix_complex_ptr(U, 0, 0);
  int i, j, k;

  /* Implement correlations between density parameters. This requires that under
   * all circumstances the scaling of the matter density is performed _after_
   * calling set_oscillation_parameters! At present, this works only with
   * the hybrid minimizer (GLB_MIN_POWELL)! */
//  for (j=0; j < glb_num_of_exps; j++)
//    if (density_corr[j] != j)
//      glbSetDensityParams(p, glbGetDensityParams(p, density_corr[j]), j);

  /* Copy standard oscillation parameters */
  th12  = glbGetOscParams(p, GLB_THETA_12);
  th13  = glbGetOscParams(p, GLB_THETA_13);
  th23  = glbGetOscParams(p, GLB_THETA_23);
  delta = glbGetOscParams(p, GLB_DELTA_CP);
  mq[0] = glbGetOscParams(p, GLB_DM_31);
  mq[1] = glbGetOscParams(p, GLB_DM_31) + glbGetOscParams(p, GLB_DM_21);
  mq[2] = glbGetOscParams(p, GLB_DM_31) + glbGetOscParams(p, GLB_DM_31);

  /* Copy non-standard parameters */
  k = ABS_EPS_S_EE;
  for (i=0; i < GLB_NU_FLAVOURS; i++)
  {
    for (j=0; j < GLB_NU_FLAVOURS; j++)
    {
      epsilon_s_plus_1[i][j] = glbGetOscParams(p,k) * cexp(I*glbGetOscParams(p,k+1));
      k += 2;
    }
    epsilon_s_plus_1[i][i] += 1.0;
  }

  k = EPS_M_EE;
  for (i=0; i < GLB_NU_FLAVOURS; i++)
  {
    epsilon_m[i][i] = glbGetOscParams(p,k);
    k++;
    for (j=i+1; j < GLB_NU_FLAVOURS; j++)
    {
      epsilon_m[i][j] = glbGetOscParams(p,k) * cexp(I*glbGetOscParams(p,k+1));
      epsilon_m[j][i] = conj(epsilon_m[i][j]);
      k += 2;
    }
  }
  
  k = ABS_EPS_D_EE;
  for (i=0; i < GLB_NU_FLAVOURS; i++)
  {
    for (j=0; j < GLB_NU_FLAVOURS; j++)
    {
      epsilon_d_plus_1[i][j] = glbGetOscParams(p,k) * cexp(I*glbGetOscParams(p,k+1));
      k += 2;
    }
    epsilon_d_plus_1[i][i] += 1.0;
  }

  /* Compute vacuum mixing matrix */
  _U[0][0] = cos(th12)*cos(th13);
  _U[0][1] = sin(th12)*cos(th13);
  _U[0][2] = sin(th13) * cexp(-I * delta);

  _U[1][0] = -sin(th12)*cos(th23) - cos(th12)*sin(th23)*sin(th13) * cexp(I*delta);
  _U[1][1] =  cos(th12)*cos(th23) - sin(th12)*sin(th23)*sin(th13) * cexp(I*delta);
  _U[1][2] =  sin(th23)*cos(th13);

  _U[2][0] =  sin(th12)*sin(th23) - cos(th12)*cos(th23)*sin(th13) * cexp(I*delta);
  _U[2][1] = -cos(th12)*sin(th23) - sin(th12)*cos(th23)*sin(th13) * cexp(I*delta);
  _U[2][2] =  cos(th23)*cos(th13);
 
  /* Calculate energy independent matrix H0 * E */
  gsl_matrix_complex_set_zero(H0_template);
  gsl_matrix_complex_set_zero(H);
  for (i=0; i < GLB_NU_FLAVOURS; i++)
    gsl_matrix_complex_set(H0_template, i, i, gsl_complex_rect(0.5*mq[i], 0.0));

  gsl_matrix_complex *T = gsl_matrix_complex_alloc(GLB_NU_FLAVOURS, GLB_NU_FLAVOURS);
  gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, GSL_COMPLEX_ONE, H0_template, U, /* T=H0.U^\dagger */
                 GSL_COMPLEX_ZERO, T);
  gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, U, T,             /* H0=U.T */
                 GSL_COMPLEX_ZERO, H0_template);
  gsl_matrix_complex_free(T);
  return 0;
}


/***************************************************************************
 * Function nsi_get_oscillation_parameters                                 *
 ***************************************************************************
 * Returns the current set of oscillation parameters.                      *
 ***************************************************************************/
int nsi_get_oscillation_parameters(glb_params p, void *user_data)
{
  int i, j, k;
  glbDefineParams(p, th12, th13, th23, delta, mq[1] - mq[0], mq[2] - mq[0]);
  
  /* Copy non-standard parameters */
  k = ABS_EPS_S_EE;
  for (i=0; i < GLB_NU_FLAVOURS; i++)
    for (j=0; j < GLB_NU_FLAVOURS; j++)
    {
      if (i == j)
      {
        glbSetOscParams(p, cabs(epsilon_s_plus_1[i][j] - 1.0), k);
        glbSetOscParams(p, carg(epsilon_s_plus_1[i][j] - 1.0), k+1);
      }
      else
      {
        glbSetOscParams(p, cabs(epsilon_s_plus_1[i][j]), k);
        glbSetOscParams(p, carg(epsilon_s_plus_1[i][j]), k+1);
      }
      k += 2;
    }

  k = EPS_M_EE;
  for (i=0; i < GLB_NU_FLAVOURS; i++)
  {
    glbSetOscParams(p, epsilon_m[i][i], k);
    k++;
    for (j=i+1; j < GLB_NU_FLAVOURS; j++)
    {
      glbSetOscParams(p, cabs(epsilon_m[i][j]), k);
      glbSetOscParams(p, carg(epsilon_m[i][j]), k+1);
      k += 2;
    }
  }
  
  k = ABS_EPS_D_EE;
  for (i=0; i < GLB_NU_FLAVOURS; i++)
    for (j=0; j < GLB_NU_FLAVOURS; j++)
    {
      if (i == j)
      {
        glbSetOscParams(p, cabs(epsilon_d_plus_1[i][j] - 1.0), k);
        glbSetOscParams(p, carg(epsilon_d_plus_1[i][j] - 1.0), k+1);
      }
      else
      {
        glbSetOscParams(p, cabs(epsilon_d_plus_1[i][j]), k);
        glbSetOscParams(p, carg(epsilon_d_plus_1[i][j]), k+1);
      }
      k += 2;
    }

  return 0;
}


/***************************************************************************
 * Function nsi_hamiltonian_cd                                             *
 ***************************************************************************
 * Calculates the Hamiltonian for neutrinos (cp_sign=1) or antineutrinos   *
 * (cp_sign=-1) with energy E, propagating in matter of density V          *
 * (> 0 even for antineutrinos) and stores the result in H.                *
 ***************************************************************************/
int nsi_hamiltonian_cd(double E, double V, int cp_sign)
{
  double inv_E = 1.0 / E;
  double complex (*_H)[GLB_NU_FLAVOURS]
    = (double complex (*)[GLB_NU_FLAVOURS]) gsl_matrix_complex_ptr(H, 0, 0);
  double complex (*_H0_template)[GLB_NU_FLAVOURS]
    = (double complex (*)[GLB_NU_FLAVOURS]) gsl_matrix_complex_ptr(H0_template, 0, 0);
  int i, j;

  if (cp_sign > 0)
  {
    for (i=0; i < GLB_NU_FLAVOURS; i++)
      for (j=0; j < GLB_NU_FLAVOURS; j++)
        _H[i][j] = _H0_template[i][j] * inv_E  +  V*epsilon_m[i][j];
  }
  else
  {
    for (i=0; i < GLB_NU_FLAVOURS; i++)
      for (j=0; j < GLB_NU_FLAVOURS; j++)
        _H[i][j] = conj(_H0_template[i][j] * inv_E  -  V*epsilon_m[i][j]);
                                                /* delta_CP -> -delta_CP */
  }
 
  _H[0][0] = _H[0][0] + cp_sign*V;

// _H[i][j] = _H[i][j] + epsilon_m[i][j]
// _H[j][i] = _H[j][i] + conj(epsilon_m[i][j]);
//                                   COMPLEX!
//  for anti-neutrinos:
// _H[i][j] = _H[i][j] - conj(epsilon_m[i][j]);
// _H[j][i] = _H[j][i] - epsilon_m[i][j];
  return 0;
}


/***************************************************************************
 * Function nsi_S_matrix_cd                                                *
 ***************************************************************************
 * Calculates the S matrix for neutrino oscillations in matter of constant *
 * density using a fast eigenvalue solver optimized to 3x3 matrices.       *
 ***************************************************************************
 * Parameters:                                                             *
 *   E: Neutrino energy                                                    *
 *   L: Baseline                                                           *
 *   V: Matter potential (must be > 0 even for antineutrinos)              *
 *   cp_sign: +1 for neutrinos, -1 for antineutrinos                       *
 ***************************************************************************/
int nsi_S_matrix_cd(double E, double L, double V, int cp_sign)
{
  /* Introduce some abbreviations */
  double complex (*_S)[3]  = (double complex (*)[3]) gsl_matrix_complex_ptr(S,0,0);
  double complex (*_Q)[3]  = (double complex (*)[3]) gsl_matrix_complex_ptr(Q,0,0);
  double complex (*_T0)[3] = (double complex (*)[3]) gsl_matrix_complex_ptr(T0,0,0);
  double *_lambda = gsl_vector_ptr(lambda,0);
  int status;
  int i, j, k;
  
  if (fabs(V) < V_THRESHOLD)                       /* Vacuum */
  {
    /* Use vacuum mixing angles and masses */
    double inv_E = 0.5/E;
    for (i=0; i < GLB_NU_FLAVOURS; i++)
      _lambda[i] = mq[i] * inv_E;

    if (cp_sign > 0)
      gsl_matrix_complex_memcpy(Q, U);
    else
    {
      double complex (*_U)[3]  = (double complex (*)[3]) gsl_matrix_complex_ptr(U,0,0);
      for (i=0; i < GLB_NU_FLAVOURS; i++)
        for (j=0; j < GLB_NU_FLAVOURS; j++)
          _Q[i][j] = conj(_U[i][j]);
    }
  }
  else                                             /* Matter */
  {
    double complex (*_H)[3] = (double complex (*)[3]) gsl_matrix_complex_ptr(H,0,0);
    
    /* Calculate neutrino Hamiltonian */
    if ((status=nsi_hamiltonian_cd(E, V, cp_sign)) != 0)
      return status;
    
    /* Calculate eigenvalues of Hamiltonian */
    if ((status=zheevh3(_H, _Q, _lambda)) != 0)
      return status;
  }

  /* Calculate S-Matrix in mass basis in matter ... */
  double phase;
  gsl_matrix_complex_set_zero(S);
  for (i=0; i < GLB_NU_FLAVOURS; i++)
  {
    phase    = -L * _lambda[i];
    _S[i][i] = cos(phase) + I*sin(phase);
  } 
  
  /* ... and transform it to the flavour basis */
  gsl_matrix_complex_set_zero(T0);
  double complex *p = &_T0[0][0];
  for (i=0; i < GLB_NU_FLAVOURS; i++)              /* T0 = S.Q^\dagger */
    for (j=0; j < GLB_NU_FLAVOURS; j++)
    {
      for (int k=0; k < GLB_NU_FLAVOURS; k++)
      {
        *p += ( creal(_S[i][k])*creal(_Q[j][k])+cimag(_S[i][k])*cimag(_Q[j][k]) )
                + I * ( cimag(_S[i][k])*creal(_Q[j][k])-creal(_S[i][k])*cimag(_Q[j][k]) );
      }
      p++;
    }
  gsl_matrix_complex_set_zero(S);
  p = &_S[0][0];
  for (i=0; i < GLB_NU_FLAVOURS; i++)              /* S = Q.T0 */
    for (j=0; j < GLB_NU_FLAVOURS; j++)
    {
      for (k=0; k < GLB_NU_FLAVOURS; k++)
      {
        *p += ( creal(_Q[i][k])*creal(_T0[k][j])-cimag(_Q[i][k])*cimag(_T0[k][j]) )
                + I * ( cimag(_Q[i][k])*creal(_T0[k][j])+creal(_Q[i][k])*cimag(_T0[k][j]) );
      }
      p++;
    }

  /* Incorporate non-standard interactions in the source and in the detector */
  if (cp_sign > 0)
  {
    gsl_matrix_complex_set_zero(T0);
    for (i=0; i < GLB_NU_FLAVOURS; i++)            /* T0 = S.(1+epsilon_s) */
      for (j=0; j < GLB_NU_FLAVOURS; j++)
        for (k=0; k < GLB_NU_FLAVOURS; k++)
          _T0[i][j] += _S[i][k] * epsilon_s_plus_1[k][j];
    gsl_matrix_complex_set_zero(S);
    for (i=0; i < GLB_NU_FLAVOURS; i++)            /* S = (1+epsilon_d).T0 */
      for (j=0; j < GLB_NU_FLAVOURS; j++)
        for (k=0; k < GLB_NU_FLAVOURS; k++)
          _S[i][j] += epsilon_d_plus_1[i][k] * _T0[k][j];
  }
  else
  {
    gsl_matrix_complex_set_zero(T0);
    for (i=0; i < GLB_NU_FLAVOURS; i++)            /* T0 = S.conj(1+epsilon_s) */
      for (j=0; j < GLB_NU_FLAVOURS; j++)
        for (k=0; k < GLB_NU_FLAVOURS; k++)
          _T0[i][j] += _S[i][k] * conj(epsilon_s_plus_1[k][j]);
    gsl_matrix_complex_set_zero(S);
    for (i=0; i < GLB_NU_FLAVOURS; i++)            /* S = conj(1+epsilon_d).T0 */
      for (j=0; j < GLB_NU_FLAVOURS; j++)
        for (k=0; k < GLB_NU_FLAVOURS; k++)
          _S[i][j] += conj(epsilon_d_plus_1[i][k]) * _T0[k][j];
  }

// S --> epsilon_d_plus_1 . S . epsilon_s_plus_1
// for anti-nu: S --> epsilon_d_plus_1^* . S . epsilon_s_plus_1^*

  return 0;
}


/***************************************************************************
 * Function nsi_filtered_probability_matrix_cd                             *
 ***************************************************************************
 * Calculates the probability matrix for neutrino oscillations in matter   *
 * of constant density, including a low pass filter to suppress aliasing   *
 * due to very fast oscillations.                                          *
 ***************************************************************************
 * Parameters:                                                             *
 *   P: Storage buffer for the probability matrix                          *
 *   E: Neutrino energy                                                    *
 *   L: Baseline                                                           *
 *   V: Matter potential (must be > 0 even for antineutrinos)              *
 *   sigma: Width of Gaussian filter                                       *
 *   cp_sign: +1 for neutrinos, -1 for antineutrinos                       *
 ***************************************************************************/
int nsi_filtered_probability_matrix_cd(double P[3][3], double E, double L, double V,
                                       double sigma, int cp_sign)
{
  /* Introduce some abbreviations */
  double complex (*_Q)[3]  = (double complex (*)[3]) gsl_matrix_complex_ptr(Q,0,0);
  double complex (*_T0)[3] = (double complex (*)[3]) gsl_matrix_complex_ptr(T0,0,0);
  double complex (*_Q1)[3] = (double complex (*)[3]) gsl_matrix_complex_ptr(Q1,0,0);
  double complex (*_Q2)[3] = (double complex (*)[3]) gsl_matrix_complex_ptr(Q2,0,0);
  double *_lambda = gsl_vector_ptr(lambda,0);
  int status;
  int i, j, k, l;
 
  /* Vacuum: Use vacuum mixing angles and masses */
  if (fabs(V) < V_THRESHOLD)
  {
    double inv_E = 0.5/E;
    for (i=0; i < GLB_NU_FLAVOURS; i++)
      _lambda[i] = mq[i] * inv_E;

    if (cp_sign > 0)
      gsl_matrix_complex_memcpy(Q, U);
    else
    {
      double complex (*_U)[3]  = (double complex (*)[3]) gsl_matrix_complex_ptr(U,0,0);
      for (i=0; i < GLB_NU_FLAVOURS; i++)
        for (j=0; j < GLB_NU_FLAVOURS; j++)
          _Q[i][j] = conj(_U[i][j]);
    }
  }

  /* Matter: Rediagonalize Hamiltonian */
  else
  {
    double complex (*_H)[3] = (double complex (*)[3]) gsl_matrix_complex_ptr(H,0,0);
    
    /* Calculate neutrino Hamiltonian */
    if ((status=nsi_hamiltonian_cd(E, V, cp_sign)) != 0)
      return status;
    
    /* Calculate eigenvalues and eigenvectors of Hamiltonian */
    if ((status=zheevh3(_H, _Q, _lambda)) != 0)
      return status;
  }

  /* Define Q_1^\dag = Q^\dag . (1 + \eps^s) and Q_2 = (1 + \eps^d) . Q */
  /* (for anti-neutrinos: \eps^{s,d} -> (\eps^{s,d})^*             */
  gsl_matrix_complex_set_zero(Q1);
  gsl_matrix_complex_set_zero(Q2);
  if (cp_sign > 0)
  {
    for (i=0; i < GLB_NU_FLAVOURS; i++)
      for (j=0; j < GLB_NU_FLAVOURS; j++)
        for (k=0; k < GLB_NU_FLAVOURS; k++)
          _Q1[i][j] += conj(epsilon_s_plus_1[k][i]) * _Q[k][j];
    for (i=0; i < GLB_NU_FLAVOURS; i++)
      for (j=0; j < GLB_NU_FLAVOURS; j++)
        for (k=0; k < GLB_NU_FLAVOURS; k++)
          _Q2[i][j] += epsilon_d_plus_1[i][k] * _Q[k][j];
  }
  else
  {
    for (i=0; i < GLB_NU_FLAVOURS; i++)
      for (j=0; j < GLB_NU_FLAVOURS; j++)
        for (k=0; k < GLB_NU_FLAVOURS; k++)
          _Q1[i][j] += epsilon_s_plus_1[k][i] * _Q[k][j];
    for (i=0; i < GLB_NU_FLAVOURS; i++)
      for (j=0; j < GLB_NU_FLAVOURS; j++)
        for (k=0; k < GLB_NU_FLAVOURS; k++)
          _Q2[i][j] += conj(epsilon_d_plus_1[i][k]) * _Q[k][j];
  }
        

  /* Calculate probability matrix (see GLoBES manual for a discussion of the algorithm) */
  double phase, filter_factor;
  double t = -0.5/1.0e-18 * SQR(sigma) / SQR(E);
  gsl_matrix_complex_set_zero(T0);
  for (i=0; i < GLB_NU_FLAVOURS; i++)
    for (j=i+1; j < GLB_NU_FLAVOURS; j++)
    {
      phase         = -L * (_lambda[i] - _lambda[j]);
      filter_factor = exp(t * SQR(phase));
      _T0[i][j]     = filter_factor * (cos(phase) + I*sin(phase));
    }

  for (k=0; k < GLB_NU_FLAVOURS; k++)
    for (l=0; l < GLB_NU_FLAVOURS; l++)
    {
      P[k][l] = 0.0;
      for (i=0; i < GLB_NU_FLAVOURS; i++)
      {
        complex t = conj(_Q1[k][i]) * _Q2[l][i];
        for (j=i+1; j < GLB_NU_FLAVOURS; j++)
          P[k][l] += 2.0 * creal(_Q1[k][j] * conj(_Q2[l][j]) * t * _T0[i][j]);
        P[k][l] += SQR_ABS(_Q1[k][i]) * SQR_ABS(_Q2[l][i]);
      }
    }
    
  return 0;
}


/***************************************************************************
 * Function nsi_probability_matrix                                         *
 ***************************************************************************
 * Calculates the neutrino oscillation probability matrix.                 *
 ***************************************************************************
 * Parameters:                                                             *
 *   P:       Buffer for the storage of the matrix                         *
 *   cp_sign: +1 for neutrinos, -1 for antineutrinos                       *
 *   E:       Neutrino energy (in GeV)                                     *
 *   psteps:  Number of layers in the matter density profile               *
 *   length:  Lengths of the layers in the matter density profile in km    *
 *   density: The matter densities in g/cm^3                               *
 *   filter_sigma: Width of low-pass filter or <0 for no filter            *
 *   user_data: Unused here, should be NULL                                *
 ***************************************************************************/
int nsi_probability_matrix(double P[3][3], int cp_sign, double E,
    int psteps, const double *length, const double *density,
    double filter_sigma, void *user_data)
{
  int status;
  int i, j;

  /* Convert energy to eV */
  E *= 1.0e9;
  
  if (filter_sigma > 0.0)                     /* With low-pass filter */
  {
    if (psteps == 1)
      nsi_filtered_probability_matrix_cd(P, E, GLB_KM_TO_EV(length[0]),
                density[0]*GLB_V_FACTOR*GLB_Ne_MANTLE, filter_sigma, cp_sign);
    else
      return -1;
  }
  else                                        /* Without low-pass filter */
  {
    if (psteps > 1)
    {
      gsl_matrix_complex_set_identity(S1);                                 /* S1 = 1 */
      for (i=0; i < psteps; i++)
      {
        status = nsi_S_matrix_cd(E, GLB_KM_TO_EV(length[i]), density[i]*GLB_V_FACTOR*GLB_Ne_MANTLE, cp_sign);
        if (status != 0)
          return status;
        gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, S, S1, /* T0 = S.S1 */
                       GSL_COMPLEX_ZERO, T0);
        gsl_matrix_complex_memcpy(S1, T0);                                 /* S1 = T0 */
      } 
      gsl_matrix_complex_memcpy(S, S1);                                    /* S = S1 */
    }
    else
    {
      status = nsi_S_matrix_cd(E, GLB_KM_TO_EV(length[0]), density[0]*GLB_V_FACTOR*GLB_Ne_MANTLE, cp_sign);
      if (status != 0)
        return status;
    }

    double complex (*_S)[3] = (double complex (*)[3]) gsl_matrix_complex_ptr(S,0,0);
    for (i=0; i < GLB_NU_FLAVOURS; i++)
      for (j=0; j < GLB_NU_FLAVOURS; j++)
        P[j][i] = SQR_ABS(_S[i][j]);
  }

  return 0;
}


