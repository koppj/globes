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



/* These routines compute the smear matrix in the same way as CalcSmearA
 * in GlobesMESbeta.m. CalcSmearB is also implemented, but it seems
 * that the algorithm can be improved.
 *
 * Patrick Huber (c) 2004
 */


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* That's new -- using the GNU Scientific Library */
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_min.h>


#include "glb_error.h"
#include "glb_types.h"
#include "glb_smear.h"


#define TOL 1.0E-7

/* An example for a pretty general sigma(E,params)
   sigma(E)/E = in[0] + in[1]/sqrt(E) + in[2]/E */

static double my_sigma(double x, double *in)
{
  return in[0]*x + in[1]*sqrt(x) + in[2];
}


/* Prototypes for private things */
void glb_init_bin_data(glb_smear *in);



/* First some functions for checking user input.
 *
 *
 * The following call sequence is our safety net:
 *   glb_smear *s;  
 *   s=glb_smear_alloc();
 *    user input (via parser)
 *   glb_default_smear(s);
 *   SetupSmearMatixX(s);
 *    now it is safe to use function like BinBoundary ...
 *  SmearMatrixX(s,...);
 * Once everthing is finished
 *  glb_smear_free(s);
 * Other call sequences may lead to undefined behaviour!
 */


/* These are called by init_glb_smear and glb_default_smear */
glb_option_type *glb_option_type_alloc()
{
  glb_option_type *in;
  in=(glb_option_type *) glb_malloc(sizeof(glb_option_type)); 
  in->corr_fac=-1;
  in->confidence_level=-1;
  in->offset=-1;
  in->low_bound=-1;
  in->up_bound=-1;
  return in;
}

glb_option_type *glb_option_type_reset(glb_option_type *in)
{ 
  in->corr_fac=-1;
  in->confidence_level=-1;
  in->offset=-1;
  in->low_bound=-1;
  in->up_bound=-1;
  return in;
}

void glb_option_type_free(glb_option_type *in)
{    
    glb_free(in);
}

int glb_default_option_type(glb_option_type *in, const glb_smear *head)
{
  int s,d;
  s=0;
  d=0;
  if(in->corr_fac==-1) {in->corr_fac=1;d=2;}
  if(in->corr_fac<=0) {glb_error("Correction factor has to be larger"
				 " than 0!");s=1;}
  if(in->confidence_level==-1) {in->confidence_level=0.999;d=2;}
  if(in->confidence_level >= 1 || in->confidence_level <= 0)
    {glb_error("The confidence level has to be inbetween 0 and 1!");s=1;}
  if(in->offset==-1)  {in->offset=0;d=2;}
  if(in->offset<0) {glb_error("Offset must not be smaller than 0!");s=1;}
  
  if(head==NULL)  
    {
      if(in->low_bound==-1) {glb_error("No low_bound given!");s=1;}
      if(in->up_bound==-1) {glb_error("No up_bound given!");s=1;}
    }
  else
    {
      if(in->low_bound==-1) {in->low_bound=head->e_min;;d=2;}
      if(in->up_bound==-1) {in->up_bound=head->e_max;d=2;}
    }
  if(in->low_bound<0) {glb_error("low_bound must be positive!");s=1;}
  if(in->up_bound<0) {glb_error("up_bound must be positive!");s=1;}
  if(in->low_bound >= in->up_bound)  
    {glb_error("low_bound must be smaller than up_bound!");s=1;}
  if(s!=0) glb_fatal("glb_option_type not properly set up!");
  return d;
}


glb_smear  *glb_smear_alloc()
{
  glb_smear *in;
  in=(glb_smear *) glb_malloc(sizeof(glb_smear));
  in->type=-1;
  in->numofbins=-1;
  in->simbins=-1;
  in->e_min=-1;
  in->e_max=-1;
  in->e_sim_min=-1;
  in->e_sim_max=-1;
  in->sigma=NULL;
  in->num_of_params=-1;
  in->sig_f=NULL;  
  in->binsize=NULL;
  in->simbinsize=NULL;
  in->bincenter=NULL;
  in->simbincenter=NULL;
  in->options=NULL;
  return in;
}

glb_smear  *glb_smear_reset(glb_smear *in)
{
  in->type=-1;
  in->numofbins=-1;
  in->simbins=-1;
  in->e_min=-1;
  in->e_max=-1;
  in->e_sim_min=-1;
  in->e_sim_max=-1;
   /* FIXME a possible memory leak -> dereferencing memory to NULL 
    * Problem: at start up sigma is not NULL, but undefined cannot be free'd
    * here !
    */
  in->sigma=NULL;
  in->num_of_params=-1;
  in->sig_f=NULL;
  
  /* FIXME a possible memory leak -> dereferencing memory to NULL */ 
  in->binsize=NULL;
  in->simbinsize=NULL;  
  in->bincenter=NULL;
  in->simbincenter=NULL;
  in->options=NULL;
  return in;
}



void glb_smear_free(glb_smear *in)
{
  if(in==NULL) return;
  if(in->binsize!=NULL) glb_free(in->binsize);
  if(in->bincenter!=NULL) glb_free(in->bincenter);

  if(in->simbinsize!=NULL) glb_free(in->simbinsize);
  if(in->simbincenter!=NULL) glb_free(in->simbincenter);
  if(in->sigma!=NULL) glb_free(in->sigma);

  if(in->options!=NULL) glb_option_type_free(in->options);
  glb_free(in);
}

int glb_default_smear(glb_smear *in)
{
  double med;
  int s,d,i;
  s=0;
  d=0;
  if(in->type==-1) {glb_error("Smearing type has not been defined!");s=1;}
  if(in->numofbins==-1) {glb_error("Number of bins has not been defined!");
  s=1;}
  if(in->numofbins<=0) {glb_error("Number of bins must be positive!");
  s=1;}

  /* simbins, simtresh and simbeam are no longer defaultized here
   * but by the SmearMatrix computation functions
   *
   *  if(in->simbins==-1) {in->simbins=in->numofbins;d=2;}
   *  if(in->simbins <=0) 
   * {glb_error("Number of simbins must be positive!");s=1;}
   */
  if(in->e_min==-1) {glb_error("No minimal energy defined!");s=1;}
  if(in->e_min<=0) {glb_error("Minimal energy must be positive!");s=1;}
  if(in->e_max==-1) {glb_error("No maximal energy defined!");s=1;}
  if(in->e_max<=in->e_min) {glb_error("Maximal energy must be larger than"
				      " the minimal energy!");s=1;}
  if(in->e_sim_min!=-1||in->e_sim_max!=-1) 
    {
      if(in->e_sim_min<=0) 
	{glb_error("Minimal sim energy must be positive!");s=1;}
      if(in->e_sim_max<=in->e_sim_min) 
	{glb_error("Maximal sim energy must belarger than the minimal"
		   " sim energy!");s=1;}
    }
  /*
   * if(in->e_sim_min==-1) {in->e_sim_min=in->e_min;d=2;}
   * if(in->e_sim_max==-1) {in->e_sim_max=in->e_max;d=2;}
   */

  if(in->num_of_params==-1) 
    {
      /*    glb_error("Number of parameters for the energy"
	      " resolution function unspecified!");
	      s=1;*/
      in->num_of_params=3;
    }

  if(in->num_of_params<0) {
    glb_error("Number of parameters for the energy"
	      " resolution function must larger"
	      " than zero!");
    s=1;
  }
  
  if(in->num_of_params==0) {in->sigma=NULL;d=2;}
  
  if(in->num_of_params>0&&in->sigma==NULL)
    {
       glb_error("No parameters for the energy"
	      " resolution function unspecified!");
       s=1;
    }
  /*if(in->sig_f==NULL) {glb_error("No energy resolution"
    " function defined!");s=1;}*/
 
 
   if(in->sig_f==NULL) {in->sig_f=&my_sigma;d=2;}
   if(in->binsize==NULL&&s==0)
    {
      /* In this case we assume that the bins are equidistant */
      in->binsize=(double *) glb_malloc(sizeof(double)*in->numofbins);
      med=(in->e_max-in->e_min)/in->numofbins;
      for(i=0;i<in->numofbins;i++) in->binsize[i]=med;
    }
  if(s==0) glb_init_bin_data(in);
  /* simbinsize will be intialized by InitSmearMatrixB */
  
  d=d+glb_default_option_type(in->options,in);
  
  if(s!=0) glb_fatal("glb_smear not properly set up!");
  return d;
}


glb_smear  *glb_copy_smear(glb_smear *dest, const glb_smear *src)
{
  int i;  
  glb_option_type *temp;
  if(dest->options==NULL) temp=glb_option_type_alloc();
  else temp=dest->options;
  dest=(glb_smear *) memmove(dest,src,sizeof(glb_smear)); 
  dest->options=temp;
  dest->options=(glb_option_type *) memmove(dest->options,src->options,
				    sizeof(glb_option_type)); 
  
 
  if(src->num_of_params <= 0) return dest;
  dest->sigma=(double *) glb_malloc(src->num_of_params*sizeof(double));
  for(i=0;i<src->num_of_params;i++) dest->sigma[i]=src->sigma[i];
  /* What about binsize etc ? */


  return dest;	  
}


/* end */

/* These functions provide data like bin boundaries and bin centers
 * in lookup tables in order to improve performance. Furthemore they
 * also give access to those data. They all are able to deal with user
 * defined bins, ie. variable width.
 */



void glb_init_bin_data(glb_smear *in)
{
  int i;
  double med;
  in->bincenter=(double *) glb_malloc(sizeof(double)*in->numofbins);
  med=in->e_min;
  for(i=0;i<in->numofbins;i++)
    {
      med += in->binsize[i];
      in->bincenter[i]=med-in->binsize[i]*0.5;
    }
}

void glb_init_sim_bin_data(glb_smear *in)
{
  int i;
  double med;
  
  in->simbincenter=(double *) glb_malloc(sizeof(double)*in->simbins);

  med=in->e_sim_min;
  for(i=0;i<in->simbins;i++)
    {
      med += in->simbinsize[i];
      in->simbincenter[i]=med-in->simbinsize[i]*0.5;
    }
  
  
}

double glb_lower_bin_boundary(int bin, const glb_smear *data)
{
  return
    data->bincenter[bin] - data->binsize[bin] * 0.5;
}

double glb_upper_bin_boundary(int bin, const glb_smear *data)
{
  return
   data->bincenter[bin] + data->binsize[bin] * 0.5;
}

double glb_bin_center(int bin, const glb_smear *data)
{
  return
    data->bincenter[bin];
}


double glb_lower_sbin_boundary(int bin, const glb_smear *data)
{
  return
    data->simbincenter[bin] - data->simbinsize[bin] * 0.5;
}

double glb_upper_sbin_boundary(int bin, const glb_smear *data)
{
  return
   data->simbincenter[bin] + data->simbinsize[bin] * 0.5;
}

double glb_sbin_center(int bin, const glb_smear *data)
{
  return
    data->simbincenter[bin];
}

/* Safe functions do boundary checking, thus they are safe but
 * should not be used in the innermost loops.
 */


static double glb_safesbin_center(int bin, const glb_smear *data)
{
  if(bin>data->simbins) return data->e_sim_max;
  if(bin<0) return data->e_sim_min;
  return
    data->simbincenter[bin];
}


static double glb_safelower_bin_boundary(int bin, const glb_smear *data)
{
  if(bin>data->numofbins) return data->e_max;
  if(bin<0) return data->e_min;
  return
    data->bincenter[bin] - data->binsize[bin] * 0.5;
}




/* end */

static double Sigma(int bin, const glb_smear *data)
{
  return data->sig_f(glb_sbin_center(bin,data),data->sigma);
}

/* It really would be nice to have an algorithm for type A matrices
 * which allows different number of simbins and bins. It is possible
 * and already implemented.
 */ 

/*
  double SmearMatrixElementA(int i, int j, const glb_smear *data)
  {
  double p1,p2;
  p1=0.5*(1+gsl_sf_erf(
  (glb_upper_bin_boundary(j,data)-glb_bin_center(i,data))/
  (sqrt(2)*Sigma(i,data))
  ));
  p2=0.5*(1+gsl_sf_erf(
  (glb_lower_bin_boundary(j,data)-glb_bin_center(i,data))/
  (sqrt(2)*Sigma(i,data))
  ));
  return p1-p2;
}
*/

/* FIXME In computing the matrix element it is important to keep
 * the indices in correct ordering. We used to compute the transpose
 * of the smearing matrix -- now it is changed according to the formulas
 * in the documentation. For a equal number of bins and simbins the difference
 * is given by the variation of sigma, which changes from sigma(E') to sigma(E)
 * in the new version.
 *
 * Furthermore the correct scaling with the binsize and simbinsize
 * has been implemented. Simbinize is multiplied as a number onto
 * the matrix element (this has to move to calc.c) and divided by the binsize
 * to cancel the factor in calc.c which is erraneously there ans should
 * be simply replaced by simbinsize and then we can get rid of all
 * binsizes here. The actual binsize is taken into account by the integration
 * over the bin width with respect to E', which is already analytically taken
 * care of and conatined in the from of the matrix elements.
 */

static double SmearMatrixElementA(int i, int j, const glb_smear *data)
  {
  double p1,p2;
  p1=0.5*(1+gsl_sf_erf(
  (glb_upper_bin_boundary(i,data)-glb_sbin_center(j,data))/
  (sqrt(2)*Sigma(j,data))
  ));
  p2=0.5*(1+gsl_sf_erf(
  (glb_lower_bin_boundary(i,data)-glb_sbin_center(j,data))/
  (sqrt(2)*Sigma(j,data))
  ));
  /* FIXME return (p1-p2)*(data->simbinsize[j]/data->binsize[i]); */
  return (p1-p2);
}

static double Round(double in)
{
  double erg,x;
  if(in<0) x= -in;
  else x=in;
  if( fabs(x-floor(x))<0.5) erg=floor(x); 
     else erg=ceil(x);
  if(in<0) erg= -erg;
  return erg;

}

static void SetupSmearMatrixA(glb_smear *test,const struct experiment *head)
{
  int i;
  double med;
  if(test->e_sim_min==-1) 
    {
      if(head->simtresh==-1)
	test->e_sim_min=head->emin;
      else
	test->e_sim_min=head->simtresh;
    }

  if(test->e_sim_max==-1) 
    {    
      if(head->simbeam==-1)
	test->e_sim_max=head->emax;
      else
	test->e_sim_max=head->simbeam;
    }

  if(test->simbins==-1) 
    {
      if(head->simbins==-1)
	test->simbins=head->numofbins;
      else
	test->simbins=head->simbins;
    }

  if(test->simbinsize==NULL)
    {
      test->simbinsize=(double *) glb_malloc(sizeof(double)*test->simbins);
      med=(test->e_sim_max-test->e_sim_min)/test->simbins;
      for(i=0;i<test->simbins;i++) test->simbinsize[i]=med;
    }
 
  glb_init_sim_bin_data(test);
  return;
}


void glb_set_up_smear_data(glb_smear *test,const struct experiment *head)
{
  int i;
  double med;
  if(test->e_min==-1) test->e_min=head->emin;
  
  if(test->e_max==-1) test->e_max=head->emax;

  if(test->numofbins==-1) test->numofbins=head->numofbins;
     
  if(test->binsize==NULL)
    {
      test->binsize=(double *) glb_malloc(sizeof(double)*test->numofbins);
      med=(test->e_max-test->e_min)/test->numofbins;
      for(i=0;i<test->numofbins;i++) test->binsize[i]=med;
    }
 
  glb_init_bin_data(test);
  
  

  if(test->e_sim_min==-1) 
    {
      if(head->simtresh==-1)
	test->e_sim_min=head->emin;
      else
	test->e_sim_min=head->simtresh;
    }

  if(test->e_sim_max==-1) 
    {    
      if(head->simbeam==-1)
	test->e_sim_max=head->emax;
      else
	test->e_sim_max=head->simbeam;
    }

  if(test->simbins==-1) 
    {
      if(head->simbins==-1)
	test->simbins=head->numofbins;
      else
	test->simbins=head->simbins;
    }

  if(test->simbinsize==NULL)
    {
      test->simbinsize=(double *) glb_malloc(sizeof(double)*test->simbins);
      med=(test->e_sim_max-test->e_sim_min)/test->simbins;
      for(i=0;i<test->simbins;i++) test->simbinsize[i]=med;
    }
 
  glb_init_sim_bin_data(test);
  
  test->sig_f=NULL;
  test->num_of_params=0;
  test->type=0;  

  return;
}



static double** SmearMatrixA(glb_smear *data, int **lowrange, int **uprange,
		      const struct experiment *head)
{
  int i,j,nonzero,zero,prev;
  double erg;
  double **out=NULL;

  int *low,*up;

  SetupSmearMatrixA(data,head);
  if(data->type != GLB_TYPE_A)
    {
      glb_error("Type mismatch in SmearMatrixA"); 
      return NULL;
    }

  out=(double**) glb_malloc(sizeof(double*) * data->numofbins);
  up=(int*) glb_malloc(sizeof(int) * data->numofbins);
  low=(int*) glb_malloc(sizeof(int) * data->numofbins);
  
  for(i=0;i<data->numofbins;i++)
    {
      nonzero=0;
      zero=0;
      prev=0;
      out[i]=NULL;
      for(j=0;j<data->simbins;j++) 
	{
	  
	  erg=Round(1/TOL*SmearMatrixElementA(i,j,data))*TOL;
	  if(erg!=0)
	    { 
	      prev=1;
	      nonzero++;
	      out[i]=(double*) glb_realloc(out[i],sizeof(double)*nonzero);
	      out[i][nonzero-1]=erg;
	    }
	  else
	    {
	      if(prev==0) zero++;
	    }
	
	}

      low[i]=zero;
      up[i]=zero+nonzero-1;
     
      
     } 
  *lowrange=&low[0];
  *uprange=&up[0];
  return out;
}


/* Here comes the stuff for type B matrices. 
 * Clearly a somewhat more transparent design is needed and
 * the accuracy obtained with this algorithm is probably 
 * not really needed. Still the basic concept of the algorithm
 * is pretty nice and very flexible. Perhaps by waiving the
 * Filter feature and allowing an adaptive simbin spacing one
 * can optimize the algorithm and make it the default choice.
 */

static double chooz(double x, double *in)
{
  if(x<=0.0018) return 0.00005;
  else return 0.05*sqrt((x-0.0018+0.001)*1000)/1000.0;
}

/* This would die if the Filter feature is abandoned */
static double truncated_sigma(double x, const glb_smear *in)
{
  double sigma=in->sig_f(in->e_min,in->sigma)/(in->options->corr_fac);
  return GSL_MAX(
		 1E-5,
		 sqrt(fabs(pow(in->sig_f(x,in->sigma),2)-sigma*sigma))
		 );
}


/* This reflects the API for the GSL */
typedef struct
{
  glb_smear *in;
  int bin;
  double root;
  double intermediate;
} glb_strange_type;

/* Since the second argument of the integrand is a pointer to void ...*/
static double BinKernel(double y, void *bin )
{
  double www,bw,norm;
  int i;
  glb_strange_type temp;
  glb_smear *in;
  temp=*(glb_strange_type *) bin;
  i=temp.bin;
  in=temp.in;
  www=glb_bin_center(i,in);
  /* FIXME Here should be a factor 1/2 
   * This problem is now understood -- the algorithm for type B matrices
   * is nothing more than to evaluate in a complicated manner the summation
   * boundaries for type A matrices. The reason is: int Si(y) = Pi.
   * Thus the complicated type B stuff is superflous, besides its nice
   * filter feature. 
   *
   * PS.: The version below, where bw is the *half* binwidth is correct
   */
  bw=(glb_upper_bin_boundary(i,in)-glb_lower_bin_boundary(i,in))/2.0;
  norm=2.0*bw; 

 
  return
    (gsl_sf_erf((bw + www - y)/(sqrt(2.0)*truncated_sigma(y,in))) + 
     
     gsl_sf_erf((bw - www + y)/(sqrt(2.0)*truncated_sigma(y,in))))/(2.0*norm)
    ;
 
}


static double BinKernelM(double y, void *bin )
{
  return -BinKernel(y,bin);
}


static double Trueglb_bin_center(void *bin)
{
  /* That was a nice try -- but it is useless */
  
  int i;
  double a,b;
  glb_strange_type temp;
  glb_smear *in;
  int status;
  int iter = 0, max_iter = 1000;
  const gsl_min_fminimizer_type *T;
  gsl_min_fminimizer *s;
  double m=0;
  

  gsl_function F;
  temp=*(glb_strange_type *) bin;
  i=temp.bin;
  in=temp.in;    
  F.function = &BinKernelM;
  F.params = bin ;
  T = gsl_min_fminimizer_brent;
  s = gsl_min_fminimizer_alloc (T);

  a=glb_lower_bin_boundary(i,in);
  b=glb_upper_bin_boundary(i,in);
  m=glb_bin_center(i,in);
  /*
    gsl_min_fminimizer_set (s, &F, m, a, b);
    
    do
    {
    iter++;
    status = gsl_min_fminimizer_iterate (s);
    
    m = gsl_min_fminimizer_x_minimum (s);
    a = gsl_min_fminimizer_x_lower (s);
    b = gsl_min_fminimizer_x_upper (s);
    
    status 
    = gsl_min_test_interval (a, b, TOL, TOL);
    }
    while (status == GSL_CONTINUE && iter < max_iter);
  */
  gsl_min_fminimizer_free(s);
  return m;
}


/* ...strange but it works. */
static double BKi(double y, void *in)
{
  double x_l, x_h;
  double result, error;
  glb_strange_type temp;
  
  gsl_function F;
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
  temp=*(glb_strange_type *) in;  
  //temp.in=in;
  //temp.bin=i;
  F.function=&BinKernel;
  F.params=&temp;
  x_l=glb_bin_center(temp.bin,temp.in);
  x_l=Trueglb_bin_center(in);
  x_h=y;
 
    gsl_integration_qag(&F,x_l,x_h,TOL,TOL,1000,1,w,&result,&error);

  gsl_integration_workspace_free(w);
  return result;
}

static double  BinKernelR(double y, void *in)
{
  glb_strange_type temp;
  temp=*(glb_strange_type *) in;  
  return BinKernel(y,in)-temp.root;
}

static double BMK(double y, void *in)
{
  int iter,status;
  double result, error,bw,f_lo,f_hi;
  double x_hi,x_lo,erg;
  glb_strange_type temp;
  glb_strange_type *pt;
  gsl_function F;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  
  iter=0;
  T=gsl_root_fsolver_brent;
  s=gsl_root_fsolver_alloc(T);
  temp=*(glb_strange_type *) in;
  temp.root=BinKernel(y,in);
  
  F.function=&BinKernelR;
  F.params=&temp;
  /* Factor 2 error FIXME */
  bw=glb_upper_bin_boundary(temp.bin,temp.in)-glb_lower_bin_boundary(temp.bin,temp.in);

  x_lo=Trueglb_bin_center(in);
  /* heuristic */
  x_hi=glb_bin_center(temp.bin,temp.in)+bw+
    3*truncated_sigma(glb_bin_center(temp.bin,temp.in),temp.in);
 
  /* Here we make sure that x_lo and x_hi really bracket the zero */
  f_lo=BinKernelR(x_lo,&temp);
  f_hi=BinKernelR(x_hi,&temp);
  
  /* Maybe one of the interval boundaries is the root */
  if(fabs(f_hi)<TOL) return x_hi;
  if(fabs(f_lo)<TOL) return x_lo;

  while((f_lo<0&&f_hi<0)||(f_lo>0&&f_hi>0))
    {
      x_hi = x_hi + truncated_sigma(glb_bin_center(temp.bin,temp.in),temp.in); 
      if(x_hi>temp.in->options->up_bound) 
	{
	  glb_fatal("Conflict with up_bound");
	  return -1; 
	  
	}
      f_hi=BinKernelR(x_hi,&temp);
    }

 
  
  gsl_root_fsolver_set(s,&F,x_lo,x_hi);
  do
    {
      iter++;
      status = gsl_root_fsolver_iterate (s);   
      result = gsl_root_fsolver_root (s);
      x_lo = gsl_root_fsolver_x_lower (s);
      x_hi = gsl_root_fsolver_x_upper (s);  
     
      status = gsl_root_test_interval (x_lo, x_hi,TOL,TOL);

    }
  while (status == GSL_CONTINUE && iter < 1000);
  gsl_root_fsolver_free(s);
  erg= BKi(result,in)-BKi(y,in);
  pt=(glb_strange_type *) in;
  pt->intermediate=result; 
  return erg;
}

static double BMKR(double y, void *in)
{
  glb_strange_type temp;
  temp = *(glb_strange_type *) in;  
  return BMK(y,in)-temp.root;
}

static double BinRange(double cl, double *low, double *up, double *area, void *in)
{
  int iter,status;
  double result, error,bw;
  double x_hi,x_lo,f_lo,f_hi;
  glb_strange_type temp;
  gsl_function F;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  
  iter=0;
  T=gsl_root_fsolver_brent;
  s=gsl_root_fsolver_alloc(T);
  temp=*(glb_strange_type *) in;
  
  

  temp.root=cl;

  F.function=&BMKR;
  F.params=&temp;
  /* Factor 2 error FIXME */
  bw=glb_upper_bin_boundary(temp.bin,temp.in)-glb_lower_bin_boundary(temp.bin,temp.in);
  x_hi=Trueglb_bin_center(in);
  x_lo=glb_bin_center(temp.bin,temp.in)-bw
    -3*truncated_sigma(glb_bin_center(temp.bin,temp.in),temp.in);
  /* Here we make sure that x_lo and x_hi really bracket the zero */
  f_lo=BMKR(x_lo,&temp);
  f_hi=BMKR(x_hi,&temp);
  
  /* Maybe one of the interval boundaries is the root */
  if(fabs(f_hi)<TOL) return x_hi;
  if(fabs(f_lo)<TOL) return x_lo;

  while((f_lo<0&&f_hi<0)||(f_lo>0&&f_hi>0))
    {
      x_lo= x_lo - truncated_sigma(glb_bin_center(temp.bin,temp.in),temp.in); 
     
        if(x_lo<temp.in->options->low_bound) {
	    glb_fatal("Conflict with low_bound"); 
	return -1;
	}
      f_lo=BMKR(x_lo,&temp);
    }
 
  gsl_root_fsolver_set(s,&F,x_lo,x_hi); 
  
  
  do
    {
      iter++;
      status = gsl_root_fsolver_iterate (s);   
      result = gsl_root_fsolver_root (s);
      x_lo = gsl_root_fsolver_x_lower (s);
      x_hi = gsl_root_fsolver_x_upper (s);
      status = gsl_root_test_interval (x_lo, x_hi,TOL, TOL);
    }
  while (status == GSL_CONTINUE && iter < 1000);
  gsl_root_fsolver_free(s);
  *low=result;
  *area=BMK(result,&temp);
  *up=temp.intermediate;
  return result;
}


static void SetupSmearMatrixB(glb_smear *test, const struct experiment *head)
{
  glb_strange_type t;
  int i;
  double low,up,area,erg,med;
  double sigma=test->sig_f(test->e_min,test->sigma)/(test->options->corr_fac);
  t.in=test;

  t.bin=0;
  erg=BinRange(test->options->confidence_level,&low,&up,&area, &t);
  if(test->e_sim_min==-1) 
    {
      if(head->simtresh==-1)
	test->e_sim_min=GSL_MAX(test->options->low_bound,low);
      else 
	test->e_sim_min=head->simtresh;
    }
  t.bin=test->numofbins - 1 ;
  erg=BinRange(test->options->confidence_level,&low,&up,&area, &t);
  if(test->e_sim_max==-1) 
    {
      if(head->simbeam==-1)
	test->e_sim_max=GSL_MIN(test->options->up_bound,up); 
      else
	test->e_sim_max=head->simbeam;
    }

  if(test->simbins==-1) 
    {
      if(head->simbins==-1)
	test->simbins=(int) floor((test->e_sim_max - test->e_sim_min)/sigma)+1;
      else
	test->simbins=head->simbins;
    }
  if(test->simbinsize==NULL)
    {
      test->simbinsize=(double *) glb_malloc(sizeof(double)*test->simbins);
      med=(test->e_sim_max-test->e_sim_min)/test->simbins;
      for(i=0;i<test->simbins;i++) test->simbinsize[i]=med;
    }
  glb_init_sim_bin_data(test);
  return;
}



static double SmearMatrixElementB(int i, int j, double low, 
			   double up, double area, glb_smear *data)
{
  glb_strange_type t;
  int j_min, j_max,k;
  double point,erg,s,l,u;
  t.in=data;
  t.bin=i;
  point=glb_sbin_center(j,data);

  s=glb_upper_sbin_boundary(j,data)-glb_lower_sbin_boundary(j,data); 
  
  /* Here comes a somewhat tricky part in order to implement the
   * offset feature as in CalcSmearingB. I am no sure wether I would
   * like to keep this.
   * even more tricky now since no overrun is allowed !
   */

  k=0;l=0;
  while(l<=low)
    {
      l=glb_sbin_center(k,data);
      k++;
    }
  l=GSL_MAX(data->options->low_bound,
		glb_safesbin_center(k-2-data->options->offset,data));
  u=0;k=0;
  while(u<=up)
    { 
      u=glb_sbin_center(k,data);
      k++;
    }
  u=GSL_MIN(data->options->up_bound,
	    glb_safesbin_center(k-1+data->options->offset,data));

  /* end */
  
  if(point > l && point < u)
    { 
      
      erg=s/M_PI*(
		  gsl_sf_Si((up-point)* M_PI /s)
		  -gsl_sf_Si((low-point)* M_PI /s)
		  )*BinKernel(point,&t)/area;
     
    }
  else erg=0.0;
 
  return erg;
}

static double** SmearMatrixB( glb_smear *data, int **lowrange, int **uprange, 
		       const struct experiment *head)
{
  int i,j,nonzero,zero,prev;
  double erg;
  double **out=NULL;

  glb_strange_type t;
 
  double low_b,up_b,area;
  int *low,*up;

  if(data->type != GLB_TYPE_B)
    {
      glb_error("Type mismatch in SmearMatrixB\n"); 
      return NULL;
    }
  SetupSmearMatrixB(data,head);
  t.in=data;
  
  

  out=(double**) malloc(sizeof(double*) * data->numofbins);
  up=(int*) malloc(sizeof(int) * data->numofbins);
  low=(int*) malloc(sizeof(int) * data->numofbins);
  
  for(i=0;i<data->numofbins;i++)
    {
      nonzero=0;
      zero=0;
      prev=0;
      out[i]=NULL;
      t.bin=i;
      BinRange(data->options->confidence_level,&low_b,&up_b,&area, &t);
      for(j=0;j<data->simbins;j++) 
	{
	  erg=SmearMatrixElementB(i,j,low_b,up_b,area,data);

	  if(erg!=0.0)
	    {
	      prev=1;
	      nonzero++;
	      out[i]=(double*) realloc(out[i],sizeof(double)*nonzero);
	      out[i][nonzero-1]=erg;
	    }
	  else
	    {
	      if(prev==0) zero++;
	    }
	
	}

      low[i]=zero;
      up[i]=zero+nonzero-1;
     
      
     } 
  *lowrange=&low[0];
  *uprange=&up[0];
  return out;
}


/* Type C matrices should combined the advantages of A&B types */
static double BinKernelC(double y, void *bin )
{
  double www,bw,norm;
  int i;
  glb_strange_type temp;
  glb_smear *in;
  temp=*(glb_strange_type *) bin;
  i=temp.bin;
  in=temp.in;
  www=glb_bin_center(i,in);
  /* FIXME Here should be a factor 1/2 */
  bw=(glb_upper_bin_boundary(i,in)-glb_lower_bin_boundary(i,in))/2.0;
  norm=2.0*bw; 
  /* This factor 2 is however irrelevant -- sometimes*/
 
  return
    (gsl_sf_erf((bw + www - y)/(sqrt(2.0)*in->sig_f(y,in->sigma)))+ 
     
     gsl_sf_erf((bw - www + y)/(sqrt(2.0)*in->sig_f(y,in->sigma))))/(2.0*norm)
    ;
 
}

static double SmearMatrixElementC(int i, int j, glb_smear *data)
{
  glb_strange_type t;
  int j_min, j_max,k;
  double point,erg,s,l,u;
  t.in=data;
  t.bin=i;
  point=glb_sbin_center(j,data);

  s=glb_upper_sbin_boundary(j,data)-glb_lower_sbin_boundary(j,data); 
  erg=s/M_PI*(M_PI)*BinKernelC(point,&t);     
  return erg;
}


static double** SmearMatrixC(glb_smear *data, int **lowrange, int **uprange,
		      const struct experiment *head)
{
  int i,j,nonzero,zero,prev;
  double erg;
  double **out=NULL;

  int *low,*up;

  SetupSmearMatrixA(data,head);
  if(data->type != GLB_TYPE_C)
    {
      glb_error("Type mismatch in SmearMatrixC"); 
      return NULL;
    }

  out=(double**) malloc(sizeof(double*) * data->numofbins);
  up=(int*) malloc(sizeof(int) * data->numofbins);
  low=(int*) malloc(sizeof(int) * data->numofbins);
  
  for(i=0;i<data->numofbins;i++)
    {
      nonzero=0;
      zero=0;
      prev=0;
      out[i]=NULL;
      for(j=0;j<data->simbins;j++) 
	{
	  
	  erg=Round(1/TOL*SmearMatrixElementC(i,j,data))*TOL;
	  
	  if(erg!=0)
	    { 
	      prev=1;
	      nonzero++;
	      out[i]=(double*) realloc(out[i],sizeof(double)*nonzero);
	      out[i][nonzero-1]=erg;
	    }
	  else
	    {
	      if(prev==0) zero++;
	    }
	
	}

      low[i]=zero;
      up[i]=zero+nonzero-1;
     
      
     } 
  *lowrange=&low[0];
  *uprange=&up[0];
  return out;
}


void glb_compute_smearing_matrix(double ***matrix, 
			  int **low, int **up, glb_smear *in
			   , const struct experiment *head)
{

  if(in->type==GLB_TYPE_A) *matrix=SmearMatrixA(in,low,up,head);
  if(in->type==GLB_TYPE_B) {
    glb_warning("================> filter_value & filter_state !");
    *matrix=SmearMatrixB(in,low,up,head);
  }
  if(in->type==GLB_TYPE_C) *matrix=SmearMatrixC(in,low,up,head);
  
}



#ifdef GLB_SMEAR_OWN 

/*-----------------------------------------------------------------*/
int
main (void)
{
  int i,j;
  glb_strange_type t;   
  int *low,*up;
  double **matrix;
  
  double sm[]={0.15,0,0};
  
  /* That is Reactor */  

  glb_option_type *opt;
  glb_smear *testR;
  
  /* That is needed for the new defaultizing process */
  struct experiment inl;
  inl.simtresh=-1;
  inl.simbeam=-1;
  inl.simbins=-1;
 
  testR=glb_smear_alloc();
  opt=glb_option_type_alloc();
  
  
  opt->corr_fac=1.0;
  opt->confidence_level=1-1E-3;
  opt->offset=5;
  opt->low_bound=.0015;
  opt->up_bound=0.01;

  testR->type=GLB_TYPE_B;
  testR->numofbins=62;
  testR->e_min=0.0018; 
  testR->e_max=0.008;
 

  testR->num_of_params=0;
  testR->sig_f=&chooz;
  testR->options=opt;
  
  glb_default_smear(testR);
  
  /* Here an example is defined -- Corresponding to NuFact */
  
  
  //glb_smear testNF = {GLB_TYPE_A,20,20,4,50,4,50,&sm[0],&my_sigma,NULL,NULL,opt};
  /* 
  opt->corr_fac=1.0;
  opt->confidence_level=1-1E-3;
  opt->offset=1;
  opt->low_bound=4;
  opt->up_bound=50;

  testR->type=GLB_TYPE_B;
  testR->numofbins=20;
  testR->e_min=4.0; 
  testR->e_max=50.0;
  testR->num_of_params=3;
  testR->sig_f=&my_sigma;
  testR->sigma=&sm[0];
  testR->options=opt;

  glb_default_smear(testR);*/
  /* That´s the call to our new routine */

  matrix=SmearMatrixB(testR,&low,&up,&inl);

  
  /* Just pretty printing to stdout */
     
  for(j=0;j<testR->numofbins;j++)
    {
      printf("%d %d ",low[j],up[j]);
      for(i=0;i<=up[j]-low[j];i++) printf("%10.10g ", matrix[j][i]);
      printf("\n");
    }
  
 
  glb_smear_free(testR);
 
 
  
  return 0;
}

#endif /* GLB_SMEAR_OWN */
