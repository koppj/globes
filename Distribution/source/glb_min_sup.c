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



#include <math.h>

#include <stdio.h>

#include "glb_error.h"
#include "glb_min_sup.h"

static void glb_minimizer_error(char error_text[])
{
  glb_fatal(error_text);  
}

double *glb_alloc_vec(int nl,int nh)
{
	double *v;
	v=(double *)glb_malloc((unsigned) (nh-nl+1)*sizeof(double));
	return v-nl;
}

double **glb_alloc_mat(int nrl,int nrh, int ncl,int nch)
{
	int i;
	double **m;

	m=(double **) glb_malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
	m -= nrl;
	for(i=nrl;i<=nrh;i++) 
	  {
		m[i]=(double *) glb_malloc((unsigned) 
					   (nch-ncl+1)*sizeof(double));
		m[i] -= ncl;
	  }
	return m;
}

void glb_free_vec(double *v,int nl,int nh)
{
	glb_free((char*) (v+nl));
}

void glb_free_mat(double **m,int nrl,int nrh,int ncl,int nch)
{
	int i;
	for(i=nrh;i>=nrl;i--) glb_free((char*) (m[i]+ncl));
	glb_free((char*) (m+nrl));
}




#define TOL 2.0e-4


typedef struct {
  int ncom;
  double *pcom;
  double *xicom;
  double (*nrfunc)(double*);
} glb_min_data;


static double one_dim_projection(double x,glb_min_data *in)
{
	int j;
	double f,*xt;

	xt=glb_alloc_vec(1,in->ncom);
	for (j=1;j<=in->ncom;j++) xt[j]=in->pcom[j]+x*in->xicom[j];
	f=in->nrfunc(xt);
	glb_free_vec(xt,1,in->ncom);
	return f;
}



#define ITMAX 100
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define SIGN(a,b) ((b) > 0.0 ? fabs(a) : -fabs(a))
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

static double glb_brent_min(double ax,double bx,double cx,
			    double (*f)(double,glb_min_data*),
			    double tol,double *xmin,glb_min_data *in)
{
	int iter;
	double a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
	double e=0.0;

	a=((ax < cx) ? ax : cx);
	b=((ax > cx) ? ax : cx);
	x=w=v=bx;
	fw=fv=fx=(*f)(x,in);
	for (iter=1;iter<=ITMAX;iter++) {
		xm=0.5*(a+b);
		tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
		if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
			*xmin=x;
			return fx;
		}
		if (fabs(e) > tol1) {
			r=(x-w)*(fx-fv);
			q=(x-v)*(fx-fw);
			p=(x-v)*q-(x-w)*r;
			q=2.0*(q-r);
			if (q > 0.0) p = -p;
			q=fabs(q);
			etemp=e;
			e=d;
			if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
				d=CGOLD*(e=(x >= xm ? a-x : b-x));
			else {
				d=p/q;
				u=x+d;
				if (u-a < tol2 || b-u < tol2)
					d=SIGN(tol1,xm-x);
			}
		} else {
			d=CGOLD*(e=(x >= xm ? a-x : b-x));
		}
		u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
		fu=(*f)(u,in);
		if (fu <= fx) {
			if (u >= x) a=x; else b=x;
			SHFT(v,w,x,u)
			SHFT(fv,fw,fx,fu)
		} else {
			if (u < x) a=u; else b=u;
			if (fu <= fw || w == x) {
				v=w;
				w=u;
				fv=fw;
				fw=fu;
			} else if (fu <= fv || v == x || v == w) {
				v=u;
				fv=fu;
			}
		}
	}
	glb_minimizer_error("glb_brent_min");
	*xmin=x;
	return fx;
}

#undef ITMAX
#undef CGOLD
#undef ZEPS
#undef SIGN




#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define SIGN(a,b) ((b) > 0.0 ? fabs(a) : -fabs(a))
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

static void bracket(double *ax,double *bx,double *cx,double *fa,double *fb,
	    double *fc, double (*func)(double,glb_min_data*),glb_min_data *in)
{
	double ulim,u,r,q,fu,dum;

	*fa=(*func)(*ax,in);
	*fb=(*func)(*bx,in);
	if (*fb > *fa) {
		SHFT(dum,*ax,*bx,dum)
		SHFT(dum,*fb,*fa,dum)
	}
	*cx=(*bx)+GOLD*(*bx-*ax);
	*fc=(*func)(*cx,in);
	while (*fb > *fc) {
		r=(*bx-*ax)*(*fb-*fc);
		q=(*bx-*cx)*(*fb-*fa);
		u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
			(2.0*SIGN(MAX(fabs(q-r),TINY),q-r));
		ulim=(*bx)+GLIMIT*(*cx-*bx);
		if ((*bx-u)*(u-*cx) > 0.0) {
			fu=(*func)(u,in);
			if (fu < *fc) {
				*ax=(*bx);
				*bx=u;
				*fa=(*fb);
				*fb=fu;
				return;
			} else if (fu > *fb) {
				*cx=u;
				*fc=fu;
				return;
			}
			u=(*cx)+GOLD*(*cx-*bx);
			fu=(*func)(u,in);
		} else if ((*cx-u)*(u-ulim) > 0.0) {
			fu=(*func)(u,in);
			if (fu < *fc) {
				SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
				SHFT(*fb,*fc,fu,(*func)(u,in))
			}
		} else if ((u-ulim)*(ulim-*cx) >= 0.0) {
			u=ulim;
			fu=(*func)(u,in);
		} else {
			u=(*cx)+GOLD*(*cx-*bx);
			fu=(*func)(u,in);
		}
		SHFT(*ax,*bx,*cx,u)
		SHFT(*fa,*fb,*fc,fu)
	}
}



#undef GOLD
#undef GLIMIT
#undef TINY
#undef MAX
#undef SIGN
#undef SHFT


static void line_minimization(double p[],double xi[],int n,double *fret,
			      double (*func)(double*),
			      glb_min_data *in)
{
	int j;
	double xx,xmin,fx,fb,fa,bx,ax;
	int ncom=0;
	double *pcom=0,*xicom=0,(*nrfunc)(double*);

	ncom=n;
	pcom=glb_alloc_vec(1,n);
	xicom=glb_alloc_vec(1,n);
	nrfunc=func;
	for (j=1;j<=n;j++) {
		pcom[j]=p[j];
		xicom[j]=xi[j];
	}

	in->ncom=n;
	in->pcom=pcom;
	in->xicom=xicom;
	in->nrfunc=func;
	ax=0.0;
	xx=1.0;
	bx=2.0;
	bracket(&ax,&xx,&bx,&fa,&fx,&fb,one_dim_projection,in);
	*fret=glb_brent_min(ax,xx,bx,one_dim_projection,TOL,&xmin,in);
	for (j=1;j<=n;j++) {
		xi[j] *= xmin;
		p[j] += xi[j];
	}
	glb_free_vec(xicom,1,n);
	glb_free_vec(pcom,1,n);
}







#define ITMAX 500
static double sqrarg;
#define SQR(a) (sqrarg=(a),sqrarg*sqrarg)

void glb_powell(double p[],double **xi,int n,
		double ftol,int *iter,double *fret,
		double (*func)(double*))
{
  int i,ibig,j;
  double t,fptt,fp,del;
  double *pt,*ptt,*xit;
  glb_min_data in;
  pt=glb_alloc_vec(1,n);
  ptt=glb_alloc_vec(1,n);
  xit=glb_alloc_vec(1,n);
	*fret=(*func)(p);
	for (j=1;j<=n;j++) pt[j]=p[j];
	for (*iter=1;;(*iter)++) {
		fp=(*fret);
		ibig=0;
		del=0.0;
		for (i=1;i<=n;i++) {
			for (j=1;j<=n;j++) xit[j]=xi[j][i];
			fptt=(*fret);
			line_minimization(p,xit,n,fret,func,&in);
			if (fabs(fptt-(*fret)) > del) {
				del=fabs(fptt-(*fret));
				ibig=i;
			}
		}
		if (2.0*fabs(fp-(*fret)) <= ftol*(fabs(fp)+fabs(*fret))) {
			glb_free_vec(xit,1,n);
			glb_free_vec(ptt,1,n);
			glb_free_vec(pt,1,n);
			return;
		}
		if (*iter == ITMAX) glb_minimizer_error("empty");
		for (j=1;j<=n;j++) {
			ptt[j]=2.0*p[j]-pt[j];
			xit[j]=p[j]-pt[j];
			pt[j]=p[j];
		}
		fptt=(*func)(ptt);
		if (fptt < fp) {
			t=2.0*(fp-2.0*(*fret)+fptt)*SQR(fp-(*fret)-del)-del*SQR(fp-fptt);
			if (t < 0.0) {
				line_minimization(p,xit,n,fret,func,&in);
				for (j=1;j<=n;j++) xi[j][ibig]=xit[j];
			}
		}
	}
}


void glb_powell2(double p[],double **xi,int n,double ftol,
		 int *iter,double *fret,
		 double (*func)(double*))
{
  int i,ibig,j;
  double t,fptt,fp,del;
  double *pt,*ptt,*xit;
  glb_min_data in;
  pt=glb_alloc_vec(1,n);
  ptt=glb_alloc_vec(1,n);
  xit=glb_alloc_vec(1,n);
	*fret=(*func)(p);
	for (j=1;j<=n;j++) pt[j]=p[j];
	for (*iter=1;;(*iter)++) {
		fp=(*fret);
		ibig=0;
		del=0.0;
		for (i=1;i<=n;i++) {
			for (j=1;j<=n;j++) xit[j]=xi[j][i];
			fptt=(*fret);
			line_minimization(p,xit,n,fret,func,&in);
			if (fabs(fptt-(*fret)) > del) {
				del=fabs(fptt-(*fret));
				ibig=i;
			}
		}
		if (2.0*fabs(fp-(*fret)) <= ftol*(fabs(fp)+fabs(*fret))) {
			glb_free_vec(xit,1,n);
			glb_free_vec(ptt,1,n);
			glb_free_vec(pt,1,n);
			return;
		}
		if (*iter == ITMAX) glb_minimizer_error("empty");
		for (j=1;j<=n;j++) {
			ptt[j]=2.0*p[j]-pt[j];
			xit[j]=p[j]-pt[j];
			pt[j]=p[j];
		}
		fptt=(*func)(ptt);
		if (fptt < fp) {
			t=2.0*(fp-2.0*(*fret)+fptt)*SQR(fp-(*fret)-del)-del*SQR(fp-fptt);
			if (t < 0.0) {
				line_minimization(p,xit,n,fret,func,&in);
				for (j=1;j<=n;j++) xi[j][ibig]=xit[j];
			}
		}
	}
}



#undef ITMAX
#undef SQR
#undef TOL
