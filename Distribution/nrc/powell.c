#include <math.h>

#define ITMAX 500
static double sqrarg;
#define SQR(a) (sqrarg=(a),sqrarg*sqrarg)

void powell(p,xi,n,ftol,iter,fret,func)
double p[],**xi,ftol,*fret,(*func)();
int n,*iter;
{
	int i,ibig,j;
	double t,fptt,fp,del;
	double *pt,*ptt,*xit,*glb_alloc_vec();
	void linmin(),glb_minimizer_error(),glb_free_vec();

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
			linmin(p,xit,n,fret,func);
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
		if (*iter == ITMAX) glb_minimizer_error("Too many iterations in routine POWELL");
		for (j=1;j<=n;j++) {
			ptt[j]=2.0*p[j]-pt[j];
			xit[j]=p[j]-pt[j];
			pt[j]=p[j];
		}
		fptt=(*func)(ptt);
		if (fptt < fp) {
			t=2.0*(fp-2.0*(*fret)+fptt)*SQR(fp-(*fret)-del)-del*SQR(fp-fptt);
			if (t < 0.0) {
				linmin(p,xit,n,fret,func);
				for (j=1;j<=n;j++) xi[j][ibig]=xit[j];
			}
		}
	}
}

#undef ITMAX
#undef SQR
