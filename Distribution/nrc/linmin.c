#define TOL 2.0e-4

int ncom=0;	/* defining declarations */
double *pcom=0,*xicom=0,(*nrfunc)();

void linmin(p,xi,n,fret,func)
double p[],xi[],*fret,(*func)();
int n;
{
	int j;
	double xx,xmin,fx,fb,fa,bx,ax;
	double brent(),f1dim(),*glb_alloc_vec();
	void mnbrak(),glb_free_vec();

	ncom=n;
	pcom=glb_alloc_vec(1,n);
	xicom=glb_alloc_vec(1,n);
	nrfunc=func;
	for (j=1;j<=n;j++) {
		pcom[j]=p[j];
		xicom[j]=xi[j];
	}
	ax=0.0;
	xx=1.0;
	bx=2.0;
	mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,f1dim);
	*fret=brent(ax,xx,bx,f1dim,TOL,&xmin);
	for (j=1;j<=n;j++) {
		xi[j] *= xmin;
		p[j] += xi[j];
	}
	glb_free_vec(xicom,1,n);
	glb_free_vec(pcom,1,n);
}

#undef TOL
