#define TOL 2.0e-4

int ncom2=0;	/* defining declarations */
double *pcom2=0,*xicom2=0,(*nrfunc2)();

double f1dim2(x)
double x;
{
	int j;
	double f,*xt,*glb_alloc_vec();
	void glb_free_vec();

	xt=glb_alloc_vec(1,ncom2);
	for (j=1;j<=ncom2;j++) xt[j]=pcom2[j]+x*xicom2[j];
	f=(*nrfunc2)(xt);
	glb_free_vec(xt,1,ncom2);
	return f;
}

void linmin2(p,xi,n,fret,func)
double p[],xi[],*fret,(*func)();
int n;
{
	int j;
	double xx,xmin,fx,fb,fa,bx,ax;
	double brent(),f1dim2(),*glb_alloc_vec();
	void mnbrak(),glb_free_vec();

	ncom2=n;
	pcom2=glb_alloc_vec(1,n);
	xicom2=glb_alloc_vec(1,n);
	nrfunc2=func;
	for (j=1;j<=n;j++) {
		pcom2[j]=p[j];
		xicom2[j]=xi[j];
	}
	ax=0.0;
	xx=1.0;
	bx=2.0;
	mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,f1dim2);
	*fret=brent(ax,xx,bx,f1dim2,TOL,&xmin);
	for (j=1;j<=n;j++) {
		xi[j] *= xmin;
		p[j] += xi[j];
	}
	glb_free_vec(xicom2,1,n);
	glb_free_vec(pcom2,1,n);
}

#undef TOL
