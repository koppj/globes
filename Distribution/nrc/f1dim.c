extern int ncom;	/* defined in LINMIN */
extern double *pcom,*xicom,(*nrfunc)();

double f1dim(x)
double x;
{
	int j;
	double f,*xt,*glb_alloc_vec();
	void glb_free_vec();

	xt=glb_alloc_vec(1,ncom);
	for (j=1;j<=ncom;j++) xt[j]=pcom[j]+x*xicom[j];
	f=(*nrfunc)(xt);
	glb_free_vec(xt,1,ncom);
	return f;
}
