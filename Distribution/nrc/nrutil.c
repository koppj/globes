
#include <malloc.h>
#include <stdio.h>

void glb_minimizer_error(char error_text[])
{
  glb_fatal("Error in minimizer");  
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

