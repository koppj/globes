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




#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex>



// compatability with gcc3.0 and higher
using namespace std;


extern "C"{ 
#undef GLS_ERROR_H
#include "glb_error.h"
}
//interestingly enough this does *not* work 

#define FLOAT double

//Phase factor filtering for constant densities

static double filter=1E9;
static int filter_state=0;


// interface to zgeev_ of netlib 
extern "C" void zgeev_(char* jobvl, char* jobvr, int* n,
 complex<FLOAT>* a, int* lda, complex<FLOAT>* w,
 complex<FLOAT>* vl, int* ldvl, complex<FLOAT>* vr, int* ldvr,
 complex<FLOAT>* work, int* lwork, FLOAT* rwork, int* info); 



// fix the libf2c problem of the SUSE distribution
int MAIN__;



// complex mixing matrix in vacuum
static complex<FLOAT> u[3][3];

// complex mixing matrix in matter
static complex<FLOAT> umat[3][3];

// complex matrix to be diagonalized
static complex<FLOAT> md[3][3];

// squared masses in vacuum
static double mq[3];



// vacuum parameters
static double t12,t23,t13,d1;



// Matter density in g/cm^3
static double matter_density;

// Neutrino energy in MeV
static double neutrino_energy;

// Anti-neutrinos (-1) Neutrinos(+1)
static int anti_neutrinos = 1;

//solar mixing angle which may be changed but will not be varied
static double solar_mixing;

// working arrays and variables for dsyev
#define LWORK 100

static char jobvl='N';
static char jobvr='V';
static int n=3;
static complex<FLOAT> a[3][3];
static complex<FLOAT> w[3];
static complex<FLOAT> vl[3][3];
static complex<FLOAT> vr[3][3];
static complex<FLOAT> work[LWORK];
static FLOAT rwork[6];
static int lda=3;
static int ldvl=3;
static int ldvr=3;
static int lwork=LWORK;
static int info;

// Matter profile storage
static double* length=NULL;    // Pointer on array of length values of the matter profile
static double* density=NULL;   // Pointer on array of density values of the matter profile
static int psteps=0;        // number of layers in the matter profile





// -----------------------------------------------------------
// Diagonalization of the hamiltonian
// (based on LAPACK routine cgeev) 
// -----------------------------------------------------------
static int ComputeMixingMatrix()
{
	int i,j,k;


	// U*diag(mq1,mq2,mq3,mq4)*dagger(U) + diag(AMSW,0,0) + NSI
	for (i=0; i<3; i++)
	{
		for (j=0; j<3; j++)
		{
			md[i][j] = complex<FLOAT>(0,0);
			for (k=0; k<3; k++)
				md[i][j] +=	u[i][k]*mq[k]*conj(u[j][k]);
		}
	}
	 
	md[0][0] += 0.1*1.5*anti_neutrinos*(10E-4)*neutrino_energy*0.5*matter_density;
		  
		zgeev_(&jobvl,&jobvr,&n,md[0],&lda,w,vl[0],&ldvl,vr[0],&ldvr,work,&lwork,rwork,&info);

	return info;
}


// -----------------------------------------------------------
// Propagation of the neutrino flavour state vector through a 
// layer of matter
// -----------------------------------------------------------
// psi         neutrino state vector in flavour basis
// panti       neutrino or antineutrino (1,-1)
// length      layer length in km
// rho         matter density in g/cm3
// en          neutrino energy in GeV
// -----------------------------------------------------------
static void PropagateInMatter(complex<FLOAT>* psi, int panti, FLOAT length, FLOAT en, FLOAT rho)
{
	complex<FLOAT> newpsi[3];
	FLOAT kappa = 2*1.2669*length/en;

	neutrino_energy = en;
	matter_density = rho;
	int info = ComputeMixingMatrix();
	
	if (info!=0)
	{
		fprintf(stderr,"dsyev error\n");
		exit(1);
	}
	
	
	if (panti==1)
	{
		// Neutrinos
		for (int i=0; i<3; i++)
		{
			newpsi[i] = complex<FLOAT>(0,0);
			for (int j=0; j<3; j++)
			{
				for (int k=0; k<3; k++)
				{
					newpsi[i] +=
					conj(vr[k][i])*exp(complex<FLOAT>(0,-1)*kappa*w[k])*vr[k][j]*psi[j];
				}
			}	 
		}
	}
	else
	{
		// Antineutrinos
		for (int i=0; i<3; i++)
		{
			newpsi[i] = complex<FLOAT>(0,0);
			for (int j=0; j<3; j++)
			{
				for (int k=0; k<3; k++)
				{
					newpsi[i] +=
					vr[k][i]*exp(complex<FLOAT>(0,-1)*kappa*w[k])*conj(vr[k][j])*psi[j];
				}
			}	 
		}
	}
	
	
	for (int i=0; i<3; i++) psi[i] = newpsi[i];
}




// -----------------------------------------------------------
// Oscillation probability in vacuum
// -----------------------------------------------------------
// pl        initial flavour
// pm        final flavour
// panti     neutrino or antineutrino (1,-1)
// pen       neutrino energy in GeV
// plength   layer length in km
// -----------------------------------------------------------
extern "C" FLOAT glbVacuumProbability(int pl, int pm, int panti,double pen, double plength)
{	
	complex<FLOAT> psi0[3];
	complex<FLOAT> psi1[3];
	complex<FLOAT> newpsi[3];
	complex<FLOAT> amplitude;
	FLOAT probability;

	// Ausgangs und Endzustand (Flavour!!)
	for (int i=0; i<3; i++)
	{
		psi0[i] = complex<FLOAT>(0,0);
		psi1[i] = complex<FLOAT>(0,0);
	}
	psi0[pl-1]=complex<FLOAT>(1,0);
	psi1[pm-1]=complex<FLOAT>(1,0);


	// Propagation
	FLOAT kappa = 2*1.2669*plength/pen;
	
	if (panti==1)
	{
		// Neutrinos
		for (int i=0; i<3; i++)
		{
			newpsi[i] = complex<FLOAT>(0,0);
			for (int j=0; j<3; j++)
			{
				for (int k=0; k<3; k++)
				{
					newpsi[i] += u[i][k]*exp(complex<FLOAT>(0,-1)*kappa*mq[k])*conj(u[j][k])*psi0[j];
				
				}
			}	 
		}
	}
	else
	{
		// Antineutrinos
		for (int i=0; i<3; i++)
		{
			newpsi[i] = complex<FLOAT>(0,0);
			for (int j=0; j<3; j++)
			{
				for (int k=0; k<3; k++)
				{
					newpsi[i] += conj(u[i][k])*exp(complex<FLOAT>(0,-1)*kappa*mq[k])*u[j][k]*psi0[j];
				}
			}	 
		}
	}
	
	amplitude=complex<FLOAT>(0,0);
	for (int i=0; i<3; i++) amplitude += psi1[i]*newpsi[i];
	probability = real(amplitude*conj(amplitude));
	return probability;
}

extern "C" void glb_set_filter(double x)
{
  filter = x;
}




extern "C" FLOAT glb_filtered_vac_probability(int pl, int pm, int panti,double pen, double plength)
{	
	complex<FLOAT> psi0[3];
	complex<FLOAT> psi1[3];
	complex<FLOAT> newpsi[3];
	complex<FLOAT> amplitude;
	complex<FLOAT> aux[3][3];
	FLOAT probability=0;


	// Ausgangs und Endzustand (Flavour!!)
	for (int i=0; i<3; i++)
	{
	  for(int j=0;j<3;j++)
	    {
	      aux[i][j]=complex<FLOAT>(0,0);
	    }
		psi0[i] = complex<FLOAT>(0,0);
		psi1[i] = complex<FLOAT>(0,0);
	}
	psi0[pl-1]=complex<FLOAT>(1,0);
	psi1[pm-1]=complex<FLOAT>(1,0);


	// Propagation
	FLOAT kappa = 2*1.2669*plength/pen;
	
	if (panti==1)
	{
		// Neutrinos
		for (int i=0; i<3; i++)
		{
			newpsi[i] = complex<FLOAT>(0,0);
			for (int j=0; j<3; j++)
			{
				for (int k=0; k<3; k++)
				{
				  aux[i][k] += u[i][k]*exp(complex<FLOAT>(0,-1)*kappa*mq[k])*conj(u[j][k])*psi0[j];
				}
			}	 
		}
	}
	else
	{
		// Antineutrinos
		for (int i=0; i<3; i++)
		{
			newpsi[i] = complex<FLOAT>(0,0);
			for (int j=0; j<3; j++)
			{
				for (int k=0; k<3; k++)
				{
				 aux[i][k] += conj(u[i][k])*exp(complex<FLOAT>(0,-1)*kappa*mq[k])*u[j][k]*psi0[j];
				}
			}	 
		}
	}
	
	amplitude=complex<FLOAT>(0,0);

	   for(int i=0;i<3;i++)
	      {
		for(int i1=0;i1<3;i1++)
		  {
		    for(int k=0;k<3;k++)
		      {
			for(int k1=0;k1<3;k1++)
			  {
			    probability += real(psi1[i]*aux[i][k]*conj(aux[i1][k1]*psi1[i1]))*exp(-(mq[k]-mq[k1])*kappa*(mq[k]-mq[k1])*kappa/(pen*filter)/(pen*filter));
			  
			  }
		      }
		  }
	      }
	return probability;
}

extern "C" FLOAT glb_filtered_const_density_probability(int pl, int pm, int panti,double pen)
{	
	complex<FLOAT> psi0[3];
	complex<FLOAT> psi1[3];
	complex<FLOAT> newpsi[3];
	complex<FLOAT> amplitude;
	complex<FLOAT> aux[3][3];
	double plength;
	FLOAT probability=0;
	anti_neutrinos=panti;
	neutrino_energy = pen;
	matter_density = density[0];
	plength=length[0];
	int info = ComputeMixingMatrix();
	
	if (info!=0)
	{
		fprintf(stderr,"dsyev error\n");
		exit(1);
	}

	// Ausgangs und Endzustand (Flavour!!)
	for (int i=0; i<3; i++)
	{
	  for(int j=0;j<3;j++)
	    {
	      aux[i][j]=complex<FLOAT>(0,0);
	    }
		psi0[i] = complex<FLOAT>(0,0);
		psi1[i] = complex<FLOAT>(0,0);
	}
	psi0[pl-1]=complex<FLOAT>(1,0);
	psi1[pm-1]=complex<FLOAT>(1,0);


	// Propagation
	FLOAT kappa = 2*1.2669*plength/pen;
	
	if (panti==1)
	{
		// Neutrinos
		for (int i=0; i<3; i++)
		{
			newpsi[i] = complex<FLOAT>(0,0);
			for (int j=0; j<3; j++)
			{
				for (int k=0; k<3; k++)
				{
				  aux[i][k] += conj(vr[k][i])*exp(complex<FLOAT>(0,-1)*kappa*w[k])*vr[k][j]*psi0[j];
				}
			}	 
		}
	}
	else
	{
		// Antineutrinos
		for (int i=0; i<3; i++)
		{
			newpsi[i] = complex<FLOAT>(0,0);
			for (int j=0; j<3; j++)
			{
				for (int k=0; k<3; k++)
				{
				 aux[i][k] += vr[k][i]*exp(complex<FLOAT>(0,-1)*kappa*w[k])*conj(vr[k][j])*psi0[j];
				}
			}	 
		}
	}
	
	amplitude=complex<FLOAT>(0,0);

	   for(int i=0;i<3;i++)
	      {
		for(int i1=0;i1<3;i1++)
		  {
		    for(int k=0;k<3;k++)
		      {
			for(int k1=0;k1<3;k1++)
			  {
			    probability += real(psi1[i]*aux[i][k]*conj(aux[i1][k1]*psi1[i1]))*exp(-abs((w[k]-w[k1]))*kappa*abs((w[k]-w[k1]))*kappa/(filter*pen)/(filter*pen));
			  
			  }
		      }
		  }
	      }
	return probability;
}

static void FCDglb_probability_matrix(FLOAT prob[3][3], int panti, double pen)
{	
	complex<FLOAT> psi0[3];
	complex<FLOAT> psi1[3];
	complex<FLOAT> newpsi[3];

	complex<FLOAT> aux[3][3];
	double plength;
	FLOAT probability=0;
	anti_neutrinos=panti;
	neutrino_energy = pen;
	matter_density = density[0];
	plength=length[0];
	int info = ComputeMixingMatrix();
	
	if (info!=0)
	{
		fprintf(stderr,"dsyev error\n");
		exit(1);
	}
	for(int pl=0;pl<3;pl++)
	  {
	    for(int pm=0;pm<3;pm++)
	      {
		probability=0;
		// Ausgangs und Endzustand (Flavour!!)
		for (int i=0; i<3; i++)
		  {
		    for(int j=0;j<3;j++)
		      {
			aux[i][j]=complex<FLOAT>(0,0);
		      }
		    psi0[i] = complex<FLOAT>(0,0);
		    psi1[i] = complex<FLOAT>(0,0);
		  }
		psi0[pl]=complex<FLOAT>(1,0);
		psi1[pm]=complex<FLOAT>(1,0);
		
		
		// Propagation
		FLOAT kappa = 2*1.2669*plength/pen;
		
		if (panti==1)
		  {
		    // Neutrinos
		    for (int i=0; i<3; i++)
		      {
			newpsi[i] = complex<FLOAT>(0,0);
			for (int j=0; j<3; j++)
			  {
			    for (int k=0; k<3; k++)
			      {
				aux[i][k] += conj(vr[k][i])*exp(complex<FLOAT>(0,-1)*kappa*w[k])*vr[k][j]*psi0[j];
			      }
			  }	 
		      }
		  }
		else
		  {
		    // Antineutrinos
		    for (int i=0; i<3; i++)
		      {
			newpsi[i] = complex<FLOAT>(0,0);
			for (int j=0; j<3; j++)
			  {
			    for (int k=0; k<3; k++)
			      {
				aux[i][k] += vr[k][i]*exp(complex<FLOAT>(0,-1)*kappa*w[k])*conj(vr[k][j])*psi0[j];
			      }
			  }	 
		      }
		  }
		
		
		
		for(int i=0;i<3;i++)
		  {
		    for(int i1=0;i1<3;i1++)
		      {
			for(int k=0;k<3;k++)
			  {
			    for(int k1=0;k1<3;k1++)
			      {
				probability += real(psi1[i]*aux[i][k]*conj(aux[i1][k1]*psi1[i1]))*exp(-abs((w[k]-w[k1]))*kappa*abs((w[k]-w[k1]))*kappa/(pen*filter)/(pen*filter));
				
			      }
			  }
		      }
		  }
		//here
		prob[pl][pm]=probability;
	      }
	  }
       
}

static void FVACglb_probability_matrix(FLOAT prob[3][3], int panti, double pen)
{	
	complex<FLOAT> psi0[3];
	complex<FLOAT> psi1[3];
	complex<FLOAT> newpsi[3];

	complex<FLOAT> aux[3][3];
	double plength;
	FLOAT probability=0;
	anti_neutrinos=panti;
	neutrino_energy = pen;
	matter_density = density[0];
	plength=length[0];

	for(int pl=0;pl<3;pl++)
	  {
	    for(int pm=0;pm<3;pm++)
	      {
		probability=0;
		// Ausgangs und Endzustand (Flavour!!)
		for (int i=0; i<3; i++)
		  {
		    for(int j=0;j<3;j++)
		      {
			aux[i][j]=complex<FLOAT>(0,0);
		      }
		    psi0[i] = complex<FLOAT>(0,0);
		    psi1[i] = complex<FLOAT>(0,0);
		  }
		psi0[pl]=complex<FLOAT>(1,0);
		psi1[pm]=complex<FLOAT>(1,0);
		
		
		// Propagation
		FLOAT kappa = 2*1.2669*plength/pen;
		
		if (panti==1)
		  {
		    // Neutrinos
		    for (int i=0; i<3; i++)
		      {
			newpsi[i] = complex<FLOAT>(0,0);
			for (int j=0; j<3; j++)
			  {
			    for (int k=0; k<3; k++)
			      {
				aux[i][k] += u[i][k]*exp(complex<FLOAT>(0,-1)*kappa*mq[k])*conj(u[j][k])*psi0[j];
			      }
			  }	 
		      }
		  }
		else
		  {
		    // Antineutrinos
		    for (int i=0; i<3; i++)
		      {
			newpsi[i] = complex<FLOAT>(0,0);
			for (int j=0; j<3; j++)
			  {
			    for (int k=0; k<3; k++)
			      {
				aux[i][k] += conj(u[i][k])*exp(complex<FLOAT>(0,-1)*kappa*mq[k])*u[j][k]*psi0[j];
			      }
			  }	 
		      }
		  }
		
		
		
		for(int i=0;i<3;i++)
		  {
		    for(int i1=0;i1<3;i1++)
		      {
			for(int k=0;k<3;k++)
			  {
			    for(int k1=0;k1<3;k1++)
			      {
				probability += real(psi1[i]*aux[i][k]*conj(aux[i1][k1]*psi1[i1]))*exp(-abs((mq[k]-mq[k1]))*kappa*abs((mq[k]-mq[k1]))*kappa/(pen*filter)/(pen*filter));
				
			      }
			  }
		      }
		  }
		//here
		prob[pl][pm]=probability;
	      }
	  }
       
}



// -----------------------------------------------------------
// Set the fundamental mixing parameters and precompute 
// mixing matrix
// -----------------------------------------------------------
// dm21    mass square difference dm^2_21
// dm31    mass square difference dm^2_31
// pt12    mixing angle theta_12
// pt12    mixing angle theta_12
// pt13    mixing angle theta_13
// pt23	   mixing angle theta_23
// pdelta1 CP phase
// -----------------------------------------------------------
extern "C" void glb_set_parameters(double dm21, double dm31, double pt12, double pt13, double pt23, double pdelta1)
{
	mq[0] = abs(dm31);
	mq[1] = abs(dm31)+dm21;
	mq[2] = abs(dm31)+dm31;

	t12 = pt12;
	t13 = pt13;
	t23 = pt23;
	d1 = pdelta1;

	

	u[0][0] = complex<FLOAT>(cos(t12)*cos(t13),0);
	u[0][1] = complex<FLOAT>(sin(t12)*cos(t13),0);
	u[0][2] = complex<FLOAT>(sin(t13)*cos(d1),-sin(t13)*sin(d1));

	u[1][0] = complex<FLOAT>(-sin(t12)*cos(t23)-cos(t12)*sin(t23)*sin(t13)*cos(d1),-cos(t12)*sin(t23)*sin(t13)*sin(d1));
	u[1][1] = complex<FLOAT>(cos(t12)*cos(t23)-sin(t12)*sin(t23)*sin(t13)*cos(d1),-sin(t12)*sin(t23)*sin(t13)*sin(d1));
	u[1][2] = complex<FLOAT>(sin(t23)*cos(t13),0);

	u[2][0] = complex<FLOAT>(sin(t12)*sin(t23)-cos(t12)*cos(t23)*sin(t13)*cos(d1),-cos(t12)*cos(t23)*sin(t13)*sin(d1));
	u[2][1] = complex<FLOAT>(-cos(t12)*sin(t23)-sin(t12)*cos(t23)*sin(t13)*cos(d1),-sin(t12)*cos(t23)*sin(t13)*sin(d1));
	u[2][2] = complex<FLOAT>(cos(t23)*cos(t13),0);
	
}




// -----------------------------------------------------------
// Oscillation probability through predefined matter profile
// -----------------------------------------------------------
// pl        initial flavour
// pm        final flavour
// panti     neutrino or antineutrino (1,-1)
// pen       neutrino energy in GeV
// -----------------------------------------------------------
extern "C" FLOAT glbProfileProbability(int pl, int pm, int panti,
double pen)
{	

	complex<FLOAT> psi0[3];
	complex<FLOAT> psi1[3];
	complex<FLOAT> amplitude;
	FLOAT probability;
	
	anti_neutrinos = panti;

	// Ausgangs und Endzustand (Flavour!!)
	for (int i=0; i<3; i++)
	{
		psi0[i] = complex<FLOAT>(0,0);
		psi1[i] = complex<FLOAT>(0,0);
	}
	psi0[pl-1]=complex<FLOAT>(1,0);
	psi1[pm-1]=complex<FLOAT>(1,0);
			
	for (int i=0; i<psteps; i++)	
	   PropagateInMatter(psi0,panti,length[i],pen,density[i]);
	
	amplitude=complex<FLOAT>(0,0);
	for (int i=0; i<3; i++) amplitude += psi1[i]*psi0[i];
	probability = real(amplitude*conj(amplitude));
	
	return probability;
}



// -----------------------------------------------------------
// Probability matrix through predefined matter profile
// -----------------------------------------------------------
// panti     neutrino or antineutrino (1,-1)
// pen       neutrino energy in GeV
// -----------------------------------------------------------
// -----------------------------------------------------------
// Calculation of the amplitude matrix through a 
// layer of matter
// -----------------------------------------------------------
// panti       neutrino or antineutrino (1,-1)
// length      layer length in km
// rho         matter density in g/cm3
// en          neutrino energy in GeV
// -----------------------------------------------------------
static void MultiplyAmplitudeMatrix(complex<FLOAT> amp[3][3], int panti, FLOAT length, FLOAT en, FLOAT rho)
{
	complex<FLOAT> newamp[3][3];

	FLOAT kappa = 2*1.2669*length/en;

	neutrino_energy = en;
	matter_density = rho;
	anti_neutrinos=panti; // this is a bug fix
	int info = ComputeMixingMatrix();
	
	if (info!=0)
	{
		fprintf(stderr,"dsyev error\n");
		exit(1);
	}
	
	
	for (int i=0; i<3; i++)
	{
		for (int j=0; j<3; j++)
		{
			newamp[i][j] = complex<FLOAT>(0,0);
			for (int k=0; k<3; k++)
			{
				for (int l=0; l<3; l++)
				{
					if (panti==1)
						newamp[i][j] += conj(vr[k][i])*exp(complex<FLOAT>(0,-1)*kappa*w[k])*vr[k][l]*amp[l][j];
					else
						newamp[i][j] += vr[k][i]*exp(complex<FLOAT>(0,-1)*kappa*w[k])*conj(vr[k][l])*amp[l][j];
				}
			}
		}
	}
	for (int i=0; i<3; i++)
	{
		for (int j=0; j<3; j++)
		{
			amp[i][j] = newamp[i][j];
		}
	}
}

extern "C" void glb_switch_filter(int on)
{
  filter_state=on;
}

extern "C" void glb_probability_matrix(FLOAT prob[3][3], int panti, double pen)
{	
	complex<FLOAT> amp[3][3];
	if(filter_state==1&&psteps==1)
	  {
	    if(density[0]<0.001)
	      {
		FVACglb_probability_matrix(prob,panti,pen);
	      }
	    else
	      {	      
		FCDglb_probability_matrix(prob,panti,pen);
	      }
	  }
	else
	  {
	    // Unity matrix as starting point
	    for (int i=0; i<3; i++)
	      for (int j=0; i<3; i++)
		amp[i][j] = complex<FLOAT>(0,0);
	    for (int i=0; i<3; i++)
	      amp[i][i]=complex<FLOAT>(1,0);
	    
	    for (int i=0; i<psteps; i++)
	      MultiplyAmplitudeMatrix(amp,panti,length[i],pen,density[i]);
	    
	    // amplitudes -> probabilities
	    for (int i=0; i<3; i++)
	      for (int j=0; j<3; j++)
		// here is the origin of the CP-sign bug
		// it is fixed by replacing prob[i][j] with prob[j][i]
		prob[j][i] =  real(amp[i][j]*conj(amp[i][j]));
	  }
}



extern "C" void glb_set_c_vacuum_parameters(double pt12, double pt13, double pt23, double pdelta1)
{
  solar_mixing = pt12; // this allows to set th12 to any value
	t12 = solar_mixing;
	t13 = pt13;
	t23 = pt23;
	d1 = pdelta1;

	

	u[0][0] = complex<FLOAT>(cos(t12)*cos(t13),0);
	u[0][1] = complex<FLOAT>(sin(t12)*cos(t13),0);
	u[0][2] = complex<FLOAT>(sin(t13)*cos(d1),-sin(t13)*sin(d1));

	u[1][0] = complex<FLOAT>(-sin(t12)*cos(t23)-cos(t12)*sin(t23)*sin(t13)*cos(d1),-cos(t12)*sin(t23)*sin(t13)*sin(d1));
	u[1][1] = complex<FLOAT>(cos(t12)*cos(t23)-sin(t12)*sin(t23)*sin(t13)*cos(d1),-sin(t12)*sin(t23)*sin(t13)*sin(d1));
	u[1][2] = complex<FLOAT>(sin(t23)*cos(t13),0);

	u[2][0] = complex<FLOAT>(sin(t12)*sin(t23)-cos(t12)*cos(t23)*sin(t13)*cos(d1),-cos(t12)*cos(t23)*sin(t13)*sin(d1));
	u[2][1] = complex<FLOAT>(-cos(t12)*sin(t23)-sin(t12)*cos(t23)*sin(t13)*cos(d1),-sin(t12)*cos(t23)*sin(t13)*sin(d1));
	u[2][2] = complex<FLOAT>(cos(t23)*cos(t13),0);
	
}


extern "C" double* glb_get_vacuum_parameters()
{
  double* erg;
  erg=(double*) glb_malloc(4*sizeof(double));
  erg[0]=t12;
  erg[1]=t13;
  erg[2]=t23;
  erg[3]=d1;
  return erg;
}

extern "C" void glb_set_c_squared_masses(double mq1, double mq2, double mq3)
{
	mq[0] = mq1;
	mq[1] = mq2;
	mq[2] = mq3;
	
} 

extern "C" double* glb_get_squared_masses()
{
  double* erg;
  erg=(double*) glb_malloc(2*sizeof(double));
  erg[0]=mq[1];
  erg[1]=mq[2]-mq[1];
  return erg;
	
} 

/* made length, density and psteps local, thus I have to provide set and
 * get functions */

extern "C" double *glb_get_density_ptr()
{
  return &density[0];
}

extern "C" double *glb_get_length_ptr()
{
  return &length[0];
}

extern "C" int glb_get_psteps()
{
  return psteps;
}

extern "C" void glb_set_density_ptr(double *d)
{
  density=d;
}

extern "C" void glb_set_length_ptr(double *l)
{
  length=l;
}

extern "C" void glb_set_psteps(int p)
{
  psteps=p;
}




