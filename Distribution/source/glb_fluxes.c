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



/* Units and Normalization */

/* Energies GeV */
/* Length km */
/* dm^2 eV^2 */
/* densitie g/cm^3 */
/* X-section 10^-38 cm^2 */
/* Targetmass kiloton (10^6 kg) */
/* Fluxes are normalized to a baseline of 1km */
/* stored_muons is for NuFact only */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "glb_types.h"
#include "glb_error.h"
#include "glb_multiex.h"

#define PI 3.1415
#define NORM 6.02204*1E-12 //10^-38/u (atomic mass unit)
 
static double mmu = 0.1056;           // muon mass in GeV
static double beam_pol = 0;           // beam polaristaion
static double beam_en = 20;           // beam energy in GeV
static double beta,rgamma;             // relativistic quantities
static double stored_muons = 2.0E20;  // stored muons per year 
// or protons on target
// care must be taken with this, it is a different number for each experiment
static double power;
static int    channel = 1;            // channel
static int NFLUX = 500;
static double xsec[32][1001][7];

void glb_flux_loader(glb_flux *data, int number, int polarity);


//this loads the fluxes for all non-NuFact beams
// the files must be of the form
// Energy[GeV] Flux[?]
// You should check the normalisation !!!
// the file has 501 lines


//this part is by Martin and for the Nufact only

void glb_set_max_energy(double en)
{
  beam_en = en;
  rgamma = en/mmu;
  beta = -sqrt(1-1/rgamma/rgamma);
}

double glb_get_max_energy()
{
  return beam_en;
}

static void set_beam_polarisation(double pol)
{
  beam_pol = pol;
}

static double check_beam_polarisation()
{
  return beam_pol;
}

static void set_stored_muons(double n)
{
  stored_muons = n;
}

static double check_stored_muons()
{
  return stored_muons;
}


static void set_power(double n)
{
  power = n;
}

static double check_power()
{
  return power;
}


static void set_channel(double ch)
{
  channel = ch;
}

static double flux_mu(double en, double baseline, double cost_lab)
{
  double x_lab, cost, x, xmax;
  baseline *=1000;

  x_lab = 2*en/mmu;
  cost = (beta+cost_lab)/(1+beta*cost_lab);
  x = x_lab/(rgamma*(1-beta*cost_lab));
  xmax = rgamma*(1-beta*cost_lab);

  if (x_lab < xmax)
    return rgamma*(1 - beta*cost) /baseline/baseline/mmu * stored_muons * x*x/PI * ((3-2*x)-channel*(1-2*x)*beam_pol*cost);
  else
    return 0.0;
}

static double flux_e(double en, double baseline, double cost_lab)
{
  double x_lab, cost, x, xmax;
  baseline *=1000;

  x_lab = 2*en/mmu;
  cost = (beta+cost_lab)/(1+beta*cost_lab);
  x = x_lab/(rgamma*(1-beta*cost_lab));
  xmax = rgamma*(1-beta*cost_lab);

  if (x_lab< xmax)
    return rgamma*(1 - beta*cost) /baseline/baseline/mmu * stored_muons* 6*x*x/PI * ((1-x)-channel*(1-x)*beam_pol*cost);
  else
    return 0.0;
}

static double fastflux_mu(double en, double baseline)
{
  double x_lab, x, xmax;
  baseline *=1000;

  x_lab = 2*en/mmu;
  x = x_lab/(rgamma*(1-beta));
  xmax = rgamma*(1-beta);

  if (x_lab < xmax)
    return rgamma*(1 - beta) /baseline/baseline/mmu * x*x/PI * (3-2*x);
  else
    return 0.0;
}

static double fastflux_e(double en, double baseline)
{
  double x_lab, x, xmax;
  baseline *=1000;

  x_lab = 2*en/mmu;
  x = x_lab/(rgamma*(1-beta));
  xmax = rgamma*(1-beta);

  if (x_lab < xmax)
    return rgamma*(1 - beta) /baseline/baseline/mmu *  6*x*x/PI * (1-x);
  else
    return 0.0;
}
//here ends the NuFact part

static double nufact_flux(double en, double baseline, int polarity, int l, int anti)
{
  double ergebnis;
  ergebnis=0.0;
  if (l==1 && anti == 1 && polarity == 1)
    {
      ergebnis=fastflux_e(en,baseline);
    }
  if (l==2 && anti == -1 && polarity == 1)
    {
      ergebnis=fastflux_mu(en,baseline);
    }
  if (l==1 && anti == -1 && polarity == -1)
    {
      ergebnis=fastflux_e(en,baseline);
    }
  if (l==2 && anti == 1 && polarity == -1)
    {
      ergebnis=fastflux_mu(en,baseline);
    }


  return ergebnis;
}



static void nufact_flux_setup(glb_flux *data, int id,int pl)
{
  /* the id argument is unused */
  double smb,pb,de;
  int i;
  de=data->parent_energy/501.0;
  smb=stored_muons;
  pb=power;
  stored_muons=data->stored_muons;
  power=data->target_power;
  glb_set_max_energy(data->parent_energy);
  data->flux_storage=glb_alloc_flux_storage(501);
  for(i=0;i<501;i++)
    {
      data->flux_storage[i][0]=i*de;
      data->flux_storage[i][1]=nufact_flux(i*de,1.0,pl,1,+1); 
      data->flux_storage[i][2]=nufact_flux(i*de,1.0,pl,2,+1); 
      data->flux_storage[i][3]=nufact_flux(i*de,1.0,pl,3,+1); 
      data->flux_storage[i][4]=nufact_flux(i*de,1.0,pl,1,-1); 
      data->flux_storage[i][5]=nufact_flux(i*de,1.0,pl,2,-1); 
      data->flux_storage[i][6]=nufact_flux(i*de,1.0,pl,3,-1); 
    }
  power=pb;
  stored_muons=smb;
}

void glb_init_fluxtables(glb_flux *data,int pos)
{
  if(data->builtin==1) nufact_flux_setup(data,pos,+1);
  if(data->builtin==2) nufact_flux_setup(data,pos,-1);
  if(data->builtin==0) glb_flux_loader(data,pos,+1);
  
}



// Here we settle the normalization issue for all beam types


static double norm2(int type)
{
  double erg;
  switch (type)
    {
    case 1: erg=NORM*100; //NuFact
      break;
    case 2: erg=NORM*100; //NuFact
      break;   
    case 3: erg=NORM*295*295*9.92033*1E6;//NUMI Beam (9@712) 
      break;
    case 4: erg=NORM*295*295*9.92033*1E6; //NUMI Beam (9@712) 
      break;
    case 5: erg=NORM*295*295*9.92033*1E6; //NUMI Beam (12@712) 
      break;
    default: erg=NORM*295*295*9.92033*1E6;
      break;
    }
  return erg;
}


 
double glb_flux_calc(double en, double baseline, 
	    int type, int l, int anti, const glb_flux *data)
{
   /* the type argument is unused */
  double incre;
  double ergebnis;
  int lowerind;
  int higherind;
  int part=1;
  if (l==1 && anti == 1)
    {
      part = 1;
    }
 if (l==2 && anti == 1)
    {
      part = 2;
    }
 if (l==3 && anti == 1)
    {
      part = 3;
    }
 if (l==1 && anti == -1)
    {
      part = 4;
    }
 if (l==2 && anti == -1)
    {
      part = 5;
    }
 if (l==3 && anti == -1)
    {
      part = 6;
    }
 
 
 incre = (data->flux_storage[500][0]-data->flux_storage[0][0])/501.0;
 lowerind = floor((en-data->flux_storage[0][0])/incre);
 higherind = lowerind + 1;

 if (lowerind<0||higherind>500) 
   {
     ergebnis=0.0;
   }
 else
   {
     ergebnis=((data->flux_storage[higherind][part]-
		data->flux_storage[lowerind][part])*(((en-data->flux_storage[0][0])
						  /incre)-lowerind)
	       +data->flux_storage[lowerind][part])
       /baseline/baseline;
   }     
 ergebnis = ergebnis*norm2(data->builtin) 
   * data->time * data->target_power * data->stored_muons * data->norm;
 return ergebnis;
}

 

/* Here begins the part where the cross-sections are read from 
 * file and made accessible to the program */


/***********************************************************
 * Flux-Loader                                             *
 ***********************************************************/

void glb_flux_loader(glb_flux *data, int ident, int polarity)
{
   /* the polarity  argument is unused */
  int i,number;

 
  FILE *fp = glb_fopen(data->file_name,"r");
  number=ident;
 
  
  data->flux_storage=glb_alloc_flux_storage(501);
      if (fp!=NULL)
	{
	  for (i=0; i<=NFLUX; i++) fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf",&data->flux_storage[i][0],&data->flux_storage[i][1],&data->flux_storage[i][2], &data->flux_storage[i][3],&data->flux_storage[i][4],&data->flux_storage[i][5],&data->flux_storage[i][6]);
	  fclose(fp);
	}
      else
	{
	  fprintf(stderr,"Error: Could not open %s\n",data->file_name); 
	} 
}

void glb_X_section_loader(char filename[], int ident)
{
  int i,number;
  FILE* fp;
  number=ident;
  if (number>31) number=31;
  fp = glb_fopen(filename, "r");
  if (fp!=NULL)
  {
    for (i=0; i<=1000; i++) fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf",&xsec[number][i][0],&xsec[number][i][1],&xsec[number][i][2], &xsec[number][i][3],&xsec[number][i][4],&xsec[number][i][5],&xsec[number][i][6]);
  }
  else
  {
    fprintf(stderr,"Error: Could not open %s\n",filename);
  }

  fclose(fp);
}

double glbXSection(int ident, double enl, int l, int anti)
{
  double incre;
  double ergebnis;
  double en;
  int lowerind;
  int higherind;
  int part;

  en = log10(enl);
  part = (int) l + 3*(1-anti)/2;

  incre = (xsec[ident][1000][0]-xsec[ident][0][0])/1000.0;
  lowerind = floor((en-xsec[ident][0][0])/incre); 
  higherind = lowerind + 1;

  if (lowerind<0 || higherind>1000) ergebnis=0.0; 
  else
  {
    ergebnis=((xsec[ident][higherind][part]-xsec[ident][lowerind][part])*(((en-xsec[ident][0][0])/incre)-lowerind)+xsec[ident][lowerind][part]);
  }
  return ergebnis*enl;
}

