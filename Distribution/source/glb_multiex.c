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
#include <string.h>
#include <math.h>
#include <globes/globes.h>

#include "glb_probability.h"
#include "glb_fluxes.h"
#include "glb_rate_engine.h"
#include "nrutil.h"
#include "glb_error.h"

#include "glb_minimize.h"
#include "glb_smear.h"
#include "glb_types.h"
#include "glb_multiex.h"
#include "glb_version.h"

#define FLOAT double

#define PI 3.1415


static struct experiment MInitMemory0(struct experiment in);

/* here come the glb_flux functions */

/* Dynamically allocating flux storage */
double** glb_alloc_flux_storage(size_t lines)
{
  double **temp;
  size_t k;
  temp = (double **) glb_malloc(sizeof(double *)*(lines+1));
  for(k=0;k<lines;k++)
    {
      temp[k]=(double *) glb_malloc(sizeof(double)*7);
    }
  temp[lines]=NULL;
  return temp;  
}
 
void glb_free_flux_storage(double **stale)
{
  size_t i;
  if(stale!=NULL){
  for(i=0;stale[i]!=NULL;i++) glb_free(stale[i]);
  glb_free(stale);
  stale=NULL;
  }
}

glb_flux *glb_flux_alloc()
{
  glb_flux *temp;
  temp=(glb_flux *) glb_malloc(sizeof(glb_flux));
  temp->builtin=-1;
  temp->file_name=NULL;
  temp->time=-1;
  temp->target_power=-1;
  temp->stored_muons=-1;
  temp->parent_energy=-1;
  temp->norm=-1;
  temp->flux_storage=NULL;
  return temp;
}

glb_flux *glb_flux_reset(glb_flux *temp)
{
  temp->builtin=-1;
  /* FIXME memory leak */
  temp->file_name=NULL;
  temp->time=-1;
  temp->target_power=-1;
  temp->stored_muons=-1;
  temp->parent_energy=-1;
  temp->norm=-1; 
  /* FIXME serious memory leak */
  glb_free_flux_storage(temp->flux_storage);
  return temp;
}



void glb_flux_free(glb_flux *stale)
{
  if(stale!=NULL)
    {
      if(stale->file_name!=NULL) glb_free(stale->file_name);
      if(stale->flux_storage!=NULL) glb_free_flux_storage(stale->flux_storage);
      glb_free(stale);
    }
}

int glb_default_flux(glb_flux *in)
{
  int s=0;
  if(in->builtin==-1) in->builtin=0;
  if(in->time==-1) in->time=1;
  if(in->target_power==-1) in->target_power=1;
  if(in->stored_muons==-1) in->stored_muons=1;
  if(in->norm==-1) in->norm=1;  
  if(in->builtin!=0)
    {
      if(in->parent_energy==-1) {glb_error("No parent energy defined!");s=1;}
    }
  if(in->builtin==0)
    {
      if(in->file_name==NULL) {glb_error("No flux file specified!");s=1;}
    }
  if(s!=0) glb_error("glb_flux not properly setup");
  return s;
}

glb_flux  *cpy_glb_flux(glb_flux *dest, const glb_flux *src)
{
  dest=(glb_flux *) memcpy(dest,src,sizeof(glb_flux)); 
  if(src->file_name!=NULL)
    {
      dest->file_name=(char *) strdup(src->file_name);
      if(dest->file_name==NULL) glb_fatal("Erros in strdup");
    }
  /* FIXME flux_storage is neither copied nor freed */ 
  return dest;	  
}


/* here come the glb_xsec functions */

/* Dynamically allocating flux storage */
double** glb_alloc_xsec_storage(size_t lines)
{
  double **temp;
  size_t k;
  temp = (double **) glb_malloc(sizeof(double *)*(lines+1));
  for(k=0;k<lines;k++)
    {
      temp[k]=(double *) glb_malloc(sizeof(double)*7);
    }
  temp[lines]=NULL;
  return temp;  
}
 
void glb_free_xsec_storage(double **stale)
{
  size_t i;
  if(stale!=NULL){
  for(i=0;stale[i]!=NULL;i++) glb_free(stale[i]);
  glb_free(stale);
  stale=NULL;
  }
}

glb_xsec *glb_xsec_alloc()
{
  glb_xsec *temp;
  temp=(glb_xsec *) glb_malloc(sizeof(glb_xsec));
  temp->builtin=-1;
  temp->file_name=NULL;
  temp->xsec_storage=NULL;
  return temp;
}

glb_xsec *glb_xsec_reset(glb_xsec *temp)
{
  temp->builtin=-1;
  /* FIXME memory leak */
  temp->file_name=NULL;
  glb_free_xsec_storage(temp->xsec_storage);
  return temp;
}



void glb_xsec_free(glb_xsec *stale)
{
  if(stale!=NULL)
    {
      if(stale->file_name!=NULL) glb_free(stale->file_name);
      if(stale->xsec_storage!=NULL) glb_free_xsec_storage(stale->xsec_storage);
      glb_free(stale);
    }
}

int glb_default_xsec(glb_xsec *in)
{
  int s=0;
  if(in->builtin==-1) in->builtin=0;
  /* if(in->xsec_storage==NULL) {glb_error("No storage for X-section allocated");
     s=-1;}*/
  if(s!=0) glb_error("glb_xsec not properly setup");
  return s;
}

glb_xsec  *cpy_glb_xsec(glb_xsec *dest, const glb_xsec *src)
{
  size_t i;
  dest=(glb_xsec *) memcpy(dest,src,sizeof(glb_xsec)); 
  if(src->file_name!=NULL)
    {
      dest->file_name=(char *) strdup(src->file_name);
      if(dest->file_name==NULL) glb_fatal("Error in strdup");
    }
  if(src->xsec_storage!=NULL)
    {
      dest->xsec_storage=glb_alloc_xsec_storage(1001);
      for(i=0;i<1001;i++) dest->xsec_storage[i]=src->xsec_storage[i];
    } 
  return dest;	  
}


//this structure contains all information needed for implementing
//the systematics for each rule and experiment. Should be roughly
//equivalent to a class in C++

/* this function definition clashes nearyl with 
 * struct systematic glb_init_systematic
 * in Mminimize.c
 */


static void init_struct_systematic(struct glb_systematic *in)
{
  in->chi_func=NULL; 
  in->dimension=-1;
  in->sp=NULL; 
  in->errors=NULL; 
  in->evalf=NULL; 
  in->info=NULL;
};



void glbInitExp(glb_exp ins)
{
  int i;
  struct experiment *in;
  in=(struct experiment *) ins;
  in->version=NULL;
  in->num_of_fluxes=-1;
  in->density_profile_type=-1;
  /* FIXME - a potential memory leak */
  for(i=0;i<32;i++) in->fluxes[i]=NULL;
  in->num_of_xsecs=-1;
  /* FIXME - a potential memory leak */
  for(i=0;i<32;i++) in->xsecs[i]=NULL;

  in->binsize=NULL;
  in->simbinsize=NULL;
  for(i=0;i<32;i++) in->errordim[i]=-1;
  for(i=0;i<32;i++) in->errordim_sys_off[i]=-1;
  for(i=0;i<32;i++) in->errordim_sys_on[i]=-1;
  in->errorfunction=-1;
  in->baseline=-1;
  in->emin=-1;
  in->emax=-1;
  in->numofbins=-1;
  in->targetmass=-1;

  in->numofchannels=-1;
  for(i=0;i<6;i++) in->listofchannels[i]=NULL;
  in->numofrules=-1;
  in->num_of_sm=-1;
  for(i=0;i<32;i++) in->lengthofrules[i]=-1;
  for(i=0;i<32;i++) in->rulescoeff[i]=NULL;
  for(i=0;i<32;i++) in->rulechannellist[i]=NULL;
  for(i=0;i<32;i++) in->signalruleerror[0][i]=-1; 
  for(i=0;i<32;i++) in->signalruleerror[1][i]=-1;
  for(i=0;i<32;i++) in->lengthofbgrules[i]=-1;
  for(i=0;i<32;i++) in->bgrulescoeff[i]=NULL;
  for(i=0;i<32;i++) in->bgrulechannellist[i]=NULL;
  for(i=0;i<32;i++) in->bgcenter[0][i]=-1;
  for(i=0;i<32;i++) in->bgcenter[1][i]=-1;
  for(i=0;i<32;i++) in->bgerror[0][i]=-1;
  for(i=0;i<32;i++) in->bgerror[1][i]=-1;
  for(i=0;i<32;i++) in->bgtcenter[0][i]=-1; 
  for(i=0;i<32;i++) in->bgtcenter[1][i]=-1;
  for(i=0;i<32;i++) in->bgterror[0][i]=-1;
  for(i=0;i<32;i++) in->bgterror[1][i]=-1;
  for(i=0;i<32;i++) {in->smear_data[i]=NULL;}
  for(i=0;i<32;i++) in->smear[i]=NULL;
  for(i=0;i<32;i++) in->lowrange[i]=NULL;
  for(i=0;i<32;i++) in->uprange[i]=NULL; 
 for(i=0;i<32;i++) in->chrb[i]=NULL; 
 for(i=0;i<32;i++) in->chra[i]=NULL; 
 in->buffer=NULL;
  in->simtresh=-1;
  in->simbeam=-1;
  in->simbins=-1; 

  in->filter_state=-1;
  in->filter_value=-1;
  in->density_center=-1;
  in->density_error=-1;
  in->psteps=-1;
  in->lengthtab=NULL;
  in->densitytab=NULL;
  in->densitybuffer=NULL;
  for(i=0;i<32;i++) in->user_pre_smearing_channel[i]=NULL;
  for(i=0;i<32;i++) in->user_post_smearing_channel[i]=NULL;
  for(i=0;i<32;i++) in->user_pre_smearing_background[i]=NULL;
  for(i=0;i<32;i++) in->user_post_smearing_background[i]=NULL;
  for(i=0;i<32;i++) in->energy_window[i][0]=-1;
  for(i=0;i<32;i++) in->energy_window[i][1]=-1;
  in->chirate=NULL;
  for(i=0;i<32;i++)
    {
      in->SignalRates[i]=NULL;
      in->BackgroundRates[i]=NULL;
      in->rates0[i]=NULL;
      in->rates1[i]=NULL;
      in->rates1T[i]=NULL;
      in->rates1BG[i]=NULL;
      in->rates1BGT[i]=NULL;
      in->ratevec[i]=NULL;
      in->ratevecBG[i]=NULL;
    }
  in->energy_tab=NULL;
  for(i=0;i<32;i++) in->no_osc_background[i]=NULL;
  for(i=0;i<32;i++) init_struct_systematic(&(in->sys[i])); 
  

}


static void my_free(void *ptr)
{
  if(ptr==NULL) return;
  glb_free(ptr);
  ptr=NULL;
}

void glbFreeExp(glb_exp ins)
{
  int i,j;
  struct experiment *in;
  in=(struct experiment *) ins;
  
  glb_free(in->version);
  for(i=0;i<32;i++)glb_flux_free(in->fluxes[i]);
  for(i=0;i<32;i++)glb_xsec_free(in->xsecs[i]);

  for(i=0;i<6;i++) my_free(in->listofchannels[i]);
  my_free(in->binsize);
  my_free(in->simbinsize);

 
  for(i=0;i<in->numofrules;i++)
    {
      my_free(in->rulescoeff[i]);
      my_free(in->rulechannellist[i]);
      my_free(in->bgrulescoeff[i]);
      my_free(in->bgrulechannellist[i]);
    }
 

  
  for(i=0;i<32;i++) glb_smear_free(in->smear_data[i]);
  for(i=0;i<in->num_of_sm;i++)
    {
      my_free(in->lowrange[i]);
      my_free(in->uprange[i]);
      for(j=0;j<in->numofbins;j++) my_free(in->smear[i][j]);
      my_free(in->smear[i]);
    }
  for(i=0;i<in->numofchannels;i++) 
    {
      my_free(in->chrb[i]);
      my_free(in->chra[i]);
      my_free(in->user_pre_smearing_channel[i]);
      my_free(in->user_post_smearing_channel[i]);
      my_free(in->user_pre_smearing_background[i]);
      my_free(in->user_post_smearing_background[i]);
    } 
  my_free(in->buffer);
  my_free(in->lengthtab);
  my_free(in->densitytab);
  my_free(in->densitybuffer);
  my_free(in->chirate);
  
  for(i=0;i<in->numofrules;i++)
    {
      my_free(in->SignalRates[i]);
      my_free(in->BackgroundRates[i]);
      my_free(in->rates0[i]);
      my_free(in->rates1[i]);
      my_free(in->rates1T[i]);
      my_free(in->rates1BG[i]);
      my_free(in->rates1BGT[i]);
      my_free(in->ratevec[i]);
      my_free(in->ratevecBG[i]);
    }
  my_free(in->energy_tab);
  my_free(in);

}



	



glb_exp glbAllocExp()
{
  struct experiment *in; 
  in=(struct experiment *) glb_malloc(sizeof(struct experiment));
  glbInitExp((glb_exp) in);
  return (glb_exp) in;
}


static int setup_density_profile(glb_exp ins)
{
  int s=0;
  double *lb,*db;
  size_t l;
  double ttl=0;
  int i;
  struct experiment *in;
  in=(struct experiment *) ins;
  if(ins->density_profile_type==-1) 
    {glb_error("No profile type specified");s=-1;}
  if(ins->density_profile_type==1)
    {
      glb_free(ins->lengthtab);
      glb_free(ins->densitytab);
      ins->psteps=1;
      if(ins->baseline > 0)
	{
	  s+=glbAverageDensityProfile(ins->baseline,&lb,&db);
	  ins->densitytab=db;
	  ins->lengthtab=lb;
	}
      else
	{
	  glb_error("Baseline must be a positive number");
	  s=-1;
	}
    }

  if(ins->density_profile_type==2)
    {
      glb_free(ins->lengthtab);
      glb_free(ins->densitytab);
     
      if(ins->baseline > 0)
	{
	  if(ins->psteps > 0)
	    {
	      l=(size_t) ins->psteps;
	      s+=glbStaceyProfile(ins->baseline,l,&lb,&db);
	      ins->densitytab=db;
	      ins->lengthtab=lb;
	    }
	  else
	    {
	      glb_error("Densitysteps must be a positive number");
	      s=-1;
	    }
	      
	}
      else
	{
	  glb_error("Baseline must be a positive number");
	  s=-1;
	}
    }
  
  if(ins->density_profile_type==3)
    {
      if(ins->lengthtab!=NULL && ins->psteps > 0)
      {
	for(i=0;i<ins->psteps;i++) ttl+=ins->lengthtab[i];
      }
      if(ttl!=ins->baseline&&ins->baseline!=-1) 
	glb_warning("Override baseline with sum of lengthtab");
      ins->baseline=ttl;
    }
  
  
  return s;
}

int glbDefaultExp(glb_exp ins)
{
  int i,k,status,def,ct;
  
  
  struct experiment *in;
  in=(struct experiment *) ins;
 
  
  
  status=0;
  def=0;
  status+=setup_density_profile(ins);
  if(in->version==NULL) 
    {
      glb_warning("Missing version in AEDL file");
      def=-1;
      in->version=(char *) strdup(glb_release_version);
    }
 
  if(glbTestReleaseVersion(in->version)<0) {
    glb_warning("AEDL file has a more recent version number than the"
		" installed globes package");}
   

  if(in->num_of_xsecs<1)  {glb_error("No X-section selected!");status=-1;}
  if(in->num_of_xsecs>31)  {glb_error("To many X-sections!");status=-1;}
  if(in->num_of_fluxes<1)  {glb_error("No flux selected!");status=-1;}
  if(in->num_of_fluxes>31)  {glb_error("To many fluxes!");status=-1;}
  if(in->num_of_fluxes>0&&in->num_of_fluxes<32)
    {
      for(i=0;i<in->num_of_fluxes;i++)
	{
	  if(in->fluxes[i]==NULL) {glb_error("Flux specs missing");status=-1;}
	  else
	    {
	      if(glb_default_flux(in->fluxes[i])==0)
		glb_init_fluxtables(in->fluxes[i],i);
	      else
		status=-1;
	    }
	}
    }

  if(in->num_of_xsecs>0&&in->num_of_xsecs<32)
    {
      for(i=0;i<in->num_of_xsecs;i++)
	{
	  if(in->xsecs[i]==NULL) {glb_error("X-section specs missing")
				    ;status=-1;}
	  else
	    {
	      if(glb_default_xsec(in->xsecs[i])==0)
		glb_init_xsectables(in->xsecs[i]);
	      else
		status=-1;
	    }
	}
    }
 
  for(i=0;i<in->numofrules;i++) 
    {
      


      if(in->errordim_sys_off[i]==-1){
	in->errordim_sys_off[i]=2;
	// This is a safe setting, i.e. no systematics are
	// considered
	def=-1;}

      if(in->errordim_sys_on[i]==-1){
	in->errordim_sys_on[i]=2;
	// This is a safe setting, i.e. no systematics are
	// considered
	def=-1;}

      in->errordim[i]=in->errordim_sys_on[i];
     
    }
  if(in->errorfunction==-1){in->errorfunction=0;def=-1;}
  if(in->baseline==-1){glb_error("No baseline specified!");status=-1;}
  if(in->emin==-1){glb_error("No emin specified!");status=-1;}
  if(in->emax==-1){glb_error("No emax specified!");status=-1;}
  if(in->numofbins==-1){glb_error("numofbins is not set!");status=-1;}
  if(in->numofbins<=0) { glb_error("Too few bins defined!");status=-1;}

  /* It's okay if they are NULL or anything else ;-) */
  if(in->binsize==NULL) {in->binsize=NULL;}
  if(in->simbinsize==NULL) {in->simbinsize=NULL;}
  
  if(in->targetmass==-1){glb_error("No targetmass specified!");status=-1;}
 
  if(in->numofchannels==-1){
    glb_error("numofchannels not specified!");status=-1;}
  
  for(i=0;i<6;i++){
    if(in->listofchannels[i]==NULL){glb_error("listofchannels not specified!");
				      status=-1;}}
  if(in->numofrules==-1){glb_error("numofrules not specified!");status=-1;}
  for(i=0;i<in->numofrules;i++)
    {
      if(in->lengthofrules[i]==-1){glb_error("No targetmass specified!");
      status=-1;}
      if(in->rulescoeff[i]==NULL){glb_error("No rulescoeff specified!");
      status=-1;}
      if(in->rulechannellist[i]==NULL){glb_error("No rulechannellist" 
						  " specified!");status=-1;}
      if(in->signalruleerror[0][i]==-1){glb_error("No signalruleerror[0]"
						   " specified!");status=-1;} 
      if(in->signalruleerror[1][i]==-1){glb_error("No signalruleerros[1]"
						   " specified!");status=-1;}
      if(in->lengthofbgrules[i]==-1){glb_error("No lengthofbgrules"
						" specified!");status=-1;}
      if(in->bgrulescoeff[i]==NULL){glb_error("No bgrulescoeff specified!");
      status=-1;}
      if(in->bgrulechannellist[i]==NULL){glb_error("No bgruloechannellist"
						    " specified!");status=-1;}
      if(in->bgcenter[0][i]==-1){glb_error("No bgcenter[0] specified!");
      status=-1;}
      if(in->bgcenter[1][i]==-1){glb_error("No bgcenter[1] specified!");
      status=-1;}
      if(in->bgerror[0][i]==-1){glb_error("No bgerror[0] specified!");
      status=-1;}
      if(in->bgerror[1][i]==-1){glb_error("No bgerror[1] specified!");
      status=-1;}
      /* if(in->bgtcenter[0][i]==-1){glb_error("No bgtcenter[0] specified!");
      status=-1;} 
      if(in->bgtcenter[1][i]==-1){glb_error("No bgtcenter[1] specified!");
      status=-1;}
      if(in->bgterror[0][i]==-1){glb_error("No bgterror[0] specified!");
      status=-1;}
      if(in->bgterror[1][i]==-1){glb_error("No bgterror[1] specified!");
      status=-1;}*/
    }
 

  if(in->filter_state==-1){in->filter_state=1;def=-1;}
  if(in->filter_value==-1){in->filter_value=0;def=-1;}
  
  if(in->num_of_sm==-1){
    glb_error("No smearing data specified!");status=-1;}
  
  for(i=0;i<in->num_of_sm;i++) 
    {
      
      if(in->smear[i]!=NULL)
	{
	  glb_set_up_smear_data(in->smear_data[i],in);
	  glb_default_smear(in->smear_data[i],in);
	}
      /* Here the smear_data is used to produce a smear matrix */
      if(in->smear[i]==NULL)
	{
	  /* Right now we do not support different binnings
	   * in the reconstructed energy 
	   */
	  in->smear_data[i]->numofbins=in->numofbins;
	  in->smear_data[i]->e_min=in->emin;
	  in->smear_data[i]->e_max=in->emax;


	  if(glb_default_smear(in->smear_data[i],in)==1) def=-1;  

	  glb_compute_smearing_matrix(&in->smear[i], 
				&in->lowrange[i],&in->uprange[i],
				in->smear_data[i],in);
	
	  in->simtresh=in->smear_data[i]->e_sim_min;
	  in->simbeam=in->smear_data[i]->e_sim_max;
	  in->simbins=in->smear_data[i]->simbins;
	
 

      }
      
    

      if(in->smear[i]==NULL){glb_error("No smear matrix defined!");
      status=-1;}
      if(in->lowrange[i]==NULL){glb_error("No lowrange defined!");status=-1;}
      if(in->uprange[i]==NULL){glb_error("No uprange defined!");status=-1;}

    }
      if(in->simtresh==-1){glb_error("No simtresh defined!");status=-1;}
      if(in->simbeam==-1){glb_error("No simbeam defined!");status=-1;}  
      if(in->simbins==-1){glb_error("No simbins defined!");status=-1;}
      if(in->simbins<in->numofbins){glb_error("Less sampling points than"
					      " bins");status=-1;} 
  //---------------------------------------------------------

    *in=MInitMemory0(*in);

  //-------------------------------------------------------

  
  if(in->density_center==-1){in->density_center=1;def=-1;}
  if(in->density_error==-1){in->density_error=0.05;def=-1;}
  if(in->psteps==-1){glb_error("psteps not defined!");status=-1;}
  if(in->lengthtab==NULL){glb_error("lengthtab not allocated!");status=-1;}
  if(in->densitytab==NULL){glb_error("densitytab not allocated!");status=-1;}
  if(in->densitybuffer==NULL){glb_error("densitybuffer not allocated!");
  status=-1;}
  // here we use lot of defaults -- and its not that clear how many
  // there are numofrules ? numofchannels ?
  for(i=0;i<in->numofchannels;i++)
    {
      if(in->simbins<=0) glb_error("Too few simbins defined!");
      else
	{

	  if(in->user_pre_smearing_channel[i]==NULL)
	    {
	      in->user_pre_smearing_channel[i]=(double*) 
		glb_malloc(sizeof(double)* in->simbins);
	      for(k=0;k<in->simbins;k++) 
		in->user_pre_smearing_channel[i][k]=1;
	      def=-1; 
	    }
	  else
	    {
	      for(ct=0;in->user_pre_smearing_channel[i][ct]!=-1;ct++) ct=ct;
	      if(ct!=in->simbins) {glb_error("user_pre_smearing_channel"
					     " has not simbins elements");
	      status=-1;}
	    }
	    
	  if(in->user_pre_smearing_background[i]==NULL)
	    {
	      in->user_pre_smearing_background[i]=(double*) 
		glb_malloc(sizeof(double)*in->simbins);
	      for(k=0;k<in->simbins;k++) 
		{in->user_pre_smearing_background[i][k]=0;}
	      def=-1; 
	    }
	  else
	    {
	      for(ct=0;in->user_pre_smearing_background[i][ct]!=-1;ct++) ct=ct;
	      if(ct!=in->simbins) {glb_error("user_pre_smearing_background"
					     " has not simbins elements");
	      status=-1;}
	    }
	    
	  if(in->user_post_smearing_channel[i]==NULL)
	    {
	      in->user_post_smearing_channel[i]=(double*) 
		glb_malloc(sizeof(double)*in->numofbins);
	      for(k=0;k<in->numofbins;k++) 
		in->user_post_smearing_channel[i][k]=1;
	      def=-1; 
	    }
	  else
	    {
	      for(ct=0;in->user_post_smearing_channel[i][ct]!=-1;ct++) ct=ct;
	      if(ct!=in->numofbins) {glb_error("user_post_smearing_channel"
					     " has not numofbins elements");
	      status=-1;}
	    }
	  if(in->user_post_smearing_background[i]==NULL)
	    {
	      in->user_post_smearing_background[i]=(double*) 
		glb_malloc(sizeof(double)*in->numofbins);
	      for(k=0;k<in->numofbins;k++) 
		{in->user_post_smearing_background[i][k]=0;}
	      def=-1; 
	    }
	  else
	    {
	      for(ct=0;in->user_post_smearing_background[i][ct]!=-1;ct++) 
		ct=ct;
	      if(ct!=in->numofbins) {glb_error("user_post_smearing_background"
					     " has not numofbins elements");
	      status=-1;}
	    }
	}
      if(in->energy_window[i][0]==-1){in->energy_window[i][0]=in->emin;def=-1;}
      if(in->energy_window[i][1]==-1){in->energy_window[i][1]=in->emax;def=-1;}
    }
  if(in->chirate==NULL){glb_error("chirate not allocated!");
  status=-1;}
  for(i=0;i<in->numofrules;i++)
    {
      if(in->SignalRates[i]==NULL||
	 in->BackgroundRates[i]==NULL||
	 in->rates0[i]==NULL||
	 in->rates1[i]==NULL||
	 in->rates1T[i]==NULL||
	 in->rates1BG[i]==NULL||
	 in->rates1BGT[i]==NULL||
	 in->ratevec[i]==NULL||
	 in->ratevecBG[i]==NULL){glb_error("No memory for ratevectors"
					   " allocated!");status=-1;};

      if(in->chrb[i]==NULL||in->chra[i]==NULL
	 ||in->buffer==NULL){glb_error("No memory for convolution"
					   " allocated!");status=-1;};
    }
  if(in->energy_tab==NULL){glb_error("energy_tab not allocated!");
  status=-1;}
  // okay thats missing
  /*
    for(i=0;i<32;i++) in->no_osc_background[i]=NULL;
    for(i=0;i<32;i++) init_struct_systematic(&(in->sys[i]));
  */
  // issue a final report
  if(status!=0) glb_fatal("Incompletely defined experiment!");
  if(def!=0) glb_warning("Incompletely defined experiment!"
			 " Defaults are going to be used!");
  return def+status;
};




//#include "StructExp.c"




// function setting glb_num_of_exps

static void glbSetNumberOfExperiments(int in)
{
  glb_num_of_exps = in;
  return;
}


char* glbSetSys(int select, int rule, int experiment)
{
  glb_experiment_list[experiment]->sys[rule]=sys_list[select];
  return sys_list[select].info;
}

void glbSetSysErrors(double *val, int rule, int experiment)
{
  int i;
  for(i=0;i<(glb_experiment_list[experiment]->sys[rule]).dimension;i++)
    {
      (glb_experiment_list[experiment]->sys[rule]).errors[i]=val[i];
    }
}

void glbSetSysCenter(double *val, int rule, int experiment)
{
  int i;
  for(i=0;i<(glb_experiment_list[experiment]->sys[rule]).dimension;i++)
    {
      (glb_experiment_list[experiment]->sys[rule]).sp[i]=val[i];
    }
}

// this function moves the pointers to the ratevectors to the correct address
// allocated with MInitFile().
// this is crucial for multiple experiment support
// without reprogramming everything !
// here we pass a pointer for efficiency since this function is used in every step

static void MMovePointers(struct experiment *in)
{
  int k;
  for (k=0;k<in->numofrules;k++) 
    { 
      glb_calc_rates_0[k] = in->rates0[k];
      glb_calc_rates_1[k] = in->rates1[k];
      glb_calc_rates_1T[k] = in->rates1T[k];
      glb_calc_rates_1BG[k] = in->rates1BG[k];
      glb_calc_rates_1BGT[k] = in->rates1BGT[k];
      glb_calc_ratevec[k] = in->ratevec[k];
      glb_calc_glb_calc_ratevec_bg[k] = in->ratevecBG[k]; 
      glb_calc_signal_rates[k] = in->SignalRates[k]; 
      glb_calc_bg_rates[k] = in->BackgroundRates[k];
      /* FIXME -- wrong number */ 
      sys_calc[k]=in->sys[k];    
    } 
  for(k=0;k<in->num_of_sm;k++) glb_calc_smear_data[k]=in->smear_data[k];
  for(k=0;k<in->numofchannels;k++)
    {
      glb_calc_chrb[k]=in->chrb[k];
      glb_calc_chra[k]=in->chra[k];
      glb_calc_user_pre_sm_channel[k]=in->user_pre_smearing_channel[k];
      glb_calc_user_post_sm_channel[k]= in->user_post_smearing_channel[k];
      glb_calc_user_pre_sm_background[k]=in->user_pre_smearing_background[k];
      glb_calc_user_post_sm_background[k]=in-> user_post_smearing_background[k];
    }
  
  for(k=0;k<in->num_of_fluxes;k++) glb_calc_fluxes[k]=in->fluxes[k];
  for(k=0;k<in->num_of_xsecs;k++) glb_calc_xsecs[k]=in->xsecs[k];
  glb_chirate = in->chirate;
  glb_calc_energy_tab = in->energy_tab;
  glb_calc_buffer=in->buffer;
  return;
}

// memory allocation for all the ratevecvtors

static struct experiment MInitMemory0(struct experiment in)
{
 struct experiment out;
  int k;
  int len,l2;
  out=in;
  len=out.numofbins;
  if(out.numofbins>=out.simbins) l2=out.numofbins;
  else l2=out.simbins;
  if(len<=0) 
    {glb_error("Too few bins defined!");}
  else
    {
      for (k=0;k<out.numofrules;k++) 
	{ 
	  out.rates0[k] =  (double*) glb_malloc( len*sizeof(double));
	  out.rates1[k] = (double*) glb_malloc( len*sizeof(double));
	  out.rates1T[k] = (double*) glb_malloc( len*sizeof(double));
	  out.rates1BG[k] = (double*) glb_malloc( len*sizeof(double));
	  out.rates1BGT[k] = (double*) glb_malloc( len*sizeof(double));
	  out.ratevec[k] = (double*) glb_malloc( len*sizeof(double));
	  out.ratevecBG[k] = (double*) glb_malloc( len*sizeof(double));  
	  if(out.simbins<=0) glb_error("Too few simbins defined!");
	  else
	    {
	      out.SignalRates[k] = (double*) glb_malloc(len*sizeof(double));
	      out.BackgroundRates[k] = (double*) glb_malloc(len*sizeof(double)); 
	    }
	}
      out.buffer=(double*) glb_malloc(l2*sizeof(double));
      out.energy_tab= (double*) glb_malloc(len*sizeof(double)); 
      out.chirate=(double*) glb_malloc(len*sizeof(double));
    }
  for(k=0;k<out.numofchannels;k++)
    {
      if(out.simbins<=0) glb_error("Too few simbins defined!");
      else out.chrb[k] = (double*) glb_malloc(out.simbins*sizeof(double)); 
      if(out.numofbins<=0) glb_error("Too few bins defined!"); 
      else out.chra[k] = (double*) glb_malloc(out.numofbins*sizeof(double));
    } 
  if(out.psteps<=0) glb_error("Too few density steps defined!");
  else if (out.densitybuffer==NULL) out.densitybuffer=(double*) 
				      glb_malloc(out.psteps*sizeof(double));
  return out;
}






static void glb_set_profile_scaling_sec(struct experiment *in)
{
  int k;
  glb_set_length_ptr(in->lengthtab);
  glb_set_psteps(in->psteps);

  for(k=0;k<in->psteps;k++)
    {
      in->densitybuffer[k]=(double) 1.0*in->densitytab[k];
    }
  glb_set_density_ptr(in->densitybuffer);

}


// This function actually sets all switches and parameters to the values specified
// in an struct experiment object and calls MMovePointers().
// here we pass a pointer for efficiency since this function is used in every step 

void glbSetExperiment(glb_exp in)
{
  int i;
  struct experiment *s;
  s=(struct experiment *) in;
  //  set_type(s->beamtype);
  glb_set_error_function(s->errorfunction);
 

  //analog of SetBase

 


 
 

  // pointer moving for density profile
  glb_set_psteps(s->psteps);
  glb_set_length_ptr(s->lengthtab);
  glb_set_baseline(s->baseline);
  //density=in->densitytab;
  glb_set_profile_scaling_sec(s);
  //analog of SetEnergy

  glb_set_energy_treshold(s->emin);
  glb_set_max_energy(s->emax);  
  glb_set_number_of_bins(s->numofbins);
  
 
  //analog of SetLuminosity
  
  glb_set_target_mass(s->targetmass);

  
  glb_set_num_of_channels(s->numofchannels);
  for (i=0;i<s->numofchannels;i++) glb_set_channel(i,s->listofchannels[0][i],
					       s->listofchannels[1][i],
					       s->listofchannels[2][i],
					       s->listofchannels[3][i],
					       s->listofchannels[4][i],
					       s->listofchannels[5][i]
					       );
 
  // analog of SetRules

 
  glb_set_number_of_rules(s->numofrules);
  for (i=0;i<s->numofrules;i++)
    {
      glb_set_rule(i,s->lengthofrules[i],
		   (int*) s->rulechannellist[i],s->rulescoeff[i]);
      glb_set_signal_errors(i,
			    (double) s->signalruleerror[0][i],
			    (double) s->signalruleerror[1][i]);
      glb_set_bg_rule(i,s->lengthofbgrules[i],
			      s->bgrulechannellist[i],s->bgrulescoeff[i]);
      glb_set_bg_center(i,
			(double) s->bgcenter[0][i],(double) s->bgcenter[1][i]);
      glb_set_bg_errors(i,(double) s->bgerror[0][i],(double) s->bgerror[1][i]);
      glb_calc_set_tresh_center(i,(double) s->bgtcenter[0][i],
				(double) s->bgtcenter[1][i]);
      glb_calc_set_tresh_errors(i,(double) s->bgterror[0][i],
				(double) s->bgterror[1][i]);
      glb_set_errordim(s->errordim[i],i);  
      
      glb_calc_smearing[i]=s->smear[i]; 
       
      //changed on 14 Feb 03
      glb_calc_uprange[i]=s->uprange[i];
      glb_calc_lowrange[i]=s->lowrange[i];
      // this is preliminary (0 should be replaced by i)
      glb_calc_simbins=s->simbins;
      glb_calc_simtresh=s->simtresh;
      glb_calc_simbeam=s->simbeam;
      //--------------------------

      glb_calc_set_energy_window(i,s->energy_window[i][0],
				 s->energy_window[i][1]);
     
    }
  glb_switch_filter(s->filter_state);
  glb_set_filter(s->filter_value);


  
  MMovePointers(s);
  
  glb_rule_number=0;


  return;
}

// redefinition of SetRates()

void glbSetRates()
{
  int i;
  for (i=0;i<glb_num_of_exps;i++)
    {
      glbSetExperiment(glb_experiment_list[i]);
      glb_set_rates();
    }
}

// redefinition of SetNewRates()

void glbSetNewRates()
{
  int i;
  for (i=0;i<glb_num_of_exps;i++)
    {
      glbSetExperiment(glb_experiment_list[i]);
      glb_set_new_rates();
    }
}

// a bunch of comands which allow to modify
// certain elements of a struct experiment object
// in glb_experiment_list


void glb_set_profile_scaling(double scale,int i)
{
  int k;
  glb_set_length_ptr(glb_experiment_list[i]->lengthtab);
 
  glb_set_psteps(glb_experiment_list[i]->psteps);
  for(k=0;k<glb_experiment_list[i]->psteps;k++)
    {
      glb_experiment_list[i]->densitybuffer[k]=(double) 
	scale*glb_experiment_list[i]->densitytab[k];
    }
  glb_set_density_ptr(glb_experiment_list[i]->densitybuffer);

}

