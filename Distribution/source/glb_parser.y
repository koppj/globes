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



     
%{
#define YYDEBUG 1

#include <math.h>
#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <globes/globes.h>
#include "glb_smear.h"
#include "glb_multiex.h"
#include "glb_types.h"
#include "glb_error.h"
#include "glb_lexer.h"  
#include "glb_parser_type.h"
#include "glb_fluxes.h"



#define TWICE 100  

#define UNTYPE -1
#define INT 0
#define DOUBLE 1
#define INT_LIST 2
#define DOUBLE_LIST 3
#define DOUBLE_LIST_INDEXED 4
#define INT_LIST_INDEXED 5
#define ENERGY_MATRIX 6
#define INT_INDEXED 7
#define INT_INDEXED_PAIR 8
#define DOUBLE_INDEXED 9
#define DOUBLE_INDEXED_PAIR 10
#define DOUBLE_INDEXED_PAIR_INV 11
#define CHAR 12
#define FUN 13
#define DOUBLE_LIST_INDEXED_SL 14
#define DOUBLE_LIST_INDEXED_BL 15
#define COUNTER 16
  
  static int exp_count=1;
  static int energy_len;
  static int loc_count; 
  static int energy_count=-1;
  static int cross_count=-1;
  static int flux_count=-1;
  static struct experiment buff;
  static struct experiment buff_list[32];
  static glb_smear ibf;
  static glb_option_type opt;
  static glb_flux flt;
  static glb_xsec xsc;
  static char *context;




  typedef struct
  {
    char *token; // the string which is used as lhs in the ini-file
    int scalar; // data type flag 
    double rl;
    double ru; // allowed range
    void *ptr; // this is a pointer the corresponding part of the exp structure
    void *len; /* this is a pointer which points to the length of the vector
		*  in a struct experiment. Example: if we parse densitytab, 
		*  this things points
		*  to psteps 
		*/
    
    char *ctx; /* here the type of the environment is given, e.g. rule
		* thus the corresponding token is matched onyl within
		* a rule type environment 
		*/
  } glb_parser_decl;
  

  static glb_parser_decl token_list[]={
    {"$version",CHAR,0,1E8,&buff.version,NULL,"global"},
    {"$parent_energy",DOUBLE,0,1E8,&buff.emax,NULL,"global"},
    {"$target_mass",DOUBLE,0,1E8,&buff.targetmass,NULL,"global"},
    
    {"$baseline",DOUBLE,0,2*GLB_EARTH_RADIUS,&buff.baseline,NULL,"global"},
    {"$profiletype",INT,1,3,&buff.density_profile_type,NULL,"global"},
    {"$binsize",DOUBLE_LIST,0,1E8,&buff.binsize,&buff.numofbins,"global"}, 
    {"$simbinsize",DOUBLE_LIST,0,1E8,&buff.simbinsize,
     &buff.simbins,"global"}, 

   
   {"$numofbins",COUNTER,0,500,&buff.numofbins,NULL,"global"},
   {"$densitytab",DOUBLE_LIST,0,100,&buff.densitytab,&buff.psteps,"global"},
   {"$lengthtab",DOUBLE_LIST,0,13000,&buff.lengthtab,&buff.psteps,"global"},
   {"rulechannellist",INT_LIST_INDEXED,0,10E8,&buff.rulechannellist[0],
    &buff.lengthofrules[0],"rule"},
   {"rulescoeff",DOUBLE_LIST_INDEXED,0,10E8,&buff.rulescoeff[0],
    &buff.lengthofrules[0],"rule"},
   
   {"bgrulechannellist",INT_LIST_INDEXED,0,10E8,&buff.bgrulechannellist[0],
    &buff.lengthofbgrules[0],"rule"},
   {"bgrulescoeff",DOUBLE_LIST_INDEXED,0,10E8,&buff.bgrulescoeff[0],
    &buff.lengthofbgrules[0],"rule"},
   {"rule",UNTYPE,-1,1,NULL,&buff.numofrules,"global"},
   {"@signal",UNTYPE,-1,1,NULL,&loc_count,"rule"},
   {"@background",UNTYPE,-1,1,NULL,&loc_count,"rule"},

    {"@pre_smearing_efficiencies",DOUBLE_LIST_INDEXED,0,1E8,
     &buff.user_pre_smearing_channel[0],&loc_count,"channel"},
   
    {"@pre_smearing_background",DOUBLE_LIST_INDEXED,0,1E8,
     &buff.user_pre_smearing_background[0],&loc_count,"channel"},

    {"@post_smearing_efficiencies",DOUBLE_LIST_INDEXED,0,1E8,
     &buff.user_post_smearing_channel[0],&loc_count,"channel"},
   
    {"@post_smearing_background",DOUBLE_LIST_INDEXED,0,1E8,
     &buff.user_post_smearing_background[0],&loc_count,"channel"},
    
   
   {"@errordim_sys_on",INT_INDEXED,0,20
    ,&buff.errordim_sys_on[0],&loc_count,"rule"},
    {"@errordim_sys_off",INT_INDEXED,0,20
    ,&buff.errordim_sys_off[0],&loc_count,"rule"},
   {"channel",UNTYPE,-1,1,NULL,&buff.numofchannels,"global"},
   {"@channel",INT_LIST_INDEXED,-1,1,&buff.listofchannels[0],
    &buff.numofchannels,"channel"},
   {"energy",UNTYPE,-1,1,NULL,&buff.num_of_sm,"global"},
   {"@energy",ENERGY_MATRIX,-1,1,&buff.smear[0],&loc_count,"energy"},

   {"@signalerror",DOUBLE_INDEXED_PAIR
    ,0,100,&buff.signalruleerror[0],&loc_count,"rule"},
   {"@backgrounderror",DOUBLE_INDEXED_PAIR,
    0,100,&buff.bgerror[0],&loc_count,"rule"},
   {"@backgroundcenter",DOUBLE_INDEXED_PAIR
    ,0,100,&buff.bgcenter[0],&loc_count,"rule"},
   {"@treshold_setttings" ,DOUBLE_INDEXED_PAIR,
    0,100,&buff.bgtcenter[0],&loc_count,"rule"},
   {"@treshold_error" ,DOUBLE_INDEXED_PAIR,
    0,100,&buff.bgterror[0],&loc_count,"rule"}, 

   {"@energy_window" ,DOUBLE_INDEXED_PAIR_INV,0,100,&buff.energy_window[0],
    &loc_count,"rule"},
   {"$densitysteps",COUNTER,1,200,&buff.psteps,NULL,"global"},
  
   {"$errorfunction",INT,0,1,&buff.errorfunction,NULL,"global"},

   {"$filter_state",INT,0,1,&buff.filter_state,NULL,"global"},
   {"$filter_value",DOUBLE,0,1E8,&buff.filter_value,NULL,"global"},
   {"$simbins" ,COUNTER,0,500,&buff.simbins,NULL,"global"}, 
   {"$simtresh" ,DOUBLE,0,1E8,&buff.simtresh,NULL,"global"},
   {"$simbeam" ,DOUBLE,0,1E8,&buff.simbeam,NULL,"global"},


   {"$emin",DOUBLE,0,1E8,&buff.emin,NULL,"global"}, 
   {"$emax",DOUBLE,0,1E8,&buff.emax,NULL,"global"},
   {"cross",UNTYPE,0,20,NULL,&buff.num_of_xsecs,"global"},
   {"@cross_file",CHAR,0,20,&xsc.file_name,NULL,"cross"},


   {"flux",UNTYPE,0,20,NULL,&buff.num_of_fluxes,"global"},
   {"@flux_file",CHAR,0,20,&flt.file_name,NULL,"flux"},

    {"@builtin",INT,0,2,&flt.builtin,NULL,"flux"},
    {"@time",DOUBLE,0,200,&flt.time,NULL,"flux"},
    {"@power",DOUBLE,0,1E8,&flt.target_power,NULL,"flux"},
    {"@stored_muons",DOUBLE,0,1E30,&flt.stored_muons,NULL,"flux"},
    {"@parent_energy",DOUBLE,0,1E8,&flt.parent_energy,NULL,"flux"},
    {"@norm",DOUBLE,0,1E30,&flt.norm,NULL,"flux"},
    
 

   {"@type" ,INT,0,3,&ibf.type,&loc_count,"energy"},
   {"@sigma_e" ,DOUBLE_LIST,0,10,&ibf.sigma,&ibf.num_of_params,"energy"},
   {"@sigma_function" ,FUN,0,1000,&ibf.sig_f,NULL,"energy"},
    
  

   {NULL,UNTYPE,0,0,NULL,NULL,"global"}

  
};



 
static void grp_end(char* name)
   {    
     
     if(strncmp(name,"energy",6)==0 )
       {
	 if(buff.num_of_sm-1 >= 0)   
	   {
	     ibf.options=glb_option_type_alloc();
	     ibf.options=(glb_option_type *) memmove(ibf.options,&opt,
					     sizeof(glb_option_type));
	     
	     if(buff.smear_data[buff.num_of_sm-1]==NULL)
	       buff.smear_data[buff.num_of_sm-1]=glb_smear_alloc();
	     buff.smear_data[buff.num_of_sm-1]=
	       glb_copy_smear(buff.smear_data[buff.num_of_sm-1],&ibf);
	     glb_option_type_free(ibf.options);
	     glb_option_type_reset(&opt);
	     if(ibf.sigma!=NULL) glb_free(ibf.sigma);
	     glb_smear_reset(&ibf);
	     
	   }

       } 

     if(strncmp(name,"flux",4)==0 )
       {
	 if(buff.num_of_fluxes > 0)   
	   {
	     if(buff.fluxes[buff.num_of_fluxes-1]==NULL)
	       buff.fluxes[buff.num_of_fluxes-1]=glb_flux_alloc();
	     buff.fluxes[buff.num_of_fluxes-1]=
	       cpy_glb_flux(buff.fluxes[buff.num_of_fluxes-1],&flt);
	    glb_flux_reset(&flt);
	   }
       }  

     if(strncmp(name,"cross",4)==0 )
       {
	 if(buff.num_of_xsecs > 0)   
	   {
	     if(buff.xsecs[buff.num_of_xsecs-1]==NULL)
	       buff.xsecs[buff.num_of_xsecs-1]=glb_xsec_alloc();
	     buff.xsecs[buff.num_of_xsecs-1]=
	       cpy_glb_xsec(buff.xsecs[buff.num_of_xsecs-1],&xsc);
	     glb_xsec_reset(&xsc);
	   }
       }  

     glb_free(context);
     context = (char *) strdup("global");
     
   }
  
static int set_channel_data(int x[6],int loc_count)
   {
     // I am sorry for this -- it's a kludge
     int i;

     for(i=0;i<6;i++) {
       buff.listofchannels[i]=
			(int*) glb_realloc(buff.listofchannels[i]
					,sizeof(int)*loc_count);

       buff.listofchannels[i][loc_count-1]=x[i];
     }
     

     return 0;
   } 


static int step_counter(char *name)
{
  int i=0;
  int *ibf;
  ibf=NULL;
  for(i=0;token_list[i].token!=NULL;i++)
    {
      if(strncmp(name,token_list[i].token,
		 strlen(token_list[i].token))==0 )
	{
	 
		     ibf=(int*) token_list[i].len;
		     if(*ibf==-1) *ibf=1;// first time encounter
		     else (*ibf)++;

	}
    }
  if(ibf!=NULL) return *ibf; // otherwise a SEGFAULT occurs when there
  // is no matching name
  return 0;
}

static int set_fnct(char *name,void *in)
{
  int i;
  sigfun *dbf;
  
 
  for(i=0;token_list[i].token!=NULL;i++)
    {
      if(strncmp(name,token_list[i].token,
		 strlen(token_list[i].token))==0&&
	 strncmp(context,token_list[i].ctx,
		 strlen(token_list[i].ctx))==0 )
	{
	     if(token_list[i].scalar==FUN) //double
	       {
		 dbf=(sigfun *) token_list[i].ptr;
		 *dbf=(sigfun) in;	 
	      	 return 0;
	       }
	     else       
	       {
		 fprintf(stderr,"Error: Value for %s out of range\n",
			 token_list[i].token);
		 return 2;
	       }
	}
      
      
    }   
  
  
  return 1;   
}


static int set_exp(char *name,double value,int scalar)
{
  int i;
  double *dbf;
  int *ibf;
  int *xibf;
  for(i=0;token_list[i].token!=NULL;i++)
    {
      if(strncmp(name,token_list[i].token,
		 strlen(token_list[i].token))==0 &&
	 strncmp(context,token_list[i].ctx,
		 strlen(token_list[i].ctx))==0
	 )
	{
	     if(token_list[i].scalar==DOUBLE) //double
	       {
		 if(value >= token_list[i].rl && value <= token_list[i].ru)
		   {
		     dbf=(double*) token_list[i].ptr;
		     *dbf=value;
		     return 0;
		   }
		 else       
		   {
		     fprintf(stderr,"Error: Value for %s out of range\n",
			     token_list[i].token);
		     return 2;
		   }
	       }
	     if(token_list[i].scalar==INT) //int
	       {
		 if(value >= token_list[i].rl && value <= token_list[i].ru)
		   {
		     ibf=(int*) token_list[i].ptr;
		     *ibf=value;
		     return 0;
		   }
		 else       
		   {
		     fprintf(stderr,"Error: Value for %s out of range\n",
			     token_list[i].token);
		     return 2;
		   }
	       }

	     if(token_list[i].scalar==COUNTER) //counter
	       {
		 if(value >= token_list[i].rl && value <= token_list[i].ru)
		   {
		     
		     ibf=(int*) token_list[i].ptr;
		     if(!((*ibf == -1) || (*ibf == (int) value))) {
		       glb_warning("Given length does not"
				   " match actual length");
		       return 2;

		     }
		     *ibf=value;
		     return 0;
		   }
		 else       
		   {
		     fprintf(stderr,"Error: Value for %s out of range\n",
			     token_list[i].token);
		     return 2;
		   }
	       }

	     if(token_list[i].scalar==INT_INDEXED) //int
	       {
		 if(value >= token_list[i].rl && value <= token_list[i].ru)
		   {
		     xibf=(int*) token_list[i].ptr;
		     xibf[loc_count-1]=(int) value;
		     return 0;
		   }
		 else       
		   {
		     fprintf(stderr,"Error: Value for %s out of range\n",
			     token_list[i].token);
		     return 2;
		   }
	       }
	     
	     if(token_list[i].scalar==DOUBLE_INDEXED) //int
	       {
		 if(value >= token_list[i].rl && value <= token_list[i].ru)
		   {
		     dbf=(double*) token_list[i].ptr;
		     dbf[loc_count-1]=(double) value;
		     return 0;
		   }
		 else       
		   {
		     fprintf(stderr,"Error: Value for %s out of range\n",
			     token_list[i].token);
		     return 2;
		   }
	       }

	  
	   }   
       }

   return 1;   
}

static int set_pair(char *name,double value,double value2,int scalar)
{
  int i;
  double *dbf;
  int *ibf;
  for(i=0;token_list[i].token!=NULL;i++)
    {
      if(strncmp(name,token_list[i].token,
		 strlen(token_list[i].token))==0&&
	 strncmp(context,token_list[i].ctx,
		 strlen(token_list[i].ctx))==0 )
	{
	     
	     if(token_list[i].scalar==INT_INDEXED_PAIR) //int
	       {
		 if(value >= token_list[i].rl && value <= token_list[i].ru)
		   {

		     ibf=(int*) token_list[i].ptr;
		     ibf[(loc_count-1)+0*32]=(int) value;
		     return 0;
		   }
		 else       
		   {
		     fprintf(stderr,"Error: Value for %s out of range\n",
			     token_list[i].token);
		     return 2;
		   }
	       }
	     
	     if(token_list[i].scalar==DOUBLE_INDEXED_PAIR) //int
	       {
		 if(value >= token_list[i].rl && value <= token_list[i].ru)
		   {
		    
		     dbf=(double*) token_list[i].ptr;
		     dbf[(loc_count-1)+0*32]=(double) value;
		     dbf[(loc_count-1)+1*32]=(double) value2;
		     return 0;
		   }
		 else       
		   {
		     fprintf(stderr,"Error: Value for %s out of range\n",
			     token_list[i].token);
		     return 2;
		   }
	       }

	     if(token_list[i].scalar==DOUBLE_INDEXED_PAIR_INV) //int
	       {
		 if(value >= token_list[i].rl && value <= token_list[i].ru)
		   {
		    
		     dbf=(double*) token_list[i].ptr;
		     dbf[(loc_count-1)*32+0]=(double) value;
		     dbf[(loc_count-1)*32+1]=(double) value2;
		     return 0;
		   }
		 else       
		   {
		     fprintf(stderr,"Error: Value for %s out of range\n",
			     token_list[i].token);
		     return 2;
		   }
	       }
	   }   
       }

   return 1;   
}

static void list_free(glb_List *stale)
{
 glb_List *ptr;
 glb_List *dummy;  
 ptr=stale;
 while(ptr != (glb_List *) NULL)
   {	
     dummy=ptr->next;
     glb_free(ptr);
     ptr=dummy;  
   }
 
}



     
static glb_List *list_cons (glb_List *tail, double newdata)
{
  glb_List *res = (glb_List*) glb_malloc(sizeof(glb_List)); 
  res->next=tail;
  res->entry=newdata;
  return res;
}


static size_t list_length (glb_List *head)
{
  size_t n;  
  for (n = 0; head; ++n)
    head = head->next;
  return n; 
}

static void showlist(glb_List *lp)
{ 
 
  if (lp){
    showlist(lp->next);             // show the tail
    printf("%.10g " , lp->entry);     // show the head
  }

}

static double list_take(glb_List *li,int k)
{
  int i;
  double erg;
  glb_List *bf;
  bf=li; 
  erg=bf->entry;
  for(i=0;i<k;i++)
    { 
      bf=bf->next;
      erg=bf->entry;
    }
  return erg;
}


static int set_exp_list(char *name,glb_List *value,int scalar)
{
  int i,k;
  double **dbf;
  int **ibf;
  int *lbf;
  int len;
  double val;
  double *list;
  int *ilist;
  for(i=0;token_list[i].token!=NULL;i++)
    {
      if(strncmp(name,token_list[i].token,strlen(token_list[i].token))==0 &&
	 strncmp(context,token_list[i].ctx,
		 strlen(token_list[i].ctx))==0)
	{
	  switch((int) token_list[i].scalar) {
	  case DOUBLE_LIST:
	    
	    //here we will have to do a lot asking asf.
	    len=list_length(value); // how long is the list
	    lbf=(int*) token_list[i].len;
	    if(*lbf==-1) *lbf=len;  // setting the length correctly in exp   
	    else if(*lbf!=len) glb_warning("Length mismatch or list"
					   " length changed");
	  
	    
	    dbf = (double**) token_list[i].ptr;
	    if(*dbf!=NULL)glb_free(*dbf);
	    list=(double*) glb_malloc(sizeof(double)*len);
	    *dbf=list; 
	    
	    
	    for(k=0;k<len;k++)
	       {
		  val=list_take(value,len-k-1);
		
		  if(val >= token_list[i].rl && val <= token_list[i].ru)
		    {
		      list[k]=val;      
		    }
		  else       
		    {
		      fprintf(stderr,"Error: Value for %s out of range\n",
			      token_list[i].token);
		      glb_free(list);
		      list=NULL;
		      return 2;
		    }
		  
	       }
	    if(scalar!=TWICE)  list_free(value);
	    return 0;
	    break;
	  
	  case DOUBLE_LIST_INDEXED: 
	    len=list_length(value); // how long is the list
	    lbf=(int*) token_list[i].len;
	    
	   
	    //  lbf[loc_count-1]=len;  // setting the length correctly in exp   
	  
	    dbf= (double**) token_list[i].ptr;
	    if(dbf[loc_count-1]!=NULL)glb_free(dbf[loc_count-1]);
	    list=(double*) glb_malloc(sizeof(double)*(len+1));
	  
	    dbf[loc_count-1]=list; 
	    list[len]=-1;
	    
	    for(k=0;k<len;k++)
	      {
		val=list_take(value,len-k-1);
		  
		if(val >= token_list[i].rl && val <= token_list[i].ru)
		  {
		    list[k]=val;

		  }
		else       
		  {
		    fprintf(stderr,"Error: In line %d: "
			    "Value for %s out of range\n",
			    glb_line_num,token_list[i].token);
		   glb_free(list);
		    return 2;
		  }
	      }
	    if(scalar!=TWICE) list_free(value);
	    return 0;
	    break;


	  case INT_LIST_INDEXED: //integer list indexed 
	    //with loc_counter
	    
	    //here we will have to do a lot asking asf.
	    len=list_length(value); // how long is the list
	    lbf=(int*) token_list[i].len; //FIXME danger !!!!
	    lbf[loc_count-1]=len;  // setting the length correctly in exp   
	    ibf= (int**) token_list[i].ptr;
	    if(ibf[loc_count-1]!=NULL)glb_free(ibf[loc_count-1]);
	    ilist=(int*) glb_malloc(sizeof(int)*len);
	   
	    ibf[loc_count-1]=ilist; 
	
	      
	    for(k=0;k<len;k++)
	      {
		val=list_take(value,len-k-1);
		
		if(val >= token_list[i].rl && val <= token_list[i].ru)
		  {
		    ilist[k]=(int) val;
		  
		  }
		else       
		  {
		    fprintf(stderr,"Error: Value for %s out of range\n",
			    token_list[i].token);
		   glb_free(ilist);
		    return 2;
		  }
	      }
	    if(scalar!=TWICE) list_free(value);
	    return 0;
	    break;	 
	  default:
	    return 1;
	    break;
	  }
	}
    }    
  return 1;    
}

static int set_exp_energy(char *name, glb_List **value)
{
  int i,k,l;
  double ***dbf;  
  int len;
  double val;
  int v1,v2;
  double **list;

  for(i=0;token_list[i].token!=NULL;i++)
    {
      if(strncmp(name,token_list[i].token,strlen(token_list[i].token))==0&&
	 strncmp(context,token_list[i].ctx,
		 strlen(token_list[i].ctx))==0 )
	{
	  switch((int) token_list[i].scalar) {
	  
	  case ENERGY_MATRIX:
	    list=(double**) glb_malloc(sizeof(double* ) * energy_len);
	    buff.lowrange[loc_count-1]=(int*) glb_malloc(energy_len*sizeof(int));
	    buff.uprange[loc_count-1]=(int*) glb_malloc(energy_len*sizeof(int));
      
	    for(l=0;l<energy_len;l++)
	      {
		len=(int) list_length(value[l]); // how long is the list
		if(len<2) {fprintf(stderr,"Error: in line %d: in @smear "
				   "number %d: sublist %d is too short!\n"
				   ,glb_line_num,loc_count,l);return 2;}
		//lbf=(int*) token_list[i].len;
	    
	   
		//lbf[loc_count-1]=len;  // setting the length correctly in exp
	
		
		
		dbf= (double***) token_list[i].ptr;
		list[l]=(double*) glb_malloc(sizeof(double)*(len-2));
		dbf[loc_count-1]=list; 

		
		v1=(int) list_take(value[l],len-0-1);
		v2=(int) list_take(value[l],len-1-1);
	      
		if(v1 >= 0 &&  v2 <= buff.simbins
		   && v2 >= v1&&v2-v1==len-3 )
		  {
		   
		    buff.lowrange[loc_count-1][l]= v1;	
		    buff.uprange[loc_count-1][l]= v2;
		    
		  }
		else
		  {
		    fprintf(stderr,"Error: In line %d: "
			    "Value for ranges in smear out of range\n",
			    glb_line_num);
		    glb_free(list[l]);
		    glb_free(buff.lowrange[loc_count-1]);
		    glb_free(buff.uprange[loc_count-1]);
		    return 2;
		  }
		
	      	for(k=0;k<len-2;k++)
		  {
		    val=list_take(value[l],len-(k+2)-1);

		    if(val >= token_list[i].rl && val <= token_list[i].ru)
		      {
			list[l][k]=val;
		
		      }
		    else       
		      {
			fprintf(stderr,"Error: In line %d: "
				"Value for %s out of range\n",
				glb_line_num,token_list[i].token);
			free(list[l]);
			return 2;
		      }
		  }
		list_free(value[l]);
	      }
	    glb_free(value);
	    return 0;
	    break;

	  }
	}
    }    
  return 1;    
}


 %}
%union {
  double  val;  /* For returning numbers.                   */
  double *dpt;  /* for rules */
  glb_List *ptr; 
  glb_List **ptrq;
  glb_symrec  *tptr;  /* For returning symbol-table pointers      */
  char *name;
  char *iname;
  int in;
  glb_namerec *nameptr;
  
}

%token <val> NUM
%token <nameptr> SFNCT
%token <tptr> VAR FNCT   /* Variable and Function */
%token <name> IDN CROSS FLUXP FLUXM
%token <name> GRP GID FNAME VERS
%token <in> SIGNAL BG GRPOPEN GRPCLOSE PM FLAVOR
%token <in> NOGLOBES CHANNEL
%token <in> RULESEP RULEMULT ENERGY
%token <nameptr> NAME RDF NDEF
%type <ptr> seq 
%type <ptrq> rule brule srule ene
%type <dpt> rulepart 
%type <val> exp 
%type <val> group
%type <in> pm
%type <nameptr> name
%type <name> cross 
%type <name> flux
%type <name> version

%expect 11
%nonassoc ','
%right '='
%left RULESEP
%left RULEMULT
%right '>'
%right '}'
%left  '{'
%left '-' '+'
%left '*' '/'
%left NEG     /* Negation--unary minus */
%right '^'    /* Exponentiation        */

%% /* Grammar rules and actions follow */

input:    /* empty */
| input line 
| NOGLOBES {YYABORT;}
| '�' {YYABORT;}
;
     
line:   '\n'
| outchannel '\n'
| outrule '\n' { /*showlist($1[0]);showlist($1[1]);*/ }
| exp '\n'  { /*printf ("\t%.10g\n", $1);*/  }
| seq '\n'  { /*showlist($1);printf("\n");*/ }
| group '\n' {}
| error '\n' {yyerrok;}
| channel '\n' {}
| ene ';' '\n' {
 set_exp_energy("@energy",$1);
}
| rule '\n' {}
| cross '\n' {/* printf("%s\n",$1);*/}
| flux '\n' {/* printf("%s\n",$1);*/}
;



exp:      NUM                { $$ = $1;                         }
| NAME                { $$ = $1->value;              }
| VAR                { $$ = $1->value.var;              }
| VAR '=' exp        { $$ = $3; $1->value.var = $3;     }
| IDN '=' exp        { if(set_exp($1,$3,0)==1) yyerror("Unknown identifier");
$$ = $3; }
| IDN '=' SFNCT        { if(set_fnct($1,$3->sf)==1) yyerror("Unknown identifier");}
| IDN '=' exp RULESEP exp  { if(set_pair($1,$3,$5,0)==1) yyerror("Unknown identifier");
$$ = $3; }
| FNCT '(' exp ')'   { $$ = (*($1->value.fnctptr))($3); }
| exp '+' exp        { $$ = $1 + $3;                    }
| exp '-' exp        { $$ = $1 - $3;                    }
| exp '*' exp        { $$ = $1 * $3;                    }
| exp '/' exp        { $$ = $1 / $3;                    }
| '-' exp  %prec NEG { $$ = -$2;                        }
| exp '^' exp        { $$ = pow ($1, $3);               }
| '(' exp ')'        { $$ = $2;                         } 
| version            { $$ = 0;}
| NDEF               {yyerror("Unknown name");YYERROR;}
;




seq:    exp  ','        {$$=list_cons(NULL,$1); }   
|    exp  ',' '\n'        {$$=list_cons(NULL,$1); }         
|  seq  exp ',' {
  glb_List *buf;
  buf=$1;
  buf=list_cons(buf,$2);
  $$=buf;   
}

|  seq  exp ',' '\n' {
  glb_List *buf;
  buf=$1;
  buf=list_cons(buf,$2);
  $$=buf;   
}
|  seq  exp  {
  glb_List *buf;
  buf=$1;
  buf=list_cons(buf,$2);
  $$=buf;   
}
| '{' seq '}'  {$$=$2; }
| '{' exp '}'  {$$=list_cons(NULL,$2); }
| IDN '=' seq        { if(set_exp_list($1,$3,3)==1) 
  yyerror("Unknown identifier");
 $$ = $3; }
;

rulepart: exp RULEMULT exp { 

  double *buf;
  buf=(double*) glb_malloc(sizeof(double)*2);
  buf[0]=$1;
  buf[1]=$3-1;
  $$=buf;
}
;

expseq: line {}
|expseq line {}
;

group: GID '(' NAME ')' 
{ if($3->value==-1) {$3->value=step_counter($1); }
  loc_count=$3->value;
  glb_free(context);
  context =(char *) strdup($1);
}
  GRPOPEN ingroup  GRPCLOSE { 
   
grp_end($1);}
| GID '(' RDF ')' GRPOPEN ingroup  GRPCLOSE {
  yyerror("Redefinition of an automatic variable");YYERROR;}  
;

ingroup: /* empty */ 
| exp {}
| rule {}
| expseq {}
| channel {}
| cross {}
| flux {}
;

version: VERS '=' FNAME {
  buff.version=strdup($3);
}

cross: CROSS '=' FNAME {
  //load_cross($3,loc_count-1);
  xsc.file_name=strdup($3);
  $$=$3;
}
;

flux: FLUXP '=' FNAME {
  //load_flux($3,loc_count-1,1);
  flt.file_name=strdup($3);
  //if(set_exp($1,$3,0)==1) yyerror("Unknown identifier");
  $$=$3;
}
;


channel: CHANNEL '=' name RULESEP pm RULESEP FLAVOR RULESEP FLAVOR RULESEP
 name RULESEP name {

  int x[6];
  x[0]=$3->value - 1 ;
  x[1]=$5;
  x[2]=$7;
  x[3]=$9;
  x[4]=(int) $11->value -1; 
  x[5]=(int) $13->value -1;


  set_channel_data(x,loc_count);



}
;

name: NAME {$$=$1}
|NDEF {yyerror("Unknown name");YYERROR;}
;
 

pm: PM {$$=$1;}
| '+' {$$=1;}
| '-' {$$=-1;}
;


ene: ENERGY '='  seq   {
  glb_List **buf;
  energy_len=1;
 
  buf=(glb_List**) glb_malloc(sizeof( glb_List* ) ); 
  buf[0]=$3;
 
  $$=buf;
}
| ene RULESEP seq 
{
  glb_List **buf;
  buf=$1;
  energy_len++;
  
  buf=(glb_List**) glb_realloc((void**) buf , sizeof( glb_List* ) * energy_len);
 
  buf[energy_len-1]=$3; 
  $$=buf; }
| ene RULESEP '\n' seq 
{
  glb_List **buf;
  buf=$1;
  energy_len++;

  buf=(glb_List**) glb_realloc((void**) buf , sizeof( glb_List* ) * energy_len);

  buf[energy_len-1]=$4; 
  $$=buf; }
;



brule: BG '=' rulepart {
  glb_List **buf;
  buf=(glb_List**) glb_malloc(sizeof(glb_List*)*2);
  buf[0]=list_cons(NULL,$3[0]);
  buf[1]=list_cons(NULL,$3[1]);
  glb_free($3);
  $$=buf;
}
| brule RULESEP rulepart {
  glb_List **buf;
  buf=$1;
  buf[0]=list_cons(buf[0],$3[0]);
  buf[1]=list_cons(buf[1],$3[1]);
  glb_free($3);
  $$=buf; 
}
;

srule: SIGNAL '=' rulepart {
  glb_List **buf;  
 
  buf=(glb_List**) glb_malloc(sizeof(glb_List*)*2);
  buf[0]=list_cons(NULL,$3[0]);
  buf[1]=list_cons(NULL,$3[1]);
  glb_free($3);
  $$=buf;
}
| srule RULESEP rulepart {
  glb_List **buf;
  buf=$1;
  buf[0]=list_cons(buf[0],$3[0]);
  buf[1]=list_cons(buf[1],$3[1]);
  glb_free($3);
  $$=buf; 
}
;

rule: brule {

int flag;  
  
  $$=$1;
  
  flag=set_exp_list("bgrulescoeff",$1[0],0);
  if(flag==1) yyerror("Unknown identifier");
  flag=set_exp_list("bgrulechannellist",$1[1],0);
  if(flag==1) yyerror("Unknown identifier");
 
  glb_free($1);
}
| srule {
  int flag;  
  $$=$1;
  flag=set_exp_list("rulescoeff",$1[0],0);
  if(flag==1) yyerror("Unknown identifier");
  flag=set_exp_list("rulechannellist",$1[1],0);
  if(flag==1) yyerror("Unknown identifier"); 

  glb_free($1);

 
}
;

outrule: NAME {
 if($1->value==-1) {$1->value=step_counter("rule");printf("loc step\n"); }
  loc_count=$1->value;
  printf("loc_count %d \n",loc_count);
  YYERROR;}
  rule { YYERROR;}
;


outchannel: NAME {
 if($1->value==-1) {$1->value=step_counter("channel");printf("loc step\n"); }
  loc_count=$1->value;
  printf("cha loc_count %d \n",loc_count); YYERROR;}
  channel { YYERROR;}
;

%%

extern glb_symrec *sym_table;

 
int
yyerror (const char *s)  /* Called by yyparse on error */
{
  if(yydebug==1) fprintf(stderr,"*****************************************\n");
  fprintf (stderr,"%s: ERROR: %s in file [%s]: in line %d\n",
	   glb_prog_name,s,glb_file_id,glb_line_num);
  if(yydebug==1) fprintf(stderr,"*****************************************\n");
  return 0;
}

int
yywarn (const char *s)  /* Called by yyparse on warning */
{
  fprintf (stderr,"%s: Warning: %s in line %d\n",
	   glb_prog_name,s,glb_line_num);
  return 0;
}


struct glb_init
{
  char *fname;
  double (*fnct)(double);
};

struct glb_init_sig
{
  char *fname;
  sigfun sf;
};
     


static struct glb_init arith_fncts[] =
  {
    {"sin",  sin},
    {"cos",  cos},
    {"atan", atan},
    {"ln",   log},
    {"exp",  exp},
    {"sqrt", sqrt},
    {NULL, NULL}
     };

static double standard_sig(double x ,double *in)
{

  return in[0]*x + in[1]*sqrt(x) + in[2];
}


static double chooz_sig(double x, double *in)
{
  if(x<=0.0018) return in[0]/1000.0;
  else return in[0]*sqrt((x-0.0018+0.001)*1000.0)/1000.0;
}

/* Number of parameters ? */
struct glb_init_sig sig_fncts[] =
  {
    {"#standard",  standard_sig},
    {"#inverse_beta", chooz_sig},
    {NULL, NULL}
  };

/* The symbol table: a chain of `struct glb_symrec'.  */
static glb_namerec *name_table = (glb_namerec *) NULL; 
static glb_symrec *sym_table = (glb_symrec *) NULL;
static glb_symrec *pre_sym_table = (glb_symrec *) NULL;

/* Put arithmetic functions in table. 
 * And all user-defined stuff.
 */
static void
init_table (void)
{
  int i;
  glb_symrec *ptr,*p;
  
  glb_namerec *sptr;
  for (i = 0; arith_fncts[i].fname != 0; i++)
    {
      ptr = glb_putsym (arith_fncts[i].fname, FNCT);
      ptr->value.fnctptr = arith_fncts[i].fnct;
    }
  ptr=pre_sym_table;
  /* Copying the contents of pre_sym_table to sym_table */
  for (i=0; ptr !=0; i++)
    {
      p=glb_putsym(ptr->name,VAR);
      p->value.var=ptr->value.var;
      ptr=ptr->next;
    }
  

  for (i = 0; sig_fncts[i].fname != 0; i++)
    {
      sptr = glb_putname (sig_fncts[i].fname,"energy",SFNCT);
      sptr->sf = sig_fncts[i].sf;
    }
}

static void
free_symtable()
{
    glb_symrec *ptr;
    glb_symrec *dummy;
    ptr=sym_table;
    while(ptr != (glb_symrec *) NULL)
      {	
	glb_free(ptr->name);
	dummy=ptr->next;
	glb_free(ptr);
	ptr=dummy;  
      }
    sym_table=NULL;
}

static void
free_presymtable()
{
    glb_symrec *ptr;
    glb_symrec *dummy;
    ptr=pre_sym_table;
    while(ptr != (glb_symrec *) NULL)
      {	
	glb_free(ptr->name);
	dummy=ptr->next;
	glb_free(ptr);
	ptr=dummy;  
      }
    pre_sym_table=NULL;
}

static void
free_nametable()
{
    glb_namerec *ptr;  
    glb_namerec *dummy;
    ptr=name_table;
    while(ptr != (glb_namerec *) NULL)
      {	
	glb_free(ptr->name);
	glb_free(ptr->context);
	dummy=ptr->next;
	glb_free(ptr);
	ptr=dummy;  
      }
    name_table=NULL;
}

glb_symrec *
glb_putsym (char *sym_name, int sym_type)
{
  glb_symrec *ptr;
  ptr = (glb_symrec *) glb_malloc (sizeof (glb_symrec));
  ptr->name = (char *) glb_malloc (strlen (sym_name) + 1);
  strcpy (ptr->name,sym_name);
  ptr->type = sym_type;
  ptr->value.var = 0; /* set value to 0 even if fctn.  */
  ptr->next = (struct glb_symrec *) sym_table;
  sym_table = ptr;
  return ptr;
}





/* Name handling */



glb_namerec *glb_getname (const char *sym_name, char* context)
{
  glb_namerec *ptr;
  for (ptr = name_table; ptr != (glb_namerec *) NULL;
       ptr = (glb_namerec *)ptr->next)
    if (strcmp (ptr->name,sym_name) == 0)
	return ptr;
  return 0;
}




glb_namerec *glb_putname (char *sym_name, char *context, int sym_type)
{
  glb_namerec *ptr;
  ptr = (glb_namerec *) glb_malloc (sizeof (glb_namerec));
  ptr->name = (char *) glb_malloc (strlen (sym_name) + 1);
  strcpy (ptr->name,sym_name);
  ptr->context = (char *) glb_malloc (strlen (context) + 1);
  strcpy (ptr->context,context);
  ptr->type = sym_type;
  ptr->value = -1; /* set value to -1 for new ones  */
  ptr->next = (struct glb_namerec *)name_table;
  name_table = ptr;
  return ptr;
}

glb_symrec *
glb_getsym (const char *sym_name)
{
  glb_symrec *ptr;
  for (ptr = sym_table; ptr != (glb_symrec *) NULL;
       ptr = (glb_symrec *)ptr->next)
    if (strcmp (ptr->name,sym_name) == 0)
      return ptr;
  return 0;
}
/* User pre-def variables */
glb_symrec *
glb_putpresym (const char *sym_name, int sym_type)
{
  glb_symrec *ptr;
  ptr = (glb_symrec *) glb_malloc (sizeof (glb_symrec));
  ptr->name = (char *) glb_malloc (strlen (sym_name) + 1);
  strcpy (ptr->name,sym_name);
  ptr->type = sym_type;
  ptr->value.var = 0; /* set value to 0 even if fctn.  */
  ptr->next = (struct glb_symrec *) pre_sym_table;
  pre_sym_table = ptr;
  return ptr;
}

glb_symrec *
glb_getpresym (const char *sym_name)
{
  glb_symrec *ptr;
  for (ptr = pre_sym_table; ptr != (glb_symrec *) NULL;
       ptr = (glb_symrec *)ptr->next)
    if (strcmp (ptr->name,sym_name) == 0)
      return ptr;
  return 0;
}


/* A new an powerful function which allows the user to define variables
 * for substitution in the AEDL files 
 */


void glbDefineAEDLVariable(const char *name, double value)
{
  glb_symrec *ptr;
  ptr=glb_getpresym(name);
  if(ptr==0) ptr = glb_putpresym (name, VAR);
  ptr->value.var = value;
  return;
}

void glbClearAEDLVariables()
{
   if(pre_sym_table!=NULL) free_presymtable(); 
}



void glb_copy_buff()
{
  buff_list[exp_count]=buff;
  exp_count++;
}

void glbReset() 
{
  glb_line_num=1;
  energy_len=1;
  energy_count=-1;
  loc_count=-1;
  flux_count=-1;
  glbInitExp(&buff);
  if(name_table!=NULL) free_nametable();
  if(sym_table!=NULL) free_symtable();
  name_table =(glb_namerec *) NULL;
  sym_table =(glb_symrec *) NULL;
  init_table ();
}

void glbResetCounters()
{ 
  flux_count=-1;
  cross_count=-1;
}

void glbResetEOF() 
{
  int i;
  exp_count=0;
  glb_line_num=1;
  energy_len=1;
  energy_count=-1;
  loc_count=-1;
  flux_count=-1;
  glbInitExp(&buff);
  for(i=0;i<32;i++)   glbInitExp(&buff_list[i]);
  if(name_table!=NULL) free_nametable();
  if(sym_table!=NULL) free_symtable();

  name_table =(glb_namerec *) NULL;
  sym_table =(glb_symrec *) NULL;
  init_table ();
}

void glb_clean_parser()
{
  if(name_table!=NULL) free_nametable();
  if(sym_table!=NULL) free_symtable();
  if(pre_sym_table!=NULL) free_presymtable();
  
}

int glbInitExperiment(char *inf,glb_exp *in, int *counter)
{
  FILE *input;
  int k,i;
  const char tch[]="%!GLoBES"; 
  char tct[11];
  struct experiment **ins;
  ins=(struct experiment **) in;
  //yydebug=1;
  context=(char *) strdup("global");
  glb_smear_reset(&ibf);
  glb_option_type_reset(&opt);
  glb_flux_reset(&flt);
  input=glb_fopen(inf,"r");
  if(input==NULL) return -1;
  /* This line produces a warning with -Wall:
   * 'warning: char format, different type arg (arg 3)'
   * I think however that understand the problem.
   */
  fscanf(input,"%8c",&tct);
  if(strncmp(tch,tct,8)!=0) {glb_warning("Not a GLoBES file!"); return 1;}
  glb_fclose(input);
  yyin=glb_fopen(inf,"r");
  if(yyin==NULL) return -1;
  glb_file_id=(char*) strdup(inf);
  glbResetEOF();
  k=yyparse ();
  glb_fclose(yyin);
  glb_free(context);
  glb_free(glb_file_id);
  if(k!=0) return -1;
  
  k=0;
  if(*counter+exp_count>32) glb_fatal("Too many experiments!");
  for(i=0;i<exp_count;i++)
    {  
      *ins[*counter+i]=buff_list[i];
      k=+glbDefaultExp(ins[*counter+i]);
    }
  (*counter)= (*counter) + exp_count;

  if(k!=0) return -2;
  
  return 0;

}