/* A Bison parser, made by GNU Bison 1.875.  */

/* Skeleton parser for Yacc-like parsing with Bison,
   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002 Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330,
   Boston, MA 02111-1307, USA.  */

/* As a special exception, when this file is copied by Bison into a
   Bison output file, you may use that output file without restriction.
   This special exception was added by the Free Software Foundation
   in version 1.24 of Bison.  */

/* Written by Richard Stallman by simplifying the original so called
   ``semantic'' parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Using locations.  */
#define YYLSP_NEEDED 0



/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     NUM = 258,
     SFNCT = 259,
     VAR = 260,
     FNCT = 261,
     IDN = 262,
     CROSS = 263,
     FLUXP = 264,
     FLUXM = 265,
     GRP = 266,
     GID = 267,
     FNAME = 268,
     VERS = 269,
     SIGNAL = 270,
     BG = 271,
     GRPOPEN = 272,
     GRPCLOSE = 273,
     PM = 274,
     FLAVOR = 275,
     NOGLOBES = 276,
     CHANNEL = 277,
     RULESEP = 278,
     RULEMULT = 279,
     ENERGY = 280,
     NAME = 281,
     RDF = 282,
     NDEF = 283,
     NEG = 284
   };
#endif
#define NUM 258
#define SFNCT 259
#define VAR 260
#define FNCT 261
#define IDN 262
#define CROSS 263
#define FLUXP 264
#define FLUXM 265
#define GRP 266
#define GID 267
#define FNAME 268
#define VERS 269
#define SIGNAL 270
#define BG 271
#define GRPOPEN 272
#define GRPCLOSE 273
#define PM 274
#define FLAVOR 275
#define NOGLOBES 276
#define CHANNEL 277
#define RULESEP 278
#define RULEMULT 279
#define ENERGY 280
#define NAME 281
#define RDF 282
#define NDEF 283
#define NEG 284




/* Copy the first part of user declarations.  */
#line 26 "glb_parser.y"

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
#define GMAX 1E100

  static int exp_count=1;
  static int energy_len;
  static int loc_count; 
  static int energy_count=-1;
  static int cross_count=-1;
  static int flux_count=-1;
  static struct glb_experiment buff;
  static struct glb_experiment buff_list[32];
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
		*  in a struct glb_experiment. Example: if we parse densitytab, 
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
    {"$parent_energy",DOUBLE,0,GMAX,&buff.emax,NULL,"global"},
    {"$target_mass",DOUBLE,0,GMAX,&buff.targetmass,NULL,"global"},
#ifdef GLB_OLD_AEDL
    {"$simbins" ,COUNTER,0,500,&buff.simbins,NULL,"global"}, 
    {"$simtresh" ,DOUBLE,0,GMAX,&buff.simtresh,NULL,"global"},
    {"$simbeam" ,DOUBLE,0,GMAX,&buff.simbeam,NULL,"global"},  
    {"$numofbins",COUNTER,0,500,&buff.numofbins,NULL,"global"},
    {"$simbinsize",DOUBLE_LIST,0,GMAX,&buff.simbinsize,
     &buff.simbins,"global"},
    {"$errorfunction",INT,0,1,&buff.errorfunction,NULL,"global"}, 
    {"@treshold_setttings" ,DOUBLE_INDEXED_PAIR,
     0,100,&buff.bgtcenter[0],&loc_count,"rule"},
    {"@treshold_error" ,DOUBLE_INDEXED_PAIR,
     0,100,&buff.bgterror[0],&loc_count,"rule"}, 
#endif /* GLB_OLD_AEDL */
    {"$sampling_points" ,COUNTER,0,500,&buff.simbins,NULL,"global"}, 
    {"$sampling_min" ,DOUBLE,0,GMAX,&buff.simtresh,NULL,"global"},
    {"$sampling_max" ,DOUBLE,0,GMAX,&buff.simbeam,NULL,"global"},
    /* FIXME -- this was not properly recognized due to the fact the first match is performed
     * by matching only the letters up to the length of the word ins this list, ie.
     * if 'bins' is before 'binsize' in this list, 'binsize' will be matched as
     * 'bins'
     */
    {"$binsize",DOUBLE_LIST,0,GMAX,&buff.binsize,
     &buff.numofbins,"global"},  
   
    {"$bins",COUNTER,0,500,&buff.numofbins,NULL,"global"},
    {"$emin",DOUBLE,0,GMAX,&buff.emin,NULL,"global"}, 
    {"$emax",DOUBLE,0,GMAX,&buff.emax,NULL,"global"},

    {"$baseline",DOUBLE,0,2*GLB_EARTH_RADIUS,&buff.baseline,NULL,"global"},
    {"$profiletype",INT,1,3,&buff.density_profile_type,NULL,"global"},
 
    {"$sampling_stepsize",DOUBLE_LIST,0,GMAX,&buff.simbinsize,
     &buff.simbins,"global"}, 

 

   
  
   {"$densitytab",DOUBLE_LIST,0,GMAX,&buff.densitytab,&buff.psteps,"global"},
   {"$lengthtab",DOUBLE_LIST,0,2*GLB_EARTH_RADIUS,
    &buff.lengthtab,&buff.psteps,"global"},
   {"rulechannellist",INT_LIST_INDEXED,0,GMAX,&buff.rulechannellist[0],
    &buff.lengthofrules[0],"rule"},
   {"rulescoeff",DOUBLE_LIST_INDEXED,0,GMAX,&buff.rulescoeff[0],
    &buff.lengthofrules[0],"rule"},
   
   {"bgrulechannellist",INT_LIST_INDEXED,0,GMAX,&buff.bgrulechannellist[0],
    &buff.lengthofbgrules[0],"rule"},
   {"bgrulescoeff",DOUBLE_LIST_INDEXED,0,GMAX,&buff.bgrulescoeff[0],
    &buff.lengthofbgrules[0],"rule"},
   {"rule",UNTYPE,-1,1,NULL,&buff.numofrules,"global"},
   {"@signal",UNTYPE,-1,1,NULL,&loc_count,"rule"},
   {"@background",UNTYPE,-1,1,NULL,&loc_count,"rule"},

    {"@pre_smearing_efficiencies",DOUBLE_LIST_INDEXED,0,GMAX,
     &buff.user_pre_smearing_channel[0],&loc_count,"channel"},
   
    {"@pre_smearing_background",DOUBLE_LIST_INDEXED,0,GMAX,
     &buff.user_pre_smearing_background[0],&loc_count,"channel"},

    {"@post_smearing_efficiencies",DOUBLE_LIST_INDEXED,0,GMAX,
     &buff.user_post_smearing_channel[0],&loc_count,"channel"},
   
    {"@post_smearing_background",DOUBLE_LIST_INDEXED,0,GMAX,
     &buff.user_post_smearing_background[0],&loc_count,"channel"},
    
   
   {"@errordim_sys_on",INT_INDEXED,0,20
    ,&buff.errordim_sys_on[0],&loc_count,"rule"},
    {"@errordim_sys_off",INT_INDEXED,0,20
    ,&buff.errordim_sys_off[0],&loc_count,"rule"},
   {"channel",UNTYPE,-1,1,NULL,&buff.numofchannels,"global"},
   {"@channel",INT_LIST_INDEXED,-1,1,&buff.listofchannels[0],
    &buff.numofchannels,"channel"},
   {"energy",UNTYPE,-1,1,NULL,&buff.num_of_sm,"global"},
   {"@energy",ENERGY_MATRIX,-1,GMAX,&buff.smear[0],&loc_count,"energy"},

   {"@signalerror",DOUBLE_INDEXED_PAIR
    ,0,100,&buff.signalruleerror[0],&loc_count,"rule"},
   {"@backgrounderror",DOUBLE_INDEXED_PAIR,
    0,100,&buff.bgerror[0],&loc_count,"rule"},
   {"@backgroundcenter",DOUBLE_INDEXED_PAIR
    ,0,100,&buff.bgcenter[0],&loc_count,"rule"},
   

   {"@energy_window" ,DOUBLE_INDEXED_PAIR_INV,0,GMAX,&buff.energy_window[0],
    &loc_count,"rule"},
   {"$densitysteps",COUNTER,1,GMAX,&buff.psteps,NULL,"global"},
  
  

   {"$filter_state",INT,GLB_OFF,GLB_ON,&buff.filter_state,NULL,"global"},
   {"$filter_value",DOUBLE,0,GMAX,&buff.filter_value,NULL,"global"},

   {"cross",UNTYPE,0,20,NULL,&buff.num_of_xsecs,"global"},
   {"@cross_file",CHAR,0,20,&xsc.file_name,NULL,"cross"},


   {"flux",UNTYPE,0,20,NULL,&buff.num_of_fluxes,"global"},
   {"@flux_file",CHAR,0,20,&flt.file_name,NULL,"flux"},

    {"@builtin",INT,1,4,&flt.builtin,NULL,"flux"},
    {"@time",DOUBLE,0,GMAX,&flt.time,NULL,"flux"},
    {"@power",DOUBLE,0,GMAX,&flt.target_power,NULL,"flux"},
    {"@stored_muons",DOUBLE,0,GMAX,&flt.stored_muons,NULL,"flux"},
    {"@parent_energy",DOUBLE,0,GMAX,&flt.parent_energy,NULL,"flux"},
    {"@norm",DOUBLE,0,GMAX,&flt.norm,NULL,"flux"},
    {"@gamma",DOUBLE,0,GMAX,&flt.gamma,NULL,"flux"},
    {"@end_point",DOUBLE,0,GMAX,&flt.end_point,NULL,"flux"},
    {"@stored_ions",DOUBLE,0,GMAX,&flt.stored_muons,NULL,"flux"},
 

   {"@type" ,INT,1,2,&ibf.type,&loc_count,"energy"},
   {"@sigma_e" ,DOUBLE_LIST,0,GMAX,&ibf.sigma,&ibf.num_of_params,"energy"},
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


 

/* Enabling traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

#if ! defined (YYSTYPE) && ! defined (YYSTYPE_IS_DECLARED)
#line 828 "glb_parser.y"
typedef union YYSTYPE {
  double  val;  /* For returning numbers.                   */
  double *dpt;  /* for rules */
  glb_List *ptr; 
  glb_List **ptrq;
  glb_symrec  *tptr;  /* For returning symbol-table pointers      */
  char *name;
  char *iname;
  int in;
  glb_namerec *nameptr;
  
} YYSTYPE;
/* Line 191 of yacc.c.  */
#line 949 "glb_parser.c"
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif



/* Copy the second part of user declarations.  */


/* Line 214 of yacc.c.  */
#line 961 "glb_parser.c"

#if ! defined (yyoverflow) || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# if YYSTACK_USE_ALLOCA
#  define YYSTACK_ALLOC alloca
# else
#  ifndef YYSTACK_USE_ALLOCA
#   if defined (alloca) || defined (_ALLOCA_H)
#    define YYSTACK_ALLOC alloca
#   else
#    ifdef __GNUC__
#     define YYSTACK_ALLOC __builtin_alloca
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's `empty if-body' warning. */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (0)
# else
#  if defined (__STDC__) || defined (__cplusplus)
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   define YYSIZE_T size_t
#  endif
#  define YYSTACK_ALLOC malloc
#  define YYSTACK_FREE free
# endif
#endif /* ! defined (yyoverflow) || YYERROR_VERBOSE */


#if (! defined (yyoverflow) \
     && (! defined (__cplusplus) \
	 || (YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  short yyss;
  YYSTYPE yyvs;
  };

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (short) + sizeof (YYSTYPE))				\
      + YYSTACK_GAP_MAXIMUM)

/* Copy COUNT objects from FROM to TO.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if 1 < __GNUC__
#   define YYCOPY(To, From, Count) \
      __builtin_memcpy (To, From, (Count) * sizeof (*(From)))
#  else
#   define YYCOPY(To, From, Count)		\
      do					\
	{					\
	  register YYSIZE_T yyi;		\
	  for (yyi = 0; yyi < (Count); yyi++)	\
	    (To)[yyi] = (From)[yyi];		\
	}					\
      while (0)
#  endif
# endif

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack)					\
    do									\
      {									\
	YYSIZE_T yynewbytes;						\
	YYCOPY (&yyptr->Stack, Stack, yysize);				\
	Stack = &yyptr->Stack;						\
	yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
	yyptr += yynewbytes / sizeof (*yyptr);				\
      }									\
    while (0)

#endif

#if defined (__STDC__) || defined (__cplusplus)
   typedef signed char yysigned_char;
#else
   typedef short yysigned_char;
#endif

/* YYFINAL -- State number of the termination state. */
#define YYFINAL  4
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   365

/* YYNTOKENS -- Number of terminals. */
#define YYNTOKENS  45
/* YYNNTS -- Number of nonterminals. */
#define YYNNTS  24
/* YYNRULES -- Number of rules. */
#define YYNRULES  77
/* YYNRULES -- Number of states. */
#define YYNSTATES  152

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   284

#define YYTRANSLATE(YYX) 						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const unsigned char yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
      41,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
      43,    44,    36,    35,    29,    34,     2,    37,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,    42,
       2,    30,    31,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,    39,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    33,     2,    32,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,    40,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    38
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const unsigned char yyprhs[] =
{
       0,     0,     3,     4,     7,     9,    11,    13,    16,    19,
      22,    25,    28,    31,    34,    38,    41,    44,    47,    49,
      51,    53,    57,    61,    65,    71,    76,    80,    84,    88,
      92,    95,    99,   103,   105,   107,   110,   114,   118,   123,
     126,   130,   134,   138,   142,   144,   147,   148,   157,   165,
     166,   168,   170,   172,   174,   176,   178,   182,   186,   190,
     204,   206,   208,   210,   212,   214,   218,   222,   227,   231,
     235,   239,   243,   245,   247,   248,   252,   253
};

/* YYRHS -- A `-1'-separated list of the rules' RHS. */
static const yysigned_char yyrhs[] =
{
      46,     0,    -1,    -1,    46,    47,    -1,    21,    -1,    40,
      -1,    41,    -1,    67,    41,    -1,    65,    41,    -1,    48,
      41,    -1,    49,    41,    -1,    52,    41,    -1,     1,    41,
      -1,    58,    41,    -1,    61,    42,    41,    -1,    64,    41,
      -1,    56,    41,    -1,    57,    41,    -1,     3,    -1,    26,
      -1,     5,    -1,     5,    30,    48,    -1,     7,    30,    48,
      -1,     7,    30,     4,    -1,     7,    30,    48,    23,    48,
      -1,     6,    43,    48,    44,    -1,    48,    35,    48,    -1,
      48,    34,    48,    -1,    48,    36,    48,    -1,    48,    37,
      48,    -1,    34,    48,    -1,    48,    39,    48,    -1,    43,
      48,    44,    -1,    55,    -1,    28,    -1,    48,    29,    -1,
      48,    29,    41,    -1,    49,    48,    29,    -1,    49,    48,
      29,    41,    -1,    49,    48,    -1,    33,    49,    32,    -1,
      33,    48,    32,    -1,     7,    30,    49,    -1,    48,    24,
      48,    -1,    47,    -1,    51,    47,    -1,    -1,    12,    43,
      26,    44,    53,    17,    54,    18,    -1,    12,    43,    27,
      44,    17,    54,    18,    -1,    -1,    48,    -1,    64,    -1,
      51,    -1,    58,    -1,    56,    -1,    57,    -1,    14,    30,
      13,    -1,     8,    30,    13,    -1,     9,    30,    13,    -1,
      22,    30,    59,    23,    60,    23,    20,    23,    20,    23,
      59,    23,    59,    -1,    26,    -1,    28,    -1,    19,    -1,
      35,    -1,    34,    -1,    25,    30,    49,    -1,    61,    23,
      49,    -1,    61,    23,    41,    49,    -1,    16,    30,    50,
      -1,    62,    23,    50,    -1,    15,    30,    50,    -1,    63,
      23,    50,    -1,    62,    -1,    63,    -1,    -1,    26,    66,
      64,    -1,    -1,    26,    68,    58,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const unsigned short yyrline[] =
{
       0,   876,   876,   877,   878,   879,   882,   883,   884,   885,
     886,   887,   888,   889,   890,   893,   894,   895,   900,   901,
     902,   903,   904,   906,   907,   909,   910,   911,   912,   913,
     914,   915,   916,   917,   918,   924,   925,   926,   933,   939,
     945,   946,   947,   952,   962,   963,   967,   966,   975,   979,
     980,   981,   982,   983,   984,   985,   988,   992,   999,  1008,
    1028,  1029,  1033,  1034,  1035,  1039,  1048,  1058,  1072,  1080,
    1090,  1099,  1109,  1122,  1136,  1136,  1145,  1145
};
#endif

#if YYDEBUG || YYERROR_VERBOSE
/* YYTNME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals. */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "NUM", "SFNCT", "VAR", "FNCT", "IDN", 
  "CROSS", "FLUXP", "FLUXM", "GRP", "GID", "FNAME", "VERS", "SIGNAL", 
  "BG", "GRPOPEN", "GRPCLOSE", "PM", "FLAVOR", "NOGLOBES", "CHANNEL", 
  "RULESEP", "RULEMULT", "ENERGY", "NAME", "RDF", "NDEF", "','", "'='", 
  "'>'", "'}'", "'{'", "'-'", "'+'", "'*'", "'/'", "NEG", "'^'", 
  "'\\247'", "'\\n'", "';'", "'('", "')'", "$accept", "input", "line", 
  "exp", "seq", "rulepart", "expseq", "group", "@1", "ingroup", "version", 
  "cross", "flux", "channel", "name", "pm", "ene", "brule", "srule", 
  "rule", "outrule", "@2", "outchannel", "@3", 0
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[YYLEX-NUM] -- Internal token number corresponding to
   token YYLEX-NUM.  */
static const unsigned short yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283,    44,
      61,    62,   125,   123,    45,    43,    42,    47,   284,    94,
     167,    10,    59,    40,    41
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const unsigned char yyr1[] =
{
       0,    45,    46,    46,    46,    46,    47,    47,    47,    47,
      47,    47,    47,    47,    47,    47,    47,    47,    48,    48,
      48,    48,    48,    48,    48,    48,    48,    48,    48,    48,
      48,    48,    48,    48,    48,    49,    49,    49,    49,    49,
      49,    49,    49,    50,    51,    51,    53,    52,    52,    54,
      54,    54,    54,    54,    54,    54,    55,    56,    57,    58,
      59,    59,    60,    60,    60,    61,    61,    61,    62,    62,
      63,    63,    64,    64,    66,    65,    68,    67
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const unsigned char yyr2[] =
{
       0,     2,     0,     2,     1,     1,     1,     2,     2,     2,
       2,     2,     2,     2,     3,     2,     2,     2,     1,     1,
       1,     3,     3,     3,     5,     4,     3,     3,     3,     3,
       2,     3,     3,     1,     1,     2,     3,     3,     4,     2,
       3,     3,     3,     3,     1,     2,     0,     8,     7,     0,
       1,     1,     1,     1,     1,     1,     3,     3,     3,    13,
       1,     1,     1,     1,     1,     3,     3,     4,     3,     3,
       3,     3,     1,     1,     0,     3,     0,     3
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const unsigned char yydefact[] =
{
       2,     4,     5,     0,     1,     0,    18,    20,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,    19,    34,
       0,     0,     6,     0,     3,     0,     0,     0,    33,     0,
       0,     0,     0,    72,    73,     0,     0,     0,    12,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,    19,     0,     0,     0,    30,     0,    35,     0,
       0,     0,     0,     0,     9,    10,    39,    11,    16,    17,
      13,     0,     0,     0,     0,    15,     8,     7,    21,     0,
      23,    22,    42,    57,    58,     0,     0,    56,     0,    70,
      68,    60,    61,     0,     0,    65,    75,    77,    41,    40,
       0,    32,    36,    27,    26,    28,    29,    31,    37,     0,
      66,    14,    69,    71,    25,     0,    46,     0,     0,     0,
      22,    38,    67,    24,     0,     0,    43,    62,    64,    63,
       0,     0,    44,    50,     0,     0,    54,    55,    53,    51,
       0,     0,    45,    48,     0,    47,     0,     0,     0,     0,
       0,    59
};

/* YYDEFGOTO[NTERM-NUM]. */
static const short yydefgoto[] =
{
      -1,     3,   132,    66,    26,    89,   134,    27,   124,   135,
      28,    29,    30,    31,    93,   130,    32,    33,    34,    35,
      36,    50,    37,    51
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -117
static const short yypact[] =
{
     -15,  -117,  -117,   133,  -117,   -37,  -117,   -22,   -24,    -6,
       4,    12,     9,    17,    39,    42,    51,    57,    61,  -117,
     235,   301,  -117,   301,  -117,   297,   253,    58,  -117,    59,
      64,    66,   -16,    69,    79,    73,    83,    84,  -117,   301,
     301,     7,    91,   107,    63,   110,   301,   301,    47,   235,
      80,   105,  -117,   311,   269,    99,    96,    82,    89,   301,
     301,   301,   301,   301,  -117,  -117,   117,  -117,  -117,  -117,
    -117,    60,   102,   301,   301,  -117,  -117,  -117,   326,   286,
    -117,    -7,   301,  -117,  -117,    93,   100,  -117,    74,  -117,
    -117,  -117,  -117,   127,   320,   301,  -117,  -117,  -117,  -117,
     285,  -117,  -117,    25,    25,    96,    96,    96,   116,   235,
     301,  -117,  -117,  -117,  -117,   301,  -117,   143,   301,    20,
      -7,  -117,   301,   326,   145,   172,   326,  -117,  -117,  -117,
     140,   172,  -117,   297,   211,   146,    59,    64,    66,    73,
     148,   147,  -117,  -117,   149,  -117,   150,   159,    47,   160,
      47,  -117
};

/* YYPGOTO[NTERM-NUM].  */
static const yysigned_char yypgoto[] =
{
    -117,  -117,    -2,    -3,   -18,   -25,  -117,  -117,  -117,    38,
    -117,  -116,   -80,   -46,   -68,  -117,  -117,  -117,  -117,   -47,
    -117,  -117,  -117,  -117
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -77
static const short yytable[] =
{
      25,    24,    54,    96,    38,    97,     1,    71,    39,   136,
       6,    80,     7,     8,     9,   136,   115,    53,    56,    40,
      57,    13,    90,    82,    41,     2,    72,    59,    60,    61,
      62,    95,    63,    52,    42,    19,    78,    79,    81,   127,
      20,    21,    43,    88,    88,   137,    94,    45,   112,   113,
      23,   137,    44,   110,   128,   129,   103,   104,   105,   106,
     107,    61,    62,     6,    63,     7,     8,     9,    94,    46,
      88,    88,    47,    91,    13,    92,   -74,   -74,   139,   138,
     149,    48,   151,   -76,   139,   138,    52,    49,    19,    85,
      86,   122,    73,    20,    21,    14,    15,   120,   118,    67,
      68,   109,    74,    23,    83,    69,    94,    70,    59,    60,
      61,    62,   123,    63,    75,   126,    59,    60,    61,    62,
      84,    63,   133,    87,    76,    77,   101,    16,   133,   100,
     102,    25,   142,     4,     5,    63,     6,   116,     7,     8,
       9,    10,    11,   111,   117,    12,   108,    13,    14,    15,
     119,    59,    60,    61,    62,    16,    63,   121,    17,    18,
     125,    19,   131,   140,   143,   145,    20,    21,   144,   141,
     147,     0,   146,     5,    22,     6,    23,     7,     8,     9,
      10,    11,   148,   150,    12,     0,    13,    14,    15,     0,
     -49,     0,     0,     0,    16,     0,     0,    17,    18,     0,
      19,     0,     0,     0,     0,    20,    21,     0,     0,     0,
       0,     0,     5,    22,     6,    23,     7,     8,     9,    10,
      11,     0,     0,    12,     0,    13,    14,    15,     0,   -52,
       0,     0,     0,    16,     0,     0,    17,    18,     6,    19,
       7,     8,     9,     0,    20,    21,     0,     0,     0,    13,
       0,     0,    22,     0,    23,     0,     6,     0,     7,     8,
      55,    52,     0,    19,     0,     0,     0,    13,    20,    21,
       0,     0,     6,     0,     7,     8,    55,     0,    23,    52,
       0,    19,     0,    13,     0,     0,     0,    21,     6,    80,
       7,     8,    55,     0,    65,    52,    23,    19,     0,    13,
       0,    99,     0,    21,     6,     0,     7,     8,    55,     0,
       0,    52,    23,    19,     0,    13,     0,     0,     0,    21,
      59,    60,    61,    62,     0,    63,    58,    52,    23,    19,
     114,    59,    60,    61,    62,    21,    63,     0,    64,     0,
      58,     0,     0,    98,    23,    59,    60,    61,    62,    58,
      63,     0,     0,     0,    59,    60,    61,    62,     0,    63,
      59,    60,    61,    62,     0,    63
};

static const short yycheck[] =
{
       3,     3,    20,    50,    41,    51,    21,    23,    30,   125,
       3,     4,     5,     6,     7,   131,    23,    20,    21,    43,
      23,    14,    47,    41,    30,    40,    42,    34,    35,    36,
      37,    49,    39,    26,    30,    28,    39,    40,    41,    19,
      33,    34,    30,    46,    47,   125,    49,    30,    73,    74,
      43,   131,    43,    71,    34,    35,    59,    60,    61,    62,
      63,    36,    37,     3,    39,     5,     6,     7,    71,    30,
      73,    74,    30,    26,    14,    28,    15,    16,   125,   125,
     148,    30,   150,    22,   131,   131,    26,    30,    28,    26,
      27,   109,    23,    33,    34,    15,    16,   100,    24,    41,
      41,    41,    23,    43,    13,    41,   109,    41,    34,    35,
      36,    37,   115,    39,    41,   118,    34,    35,    36,    37,
      13,    39,   125,    13,    41,    41,    44,    22,   131,    30,
      41,   134,   134,     0,     1,    39,     3,    44,     5,     6,
       7,     8,     9,    41,    44,    12,    29,    14,    15,    16,
      23,    34,    35,    36,    37,    22,    39,    41,    25,    26,
      17,    28,    17,    23,    18,    18,    33,    34,    20,   131,
      20,    -1,    23,     1,    41,     3,    43,     5,     6,     7,
       8,     9,    23,    23,    12,    -1,    14,    15,    16,    -1,
      18,    -1,    -1,    -1,    22,    -1,    -1,    25,    26,    -1,
      28,    -1,    -1,    -1,    -1,    33,    34,    -1,    -1,    -1,
      -1,    -1,     1,    41,     3,    43,     5,     6,     7,     8,
       9,    -1,    -1,    12,    -1,    14,    15,    16,    -1,    18,
      -1,    -1,    -1,    22,    -1,    -1,    25,    26,     3,    28,
       5,     6,     7,    -1,    33,    34,    -1,    -1,    -1,    14,
      -1,    -1,    41,    -1,    43,    -1,     3,    -1,     5,     6,
       7,    26,    -1,    28,    -1,    -1,    -1,    14,    33,    34,
      -1,    -1,     3,    -1,     5,     6,     7,    -1,    43,    26,
      -1,    28,    -1,    14,    -1,    -1,    -1,    34,     3,     4,
       5,     6,     7,    -1,    41,    26,    43,    28,    -1,    14,
      -1,    32,    -1,    34,     3,    -1,     5,     6,     7,    -1,
      -1,    26,    43,    28,    -1,    14,    -1,    -1,    -1,    34,
      34,    35,    36,    37,    -1,    39,    29,    26,    43,    28,
      44,    34,    35,    36,    37,    34,    39,    -1,    41,    -1,
      29,    -1,    -1,    32,    43,    34,    35,    36,    37,    29,
      39,    -1,    -1,    -1,    34,    35,    36,    37,    -1,    39,
      34,    35,    36,    37,    -1,    39
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const unsigned char yystos[] =
{
       0,    21,    40,    46,     0,     1,     3,     5,     6,     7,
       8,     9,    12,    14,    15,    16,    22,    25,    26,    28,
      33,    34,    41,    43,    47,    48,    49,    52,    55,    56,
      57,    58,    61,    62,    63,    64,    65,    67,    41,    30,
      43,    30,    30,    30,    43,    30,    30,    30,    30,    30,
      66,    68,    26,    48,    49,     7,    48,    48,    29,    34,
      35,    36,    37,    39,    41,    41,    48,    41,    41,    41,
      41,    23,    42,    23,    23,    41,    41,    41,    48,    48,
       4,    48,    49,    13,    13,    26,    27,    13,    48,    50,
      50,    26,    28,    59,    48,    49,    64,    58,    32,    32,
      30,    44,    41,    48,    48,    48,    48,    48,    29,    41,
      49,    41,    50,    50,    44,    23,    44,    44,    24,    23,
      48,    41,    49,    48,    53,    17,    48,    19,    34,    35,
      60,    17,    47,    48,    51,    54,    56,    57,    58,    64,
      23,    54,    47,    18,    20,    18,    23,    20,    23,    59,
      23,    59
};

#if ! defined (YYSIZE_T) && defined (__SIZE_TYPE__)
# define YYSIZE_T __SIZE_TYPE__
#endif
#if ! defined (YYSIZE_T) && defined (size_t)
# define YYSIZE_T size_t
#endif
#if ! defined (YYSIZE_T)
# if defined (__STDC__) || defined (__cplusplus)
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# endif
#endif
#if ! defined (YYSIZE_T)
# define YYSIZE_T unsigned int
#endif

#define yyerrok		(yyerrstatus = 0)
#define yyclearin	(yychar = YYEMPTY)
#define YYEMPTY		(-2)
#define YYEOF		0

#define YYACCEPT	goto yyacceptlab
#define YYABORT		goto yyabortlab
#define YYERROR		goto yyerrlab1

/* Like YYERROR except do call yyerror.  This remains here temporarily
   to ease the transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  */

#define YYFAIL		goto yyerrlab

#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)					\
do								\
  if (yychar == YYEMPTY && yylen == 1)				\
    {								\
      yychar = (Token);						\
      yylval = (Value);						\
      yytoken = YYTRANSLATE (yychar);				\
      YYPOPSTACK;						\
      goto yybackup;						\
    }								\
  else								\
    { 								\
      yyerror ("syntax error: cannot back up");\
      YYERROR;							\
    }								\
while (0)

#define YYTERROR	1
#define YYERRCODE	256

/* YYLLOC_DEFAULT -- Compute the default location (before the actions
   are run).  */

#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)         \
  Current.first_line   = Rhs[1].first_line;      \
  Current.first_column = Rhs[1].first_column;    \
  Current.last_line    = Rhs[N].last_line;       \
  Current.last_column  = Rhs[N].last_column;
#endif

/* YYLEX -- calling `yylex' with the right arguments.  */

#ifdef YYLEX_PARAM
# define YYLEX yylex (YYLEX_PARAM)
#else
# define YYLEX yylex ()
#endif

/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)			\
do {						\
  if (yydebug)					\
    YYFPRINTF Args;				\
} while (0)

# define YYDSYMPRINT(Args)			\
do {						\
  if (yydebug)					\
    yysymprint Args;				\
} while (0)

# define YYDSYMPRINTF(Title, Token, Value, Location)		\
do {								\
  if (yydebug)							\
    {								\
      YYFPRINTF (stderr, "%s ", Title);				\
      yysymprint (stderr, 					\
                  Token, Value);	\
      YYFPRINTF (stderr, "\n");					\
    }								\
} while (0)

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (cinluded).                                                   |
`------------------------------------------------------------------*/

#if defined (__STDC__) || defined (__cplusplus)
static void
yy_stack_print (short *bottom, short *top)
#else
static void
yy_stack_print (bottom, top)
    short *bottom;
    short *top;
#endif
{
  YYFPRINTF (stderr, "Stack now");
  for (/* Nothing. */; bottom <= top; ++bottom)
    YYFPRINTF (stderr, " %d", *bottom);
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)				\
do {								\
  if (yydebug)							\
    yy_stack_print ((Bottom), (Top));				\
} while (0)


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

#if defined (__STDC__) || defined (__cplusplus)
static void
yy_reduce_print (int yyrule)
#else
static void
yy_reduce_print (yyrule)
    int yyrule;
#endif
{
  int yyi;
  unsigned int yylineno = yyrline[yyrule];
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %u), ",
             yyrule - 1, yylineno);
  /* Print the symbols being reduced, and their result.  */
  for (yyi = yyprhs[yyrule]; 0 <= yyrhs[yyi]; yyi++)
    YYFPRINTF (stderr, "%s ", yytname [yyrhs[yyi]]);
  YYFPRINTF (stderr, "-> %s\n", yytname [yyr1[yyrule]]);
}

# define YY_REDUCE_PRINT(Rule)		\
do {					\
  if (yydebug)				\
    yy_reduce_print (Rule);		\
} while (0)

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YYDSYMPRINT(Args)
# define YYDSYMPRINTF(Title, Token, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef	YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   SIZE_MAX < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#if YYMAXDEPTH == 0
# undef YYMAXDEPTH
#endif

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif



#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined (__GLIBC__) && defined (_STRING_H)
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
static YYSIZE_T
#   if defined (__STDC__) || defined (__cplusplus)
yystrlen (const char *yystr)
#   else
yystrlen (yystr)
     const char *yystr;
#   endif
{
  register const char *yys = yystr;

  while (*yys++ != '\0')
    continue;

  return yys - yystr - 1;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined (__GLIBC__) && defined (_STRING_H) && defined (_GNU_SOURCE)
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
static char *
#   if defined (__STDC__) || defined (__cplusplus)
yystpcpy (char *yydest, const char *yysrc)
#   else
yystpcpy (yydest, yysrc)
     char *yydest;
     const char *yysrc;
#   endif
{
  register char *yyd = yydest;
  register const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

#endif /* !YYERROR_VERBOSE */



#if YYDEBUG
/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

#if defined (__STDC__) || defined (__cplusplus)
static void
yysymprint (FILE *yyoutput, int yytype, YYSTYPE *yyvaluep)
#else
static void
yysymprint (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE *yyvaluep;
#endif
{
  /* Pacify ``unused variable'' warnings.  */
  (void) yyvaluep;

  if (yytype < YYNTOKENS)
    {
      YYFPRINTF (yyoutput, "token %s (", yytname[yytype]);
# ifdef YYPRINT
      YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# endif
    }
  else
    YYFPRINTF (yyoutput, "nterm %s (", yytname[yytype]);

  switch (yytype)
    {
      default:
        break;
    }
  YYFPRINTF (yyoutput, ")");
}

#endif /* ! YYDEBUG */
/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

#if defined (__STDC__) || defined (__cplusplus)
static void
yydestruct (int yytype, YYSTYPE *yyvaluep)
#else
static void
yydestruct (yytype, yyvaluep)
    int yytype;
    YYSTYPE *yyvaluep;
#endif
{
  /* Pacify ``unused variable'' warnings.  */
  (void) yyvaluep;

  switch (yytype)
    {

      default:
        break;
    }
}


/* Prevent warnings from -Wmissing-prototypes.  */

#ifdef YYPARSE_PARAM
# if defined (__STDC__) || defined (__cplusplus)
int yyparse (void *YYPARSE_PARAM);
# else
int yyparse ();
# endif
#else /* ! YYPARSE_PARAM */
#if defined (__STDC__) || defined (__cplusplus)
int yyparse (void);
#else
int yyparse ();
#endif
#endif /* ! YYPARSE_PARAM */



/* The lookahead symbol.  */
int yychar;

/* The semantic value of the lookahead symbol.  */
YYSTYPE yylval;

/* Number of syntax errors so far.  */
int yynerrs;



/*----------.
| yyparse.  |
`----------*/

#ifdef YYPARSE_PARAM
# if defined (__STDC__) || defined (__cplusplus)
int yyparse (void *YYPARSE_PARAM)
# else
int yyparse (YYPARSE_PARAM)
  void *YYPARSE_PARAM;
# endif
#else /* ! YYPARSE_PARAM */
#if defined (__STDC__) || defined (__cplusplus)
int
yyparse (void)
#else
int
yyparse ()

#endif
#endif
{
  
  register int yystate;
  register int yyn;
  int yyresult;
  /* Number of tokens to shift before error messages enabled.  */
  int yyerrstatus;
  /* Lookahead token as an internal (translated) token number.  */
  int yytoken = 0;

  /* Three stacks and their tools:
     `yyss': related to states,
     `yyvs': related to semantic values,
     `yyls': related to locations.

     Refer to the stacks thru separate pointers, to allow yyoverflow
     to reallocate them elsewhere.  */

  /* The state stack.  */
  short	yyssa[YYINITDEPTH];
  short *yyss = yyssa;
  register short *yyssp;

  /* The semantic value stack.  */
  YYSTYPE yyvsa[YYINITDEPTH];
  YYSTYPE *yyvs = yyvsa;
  register YYSTYPE *yyvsp;



#define YYPOPSTACK   (yyvsp--, yyssp--)

  YYSIZE_T yystacksize = YYINITDEPTH;

  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;


  /* When reducing, the number of symbols on the RHS of the reduced
     rule.  */
  int yylen;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY;		/* Cause a token to be read.  */

  /* Initialize stack pointers.
     Waste one element of value and location stack
     so that they stay on the same level as the state stack.
     The wasted elements are never initialized.  */

  yyssp = yyss;
  yyvsp = yyvs;

  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed. so pushing a state here evens the stacks.
     */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
	/* Give user a chance to reallocate the stack. Use copies of
	   these so that the &'s don't force the real ones into
	   memory.  */
	YYSTYPE *yyvs1 = yyvs;
	short *yyss1 = yyss;


	/* Each stack pointer address is followed by the size of the
	   data in use in that stack, in bytes.  This used to be a
	   conditional around just the two extra args, but that might
	   be undefined if yyoverflow is a macro.  */
	yyoverflow ("parser stack overflow",
		    &yyss1, yysize * sizeof (*yyssp),
		    &yyvs1, yysize * sizeof (*yyvsp),

		    &yystacksize);

	yyss = yyss1;
	yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyoverflowlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
	goto yyoverflowlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
	yystacksize = YYMAXDEPTH;

      {
	short *yyss1 = yyss;
	union yyalloc *yyptr =
	  (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
	if (! yyptr)
	  goto yyoverflowlab;
	YYSTACK_RELOCATE (yyss);
	YYSTACK_RELOCATE (yyvs);

#  undef YYSTACK_RELOCATE
	if (yyss1 != yyssa)
	  YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;


      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
		  (unsigned long int) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
	YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

/* Do appropriate processing given the current state.  */
/* Read a lookahead token if we need one and don't already have one.  */
/* yyresume: */

  /* First try to decide what to do without reference to lookahead token.  */

  yyn = yypact[yystate];
  if (yyn == YYPACT_NINF)
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid lookahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = YYLEX;
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YYDSYMPRINTF ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yyn == 0 || yyn == YYTABLE_NINF)
	goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  if (yyn == YYFINAL)
    YYACCEPT;

  /* Shift the lookahead token.  */
  YYDPRINTF ((stderr, "Shifting token %s, ", yytname[yytoken]));

  /* Discard the token being shifted unless it is eof.  */
  if (yychar != YYEOF)
    yychar = YYEMPTY;

  *++yyvsp = yylval;


  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  yystate = yyn;
  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     `$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
        case 4:
#line 878 "glb_parser.y"
    {YYABORT;}
    break;

  case 5:
#line 879 "glb_parser.y"
    {YYABORT;}
    break;

  case 8:
#line 884 "glb_parser.y"
    { /*showlist($1[0]);showlist($1[1]);*/ }
    break;

  case 9:
#line 885 "glb_parser.y"
    { /*printf ("\t%.10g\n", $1);*/  }
    break;

  case 10:
#line 886 "glb_parser.y"
    { /*showlist($1);printf("\n");*/ }
    break;

  case 11:
#line 887 "glb_parser.y"
    {}
    break;

  case 12:
#line 888 "glb_parser.y"
    {yyerrok;}
    break;

  case 13:
#line 889 "glb_parser.y"
    {}
    break;

  case 14:
#line 890 "glb_parser.y"
    {
 set_exp_energy("@energy",yyvsp[-2].ptrq);
}
    break;

  case 15:
#line 893 "glb_parser.y"
    {}
    break;

  case 16:
#line 894 "glb_parser.y"
    {/* printf("%s\n",$1);*/}
    break;

  case 17:
#line 895 "glb_parser.y"
    {/* printf("%s\n",$1);*/}
    break;

  case 18:
#line 900 "glb_parser.y"
    { yyval.val = yyvsp[0].val;                         }
    break;

  case 19:
#line 901 "glb_parser.y"
    { yyval.val = yyvsp[0].nameptr->value;              }
    break;

  case 20:
#line 902 "glb_parser.y"
    { yyval.val = yyvsp[0].tptr->value.var;              }
    break;

  case 21:
#line 903 "glb_parser.y"
    { yyval.val = yyvsp[0].val; yyvsp[-2].tptr->value.var = yyvsp[0].val;     }
    break;

  case 22:
#line 904 "glb_parser.y"
    { if(set_exp(yyvsp[-2].name,yyvsp[0].val,0)==1) yyerror("Unknown identifier");
yyval.val = yyvsp[0].val; }
    break;

  case 23:
#line 906 "glb_parser.y"
    { if(set_fnct(yyvsp[-2].name,yyvsp[0].nameptr->sf)==1) yyerror("Unknown identifier");}
    break;

  case 24:
#line 907 "glb_parser.y"
    { if(set_pair(yyvsp[-4].name,yyvsp[-2].val,yyvsp[0].val,0)==1) yyerror("Unknown identifier");
yyval.val = yyvsp[-2].val; }
    break;

  case 25:
#line 909 "glb_parser.y"
    { yyval.val = (*(yyvsp[-3].tptr->value.fnctptr))(yyvsp[-1].val); }
    break;

  case 26:
#line 910 "glb_parser.y"
    { yyval.val = yyvsp[-2].val + yyvsp[0].val;                    }
    break;

  case 27:
#line 911 "glb_parser.y"
    { yyval.val = yyvsp[-2].val - yyvsp[0].val;                    }
    break;

  case 28:
#line 912 "glb_parser.y"
    { yyval.val = yyvsp[-2].val * yyvsp[0].val;                    }
    break;

  case 29:
#line 913 "glb_parser.y"
    { yyval.val = yyvsp[-2].val / yyvsp[0].val;                    }
    break;

  case 30:
#line 914 "glb_parser.y"
    { yyval.val = -yyvsp[0].val;                        }
    break;

  case 31:
#line 915 "glb_parser.y"
    { yyval.val = pow (yyvsp[-2].val, yyvsp[0].val);               }
    break;

  case 32:
#line 916 "glb_parser.y"
    { yyval.val = yyvsp[-1].val;                         }
    break;

  case 33:
#line 917 "glb_parser.y"
    { yyval.val = 0;}
    break;

  case 34:
#line 918 "glb_parser.y"
    {yyerror("Unknown name");YYERROR;}
    break;

  case 35:
#line 924 "glb_parser.y"
    {yyval.ptr=list_cons(NULL,yyvsp[-1].val); }
    break;

  case 36:
#line 925 "glb_parser.y"
    {yyval.ptr=list_cons(NULL,yyvsp[-2].val); }
    break;

  case 37:
#line 926 "glb_parser.y"
    {
  glb_List *buf;
  buf=yyvsp[-2].ptr;
  buf=list_cons(buf,yyvsp[-1].val);
  yyval.ptr=buf;   
}
    break;

  case 38:
#line 933 "glb_parser.y"
    {
  glb_List *buf;
  buf=yyvsp[-3].ptr;
  buf=list_cons(buf,yyvsp[-2].val);
  yyval.ptr=buf;   
}
    break;

  case 39:
#line 939 "glb_parser.y"
    {
  glb_List *buf;
  buf=yyvsp[-1].ptr;
  buf=list_cons(buf,yyvsp[0].val);
  yyval.ptr=buf;   
}
    break;

  case 40:
#line 945 "glb_parser.y"
    {yyval.ptr=yyvsp[-1].ptr; }
    break;

  case 41:
#line 946 "glb_parser.y"
    {yyval.ptr=list_cons(NULL,yyvsp[-1].val); }
    break;

  case 42:
#line 947 "glb_parser.y"
    { if(set_exp_list(yyvsp[-2].name,yyvsp[0].ptr,3)==1) 
  yyerror("Unknown identifier");
 yyval.ptr = yyvsp[0].ptr; }
    break;

  case 43:
#line 952 "glb_parser.y"
    { 

  double *buf;
  buf=(double*) glb_malloc(sizeof(double)*2);
  buf[0]=yyvsp[-2].val;
  buf[1]=yyvsp[0].val-1;
  yyval.dpt=buf;
}
    break;

  case 44:
#line 962 "glb_parser.y"
    {}
    break;

  case 45:
#line 963 "glb_parser.y"
    {}
    break;

  case 46:
#line 967 "glb_parser.y"
    { if(yyvsp[-1].nameptr->value==-1) {yyvsp[-1].nameptr->value=step_counter(yyvsp[-3].name); }
  loc_count=yyvsp[-1].nameptr->value;
  glb_free(context);
  context =(char *) strdup(yyvsp[-3].name);
}
    break;

  case 47:
#line 972 "glb_parser.y"
    { 
   
grp_end(context);}
    break;

  case 48:
#line 975 "glb_parser.y"
    {
  yyerror("Redefinition of an automatic variable");YYERROR;}
    break;

  case 50:
#line 980 "glb_parser.y"
    {}
    break;

  case 51:
#line 981 "glb_parser.y"
    {}
    break;

  case 52:
#line 982 "glb_parser.y"
    {}
    break;

  case 53:
#line 983 "glb_parser.y"
    {}
    break;

  case 54:
#line 984 "glb_parser.y"
    {}
    break;

  case 55:
#line 985 "glb_parser.y"
    {}
    break;

  case 56:
#line 988 "glb_parser.y"
    {
  buff.version=strdup(yyvsp[0].name);
}
    break;

  case 57:
#line 992 "glb_parser.y"
    {
  //load_cross($3,loc_count-1);
  xsc.file_name=strdup(yyvsp[0].name);
  yyval.name=yyvsp[0].name;
}
    break;

  case 58:
#line 999 "glb_parser.y"
    {
  //load_flux($3,loc_count-1,1);
  flt.file_name=strdup(yyvsp[0].name);
  //if(set_exp($1,$3,0)==1) yyerror("Unknown identifier");
  yyval.name=yyvsp[0].name;
}
    break;

  case 59:
#line 1009 "glb_parser.y"
    {

  int x[6];
  x[0]=yyvsp[-10].nameptr->value - 1 ;
  x[1]=yyvsp[-8].in;
  x[2]=yyvsp[-6].in;
  x[3]=yyvsp[-4].in;
  x[4]=(int) yyvsp[-2].nameptr->value -1; 
  x[5]=(int) yyvsp[0].nameptr->value -1;


  set_channel_data(x,loc_count);



}
    break;

  case 60:
#line 1028 "glb_parser.y"
    {yyval.nameptr=yyvsp[0].nameptr;}
    break;

  case 61:
#line 1029 "glb_parser.y"
    {yyerror("Unknown name");YYERROR;}
    break;

  case 62:
#line 1033 "glb_parser.y"
    {yyval.in=yyvsp[0].in;}
    break;

  case 63:
#line 1034 "glb_parser.y"
    {yyval.in=1;}
    break;

  case 64:
#line 1035 "glb_parser.y"
    {yyval.in=-1;}
    break;

  case 65:
#line 1039 "glb_parser.y"
    {
  glb_List **buf;
  energy_len=1;
 
  buf=(glb_List**) glb_malloc(sizeof( glb_List* ) ); 
  buf[0]=yyvsp[0].ptr;
 
  yyval.ptrq=buf;
}
    break;

  case 66:
#line 1049 "glb_parser.y"
    {
  glb_List **buf;
  buf=yyvsp[-2].ptrq;
  energy_len++;
  
  buf=(glb_List**) glb_realloc((void**) buf , sizeof( glb_List* ) * energy_len);
 
  buf[energy_len-1]=yyvsp[0].ptr; 
  yyval.ptrq=buf; }
    break;

  case 67:
#line 1059 "glb_parser.y"
    {
  glb_List **buf;
  buf=yyvsp[-3].ptrq;
  energy_len++;

  buf=(glb_List**) glb_realloc((void**) buf , sizeof( glb_List* ) * energy_len);

  buf[energy_len-1]=yyvsp[0].ptr; 
  yyval.ptrq=buf; }
    break;

  case 68:
#line 1072 "glb_parser.y"
    {
  glb_List **buf;
  buf=(glb_List**) glb_malloc(sizeof(glb_List*)*2);
  buf[0]=list_cons(NULL,yyvsp[0].dpt[0]);
  buf[1]=list_cons(NULL,yyvsp[0].dpt[1]);
  glb_free(yyvsp[0].dpt);
  yyval.ptrq=buf;
}
    break;

  case 69:
#line 1080 "glb_parser.y"
    {
  glb_List **buf;
  buf=yyvsp[-2].ptrq;
  buf[0]=list_cons(buf[0],yyvsp[0].dpt[0]);
  buf[1]=list_cons(buf[1],yyvsp[0].dpt[1]);
  glb_free(yyvsp[0].dpt);
  yyval.ptrq=buf; 
}
    break;

  case 70:
#line 1090 "glb_parser.y"
    {
  glb_List **buf;  
 
  buf=(glb_List**) glb_malloc(sizeof(glb_List*)*2);
  buf[0]=list_cons(NULL,yyvsp[0].dpt[0]);
  buf[1]=list_cons(NULL,yyvsp[0].dpt[1]);
  glb_free(yyvsp[0].dpt);
  yyval.ptrq=buf;
}
    break;

  case 71:
#line 1099 "glb_parser.y"
    {
  glb_List **buf;
  buf=yyvsp[-2].ptrq;
  buf[0]=list_cons(buf[0],yyvsp[0].dpt[0]);
  buf[1]=list_cons(buf[1],yyvsp[0].dpt[1]);
  glb_free(yyvsp[0].dpt);
  yyval.ptrq=buf; 
}
    break;

  case 72:
#line 1109 "glb_parser.y"
    {

int flag;  
  
  yyval.ptrq=yyvsp[0].ptrq;
  
  flag=set_exp_list("bgrulescoeff",yyvsp[0].ptrq[0],0);
  if(flag==1) yyerror("Unknown identifier");
  flag=set_exp_list("bgrulechannellist",yyvsp[0].ptrq[1],0);
  if(flag==1) yyerror("Unknown identifier");
 
  glb_free(yyvsp[0].ptrq);
}
    break;

  case 73:
#line 1122 "glb_parser.y"
    {
  int flag;  
  yyval.ptrq=yyvsp[0].ptrq;
  flag=set_exp_list("rulescoeff",yyvsp[0].ptrq[0],0);
  if(flag==1) yyerror("Unknown identifier");
  flag=set_exp_list("rulechannellist",yyvsp[0].ptrq[1],0);
  if(flag==1) yyerror("Unknown identifier"); 

  glb_free(yyvsp[0].ptrq);

 
}
    break;

  case 74:
#line 1136 "glb_parser.y"
    {
 if(yyvsp[0].nameptr->value==-1) {yyvsp[0].nameptr->value=step_counter("rule");printf("loc step\n"); }
  loc_count=yyvsp[0].nameptr->value;
  printf("loc_count %d \n",loc_count);
  YYERROR;}
    break;

  case 75:
#line 1141 "glb_parser.y"
    { YYERROR;}
    break;

  case 76:
#line 1145 "glb_parser.y"
    {
 if(yyvsp[0].nameptr->value==-1) {yyvsp[0].nameptr->value=step_counter("channel");printf("loc step\n"); }
  loc_count=yyvsp[0].nameptr->value;
  printf("cha loc_count %d \n",loc_count); YYERROR;}
    break;

  case 77:
#line 1149 "glb_parser.y"
    { YYERROR;}
    break;


    }

/* Line 991 of yacc.c.  */
#line 2508 "glb_parser.c"

  yyvsp -= yylen;
  yyssp -= yylen;


  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;


  /* Now `shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*------------------------------------.
| yyerrlab -- here on detecting error |
`------------------------------------*/
yyerrlab:
  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if YYERROR_VERBOSE
      yyn = yypact[yystate];

      if (YYPACT_NINF < yyn && yyn < YYLAST)
	{
	  YYSIZE_T yysize = 0;
	  int yytype = YYTRANSLATE (yychar);
	  char *yymsg;
	  int yyx, yycount;

	  yycount = 0;
	  /* Start YYX at -YYN if negative to avoid negative indexes in
	     YYCHECK.  */
	  for (yyx = yyn < 0 ? -yyn : 0;
	       yyx < (int) (sizeof (yytname) / sizeof (char *)); yyx++)
	    if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR)
	      yysize += yystrlen (yytname[yyx]) + 15, yycount++;
	  yysize += yystrlen ("syntax error, unexpected ") + 1;
	  yysize += yystrlen (yytname[yytype]);
	  yymsg = (char *) YYSTACK_ALLOC (yysize);
	  if (yymsg != 0)
	    {
	      char *yyp = yystpcpy (yymsg, "syntax error, unexpected ");
	      yyp = yystpcpy (yyp, yytname[yytype]);

	      if (yycount < 5)
		{
		  yycount = 0;
		  for (yyx = yyn < 0 ? -yyn : 0;
		       yyx < (int) (sizeof (yytname) / sizeof (char *));
		       yyx++)
		    if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR)
		      {
			const char *yyq = ! yycount ? ", expecting " : " or ";
			yyp = yystpcpy (yyp, yyq);
			yyp = yystpcpy (yyp, yytname[yyx]);
			yycount++;
		      }
		}
	      yyerror (yymsg);
	      YYSTACK_FREE (yymsg);
	    }
	  else
	    yyerror ("syntax error; also virtual memory exhausted");
	}
      else
#endif /* YYERROR_VERBOSE */
	yyerror ("syntax error");
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
	 error, discard it.  */

      /* Return failure if at end of input.  */
      if (yychar == YYEOF)
        {
	  /* Pop the error token.  */
          YYPOPSTACK;
	  /* Pop the rest of the stack.  */
	  while (yyss < yyssp)
	    {
	      YYDSYMPRINTF ("Error: popping", yystos[*yyssp], yyvsp, yylsp);
	      yydestruct (yystos[*yyssp], yyvsp);
	      YYPOPSTACK;
	    }
	  YYABORT;
        }

      YYDSYMPRINTF ("Error: discarding", yytoken, &yylval, &yylloc);
      yydestruct (yytoken, &yylval);
      yychar = YYEMPTY;

    }

  /* Else will try to reuse lookahead token after shifting the error
     token.  */
  goto yyerrlab2;


/*----------------------------------------------------.
| yyerrlab1 -- error raised explicitly by an action.  |
`----------------------------------------------------*/
yyerrlab1:

  /* Suppress GCC warning that yyerrlab1 is unused when no action
     invokes YYERROR.  */
#if defined (__GNUC_MINOR__) && 2093 <= (__GNUC__ * 1000 + __GNUC_MINOR__) \
    && !defined __cplusplus
  __attribute__ ((__unused__))
#endif


  goto yyerrlab2;


/*---------------------------------------------------------------.
| yyerrlab2 -- pop states until the error token can be shifted.  |
`---------------------------------------------------------------*/
yyerrlab2:
  yyerrstatus = 3;	/* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (yyn != YYPACT_NINF)
	{
	  yyn += YYTERROR;
	  if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
	    {
	      yyn = yytable[yyn];
	      if (0 < yyn)
		break;
	    }
	}

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
	YYABORT;

      YYDSYMPRINTF ("Error: popping", yystos[*yyssp], yyvsp, yylsp);
      yydestruct (yystos[yystate], yyvsp);
      yyvsp--;
      yystate = *--yyssp;

      YY_STACK_PRINT (yyss, yyssp);
    }

  if (yyn == YYFINAL)
    YYACCEPT;

  YYDPRINTF ((stderr, "Shifting error token, "));

  *++yyvsp = yylval;


  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;

#ifndef yyoverflow
/*----------------------------------------------.
| yyoverflowlab -- parser overflow comes here.  |
`----------------------------------------------*/
yyoverflowlab:
  yyerror ("parser stack overflow");
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
  return yyresult;
}


#line 1152 "glb_parser.y"


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
    {"tan",  tan},
    {"asin", asin},
    {"acos", acos},
    {"atan", atan},
    {"log",   log},
    {"exp",  exp},
    {"log10", log10},
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
/* cannot use static here, since its declared earlier as extern */
glb_symrec *sym_table = (glb_symrec *) NULL;
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

static glb_naming *glb_putnames (char *sym_name, char *context, int value,
				 glb_naming *in)
{
  glb_naming *ptr;
  ptr = (glb_naming *) glb_malloc (sizeof (glb_naming));
  ptr->name = (char *) glb_malloc (strlen (sym_name) + 1);
  strcpy (ptr->name,sym_name);
  ptr->context = (char *) glb_malloc (strlen (context) + 1);
  strcpy (ptr->context,context);
  ptr->value = value; /* set value to -1 for new ones  */
  ptr->next = (struct glb_naming *) in;
  //in = ptr;
  return ptr;
}

static glb_naming *copy_names (glb_naming *in)
{
  
  glb_namerec *ptr;
  for (ptr = name_table; ptr != (glb_namerec *) NULL;
       ptr = (glb_namerec *)ptr->next)
    {
      in=glb_putnames(ptr->name,ptr->context,ptr->value,in);
    }
  return in;
}



glb_namerec *glb_getname (const char *sym_name, char* context)
{
  glb_namerec *ptr;
  for (ptr = name_table; ptr != (glb_namerec *) NULL;
       ptr = (glb_namerec *)ptr->next)
    if (strcmp (ptr->name,sym_name) == 0 )
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
  /* I am not sure how well this assigment really works */
  buff.names=copy_names(buff.names);
  buff_list[exp_count]=buff; 
  exp_count++;
}

void glbReset() 
{
  glb_line_num=0;
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
  glb_line_num=0;
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
  struct glb_experiment **ins;
  ins=(struct glb_experiment **) in;
  // yydebug=1;
  context=(char *) strdup("global");
  glb_smear_reset(&ibf);
  glb_option_type_reset(&opt);
  glb_flux_reset(&flt);
  input=glb_fopen(inf,"r");
  if(input==NULL) return -2;
  /* This line produces a warning with -Wall:
   * 'warning: char format, different type arg (arg 3)'
   * I think however that understand the problem.
   */
  fscanf(input,"%8c",&tct);
  if(strncmp(tch,tct,8)!=0) {glb_error("Not a GLoBES file!"); return -2;}
  glb_fclose(input);
  yyin=glb_fopen(inf,"r");
  if(yyin==NULL) return -2;
  glb_file_id=(char*) strdup(inf);
  glbResetEOF();
  k=yyparse ();
  glb_fclose(yyin);
  glb_free(context);
  glb_free(glb_file_id);
  
  if(k!=0) return -2;
  
  k=0;
  if(*counter+exp_count>32) glb_fatal("Too many experiments!");
  for(i=0;i<exp_count;i++)
    {  
      *ins[*counter+i]=buff_list[i];
      k=+glbDefaultExp(ins[*counter+i]);
    }
  (*counter)= (*counter) + exp_count;


  if(k!=0) return -1;
  
  return 0;

}

