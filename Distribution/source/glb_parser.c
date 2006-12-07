/* A Bison parser, made by GNU Bison 2.1.  */

/* Skeleton parser for Yacc-like parsing with Bison,
   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004, 2005 Free Software Foundation, Inc.

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
   Foundation, Inc., 51 Franklin Street, Fifth Floor,
   Boston, MA 02110-1301, USA.  */

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

/* Bison version.  */
#define YYBISON_VERSION "2.1"

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
     LVAR = 260,
     VAR = 261,
     FNCT = 262,
     IDN = 263,
     CROSS = 264,
     FLUXP = 265,
     FLUXM = 266,
     SYS_ON_FUNCTION = 267,
     SYS_OFF_FUNCTION = 268,
     GRP = 269,
     GID = 270,
     FNAME = 271,
     VERS = 272,
     SIGNAL = 273,
     BG = 274,
     GRPOPEN = 275,
     GRPCLOSE = 276,
     PM = 277,
     FLAVOR = 278,
     NOGLOBES = 279,
     CHANNEL = 280,
     RULESEP = 281,
     RULEMULT = 282,
     ENERGY = 283,
     NAME = 284,
     RDF = 285,
     NDEF = 286,
     NEG = 287
   };
#endif
/* Tokens.  */
#define NUM 258
#define SFNCT 259
#define LVAR 260
#define VAR 261
#define FNCT 262
#define IDN 263
#define CROSS 264
#define FLUXP 265
#define FLUXM 266
#define SYS_ON_FUNCTION 267
#define SYS_OFF_FUNCTION 268
#define GRP 269
#define GID 270
#define FNAME 271
#define VERS 272
#define SIGNAL 273
#define BG 274
#define GRPOPEN 275
#define GRPCLOSE 276
#define PM 277
#define FLAVOR 278
#define NOGLOBES 279
#define CHANNEL 280
#define RULESEP 281
#define RULEMULT 282
#define ENERGY 283
#define NAME 284
#define RDF 285
#define NDEF 286
#define NEG 287




/* Copy the first part of user declarations.  */
#line 26 "glb_parser.y"

#define YYDEBUG 1

#include <math.h>
#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <globes/globes.h>
#include <gsl/gsl_spline.h>
#include "glb_smear.h"
#include "glb_multiex.h"
#include "glb_types.h"
#include "glb_error.h"
#include "glb_lexer.h"  
#include "glb_parser_type.h"
#include "glb_fluxes.h"
#include "glb_sys.h"



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
#define DOUBLE_PAIR 17
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
  static int errordim_sys_on=-1;         /* Temp. storage for the old numerical errordims */
  static int errordim_sys_off=-1;
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
/* JK - Has been removed
    {"$errorfunction",INT,0,1,&buff.errorfunction,NULL,"global"}, 
    {"@treshold_setttings" ,DOUBLE_INDEXED_PAIR,
     0,100,&buff.bgtcenter[0],&loc_count,"rule"},
    {"@treshold_error" ,DOUBLE_INDEXED_PAIR,
     0,100,&buff.bgterror[0],&loc_count,"rule"}, */
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
    
 
   {"@errordim_sys_on",INT,0,20,&errordim_sys_on,NULL,"rule"},
   {"@errordim_sys_off",INT,0,20,&errordim_sys_off,NULL,"rule"},
   {"@sys_on_function",CHAR,0,20,&buff.sys_on_strings[0],&loc_count,"rule"},
   {"@sys_off_function",CHAR,0,20,&buff.sys_off_strings[0],&loc_count,"rule"},

   {"channel",UNTYPE,-1,1,NULL,&buff.numofchannels,"global"},
   {"@channel",INT_LIST_INDEXED,-1,1,&buff.listofchannels[0],
    &buff.numofchannels,"channel"},
   {"energy",UNTYPE,-1,1,NULL,&buff.num_of_sm,"global"},
   {"@energy",ENERGY_MATRIX,-1,GMAX,&buff.smear[0],&loc_count,"energy"},

   {"@signalerror",DOUBLE_INDEXED_PAIR,
    0,100,&buff.signal_errors[0],&loc_count,"rule"},
   {"@backgrounderror",DOUBLE_INDEXED_PAIR,
    0,100,&buff.bg_errors[0],&loc_count,"rule"},
   {"@backgroundcenter",DOUBLE_INDEXED_PAIR,
    0,100,&buff.bg_centers[0],&loc_count,"rule"},

   {"@sys_on_errors",DOUBLE_LIST_INDEXED,0,GMAX,
    &buff.sys_on_errors[0],&loc_count,"rule"},
   {"@sys_off_errors",DOUBLE_LIST_INDEXED,0,GMAX,
    &buff.sys_off_errors[0],&loc_count,"rule"},
   

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


static void grp_start(char* name)
   {
     if(strncmp(name,"rule",4)==0 )
       {
         /* Reset variables to default values */
         errordim_sys_on  = -1.0;
         errordim_sys_off = -1.0; 
       }
   }

 
static void grp_end(char* name)
   {
     char tmp_errordim[2];    
     
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

     if(strncmp(name,"cross",5)==0 )
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

     if(strncmp(name,"rule",4)==0 )
       {
         int nr = buff.numofrules - 1;

         /* Parse old (numerical) errordims */
         if (buff.sys_on_strings[nr] == NULL  &&  errordim_sys_on >= 0)
           buff.sys_on_strings[nr] = glbConvertErrorDim(errordim_sys_on);
         if (buff.sys_off_strings[nr] == NULL  &&  errordim_sys_off >= 0)
           buff.sys_off_strings[nr] = glbConvertErrorDim(errordim_sys_off);
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
	     if(token_list[i].scalar==DOUBLE_PAIR)
	       {
		 if(value >= token_list[i].rl && value <= token_list[i].ru)
		   {
		    
		     dbf=(double*) token_list[i].ptr;
		     dbf[0]=(double) value;
		     dbf[1]=(double) value2;
		     return 0;
		   }
		 else       
		   {
		     fprintf(stderr,"Error: Value for %s out of range\n",
			     token_list[i].token);
		     return 2;
		   }
	       }
	     
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
		     dbf[(loc_count-1)*2+0]=(double) value;
		     dbf[(loc_count-1)*2+1]=(double) value2;
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


static size_t list_length (glb_List *head)
{
  size_t n;  
  for (n = 0; head; ++n)
    head = head->next;
  return n; 
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


static double glb_reverse(double x)
{
  return x;
}


static glb_List *list_cons (glb_List *tail, double newdata)
{
  glb_List *res = (glb_List*) glb_malloc(sizeof(glb_List)); 
  res->next=tail;
  res->entry=newdata;
  return res;
}

/* this functions threads a funct_t function (double f (double x))
   over a list */

static double glb_list_copy(double x)
{
  return x;
}

static glb_List *thread_list(func_t f, int reverse, int destroy ,glb_List *tail)
{
  glb_List *res, *head;
  size_t n,l;
  double nv;
  res=NULL;
  /* the problem is that threading reverses the list, hence we
     re-reverse it by default, unless someone wants to reverse himself
  */
  if(reverse==1)
    {
      head=tail;  
      for (n = 0; head; ++n)
	{
	  nv=f(head->entry);
	  res=list_cons(res,nv);
	  head = head->next;
	}
    }
  else
    {
      double *rlist,x;
      l=list_length(tail);
      rlist=(double *) malloc(sizeof(double)*l);   
      head=tail;  
      for (n = 0; head; ++n)
	{
	  rlist[n]=head->entry;
	  head = head->next;
	}

      for(n=l;n>0; n--)
	{
	  x=f(rlist[n-1]);
	  res=list_cons(res,x);
	}
      free(rlist);
    }
  if(destroy==1) list_free(tail);
  if(destroy==-1) {list_free(res);res=tail;}
  return res;
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


static glb_List *glb_interpolation(glb_List *xval,glb_List *yval,int flag,glb_List *where)
{
  glb_List *head,*res=NULL;
  size_t xl,yl,rl,i;
  double *xlist,*ylist,*rlist,x;
  gsl_interp_type type;

  if(flag==1) type=*gsl_interp_linear;
  else if (flag==2) type=*gsl_interp_cspline;
  else {glb_error("Invalid interpolation scheme flag");return NULL;}

  xl=list_length(xval);
  yl=list_length(yval);
  rl=list_length(where);
 
  if(yl!=xl) {glb_error("Xval and Yval in glb_interpolation are not of the same length"); return NULL;}

  xlist=(double*) malloc(sizeof(double)*xl);
  ylist=(double*) malloc(sizeof(double)*yl);  
  rlist=(double*) malloc(sizeof(double)*rl);

  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *spline = gsl_spline_alloc (&type,yl);

  head=xval;
  for (i = yl-1; head; i--){
    xlist[i]=head->entry;
    head = head->next;
  }
  
  head=yval;
  for (i = yl-1; head; i--){
    ylist[i]=head->entry;
    head = head->next;
  }

  gsl_spline_init(spline,xlist,ylist,yl);
  
  head=where;
  for(i=0; head; i++){
     rlist[i]=gsl_spline_eval(spline,head->entry,acc);   
    head=head->next;
 }

  for(i=rl;i>0;i--) res=list_cons(res,rlist[i-1]);
 
  free(xlist);
  free(ylist);
  free(rlist);
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);

  return res;
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

/* Enabling the token table.  */
#ifndef YYTOKEN_TABLE
# define YYTOKEN_TABLE 0
#endif

#if ! defined (YYSTYPE) && ! defined (YYSTYPE_IS_DECLARED)
#line 985 "glb_parser.y"
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
/* Line 196 of yacc.c.  */
#line 1121 "glb_parser.c"
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif



/* Copy the second part of user declarations.  */


/* Line 219 of yacc.c.  */
#line 1133 "glb_parser.c"

#if ! defined (YYSIZE_T) && defined (__SIZE_TYPE__)
# define YYSIZE_T __SIZE_TYPE__
#endif
#if ! defined (YYSIZE_T) && defined (size_t)
# define YYSIZE_T size_t
#endif
#if ! defined (YYSIZE_T) && (defined (__STDC__) || defined (__cplusplus))
# include <stddef.h> /* INFRINGES ON USER NAME SPACE */
# define YYSIZE_T size_t
#endif
#if ! defined (YYSIZE_T)
# define YYSIZE_T unsigned int
#endif

#ifndef YY_
# if YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(msgid) dgettext ("bison-runtime", msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(msgid) msgid
# endif
#endif

#if ! defined (yyoverflow) || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if defined (__STDC__) || defined (__cplusplus)
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#     define YYINCLUDED_STDLIB_H
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's `empty if-body' warning. */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (0)
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2005 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM ((YYSIZE_T) -1)
#  endif
#  ifdef __cplusplus
extern "C" {
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if (! defined (malloc) && ! defined (YYINCLUDED_STDLIB_H) \
	&& (defined (__STDC__) || defined (__cplusplus)))
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if (! defined (free) && ! defined (YYINCLUDED_STDLIB_H) \
	&& (defined (__STDC__) || defined (__cplusplus)))
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifdef __cplusplus
}
#  endif
# endif
#endif /* ! defined (yyoverflow) || YYERROR_VERBOSE */


#if (! defined (yyoverflow) \
     && (! defined (__cplusplus) \
	 || (defined (YYSTYPE_IS_TRIVIAL) && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  short int yyss;
  YYSTYPE yyvs;
  };

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (short int) + sizeof (YYSTYPE))			\
      + YYSTACK_GAP_MAXIMUM)

/* Copy COUNT objects from FROM to TO.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined (__GNUC__) && 1 < __GNUC__
#   define YYCOPY(To, From, Count) \
      __builtin_memcpy (To, From, (Count) * sizeof (*(From)))
#  else
#   define YYCOPY(To, From, Count)		\
      do					\
	{					\
	  YYSIZE_T yyi;				\
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
   typedef short int yysigned_char;
#endif

/* YYFINAL -- State number of the termination state. */
#define YYFINAL  4
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   477

/* YYNTOKENS -- Number of terminals. */
#define YYNTOKENS  52
/* YYNNTS -- Number of nonterminals. */
#define YYNNTS  24
/* YYNRULES -- Number of rules. */
#define YYNRULES  85
/* YYNRULES -- Number of states. */
#define YYNSTATES  179

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   287

#define YYTRANSLATE(YYX)						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const unsigned char yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
      44,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    50,     2,     2,     2,     2,     2,     2,
      46,    47,    39,    38,    32,    37,     2,    40,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,    45,
       2,    33,    34,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,    48,     2,    49,    42,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    36,    51,    35,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,    43,     2,     2,
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
      25,    26,    27,    28,    29,    30,    31,    41
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const unsigned short int yyprhs[] =
{
       0,     0,     3,     4,     7,     9,    11,    13,    16,    19,
      22,    25,    28,    31,    34,    38,    41,    44,    47,    49,
      51,    53,    57,    61,    65,    71,    76,    80,    84,    88,
      92,    95,    99,   103,   105,   107,   110,   114,   118,   123,
     126,   129,   133,   137,   141,   146,   150,   152,   156,   167,
     171,   173,   176,   177,   186,   194,   195,   197,   199,   201,
     203,   205,   207,   211,   215,   219,   223,   227,   241,   243,
     245,   247,   249,   251,   255,   259,   264,   268,   272,   276,
     280,   282,   284,   285,   289,   290
};

/* YYRHS -- A `-1'-separated list of the rules' RHS. */
static const yysigned_char yyrhs[] =
{
      53,     0,    -1,    -1,    53,    54,    -1,    24,    -1,    43,
      -1,    44,    -1,    74,    44,    -1,    72,    44,    -1,    55,
      44,    -1,    56,    44,    -1,    59,    44,    -1,     1,    44,
      -1,    66,    44,    -1,    69,    45,    44,    -1,    65,    44,
      -1,    63,    44,    -1,    64,    44,    -1,     3,    -1,    29,
      -1,     6,    -1,     6,    33,    55,    -1,     8,    33,    55,
      -1,     8,    33,     4,    -1,     8,    33,    55,    26,    55,
      -1,     7,    46,    55,    47,    -1,    55,    38,    55,    -1,
      55,    37,    55,    -1,    55,    39,    55,    -1,    55,    40,
      55,    -1,    37,    55,    -1,    55,    42,    55,    -1,    46,
      55,    47,    -1,    62,    -1,    31,    -1,    55,    32,    -1,
      55,    32,    44,    -1,    56,    55,    32,    -1,    56,    55,
      32,    44,    -1,    56,    55,    -1,    36,    35,    -1,    36,
      56,    35,    -1,    36,    55,    35,    -1,     8,    33,    56,
      -1,     7,    46,    56,    47,    -1,     7,    48,    49,    -1,
       5,    -1,     5,    33,    56,    -1,    50,    46,    56,    32,
      56,    32,    55,    51,    56,    47,    -1,    55,    27,    55,
      -1,    54,    -1,    58,    54,    -1,    -1,    15,    46,    29,
      47,    60,    20,    61,    21,    -1,    15,    46,    30,    47,
      20,    61,    21,    -1,    -1,    55,    -1,    65,    -1,    58,
      -1,    66,    -1,    63,    -1,    64,    -1,    17,    33,    16,
      -1,     9,    33,    16,    -1,    10,    33,    16,    -1,    12,
      33,    16,    -1,    13,    33,    16,    -1,    25,    33,    67,
      26,    68,    26,    23,    26,    23,    26,    67,    26,    67,
      -1,    29,    -1,    31,    -1,    22,    -1,    38,    -1,    37,
      -1,    28,    33,    56,    -1,    69,    26,    56,    -1,    69,
      26,    44,    56,    -1,    19,    33,    57,    -1,    70,    26,
      57,    -1,    18,    33,    57,    -1,    71,    26,    57,    -1,
      70,    -1,    71,    -1,    -1,    29,    73,    65,    -1,    -1,
      29,    75,    66,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const unsigned short int yyrline[] =
{
       0,  1033,  1033,  1034,  1035,  1036,  1039,  1040,  1041,  1042,
    1043,  1044,  1045,  1046,  1047,  1050,  1051,  1052,  1057,  1058,
    1059,  1060,  1061,  1063,  1064,  1066,  1068,  1069,  1070,  1071,
    1072,  1073,  1074,  1075,  1076,  1082,  1083,  1084,  1091,  1097,
    1103,  1104,  1105,  1106,  1110,  1111,  1112,  1113,  1114,  1117,
    1127,  1128,  1132,  1131,  1141,  1145,  1146,  1147,  1148,  1149,
    1150,  1151,  1154,  1159,  1166,  1174,  1179,  1185,  1202,  1203,
    1207,  1208,  1209,  1213,  1222,  1232,  1246,  1254,  1264,  1273,
    1283,  1296,  1310,  1310,  1319,  1319
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || YYTOKEN_TABLE
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals. */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "NUM", "SFNCT", "LVAR", "VAR", "FNCT",
  "IDN", "CROSS", "FLUXP", "FLUXM", "SYS_ON_FUNCTION", "SYS_OFF_FUNCTION",
  "GRP", "GID", "FNAME", "VERS", "SIGNAL", "BG", "GRPOPEN", "GRPCLOSE",
  "PM", "FLAVOR", "NOGLOBES", "CHANNEL", "RULESEP", "RULEMULT", "ENERGY",
  "NAME", "RDF", "NDEF", "','", "'='", "'>'", "'}'", "'{'", "'-'", "'+'",
  "'*'", "'/'", "NEG", "'^'", "'\\247'", "'\\n'", "';'", "'('", "')'",
  "'['", "']'", "'!'", "'|'", "$accept", "input", "line", "exp", "seq",
  "rulepart", "expseq", "group", "@1", "ingroup", "version", "cross",
  "flux", "rule", "channel", "name", "pm", "ene", "brule", "srule",
  "outrule", "@2", "outchannel", "@3", 0
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[YYLEX-NUM] -- Internal token number corresponding to
   token YYLEX-NUM.  */
static const unsigned short int yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283,   284,
     285,   286,    44,    61,    62,   125,   123,    45,    43,    42,
      47,   287,    94,   167,    10,    59,    40,    41,    91,    93,
      33,   124
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const unsigned char yyr1[] =
{
       0,    52,    53,    53,    53,    53,    54,    54,    54,    54,
      54,    54,    54,    54,    54,    54,    54,    54,    55,    55,
      55,    55,    55,    55,    55,    55,    55,    55,    55,    55,
      55,    55,    55,    55,    55,    56,    56,    56,    56,    56,
      56,    56,    56,    56,    56,    56,    56,    56,    56,    57,
      58,    58,    60,    59,    59,    61,    61,    61,    61,    61,
      61,    61,    62,    63,    64,    65,    65,    66,    67,    67,
      68,    68,    68,    69,    69,    69,    70,    70,    71,    71,
      65,    65,    73,    72,    75,    74
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const unsigned char yyr2[] =
{
       0,     2,     0,     2,     1,     1,     1,     2,     2,     2,
       2,     2,     2,     2,     3,     2,     2,     2,     1,     1,
       1,     3,     3,     3,     5,     4,     3,     3,     3,     3,
       2,     3,     3,     1,     1,     2,     3,     3,     4,     2,
       2,     3,     3,     3,     4,     3,     1,     3,    10,     3,
       1,     2,     0,     8,     7,     0,     1,     1,     1,     1,
       1,     1,     3,     3,     3,     3,     3,    13,     1,     1,
       1,     1,     1,     3,     3,     4,     3,     3,     3,     3,
       1,     1,     0,     3,     0,     3
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const unsigned char yydefact[] =
{
       2,     4,     5,     0,     1,     0,    18,    46,    20,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,    19,    34,     0,     0,     6,     0,     0,     3,     0,
       0,     0,    33,     0,     0,     0,     0,     0,    80,    81,
       0,     0,    12,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      19,    40,     0,     0,     0,     0,    30,     0,     0,    35,
       0,     0,     0,     0,     0,     9,    10,    39,    11,    16,
      17,    15,    13,     0,     0,     0,     0,     8,     7,     0,
      47,    21,     0,     0,    45,    23,    22,    43,    63,    64,
      65,    66,     0,     0,    62,     0,    78,    76,    68,    69,
       0,    73,    83,    85,    42,    41,     0,     0,    32,     0,
      36,    27,    26,    28,    29,    31,    37,     0,    74,    14,
      77,    79,    25,    44,     0,    52,     0,     0,     0,     0,
      22,     0,    38,    75,    24,     0,     0,    49,    70,    72,
      71,     0,     0,     0,    50,    56,     0,     0,    60,    61,
      57,    59,     0,     0,     0,    51,    54,     0,     0,    53,
       0,     0,     0,     0,     0,    48,     0,     0,    67
};

/* YYDEFGOTO[NTERM-NUM]. */
static const short int yydefgoto[] =
{
      -1,     3,   154,    77,    30,   106,   156,    31,   145,   157,
      32,    33,    34,    35,    36,   110,   151,    37,    38,    39,
      40,    58,    41,    59
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -160
static const short int yypact[] =
{
     -21,  -160,  -160,   164,  -160,   -25,  -160,     2,     6,   -39,
      14,    23,    33,    40,    45,    39,    46,    48,    53,    58,
      60,   224,  -160,    26,   422,  -160,   422,    50,  -160,   279,
      81,    65,  -160,    67,    71,    73,    75,   -18,    94,    95,
      79,    88,  -160,   307,   422,   307,    77,    99,   117,   121,
     130,   131,    29,   132,   422,   422,    19,   307,    82,   137,
    -160,  -160,   213,   355,   109,   142,   136,   407,   307,   140,
     422,   422,   422,   422,   422,  -160,  -160,   148,  -160,  -160,
    -160,  -160,  -160,   291,   150,   422,   422,  -160,  -160,   418,
     422,   435,   119,   323,  -160,  -160,   165,   422,  -160,  -160,
    -160,  -160,   116,   149,  -160,   102,  -160,  -160,  -160,  -160,
     171,   422,  -160,  -160,  -160,  -160,   422,   374,  -160,   390,
    -160,    35,    35,   136,   136,   136,   154,   307,   422,  -160,
    -160,  -160,  -160,  -160,   422,  -160,   179,   422,     8,   424,
     165,   307,  -160,   422,   435,   186,   210,   435,  -160,  -160,
    -160,   183,   406,   210,  -160,   279,   256,   191,    67,    71,
      73,    75,   198,   422,   203,  -160,  -160,   200,   -26,  -160,
     207,   307,   206,   342,    19,  -160,   208,    19,  -160
};

/* YYPGOTO[NTERM-NUM].  */
static const short int yypgoto[] =
{
    -160,  -160,    -2,    -3,   -19,   -49,  -160,  -160,  -160,    80,
    -160,  -136,   -93,   -56,   -54,  -159,  -160,  -160,  -160,  -160,
    -160,  -160,  -160,  -160
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -85
static const short int yytable[] =
{
      29,    28,   112,     1,    63,   113,   107,    45,    83,    46,
     158,    70,    71,    72,    73,   176,    74,   158,   178,    42,
      62,    66,     2,    67,    90,   171,    93,    84,    97,     6,
     148,     7,     8,     9,    10,    43,   130,   131,   111,    44,
      89,    91,    92,    16,    96,   149,   150,    47,   108,   119,
     109,   105,   105,   159,    89,    60,    48,    22,   102,   103,
     159,    61,    23,    24,   128,    89,    49,   121,   122,   123,
     124,   125,    26,    50,    72,    73,    27,    74,    51,    53,
      89,    54,   105,   105,     6,    52,    55,     8,    64,    65,
     160,    56,   161,    57,    13,    14,    68,   160,    16,   161,
      17,    18,     6,    95,     7,     8,     9,    10,   143,    78,
      60,    79,    22,   139,   140,    80,    16,    81,    24,    82,
      85,    86,   152,    87,    89,    76,    94,    26,    60,   137,
      22,   144,    88,    98,   147,    23,    24,    99,    89,    70,
      71,    72,    73,   155,    74,    26,   100,   101,   104,    27,
     155,    69,   173,    29,   165,   116,    70,    71,    72,    73,
     168,    74,    19,   135,     4,     5,   132,     6,    89,     7,
       8,     9,    10,    11,    12,   117,    13,    14,    74,    15,
     126,    16,    17,    18,   120,    70,    71,    72,    73,    19,
      74,   134,    20,    21,   129,    22,   136,   138,   142,   146,
      23,    24,    70,    71,    72,    73,   153,    74,    25,   162,
      26,     5,   166,     6,    27,     7,     8,     9,    10,    11,
      12,   167,    13,    14,   169,    15,   170,    16,    17,    18,
     172,   -55,   174,   164,   177,    19,   -82,   -82,    20,    21,
       0,    22,   -82,   -82,     0,    69,    23,    24,   114,   -84,
      70,    71,    72,    73,    25,    74,    26,     5,     0,     6,
      27,     7,     8,     9,    10,    11,    12,     0,    13,    14,
       0,    15,     0,    16,    17,    18,     0,   -58,     0,     0,
       0,    19,     0,     0,    20,    21,     0,    22,     0,     0,
       0,     0,    23,    24,     6,     0,     7,     8,     9,    10,
      25,     0,    26,     0,     0,     0,    27,     0,    16,     0,
       6,    69,     7,     8,     9,    10,    70,    71,    72,    73,
      60,    74,    22,    75,    16,     0,     6,    23,    24,     8,
      64,    65,     0,     0,     0,   127,    60,    26,    22,     0,
      16,    27,     0,    23,    24,     6,     0,     0,     8,    64,
      65,     0,    60,    26,    22,     0,     0,    27,     6,    16,
      24,     8,    64,    65,     0,     0,     0,     0,     0,    26,
     133,    60,    16,    22,     0,     0,     0,     6,    95,    24,
       8,    64,    65,     0,    60,     0,    22,     0,    26,   175,
     115,    16,    24,     6,     0,     0,     8,    64,    65,     0,
       0,    26,     0,    60,     0,    22,     0,    16,     0,     6,
       0,    24,     8,    64,    65,     0,     0,     0,     0,    60,
      26,    22,   141,    16,     0,     6,     0,    24,     8,    64,
      65,     0,     0,     0,     0,    60,    26,    22,   163,    16,
       0,     0,     0,    24,    70,    71,    72,    73,     0,    74,
      69,    60,    26,    22,   118,    70,    71,    72,    73,    24,
      74,    70,    71,    72,    73,     0,    74,     0,    26,     0,
       0,   132,    70,    71,    72,    73,     0,    74
};

static const short int yycheck[] =
{
       3,     3,    58,    24,    23,    59,    55,    46,    26,    48,
     146,    37,    38,    39,    40,   174,    42,   153,   177,    44,
      23,    24,    43,    26,    43,    51,    45,    45,    47,     3,
      22,     5,     6,     7,     8,    33,    85,    86,    57,    33,
      43,    44,    45,    17,    47,    37,    38,    33,    29,    68,
      31,    54,    55,   146,    57,    29,    33,    31,    29,    30,
     153,    35,    36,    37,    83,    68,    33,    70,    71,    72,
      73,    74,    46,    33,    39,    40,    50,    42,    33,    33,
      83,    33,    85,    86,     3,    46,    33,     6,     7,     8,
     146,    33,   146,    33,    12,    13,    46,   153,    17,   153,
      18,    19,     3,     4,     5,     6,     7,     8,   127,    44,
      29,    44,    31,   116,   117,    44,    17,    44,    37,    44,
      26,    26,   141,    44,   127,    44,    49,    46,    29,    27,
      31,   134,    44,    16,   137,    36,    37,    16,   141,    37,
      38,    39,    40,   146,    42,    46,    16,    16,    16,    50,
     153,    32,   171,   156,   156,    46,    37,    38,    39,    40,
     163,    42,    25,    47,     0,     1,    47,     3,   171,     5,
       6,     7,     8,     9,    10,    33,    12,    13,    42,    15,
      32,    17,    18,    19,    44,    37,    38,    39,    40,    25,
      42,    26,    28,    29,    44,    31,    47,    26,    44,    20,
      36,    37,    37,    38,    39,    40,    20,    42,    44,    26,
      46,     1,    21,     3,    50,     5,     6,     7,     8,     9,
      10,    23,    12,    13,    21,    15,    26,    17,    18,    19,
      23,    21,    26,   153,    26,    25,    12,    13,    28,    29,
      -1,    31,    18,    19,    -1,    32,    36,    37,    35,    25,
      37,    38,    39,    40,    44,    42,    46,     1,    -1,     3,
      50,     5,     6,     7,     8,     9,    10,    -1,    12,    13,
      -1,    15,    -1,    17,    18,    19,    -1,    21,    -1,    -1,
      -1,    25,    -1,    -1,    28,    29,    -1,    31,    -1,    -1,
      -1,    -1,    36,    37,     3,    -1,     5,     6,     7,     8,
      44,    -1,    46,    -1,    -1,    -1,    50,    -1,    17,    -1,
       3,    32,     5,     6,     7,     8,    37,    38,    39,    40,
      29,    42,    31,    44,    17,    -1,     3,    36,    37,     6,
       7,     8,    -1,    -1,    -1,    44,    29,    46,    31,    -1,
      17,    50,    -1,    36,    37,     3,    -1,    -1,     6,     7,
       8,    -1,    29,    46,    31,    -1,    -1,    50,     3,    17,
      37,     6,     7,     8,    -1,    -1,    -1,    -1,    -1,    46,
      47,    29,    17,    31,    -1,    -1,    -1,     3,     4,    37,
       6,     7,     8,    -1,    29,    -1,    31,    -1,    46,    47,
      35,    17,    37,     3,    -1,    -1,     6,     7,     8,    -1,
      -1,    46,    -1,    29,    -1,    31,    -1,    17,    -1,     3,
      -1,    37,     6,     7,     8,    -1,    -1,    -1,    -1,    29,
      46,    31,    32,    17,    -1,     3,    -1,    37,     6,     7,
       8,    -1,    -1,    -1,    -1,    29,    46,    31,    32,    17,
      -1,    -1,    -1,    37,    37,    38,    39,    40,    -1,    42,
      32,    29,    46,    31,    47,    37,    38,    39,    40,    37,
      42,    37,    38,    39,    40,    -1,    42,    -1,    46,    -1,
      -1,    47,    37,    38,    39,    40,    -1,    42
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const unsigned char yystos[] =
{
       0,    24,    43,    53,     0,     1,     3,     5,     6,     7,
       8,     9,    10,    12,    13,    15,    17,    18,    19,    25,
      28,    29,    31,    36,    37,    44,    46,    50,    54,    55,
      56,    59,    62,    63,    64,    65,    66,    69,    70,    71,
      72,    74,    44,    33,    33,    46,    48,    33,    33,    33,
      33,    33,    46,    33,    33,    33,    33,    33,    73,    75,
      29,    35,    55,    56,     7,     8,    55,    55,    46,    32,
      37,    38,    39,    40,    42,    44,    44,    55,    44,    44,
      44,    44,    44,    26,    45,    26,    26,    44,    44,    55,
      56,    55,    55,    56,    49,     4,    55,    56,    16,    16,
      16,    16,    29,    30,    16,    55,    57,    57,    29,    31,
      67,    56,    65,    66,    35,    35,    46,    33,    47,    56,
      44,    55,    55,    55,    55,    55,    32,    44,    56,    44,
      57,    57,    47,    47,    26,    47,    47,    27,    26,    55,
      55,    32,    44,    56,    55,    60,    20,    55,    22,    37,
      38,    68,    56,    20,    54,    55,    58,    61,    63,    64,
      65,    66,    26,    32,    61,    54,    21,    23,    55,    21,
      26,    51,    23,    56,    26,    47,    67,    26,    67
};

#define yyerrok		(yyerrstatus = 0)
#define yyclearin	(yychar = YYEMPTY)
#define YYEMPTY		(-2)
#define YYEOF		0

#define YYACCEPT	goto yyacceptlab
#define YYABORT		goto yyabortlab
#define YYERROR		goto yyerrorlab


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
    {								\
      yyerror (YY_("syntax error: cannot back up")); \
      YYERROR;							\
    }								\
while (0)


#define YYTERROR	1
#define YYERRCODE	256


/* YYLLOC_DEFAULT -- Set CURRENT to span from RHS[1] to RHS[N].
   If N is 0, then set CURRENT to the empty location which ends
   the previous symbol: RHS[0] (always defined).  */

#define YYRHSLOC(Rhs, K) ((Rhs)[K])
#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)				\
    do									\
      if (N)								\
	{								\
	  (Current).first_line   = YYRHSLOC (Rhs, 1).first_line;	\
	  (Current).first_column = YYRHSLOC (Rhs, 1).first_column;	\
	  (Current).last_line    = YYRHSLOC (Rhs, N).last_line;		\
	  (Current).last_column  = YYRHSLOC (Rhs, N).last_column;	\
	}								\
      else								\
	{								\
	  (Current).first_line   = (Current).last_line   =		\
	    YYRHSLOC (Rhs, 0).last_line;				\
	  (Current).first_column = (Current).last_column =		\
	    YYRHSLOC (Rhs, 0).last_column;				\
	}								\
    while (0)
#endif


/* YY_LOCATION_PRINT -- Print the location on the stream.
   This macro was not mandated originally: define only if we know
   we won't break user code: when these are the locations we know.  */

#ifndef YY_LOCATION_PRINT
# if YYLTYPE_IS_TRIVIAL
#  define YY_LOCATION_PRINT(File, Loc)			\
     fprintf (File, "%d.%d-%d.%d",			\
              (Loc).first_line, (Loc).first_column,	\
              (Loc).last_line,  (Loc).last_column)
# else
#  define YY_LOCATION_PRINT(File, Loc) ((void) 0)
# endif
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

# define YY_SYMBOL_PRINT(Title, Type, Value, Location)		\
do {								\
  if (yydebug)							\
    {								\
      YYFPRINTF (stderr, "%s ", Title);				\
      yysymprint (stderr,					\
                  Type, Value);	\
      YYFPRINTF (stderr, "\n");					\
    }								\
} while (0)

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

#if defined (__STDC__) || defined (__cplusplus)
static void
yy_stack_print (short int *bottom, short int *top)
#else
static void
yy_stack_print (bottom, top)
    short int *bottom;
    short int *top;
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
  unsigned long int yylno = yyrline[yyrule];
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu), ",
             yyrule - 1, yylno);
  /* Print the symbols being reduced, and their result.  */
  for (yyi = yyprhs[yyrule]; 0 <= yyrhs[yyi]; yyi++)
    YYFPRINTF (stderr, "%s ", yytname[yyrhs[yyi]]);
  YYFPRINTF (stderr, "-> %s\n", yytname[yyr1[yyrule]]);
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
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
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
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

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
  const char *yys = yystr;

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
  char *yyd = yydest;
  const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

# ifndef yytnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for yyerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from yytname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T
yytnamerr (char *yyres, const char *yystr)
{
  if (*yystr == '"')
    {
      size_t yyn = 0;
      char const *yyp = yystr;

      for (;;)
	switch (*++yyp)
	  {
	  case '\'':
	  case ',':
	    goto do_not_strip_quotes;

	  case '\\':
	    if (*++yyp != '\\')
	      goto do_not_strip_quotes;
	    /* Fall through.  */
	  default:
	    if (yyres)
	      yyres[yyn] = *yyp;
	    yyn++;
	    break;

	  case '"':
	    if (yyres)
	      yyres[yyn] = '\0';
	    return yyn;
	  }
    do_not_strip_quotes: ;
    }

  if (! yyres)
    return yystrlen (yystr);

  return yystpcpy (yyres, yystr) - yyres;
}
# endif

#endif /* YYERROR_VERBOSE */



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
    YYFPRINTF (yyoutput, "token %s (", yytname[yytype]);
  else
    YYFPRINTF (yyoutput, "nterm %s (", yytname[yytype]);


# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# endif
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
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep)
#else
static void
yydestruct (yymsg, yytype, yyvaluep)
    const char *yymsg;
    int yytype;
    YYSTYPE *yyvaluep;
#endif
{
  /* Pacify ``unused variable'' warnings.  */
  (void) yyvaluep;

  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

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



/* The look-ahead symbol.  */
int yychar;

/* The semantic value of the look-ahead symbol.  */
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
    ;
#endif
#endif
{
  
  int yystate;
  int yyn;
  int yyresult;
  /* Number of tokens to shift before error messages enabled.  */
  int yyerrstatus;
  /* Look-ahead token as an internal (translated) token number.  */
  int yytoken = 0;

  /* Three stacks and their tools:
     `yyss': related to states,
     `yyvs': related to semantic values,
     `yyls': related to locations.

     Refer to the stacks thru separate pointers, to allow yyoverflow
     to reallocate them elsewhere.  */

  /* The state stack.  */
  short int yyssa[YYINITDEPTH];
  short int *yyss = yyssa;
  short int *yyssp;

  /* The semantic value stack.  */
  YYSTYPE yyvsa[YYINITDEPTH];
  YYSTYPE *yyvs = yyvsa;
  YYSTYPE *yyvsp;



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
	short int *yyss1 = yyss;


	/* Each stack pointer address is followed by the size of the
	   data in use in that stack, in bytes.  This used to be a
	   conditional around just the two extra args, but that might
	   be undefined if yyoverflow is a macro.  */
	yyoverflow (YY_("memory exhausted"),
		    &yyss1, yysize * sizeof (*yyssp),
		    &yyvs1, yysize * sizeof (*yyvsp),

		    &yystacksize);

	yyss = yyss1;
	yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyexhaustedlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
	goto yyexhaustedlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
	yystacksize = YYMAXDEPTH;

      {
	short int *yyss1 = yyss;
	union yyalloc *yyptr =
	  (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
	if (! yyptr)
	  goto yyexhaustedlab;
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
/* Read a look-ahead token if we need one and don't already have one.  */
/* yyresume: */

  /* First try to decide what to do without reference to look-ahead token.  */

  yyn = yypact[yystate];
  if (yyn == YYPACT_NINF)
    goto yydefault;

  /* Not known => get a look-ahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid look-ahead symbol.  */
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
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
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

  /* Shift the look-ahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

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
#line 1035 "glb_parser.y"
    {YYABORT;}
    break;

  case 5:
#line 1036 "glb_parser.y"
    {YYABORT;}
    break;

  case 8:
#line 1041 "glb_parser.y"
    { /*showlist($1[0]);showlist($1[1]);*/ }
    break;

  case 9:
#line 1042 "glb_parser.y"
    { /*printf ("\t%.10g\n", $1);*/  }
    break;

  case 10:
#line 1043 "glb_parser.y"
    { /*showlist($1);printf("\n");*/ }
    break;

  case 11:
#line 1044 "glb_parser.y"
    {}
    break;

  case 12:
#line 1045 "glb_parser.y"
    {yyerrok;}
    break;

  case 13:
#line 1046 "glb_parser.y"
    {}
    break;

  case 14:
#line 1047 "glb_parser.y"
    {
 set_exp_energy("@energy",(yyvsp[-2].ptrq));
}
    break;

  case 15:
#line 1050 "glb_parser.y"
    {}
    break;

  case 16:
#line 1051 "glb_parser.y"
    {/* printf("%s\n",$1);*/}
    break;

  case 17:
#line 1052 "glb_parser.y"
    {/* printf("%s\n",$1);*/}
    break;

  case 18:
#line 1057 "glb_parser.y"
    { (yyval.val) = (yyvsp[0].val);                         }
    break;

  case 19:
#line 1058 "glb_parser.y"
    { (yyval.val) = (yyvsp[0].nameptr)->value;              }
    break;

  case 20:
#line 1059 "glb_parser.y"
    { (yyval.val) = (yyvsp[0].tptr)->value.var;              }
    break;

  case 21:
#line 1060 "glb_parser.y"
    { (yyval.val) = (yyvsp[0].val); (yyvsp[-2].tptr)->value.var = (yyvsp[0].val);     }
    break;

  case 22:
#line 1061 "glb_parser.y"
    { if(set_exp((yyvsp[-2].name),(yyvsp[0].val),0)==1) yyerror("Unknown identifier");
(yyval.val) = (yyvsp[0].val); }
    break;

  case 23:
#line 1063 "glb_parser.y"
    { if(set_fnct((yyvsp[-2].name),(yyvsp[0].nameptr)->sf)==1) yyerror("Unknown identifier");}
    break;

  case 24:
#line 1064 "glb_parser.y"
    { if(set_pair((yyvsp[-4].name),(yyvsp[-2].val),(yyvsp[0].val),0)==1) yyerror("Unknown identifier");
(yyval.val) = (yyvsp[-2].val); }
    break;

  case 25:
#line 1066 "glb_parser.y"
    { (yyval.val) = (*((yyvsp[-3].tptr)->value.fnctptr))((yyvsp[-1].val)); }
    break;

  case 26:
#line 1068 "glb_parser.y"
    { (yyval.val) = (yyvsp[-2].val) + (yyvsp[0].val);                    }
    break;

  case 27:
#line 1069 "glb_parser.y"
    { (yyval.val) = (yyvsp[-2].val) - (yyvsp[0].val);                    }
    break;

  case 28:
#line 1070 "glb_parser.y"
    { (yyval.val) = (yyvsp[-2].val) * (yyvsp[0].val);                    }
    break;

  case 29:
#line 1071 "glb_parser.y"
    { (yyval.val) = (yyvsp[-2].val) / (yyvsp[0].val);                    }
    break;

  case 30:
#line 1072 "glb_parser.y"
    { (yyval.val) = -(yyvsp[0].val);                        }
    break;

  case 31:
#line 1073 "glb_parser.y"
    { (yyval.val) = pow ((yyvsp[-2].val), (yyvsp[0].val));               }
    break;

  case 32:
#line 1074 "glb_parser.y"
    { (yyval.val) = (yyvsp[-1].val);                         }
    break;

  case 33:
#line 1075 "glb_parser.y"
    { (yyval.val) = 0;}
    break;

  case 34:
#line 1076 "glb_parser.y"
    {yyerror("Unknown name");YYERROR;}
    break;

  case 35:
#line 1082 "glb_parser.y"
    {(yyval.ptr)=list_cons(NULL,(yyvsp[-1].val)); }
    break;

  case 36:
#line 1083 "glb_parser.y"
    {(yyval.ptr)=list_cons(NULL,(yyvsp[-2].val)); }
    break;

  case 37:
#line 1084 "glb_parser.y"
    {
  glb_List *buf;
  buf=(yyvsp[-2].ptr);
  buf=list_cons(buf,(yyvsp[-1].val));
  (yyval.ptr)=buf;   
}
    break;

  case 38:
#line 1091 "glb_parser.y"
    {
  glb_List *buf;
  buf=(yyvsp[-3].ptr);
  buf=list_cons(buf,(yyvsp[-2].val));
  (yyval.ptr)=buf;   
}
    break;

  case 39:
#line 1097 "glb_parser.y"
    {
  glb_List *buf;
  buf=(yyvsp[-1].ptr);
  buf=list_cons(buf,(yyvsp[0].val));
  (yyval.ptr)=buf;   
}
    break;

  case 40:
#line 1103 "glb_parser.y"
    {(yyval.ptr)=NULL;}
    break;

  case 41:
#line 1104 "glb_parser.y"
    {(yyval.ptr)=(yyvsp[-1].ptr); }
    break;

  case 42:
#line 1105 "glb_parser.y"
    {(yyval.ptr)=list_cons(NULL,(yyvsp[-1].val)); }
    break;

  case 43:
#line 1106 "glb_parser.y"
    { if(set_exp_list((yyvsp[-2].name),(yyvsp[0].ptr),3)==1) 
  yyerror("Unknown identifier");
 (yyval.ptr) = (yyvsp[0].ptr); }
    break;

  case 44:
#line 1110 "glb_parser.y"
    {(yyval.ptr) = thread_list((yyvsp[-3].tptr)->value.fnctptr,(yyvsp[-3].tptr)->reverse,(yyvsp[-3].tptr)->destroy,(yyvsp[-1].ptr));}
    break;

  case 45:
#line 1111 "glb_parser.y"
    {(yyval.ptr) = (*((yyvsp[-2].tptr)->value.lfnctptr))();}
    break;

  case 46:
#line 1112 "glb_parser.y"
    { (yyval.ptr) = (yyvsp[0].tptr)->list;              }
    break;

  case 47:
#line 1113 "glb_parser.y"
    { (yyval.ptr) = (yyvsp[0].ptr); (yyvsp[-2].tptr)->list = (yyvsp[0].ptr);/* printf("%s ",$1->name);showlist($3);printf("\n");*/ }
    break;

  case 48:
#line 1114 "glb_parser.y"
    {(yyval.ptr)=glb_interpolation((yyvsp[-7].ptr),(yyvsp[-5].ptr),floor((yyvsp[-3].val)),(yyvsp[-1].ptr));}
    break;

  case 49:
#line 1117 "glb_parser.y"
    { 

  double *buf;
  buf=(double*) glb_malloc(sizeof(double)*2);
  buf[0]=(yyvsp[-2].val);
  buf[1]=(yyvsp[0].val)-1;
  (yyval.dpt)=buf;
}
    break;

  case 50:
#line 1127 "glb_parser.y"
    {}
    break;

  case 51:
#line 1128 "glb_parser.y"
    {}
    break;

  case 52:
#line 1132 "glb_parser.y"
    { if((yyvsp[-1].nameptr)->value==-1) {(yyvsp[-1].nameptr)->value=step_counter((yyvsp[-3].name)); }
  loc_count=(yyvsp[-1].nameptr)->value;
  glb_free(context);
  context =(char *) strdup((yyvsp[-3].name));
  grp_start(context);
}
    break;

  case 53:
#line 1138 "glb_parser.y"
    { 
   
grp_end(context);}
    break;

  case 54:
#line 1141 "glb_parser.y"
    {
  yyerror("Redefinition of an automatic variable");YYERROR;}
    break;

  case 56:
#line 1146 "glb_parser.y"
    {}
    break;

  case 57:
#line 1147 "glb_parser.y"
    {}
    break;

  case 58:
#line 1148 "glb_parser.y"
    {}
    break;

  case 59:
#line 1149 "glb_parser.y"
    {}
    break;

  case 60:
#line 1150 "glb_parser.y"
    {}
    break;

  case 61:
#line 1151 "glb_parser.y"
    {}
    break;

  case 62:
#line 1154 "glb_parser.y"
    {
  buff.version=strdup((yyvsp[0].name));
}
    break;

  case 63:
#line 1159 "glb_parser.y"
    {
  //load_cross($3,loc_count-1);
  xsc.file_name=strdup((yyvsp[0].name));
  (yyval.name)=(yyvsp[0].name);
}
    break;

  case 64:
#line 1166 "glb_parser.y"
    {
  //load_flux($3,loc_count-1,1);
  flt.file_name=strdup((yyvsp[0].name));
  //if(set_exp($1,$3,0)==1) yyerror("Unknown identifier");
  (yyval.name)=(yyvsp[0].name);
}
    break;

  case 65:
#line 1174 "glb_parser.y"
    {
  buff.sys_on_strings[buff.numofrules-1] = strdup((yyvsp[0].name));
}
    break;

  case 66:
#line 1179 "glb_parser.y"
    {
  buff.sys_off_strings[buff.numofrules-1] = strdup((yyvsp[0].name));
}
    break;

  case 67:
#line 1186 "glb_parser.y"
    {

  int x[6];
  x[0]=(yyvsp[-10].nameptr)->value - 1 ;
  x[1]=(yyvsp[-8].in);
  x[2]=(yyvsp[-6].in);
  x[3]=(yyvsp[-4].in);
  x[4]=(int) (yyvsp[-2].nameptr)->value -1; 
  x[5]=(int) (yyvsp[0].nameptr)->value -1;


  set_channel_data(x,loc_count);
}
    break;

  case 68:
#line 1202 "glb_parser.y"
    {(yyval.nameptr)=(yyvsp[0].nameptr);}
    break;

  case 69:
#line 1203 "glb_parser.y"
    {yyerror("Unknown name");YYERROR;}
    break;

  case 70:
#line 1207 "glb_parser.y"
    {(yyval.in)=(yyvsp[0].in);}
    break;

  case 71:
#line 1208 "glb_parser.y"
    {(yyval.in)=1;}
    break;

  case 72:
#line 1209 "glb_parser.y"
    {(yyval.in)=-1;}
    break;

  case 73:
#line 1213 "glb_parser.y"
    {
  glb_List **buf;
  energy_len=1;
 
  buf=(glb_List**) glb_malloc(sizeof( glb_List* ) ); 
  buf[0]=(yyvsp[0].ptr);
 
  (yyval.ptrq)=buf;
}
    break;

  case 74:
#line 1223 "glb_parser.y"
    {
  glb_List **buf;
  buf=(yyvsp[-2].ptrq);
  energy_len++;
  
  buf=(glb_List**) glb_realloc((void**) buf , sizeof( glb_List* ) * energy_len);
 
  buf[energy_len-1]=(yyvsp[0].ptr); 
  (yyval.ptrq)=buf; }
    break;

  case 75:
#line 1233 "glb_parser.y"
    {
  glb_List **buf;
  buf=(yyvsp[-3].ptrq);
  energy_len++;

  buf=(glb_List**) glb_realloc((void**) buf , sizeof( glb_List* ) * energy_len);

  buf[energy_len-1]=(yyvsp[0].ptr); 
  (yyval.ptrq)=buf; }
    break;

  case 76:
#line 1246 "glb_parser.y"
    {
  glb_List **buf;
  buf=(glb_List**) glb_malloc(sizeof(glb_List*)*2);
  buf[0]=list_cons(NULL,(yyvsp[0].dpt)[0]);
  buf[1]=list_cons(NULL,(yyvsp[0].dpt)[1]);
  glb_free((yyvsp[0].dpt));
  (yyval.ptrq)=buf;
}
    break;

  case 77:
#line 1254 "glb_parser.y"
    {
  glb_List **buf;
  buf=(yyvsp[-2].ptrq);
  buf[0]=list_cons(buf[0],(yyvsp[0].dpt)[0]);
  buf[1]=list_cons(buf[1],(yyvsp[0].dpt)[1]);
  glb_free((yyvsp[0].dpt));
  (yyval.ptrq)=buf; 
}
    break;

  case 78:
#line 1264 "glb_parser.y"
    {
  glb_List **buf;  
 
  buf=(glb_List**) glb_malloc(sizeof(glb_List*)*2);
  buf[0]=list_cons(NULL,(yyvsp[0].dpt)[0]);
  buf[1]=list_cons(NULL,(yyvsp[0].dpt)[1]);
  glb_free((yyvsp[0].dpt));
  (yyval.ptrq)=buf;
}
    break;

  case 79:
#line 1273 "glb_parser.y"
    {
  glb_List **buf;
  buf=(yyvsp[-2].ptrq);
  buf[0]=list_cons(buf[0],(yyvsp[0].dpt)[0]);
  buf[1]=list_cons(buf[1],(yyvsp[0].dpt)[1]);
  glb_free((yyvsp[0].dpt));
  (yyval.ptrq)=buf; 
}
    break;

  case 80:
#line 1283 "glb_parser.y"
    {

int flag;  
  
  (yyval.ptrq)=(yyvsp[0].ptrq);
  
  flag=set_exp_list("bgrulescoeff",(yyvsp[0].ptrq)[0],0);
  if(flag==1) yyerror("Unknown identifier");
  flag=set_exp_list("bgrulechannellist",(yyvsp[0].ptrq)[1],0);
  if(flag==1) yyerror("Unknown identifier");
 
  glb_free((yyvsp[0].ptrq));
}
    break;

  case 81:
#line 1296 "glb_parser.y"
    {
  int flag;  
  (yyval.ptrq)=(yyvsp[0].ptrq);
  flag=set_exp_list("rulescoeff",(yyvsp[0].ptrq)[0],0);
  if(flag==1) yyerror("Unknown identifier");
  flag=set_exp_list("rulechannellist",(yyvsp[0].ptrq)[1],0);
  if(flag==1) yyerror("Unknown identifier"); 

  glb_free((yyvsp[0].ptrq));

 
}
    break;

  case 82:
#line 1310 "glb_parser.y"
    {
 if((yyvsp[0].nameptr)->value==-1) {(yyvsp[0].nameptr)->value=step_counter("rule");printf("loc step\n"); }
  loc_count=(yyvsp[0].nameptr)->value;
  printf("loc_count %d \n",loc_count);
  YYERROR;}
    break;

  case 83:
#line 1315 "glb_parser.y"
    { YYERROR;}
    break;

  case 84:
#line 1319 "glb_parser.y"
    {
 if((yyvsp[0].nameptr)->value==-1) {(yyvsp[0].nameptr)->value=step_counter("channel");printf("loc step\n"); }
  loc_count=(yyvsp[0].nameptr)->value;
  printf("cha loc_count %d \n",loc_count); YYERROR;}
    break;

  case 85:
#line 1323 "glb_parser.y"
    { YYERROR;}
    break;


      default: break;
    }

/* Line 1126 of yacc.c.  */
#line 2873 "glb_parser.c"

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
	  int yytype = YYTRANSLATE (yychar);
	  YYSIZE_T yysize0 = yytnamerr (0, yytname[yytype]);
	  YYSIZE_T yysize = yysize0;
	  YYSIZE_T yysize1;
	  int yysize_overflow = 0;
	  char *yymsg = 0;
#	  define YYERROR_VERBOSE_ARGS_MAXIMUM 5
	  char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
	  int yyx;

#if 0
	  /* This is so xgettext sees the translatable formats that are
	     constructed on the fly.  */
	  YY_("syntax error, unexpected %s");
	  YY_("syntax error, unexpected %s, expecting %s");
	  YY_("syntax error, unexpected %s, expecting %s or %s");
	  YY_("syntax error, unexpected %s, expecting %s or %s or %s");
	  YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s");
#endif
	  char *yyfmt;
	  char const *yyf;
	  static char const yyunexpected[] = "syntax error, unexpected %s";
	  static char const yyexpecting[] = ", expecting %s";
	  static char const yyor[] = " or %s";
	  char yyformat[sizeof yyunexpected
			+ sizeof yyexpecting - 1
			+ ((YYERROR_VERBOSE_ARGS_MAXIMUM - 2)
			   * (sizeof yyor - 1))];
	  char const *yyprefix = yyexpecting;

	  /* Start YYX at -YYN if negative to avoid negative indexes in
	     YYCHECK.  */
	  int yyxbegin = yyn < 0 ? -yyn : 0;

	  /* Stay within bounds of both yycheck and yytname.  */
	  int yychecklim = YYLAST - yyn;
	  int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
	  int yycount = 1;

	  yyarg[0] = yytname[yytype];
	  yyfmt = yystpcpy (yyformat, yyunexpected);

	  for (yyx = yyxbegin; yyx < yyxend; ++yyx)
	    if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR)
	      {
		if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
		  {
		    yycount = 1;
		    yysize = yysize0;
		    yyformat[sizeof yyunexpected - 1] = '\0';
		    break;
		  }
		yyarg[yycount++] = yytname[yyx];
		yysize1 = yysize + yytnamerr (0, yytname[yyx]);
		yysize_overflow |= yysize1 < yysize;
		yysize = yysize1;
		yyfmt = yystpcpy (yyfmt, yyprefix);
		yyprefix = yyor;
	      }

	  yyf = YY_(yyformat);
	  yysize1 = yysize + yystrlen (yyf);
	  yysize_overflow |= yysize1 < yysize;
	  yysize = yysize1;

	  if (!yysize_overflow && yysize <= YYSTACK_ALLOC_MAXIMUM)
	    yymsg = (char *) YYSTACK_ALLOC (yysize);
	  if (yymsg)
	    {
	      /* Avoid sprintf, as that infringes on the user's name space.
		 Don't have undefined behavior even if the translation
		 produced a string with the wrong number of "%s"s.  */
	      char *yyp = yymsg;
	      int yyi = 0;
	      while ((*yyp = *yyf))
		{
		  if (*yyp == '%' && yyf[1] == 's' && yyi < yycount)
		    {
		      yyp += yytnamerr (yyp, yyarg[yyi++]);
		      yyf += 2;
		    }
		  else
		    {
		      yyp++;
		      yyf++;
		    }
		}
	      yyerror (yymsg);
	      YYSTACK_FREE (yymsg);
	    }
	  else
	    {
	      yyerror (YY_("syntax error"));
	      goto yyexhaustedlab;
	    }
	}
      else
#endif /* YYERROR_VERBOSE */
	yyerror (YY_("syntax error"));
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse look-ahead token after an
	 error, discard it.  */

      if (yychar <= YYEOF)
        {
	  /* Return failure if at end of input.  */
	  if (yychar == YYEOF)
	    YYABORT;
        }
      else
	{
	  yydestruct ("Error: discarding", yytoken, &yylval);
	  yychar = YYEMPTY;
	}
    }

  /* Else will try to reuse look-ahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:

  /* Pacify compilers like GCC when the user code never invokes
     YYERROR and the label yyerrorlab therefore never appears in user
     code.  */
  if (0)
     goto yyerrorlab;

yyvsp -= yylen;
  yyssp -= yylen;
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
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


      yydestruct ("Error: popping", yystos[yystate], yyvsp);
      YYPOPSTACK;
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  if (yyn == YYFINAL)
    YYACCEPT;

  *++yyvsp = yylval;


  /* Shift the error token. */
  YY_SYMBOL_PRINT ("Shifting", yystos[yyn], yyvsp, yylsp);

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
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (yychar != YYEOF && yychar != YYEMPTY)
     yydestruct ("Cleanup: discarding lookahead",
		 yytoken, &yylval);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
		  yystos[*yyssp], yyvsp);
      YYPOPSTACK;
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
  return yyresult;
}


#line 1326 "glb_parser.y"


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
  /* these two flags determine the behaviour under threading
     reverse 0 leaves the list as is (i.e reverse)
     reverse 1 reverse the list (right order)

     destroy -1 leaves the list (even the memory location) unscathed
     destroy 0  procuces a copy
     destroy 1 destroys the input list and returns a new list
  */
  int reverse;
  int destroy;
};

struct glb_init_list
{
  char *fname;
  glb_List *(*fnct)(void);
};


struct glb_init_sig
{
  char *fname;
  sigfun sf;
};
     
static double echo(double x)
{
  fprintf(stdout,"%f ",x);
  return x;
}

static double echon(double x)
{
  fprintf(stdout,"%f\n",x);
  return x;
}


static double line(double x)
{
  size_t n,i;
  n=floor(x);
  for(i=0;i<n;i++) fprintf(stdout,"\n",x);
  return x;
}

static struct glb_init arith_fncts[] =
  {
    {"sin",  sin, 0,1},
    {"cos",  cos,0,1},
    {"tan",  tan,0,1},
    {"asin", asin,0,1},
    {"acos", acos,0,1},
    {"atan", atan,0,1},
    {"log",   log,0,1},
    {"exp",  exp,0,1},
    {"log10", log10,0,1},
    {"sqrt", sqrt,0,1},
    {"reverse",glb_reverse,1,1},
    {"copy",glb_list_copy,0,0},
    {"echo",echo,0,-1},
    {"echon",echon,0,-1},
    {"line",line,0,-1},
    {NULL, NULL,0,0}
     };


//void glb_set_up_smear_data(glb_smear *test,const struct glb_experiment *head)
static glb_List *glb_bincenter(void)
{
  int i;
  glb_List *res=NULL;
  glb_smear *test;
  test=glb_smear_alloc();
  
  if(buff.numofbins<0) {glb_error("Cannot compute bincenter. Binning not set up properly."); return NULL;}

  glb_set_up_smear_data(test,&buff);

 
  for(i=0;i<test->numofbins;i++)
    {    
      res=list_cons(res,test->bincenter[i]);
    }
 
  glb_smear_free(test);
 
  return res;
}


static glb_List *glb_simbincenter(void)
{
  int i;
  glb_List *res=NULL;
  glb_smear *test;
  test=glb_smear_alloc();
  
  if(buff.simbins<0) {glb_error("Cannot compute simbincenter. Sim-binning not set up properly."); return NULL;}

  
  glb_set_up_smear_data(test,&buff);

  for(i=0;i<test->simbins;i++)
    {
      
      res=list_cons(res,test->simbincenter[i]);
    }
 
  glb_smear_free(test);
 
  
  return res;
}




static struct glb_init_list list_fncts[] =
  {
    {"bincenter",glb_bincenter},
    {"simbincenter",glb_simbincenter},
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

#define BIN_LIST 1
#define SAMPLING_LIST 2
#define DENSITY_LIST 3


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
      ptr->destroy=arith_fncts[i].destroy;
      ptr->reverse=arith_fncts[i].reverse;
 
    }
  for (i = 0; list_fncts[i].fname != 0; i++)
    {
      ptr = glb_putsym (list_fncts[i].fname, FNCT);
      ptr->value.lfnctptr = list_fncts[i].fnct;
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
	if(ptr->list!=NULL){ list_free(ptr->list);}
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
  ptr->list= NULL;
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
  glb_flux_reset(&flt);
  glb_xsec_reset(&xsc);
 
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

