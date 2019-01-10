/* A Bison parser, made by GNU Bison 3.0.4.  */

/* Bison implementation for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2015 Free Software Foundation, Inc.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "3.0.4"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1




/* Copy the first part of user declarations.  */
#line 22 "glb_parser.y" /* yacc.c:339  */

#if HAVE_CONFIG_H   /* config.h should come before any other includes */
#  include "config.h"
#endif



#define YYDEBUG 1
//int yydebug = 1; /* Uncomment this line to get debug output from the parser */
#include <math.h>
#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
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
#include "glb_minimize.h"

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
//#define INT_INDEXED_PAIR 8
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
  static struct glb_experiment *buff_list[GLB_MAX_EXP];
  static glb_smear ibf;
  static glb_option_type opt;
  static glb_flux flt;
  static glb_xsec xsc;
  static glb_nuisance nuis;
  static int errordim_sys_on=-1;         /* Temp. storage for the old numerical errordims */
  static int errordim_sys_off=-1;
  static char *context;

  int yyerror (const char *s, ...);           /* Forward declaration to suppress compiler warning */


  /* The name and symbol tables: a chain of `struct glb_symrec'.  */
  static glb_namerec *name_table = (glb_namerec *) NULL;
  glb_symrec *sym_table = (glb_symrec *) NULL;
  glb_symrec *pre_sym_table = (glb_symrec *) NULL;
          /* cannot use static here, since its declared earlier as extern */

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


  /* IMPORTANT NOTE: Token names that are substrings of other token names
     (e.g. @error and @error_list) should appear *after* that other token */
  static glb_parser_decl token_list[]={
    {"$version",CHAR,0,1E8,&buff.version,NULL,"global"},
    {"$citation",CHAR,0,1E8,&buff.citation,NULL,"global"},
    {"$parent_energy",DOUBLE,0,GMAX,&buff.emax,NULL,"global"},
    {"$oscillation_engine",CHAR,0,1E8,&buff.osc_engine.name,NULL,"global"},
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
   {"@data_flag",INT_INDEXED,0,1,&buff.data_on_off[0],&loc_count,"rule"},
   {"@data",DOUBLE_LIST_INDEXED,-GMAX,GMAX,&buff.data[0],&loc_count,"rule"},

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
   {"@sys_on_multiex_errors_sig",  ENERGY_MATRIX, 1, GLB_MAX_NUISANCE,
     &buff.sys_on_multiex_errors_sig[0], &loc_count, "rule"},
   {"@sys_on_multiex_errors_bg",   ENERGY_MATRIX, 1, GLB_MAX_NUISANCE,
     &buff.sys_on_multiex_errors_bg[0], &loc_count, "rule"},
   {"@sys_off_multiex_errors_sig", ENERGY_MATRIX, 1, GLB_MAX_NUISANCE,
     &buff.sys_off_multiex_errors_sig[0], &loc_count, "rule"},
   {"@sys_off_multiex_errors_bg",  ENERGY_MATRIX, 1, GLB_MAX_NUISANCE,
     &buff.sys_off_multiex_errors_bg[0], &loc_count, "rule"},


   {"sys", UNTYPE, 0, 20, NULL, &buff.n_nuisance, "global"},
   {"@energy_list", DOUBLE_LIST, 0, GMAX, &nuis.energy_list, &nuis.n_energies, "sys"},
   {"@error_list",  DOUBLE_LIST, 0, GMAX, &nuis.error_list,  &nuis.n_energies, "sys"},
   {"@error",       DOUBLE,      0, GMAX, &nuis.error,       NULL,             "sys"}, 
    {"@systype",       INT,      0, 20, &nuis.systype,       NULL,             "sys"},


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



    {"nuflux",UNTYPE,0,20,NULL,&buff.num_of_fluxes,"global"},
    {"@flux_file",CHAR,0,20,&flt.file_name,NULL,"nuflux"},

    {"@builtin",INT,1,4,&flt.builtin,NULL,"nuflux"},
    {"@time",DOUBLE,0,GMAX,&flt.time,NULL,"nuflux"},
    {"@power",DOUBLE,0,GMAX,&flt.target_power,NULL,"nuflux"},
    {"@stored_muons",DOUBLE,0,GMAX,&flt.stored_muons,NULL,"nuflux"},
    {"@parent_energy",DOUBLE,0,GMAX,&flt.parent_energy,NULL,"nuflux"},
    {"@norm",DOUBLE,0,GMAX,&flt.norm,NULL,"nuflux"},
    {"@gamma",DOUBLE,0,GMAX,&flt.gamma,NULL,"nuflux"},
    {"@end_point",DOUBLE,0,GMAX,&flt.end_point,NULL,"nuflux"},
    {"@stored_ions",DOUBLE,0,GMAX,&flt.stored_muons,NULL,"nuflux"},


   {"@type" ,INT,1,2,&ibf.type,&loc_count,"energy"},
   {"@sigma_e" ,DOUBLE_LIST,0,GMAX,&ibf.sigma,&ibf.num_of_params,"energy"},
   {"@sigma_function" ,FUN,0,1000,&ibf.sig_f,NULL,"energy"},



   {NULL,UNTYPE,0,0,NULL,NULL,"global"}
};


static void grp_start(char* name, int id) // id = index of the group among groups of same type
   {
     if(strncmp(name,"rule",4)==0 )
       {
         /* Reset variables to default values */
         errordim_sys_on  = -1.0;
         errordim_sys_off = -1.0;
       }
   }


static void grp_end(char* name, int id) // id = index of the group among groups of same type
   {
     if (strncmp(name,"energy",6)==0 )
       {
         if(id-1 >= 0)
           {
             ibf.options=glb_option_type_alloc();
             ibf.options=(glb_option_type *) memmove(ibf.options,&opt,
                                             sizeof(glb_option_type));

             if (buff.smear_data[id-1] != NULL)
             {
               glb_smear_free(buff.smear_data[id-1]);
               buff.smear_data[id-1] = NULL;
             }
             buff.smear_data[id-1] = glb_smear_alloc();
             buff.smear_data[id-1] = glb_copy_smear(buff.smear_data[id-1],&ibf);
             glb_option_type_free(ibf.options);
             glb_option_type_reset(&opt);
             if(ibf.sigma!=NULL) glb_free(ibf.sigma);
             glb_smear_reset(&ibf);
           }
       }

     if (strncmp(name,"flux",4)==0 )
       {
         if(id > 0)
           {
             glb_error("The 'flux' directive is deprecated (consider using 'nuflux').\n"
                       "The flux normalization may not be what you expect.\n"
                       "Please, consult the manual!");
             if(flt.builtin<=0) flt.builtin=GLB_OLD_NORM;

             if (buff.fluxes[id-1] != NULL)
             { 
               glb_free_flux(buff.fluxes[id-1]);
               buff.fluxes[id-1] = NULL;
             }
             buff.fluxes[id-1] = glb_malloc(sizeof(glb_flux));
             memset(buff.fluxes[id-1], 0, sizeof(*buff.fluxes[0]));
             glb_reset_flux(buff.fluxes[id-1]);
             if (glb_copy_flux(buff.fluxes[id-1], &flt) != GLB_SUCCESS)
               glb_error("grp_end: Error copying flux data");
             glb_reset_flux(&flt);
           }
       }


     if (strncmp(name,"nuflux",6)==0 )
       {
         if(id > 0)
           {
             if (buff.fluxes[id-1] != NULL)
             { 
               glb_free_flux(buff.fluxes[id-1]);
               buff.fluxes[id-1] = NULL;
             }
             buff.fluxes[id-1] = glb_malloc(sizeof(glb_flux));
             memset(buff.fluxes[id-1], 0, sizeof(*buff.fluxes[0]));
             glb_reset_flux(buff.fluxes[id-1]);
             if (glb_copy_flux(buff.fluxes[id-1], &flt) != GLB_SUCCESS)
               glb_error("grp_end: Error copying flux data");
             glb_reset_flux(&flt);
           }
       }

     if (strncmp(name,"cross",5)==0 )
       {
         if(id > 0)
           {
             if (buff.xsecs[id-1] != NULL)
             {
               glb_free_xsec(buff.xsecs[id-1]);
               buff.xsecs[id-1] = NULL;
             }
             buff.xsecs[id-1] = glb_malloc(sizeof(glb_xsec));
             memset(buff.xsecs[id-1], 0, sizeof(*buff.xsecs[0]));
             glb_reset_xsec(buff.xsecs[id-1]);
             if (glb_copy_xsec(buff.xsecs[id-1], &xsc) != GLB_SUCCESS)
               glb_error("grp_end: Error copying cross section data");
             glb_reset_xsec(&xsc);
           }
       }



     if (strncmp(name,"rule",4)==0 )
       {
         int nr = id - 1;

         /* Parse old (numerical) errordims */
         if (buff.sys_on_strings[nr] == NULL  &&  errordim_sys_on >= 0)
           buff.sys_on_strings[nr] = glbConvertErrorDim(errordim_sys_on);
         if (buff.sys_off_strings[nr] == NULL  &&  errordim_sys_off >= 0)
           buff.sys_off_strings[nr] = glbConvertErrorDim(errordim_sys_off);
       }

     if (strncmp(name,"sys",3) == 0)
     {
       if (nuis.name)
       {
         glb_free(nuis.name);
         nuis.name = NULL;
       }
       nuis.name = strdup(name_table->name);
       if (!buff.nuisance_params[id-1])
         buff.nuisance_params[id-1] = glb_alloc_nuisance();
       glb_nuisance *n = buff.nuisance_params[id-1];
       if (n)  memcpy(n, &nuis, sizeof(glb_nuisance));
       glbResetNuisance();
     }

     glb_free(context);
     context = (char *) strdup("global");

   }

static int set_channel_data(int x[6],int loc_count)
   {
     // I am sorry for this -- it's a kludge
     int i;

     for(i=0;i<6;i++) {
       if (loc_count >= buff.numofchannels) /* Don't realloc when parsing a re-definition */
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
/* JK 2012-06-25 Remove this restriction */
/*                     if(!((*ibf == -1) || (*ibf == (int) value))) {
                       glb_warning("Given length does not"
                                   " match actual length");
                       return 2;

                     }*/
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

             // JK - 2012-05-17 I think the following wouldn't work. It's not used,
             // so I comment it out for the moment
//             if(token_list[i].scalar==INT_INDEXED_PAIR) //int
//               {
//                 if(value >= token_list[i].rl && value <= token_list[i].ru)
//                   {
//
//                     ibf=(int*) token_list[i].ptr;
//                     ibf[(loc_count-1)+0*32]=(int) value;
//                     return 0;
//                   }
//                 else
//                   {
//                     fprintf(stderr,"Error: Value for %s out of range\n",
//                             token_list[i].token);
//                     return 2;
//                   }
//               }

             if(token_list[i].scalar==DOUBLE_INDEXED_PAIR) //int
               {
                 if(value >= token_list[i].rl && value <= token_list[i].ru)
                   {
                     if (strcmp(token_list[i].ctx, "rule") != 0)
                       fprintf(stderr, "Error: DOUBLE_INDEXED_PAIR works only inside "
                                       "a rule environment! Ignoring it.\n");
                     else
                     {
                       int delta = GLB_MAX_RULES;
                       dbf=(double*) token_list[i].ptr;
                       dbf[(loc_count-1)+0*delta]=(double) value;
                       dbf[(loc_count-1)+1*delta]=(double) value2;
                     }
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


static int set_string(char *name, char *str)
{
  int i;
  for(i=0; token_list[i].token != NULL; i++)
  {
    if (strncmp(name, token_list[i].token, strlen(token_list[i].token)) == 0)
    {
      char **p = (char **)(token_list[i].ptr);
      if (*p)  glb_free(*p);
      *p = strdup(str);
      return 0;
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
  double *xlist,*ylist,*rlist;
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

    /* PH 06/09/16 in going from GSL 1.14 to 1.15 the behavior
     changed: if head->entry is outside the range of xlist a fatal
     error is generated, i.e. no more implicit extrapolation, so we
     have to catch that by hand */

    if(head->entry<xlist[0])
      {
	glb_warning("Linear extrapolation used because x value %f is out of range %f - %f",head->entry,xlist[0],xlist[yl-1]);
	rlist[i]=ylist[0]-(xlist[0]-head->entry)*(ylist[1]-ylist[0])/(xlist[1]-xlist[0]);
      }
    else if(head->entry>xlist[yl-1])
      {
	glb_warning("Linear extrapolation used because x value %f is out of range %f - %f",head->entry,xlist[0],xlist[yl-1]);
	rlist[i]=ylist[yl-1]-(xlist[yl-1]-head->entry)*(ylist[yl-2]-ylist[yl-1])/(xlist[yl-2]-xlist[yl-1]);
      }
    else
      {   
	rlist[i]=gsl_spline_eval(spline,head->entry,acc);
      }
     head=head->next; }

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
/* JK 2012-06-26 Removed to allow for redefinitions of binning
            else if(*lbf!=len) glb_fatal("Line %d: Length mismatch or list"
                                         " length changed", glb_line_num);*/


            dbf = (double**) token_list[i].ptr;
            if(*dbf!=NULL){glb_free(*dbf);*dbf=NULL;}
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
                      *dbf=NULL;
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
            if(dbf[loc_count-1]!=NULL){glb_free(dbf[loc_count-1]);dbf[loc_count-1]=NULL;}
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
                    dbf[loc_count-1]=NULL;

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

  for(i=0; token_list[i].token != NULL; i++)
  {
    if(strncmp(name,token_list[i].token, strlen(token_list[i].token)) == 0 &&
       strncmp(context,token_list[i].ctx, strlen(token_list[i].ctx)) == 0)
    {
      switch((int) token_list[i].scalar)
      {
        case ENERGY_MATRIX:
          list = (double**) glb_malloc(sizeof(double* ) * (energy_len+1));
          buff.lowrange[loc_count-1] = (int*) glb_malloc((energy_len+1)*sizeof(int));
          buff.uprange[loc_count-1]  = (int*) glb_malloc((energy_len+1)*sizeof(int));

          /* Loop over all analysis bins */
          for(l=0; l < energy_len; l++)
          {
            len = (int) list_length(value[l]); /* how long is the list provided by the user? */
            if (len < 2) {
              fprintf(stderr, "Error: in line %d: in @smear number %d: sublist %d too short!\n",
                              glb_line_num, loc_count, l); return 2;
            }

            dbf = (double***) token_list[i].ptr;
            list[l] = (double*) glb_malloc(sizeof(double)*(len-2+1));
            dbf[loc_count-1] = list;

            v1 = (int) list_take(value[l], len-0-1); /* Sampling point range for this bin */
            v2 = (int) list_take(value[l], len-1-1);

            if(v1 >= 0  && v2 <= buff.simbins  &&
               v2 >= v1 && v2-v1 == len-3)
            {
              buff.lowrange[loc_count-1][l] = v1;
              buff.uprange[loc_count-1][l]  = v2;
            }
            else
            {
              fprintf(stderr,"Error: In line %d: Value for ranges in smear out of range\n",
                      glb_line_num);
              glb_free(list[l]);
              glb_free(buff.lowrange[loc_count-1]);
              glb_free(buff.uprange[loc_count-1]);
              return 2;
            }

            /* Loop over all sampling points contributing to the current analysis bin */
            for(k=0; k < len-2; k++)
            {
              val = list_take(value[l], len-(k+2)-1);
              if (val >= token_list[i].rl && val <= token_list[i].ru)
                list[l][k]=val;
              else
              {
                fprintf(stderr,"Error: In line %d: Value for %s out of range\n",
                        glb_line_num, token_list[i].token);
                free(list[l]);
                return 2;
              }
            } /* for (k) */
            list_free(value[l]);
            value[l] = NULL;

            list[l][len-2] = -1; /* This signals the end of the list */
          } /* for (l) */
          glb_free(value);
          value = NULL;

          list[energy_len] = NULL; /* This signals the end of the list */
          buff.lowrange[loc_count-1][energy_len] = -1;
          buff.uprange[loc_count-1][energy_len] = -1;
          return 0;
          break;
      } /* switch */
    }
  } /* for (i) */
  return 1;
}


static int set_multiex_errors(char *name, glb_List **value)
{
  int i, j, k;
  for(i=0; token_list[i].token != NULL; i++)
  {
      if(strncmp(name, token_list[i].token, strlen(token_list[i].token)) == 0  &&
         strncmp(context,token_list[i].ctx, strlen(token_list[i].ctx)) == 0)
      {
        switch((int) token_list[i].scalar)
        {
          case ENERGY_MATRIX:
            if (energy_len > GLB_MAX_RULES)
            {
              fprintf(stderr, "Error in line %d: @sys_XX_multiex_errors_YY definition too long.\n",
                      glb_line_num);
              return 2;
            }
            for (j=0; j < energy_len; j++)
            {
              int len = (int) list_length(value[j]);
              if (len > GLB_MAX_CHANNELS-1)
              {
                fprintf(stderr, "Error in line %d: Entry %d in @sys_XX_multiex_errors_YY too long.\n",
                        glb_line_num, j+1);
                return 3;
              }
              int *(*x)[GLB_MAX_CHANNELS] = (int *(*)[GLB_MAX_CHANNELS]) token_list[i].ptr;
              x[loc_count-1][j] = glb_malloc(sizeof(int) * (len+1));
              if (len > 0)
              {
                for (k=0; k < len; k++)
                {
                  double val = list_take(value[j], len-k-1);
                  if(val >= token_list[i].rl && val <= token_list[i].ru)
                    x[loc_count-1][j][k] = (int)(val+0.5 - 1);
                      /* -1 because first sys<> name has index 1 in parser */
                  else
                  {
                    fprintf(stderr, "Error in line %d: Value for %s out of range\n",
                            glb_line_num, token_list[i].token);
                    glb_free(x[loc_count-1][j]);
                    x[loc_count-1][j] = NULL;
                    return 4;
                  }
                } /* for (k) */
                x[loc_count-1][j][k] = -1; /* This signals the end of the list */
              }
              else
                x[loc_count-1][j][0] = -1;
 
              list_free(value[j]);
              value[j] = NULL;
            } /* for (j) */

            free(value);
            return 0;
            break;
        }
      }
  } /* for (i) */

  return 1;
}


 

#line 1237 "glb_parser.c" /* yacc.c:339  */

# ifndef YY_NULLPTR
#  if defined __cplusplus && 201103L <= __cplusplus
#   define YY_NULLPTR nullptr
#  else
#   define YY_NULLPTR 0
#  endif
# endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

/* In a future release of Bison, this section will be replaced
   by #include "y.tab.h".  */
#ifndef YY_YY_Y_TAB_H_INCLUDED
# define YY_YY_Y_TAB_H_INCLUDED
/* Debug traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif
#if YYDEBUG
extern int yydebug;
#endif

/* Token type.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
  enum yytokentype
  {
    NUM = 258,
    SFNCT = 259,
    BOGUS = 260,
    LVAR = 261,
    VAR = 262,
    FNCT = 263,
    IDN = 264,
    CROSS = 265,
    FLUXP = 266,
    FLUXM = 267,
    NUFLUX = 268,
    SYS_ON_FUNCTION = 269,
    SYS_OFF_FUNCTION = 270,
    SYS_MULTIEX_ERRORS = 271,
    GRP = 272,
    GID = 273,
    FNAME = 274,
    VERS = 275,
    SIGNAL = 276,
    BG = 277,
    ENERGY = 278,
    CHANNEL = 279,
    NDEF = 280,
    GRPOPEN = 281,
    GRPCLOSE = 282,
    PM = 283,
    FLAVOR = 284,
    NOGLOBES = 285,
    RULESEP = 286,
    RULEMULT = 287,
    NAME = 288,
    RDF = 289,
    ENDEXP = 290,
    ENDDET = 291,
    NEG = 292
  };
#endif
/* Tokens.  */
#define NUM 258
#define SFNCT 259
#define BOGUS 260
#define LVAR 261
#define VAR 262
#define FNCT 263
#define IDN 264
#define CROSS 265
#define FLUXP 266
#define FLUXM 267
#define NUFLUX 268
#define SYS_ON_FUNCTION 269
#define SYS_OFF_FUNCTION 270
#define SYS_MULTIEX_ERRORS 271
#define GRP 272
#define GID 273
#define FNAME 274
#define VERS 275
#define SIGNAL 276
#define BG 277
#define ENERGY 278
#define CHANNEL 279
#define NDEF 280
#define GRPOPEN 281
#define GRPCLOSE 282
#define PM 283
#define FLAVOR 284
#define NOGLOBES 285
#define RULESEP 286
#define RULEMULT 287
#define NAME 288
#define RDF 289
#define ENDEXP 290
#define ENDDET 291
#define NEG 292

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED

union YYSTYPE
{
#line 1191 "glb_parser.y" /* yacc.c:355  */

  double  val;  /* For returning numbers.                   */
  double *dpt;  /* for rules */
  glb_List *ptr;
  glb_List **ptrq;
  glb_symrec  *tptr;  /* For returning symbol-table pointers      */
  char *name;
  char *iname;
  int in;
  glb_namerec *nameptr;

#line 1363 "glb_parser.c" /* yacc.c:355  */
};

typedef union YYSTYPE YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif


extern YYSTYPE yylval;

int yyparse (void);

#endif /* !YY_YY_Y_TAB_H_INCLUDED  */

/* Copy the second part of user declarations.  */

#line 1380 "glb_parser.c" /* yacc.c:358  */

#ifdef short
# undef short
#endif

#ifdef YYTYPE_UINT8
typedef YYTYPE_UINT8 yytype_uint8;
#else
typedef unsigned char yytype_uint8;
#endif

#ifdef YYTYPE_INT8
typedef YYTYPE_INT8 yytype_int8;
#else
typedef signed char yytype_int8;
#endif

#ifdef YYTYPE_UINT16
typedef YYTYPE_UINT16 yytype_uint16;
#else
typedef unsigned short int yytype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 yytype_int16;
#else
typedef short int yytype_int16;
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif ! defined YYSIZE_T
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned int
# endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T) -1)

#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(Msgid) dgettext ("bison-runtime", Msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(Msgid) Msgid
# endif
#endif

#ifndef YY_ATTRIBUTE
# if (defined __GNUC__                                               \
      && (2 < __GNUC__ || (__GNUC__ == 2 && 96 <= __GNUC_MINOR__)))  \
     || defined __SUNPRO_C && 0x5110 <= __SUNPRO_C
#  define YY_ATTRIBUTE(Spec) __attribute__(Spec)
# else
#  define YY_ATTRIBUTE(Spec) /* empty */
# endif
#endif

#ifndef YY_ATTRIBUTE_PURE
# define YY_ATTRIBUTE_PURE   YY_ATTRIBUTE ((__pure__))
#endif

#ifndef YY_ATTRIBUTE_UNUSED
# define YY_ATTRIBUTE_UNUSED YY_ATTRIBUTE ((__unused__))
#endif

#if !defined _Noreturn \
     && (!defined __STDC_VERSION__ || __STDC_VERSION__ < 201112)
# if defined _MSC_VER && 1200 <= _MSC_VER
#  define _Noreturn __declspec (noreturn)
# else
#  define _Noreturn YY_ATTRIBUTE ((__noreturn__))
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(E) ((void) (E))
#else
# define YYUSE(E) /* empty */
#endif

#if defined __GNUC__ && 407 <= __GNUC__ * 100 + __GNUC_MINOR__
/* Suppress an incorrect diagnostic about yylval being uninitialized.  */
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN \
    _Pragma ("GCC diagnostic push") \
    _Pragma ("GCC diagnostic ignored \"-Wuninitialized\"")\
    _Pragma ("GCC diagnostic ignored \"-Wmaybe-uninitialized\"")
# define YY_IGNORE_MAYBE_UNINITIALIZED_END \
    _Pragma ("GCC diagnostic pop")
#else
# define YY_INITIAL_VALUE(Value) Value
#endif
#ifndef YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_END
#endif
#ifndef YY_INITIAL_VALUE
# define YY_INITIAL_VALUE(Value) /* Nothing. */
#endif


#if ! defined yyoverflow || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined EXIT_SUCCESS
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
      /* Use EXIT_SUCCESS as a witness for stdlib.h.  */
#     ifndef EXIT_SUCCESS
#      define EXIT_SUCCESS 0
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's 'empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (0)
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined EXIT_SUCCESS \
       && ! ((defined YYMALLOC || defined malloc) \
             && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef EXIT_SUCCESS
#    define EXIT_SUCCESS 0
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined EXIT_SUCCESS
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined EXIT_SUCCESS
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* ! defined yyoverflow || YYERROR_VERBOSE */


#if (! defined yyoverflow \
     && (! defined __cplusplus \
         || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yytype_int16 yyss_alloc;
  YYSTYPE yyvs_alloc;
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (yytype_int16) + sizeof (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

# define YYCOPY_NEEDED 1

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack_alloc, Stack)                           \
    do                                                                  \
      {                                                                 \
        YYSIZE_T yynewbytes;                                            \
        YYCOPY (&yyptr->Stack_alloc, Stack, yysize);                    \
        Stack = &yyptr->Stack_alloc;                                    \
        yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
        yyptr += yynewbytes / sizeof (*yyptr);                          \
      }                                                                 \
    while (0)

#endif

#if defined YYCOPY_NEEDED && YYCOPY_NEEDED
/* Copy COUNT objects from SRC to DST.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(Dst, Src, Count) \
      __builtin_memcpy (Dst, Src, (Count) * sizeof (*(Src)))
#  else
#   define YYCOPY(Dst, Src, Count)              \
      do                                        \
        {                                       \
          YYSIZE_T yyi;                         \
          for (yyi = 0; yyi < (Count); yyi++)   \
            (Dst)[yyi] = (Src)[yyi];            \
        }                                       \
      while (0)
#  endif
# endif
#endif /* !YYCOPY_NEEDED */

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  38
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   355

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  52
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  24
/* YYNRULES -- Number of rules.  */
#define YYNRULES  78
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  167

/* YYTRANSLATE[YYX] -- Symbol number corresponding to YYX as returned
   by yylex, with out-of-bounds checking.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   292

#define YYTRANSLATE(YYX)                                                \
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[TOKEN-NUM] -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex, without out-of-bounds checking.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
      49,    50,    44,    43,    37,    42,     2,    45,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,    51,
       2,    38,    39,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,    47,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    41,     2,    40,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,    48,     2,     2,
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
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34,
      35,    36,    46
};

#if YYDEBUG
  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,  1250,  1250,  1251,  1252,  1253,  1257,  1258,  1259,  1260,
    1261,  1265,  1266,  1267,  1268,  1269,  1274,  1278,  1283,  1289,
    1290,  1291,  1292,  1293,  1294,  1295,  1296,  1297,  1301,  1310,
    1315,  1319,  1320,  1321,  1322,  1327,  1328,  1329,  1330,  1331,
    1332,  1333,  1337,  1348,  1347,  1358,  1365,  1366,  1370,  1371,
    1372,  1373,  1374,  1375,  1376,  1377,  1381,  1391,  1401,  1413,
    1425,  1444,  1445,  1449,  1450,  1451,  1456,  1464,  1478,  1482,
    1489,  1498,  1509,  1519,  1530,  1539,  1548,  1554,  1560
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || 0
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "NUM", "SFNCT", "BOGUS", "LVAR", "VAR",
  "FNCT", "IDN", "CROSS", "FLUXP", "FLUXM", "NUFLUX", "SYS_ON_FUNCTION",
  "SYS_OFF_FUNCTION", "SYS_MULTIEX_ERRORS", "GRP", "GID", "FNAME", "VERS",
  "SIGNAL", "BG", "ENERGY", "CHANNEL", "NDEF", "GRPOPEN", "GRPCLOSE", "PM",
  "FLAVOR", "NOGLOBES", "RULESEP", "RULEMULT", "NAME", "RDF", "ENDEXP",
  "ENDDET", "','", "'='", "'>'", "'}'", "'{'", "'-'", "'+'", "'*'", "'/'",
  "NEG", "'^'", "'\\247'", "'('", "')'", "';'", "$accept", "input",
  "topleveldirective", "exp", "listcopy", "seq", "list", "rulepart",
  "group", "$@1", "ingroup", "ingroup_statement", "version", "cross",
  "flux", "nuflux", "channel", "name", "pm", "ene", "energy", "brule",
  "srule", "rule", YY_NULLPTR
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[NUM] -- (External) token number corresponding to the
   (internal) symbol number NUM (which must be that of a token).  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283,   284,
     285,   286,   287,   288,   289,   290,   291,    44,    61,    62,
     125,   123,    45,    43,    42,    47,   292,    94,   167,    40,
      41,    59
};
# endif

#define YYPACT_NINF -73

#define yypact_value_is_default(Yystate) \
  (!!((Yystate) == (-73)))

#define YYTABLE_NINF -1

#define yytable_value_is_error(Yytable_value) \
  0

  /* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
     STATE-NUM.  */
static const yytype_int16 yypact[] =
{
     216,   -73,   -32,   -28,   -22,   -30,     7,    21,   -73,   -73,
     -73,   -73,   -73,   259,   306,   -73,   306,    75,   -73,    23,
     -73,   -73,   -73,   -73,    12,   306,   106,    25,   236,     5,
      31,    20,    36,   -73,   143,   -25,    30,   280,   -73,   -73,
     306,   306,   306,   306,   306,    49,   -27,   -73,    23,   -73,
      90,    -3,    -1,    74,   -73,   -12,   -73,    35,    51,   -73,
     306,   286,   306,   -73,   306,   -73,   -73,     1,     1,    30,
      30,    30,   106,    12,   -73,   -73,    12,   -73,   -73,   306,
     -73,    76,   299,    23,    23,    99,    67,    23,    79,   -73,
     306,   -73,   151,   260,   188,    69,    80,    81,    82,    83,
      84,    85,    87,    91,    92,   -73,    23,   -73,   -73,   -73,
     -73,   -73,   -73,   -73,   107,   114,   -73,    12,   -73,   130,
     132,   133,   134,   144,    12,   306,   306,    12,     4,   306,
     306,   118,   -73,   -73,   -73,   -73,   -73,   -73,   138,    44,
     -73,   -73,   -29,   -73,   -73,   139,   -73,   -73,   -73,    12,
     306,   -73,   -19,   -73,    23,   -73,   -73,   -73,   146,   150,
     158,   152,   174,     4,   175,     4,   -73
};

  /* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
     Performed when YYTABLE does not specify something else to do.  Zero
     means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       0,    11,    38,    13,     0,     0,     0,     0,    27,     4,
      12,     9,    10,     0,     0,     5,     0,     0,     2,     7,
      41,     8,     6,    26,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    31,     0,     0,    23,     0,     1,     3,
       0,     0,     0,     0,     0,     0,     0,    39,    14,    37,
       0,     0,     0,     0,    16,    15,    34,     0,     0,    56,
       0,     0,     0,    33,     0,    32,    25,    20,    19,    21,
      22,    24,     0,     0,    18,    35,     0,    36,    28,     0,
      43,     0,     0,    29,    30,     0,     0,    17,     0,    46,
       0,    46,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    45,    48,    49,    47,    53,
      54,    55,    52,    51,    74,    75,    50,     0,    44,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,    57,    58,    59,    76,    77,    66,    78,     0,
      72,    70,    68,    62,    61,     0,    71,    73,    40,     0,
       0,    69,     0,    67,    42,    63,    65,    64,     0,     0,
       0,     0,     0,     0,     0,     0,    60
};

  /* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
     -73,   -73,   165,     0,   -73,   194,   -21,   -72,   -73,   -73,
     123,   -73,   -73,   -73,   -73,   -73,   -73,   -66,   -73,    89,
     -73,   -73,   -73,   -73
};

  /* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,    17,    18,   139,    20,    51,    21,   140,    22,    88,
      92,   108,    23,   109,   110,   111,   112,   145,   158,   138,
     113,   114,   115,   116
};

  /* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule whose
     number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_uint8 yytable[] =
{
      19,    27,   149,    47,    27,    52,    24,    56,    28,   155,
      25,    73,    64,    34,    36,    65,    37,    19,     2,    79,
      45,    46,   151,   156,   157,    48,    50,    26,    55,   143,
      40,    41,    42,    43,    64,    44,    76,   144,    57,    58,
      67,    68,    69,    70,    71,    42,    43,    75,    44,    77,
      59,    52,    56,    13,   141,    86,    29,   146,   147,    30,
      82,    55,    83,    53,    84,    40,    41,    42,    43,    60,
      44,   107,    85,   107,    61,    38,   150,    44,     1,    87,
      78,     2,     3,     4,     5,    80,    40,    41,    42,    43,
      93,    44,   106,     6,   106,     7,   131,   164,    72,   166,
       8,    81,    89,   137,    90,    91,   137,   119,    10,     1,
      11,    12,     2,     3,     4,     5,    13,    14,   120,   121,
     122,   123,   124,   125,    16,   126,     7,    62,   153,   127,
     128,     8,    40,    41,    42,    43,    62,    44,   129,    10,
      74,    40,    41,    42,    43,   130,    44,    13,    14,   132,
     154,   133,   134,   135,     1,    16,    49,     2,     3,     4,
       5,    95,    96,   136,    97,    98,    99,   100,   148,   149,
     152,     7,   101,   102,   103,   104,     8,   159,   105,   160,
      62,   162,    39,    63,    10,    40,    41,    42,    43,   161,
      44,     1,    13,    14,     2,     3,     4,     5,    95,    96,
      16,    97,    98,    99,   100,   163,   165,    35,     7,   101,
     102,   103,   104,     8,    94,   118,   142,     0,     0,     1,
       0,    10,     2,     3,     4,     5,     0,     0,     0,    13,
      14,     0,     0,     0,     6,     0,     7,    16,     0,     1,
      54,     8,     2,     3,     4,     5,     9,     0,     0,    10,
       0,    11,    12,     0,     0,     0,     7,    13,    14,     0,
       0,     8,     1,     0,    15,    16,     3,    31,    32,    10,
       0,     0,     0,     0,     0,     0,     0,    13,    14,     7,
       0,     0,     0,     0,     8,    16,     0,     0,     0,     1,
      54,     0,    10,     3,    31,    32,     0,   117,     0,    33,
       0,    14,    40,    41,    42,    43,     7,    44,    16,     1,
       0,     8,     0,     3,    31,    32,     0,     0,     0,    10,
       0,     0,    40,    41,    42,    43,     7,    44,    14,     0,
      66,     8,     0,     0,     0,    16,     0,     0,     0,    10,
       0,    40,    41,    42,    43,     0,    44,     0,    14,    74,
       0,     0,     0,     0,     0,    16
};

static const yytype_int16 yycheck[] =
{
       0,    31,    31,    24,    31,    26,    38,    28,    38,    28,
      38,    38,    37,    13,    14,    40,    16,    17,     6,    31,
       8,     9,    51,    42,    43,    25,    26,    49,    28,    25,
      42,    43,    44,    45,    37,    47,    37,    33,    33,    34,
      40,    41,    42,    43,    44,    44,    45,    50,    47,    50,
      19,    72,    73,    41,   126,    76,    49,   129,   130,    38,
      60,    61,    62,    38,    64,    42,    43,    44,    45,    49,
      47,    92,    72,    94,    38,     0,    32,    47,     3,    79,
       6,     6,     7,     8,     9,    50,    42,    43,    44,    45,
      90,    47,    92,    18,    94,    20,   117,   163,    49,   165,
      25,    50,    26,   124,    37,    26,   127,    38,    33,     3,
      35,    36,     6,     7,     8,     9,    41,    42,    38,    38,
      38,    38,    38,    38,    49,    38,    20,    37,   149,    38,
      38,    25,    42,    43,    44,    45,    37,    47,    31,    33,
      50,    42,    43,    44,    45,    31,    47,    41,    42,    19,
     150,    19,    19,    19,     3,    49,    50,     6,     7,     8,
       9,    10,    11,    19,    13,    14,    15,    16,    50,    31,
      31,    20,    21,    22,    23,    24,    25,    31,    27,    29,
      37,    29,    17,    40,    33,    42,    43,    44,    45,    31,
      47,     3,    41,    42,     6,     7,     8,     9,    10,    11,
      49,    13,    14,    15,    16,    31,    31,    13,    20,    21,
      22,    23,    24,    25,    91,    27,   127,    -1,    -1,     3,
      -1,    33,     6,     7,     8,     9,    -1,    -1,    -1,    41,
      42,    -1,    -1,    -1,    18,    -1,    20,    49,    -1,     3,
       4,    25,     6,     7,     8,     9,    30,    -1,    -1,    33,
      -1,    35,    36,    -1,    -1,    -1,    20,    41,    42,    -1,
      -1,    25,     3,    -1,    48,    49,     7,     8,     9,    33,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    41,    42,    20,
      -1,    -1,    -1,    -1,    25,    49,    -1,    -1,    -1,     3,
       4,    -1,    33,     7,     8,     9,    -1,    37,    -1,    40,
      -1,    42,    42,    43,    44,    45,    20,    47,    49,     3,
      -1,    25,    -1,     7,     8,     9,    -1,    -1,    -1,    33,
      -1,    -1,    42,    43,    44,    45,    20,    47,    42,    -1,
      50,    25,    -1,    -1,    -1,    49,    -1,    -1,    -1,    33,
      -1,    42,    43,    44,    45,    -1,    47,    -1,    42,    50,
      -1,    -1,    -1,    -1,    -1,    49
};

  /* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,     3,     6,     7,     8,     9,    18,    20,    25,    30,
      33,    35,    36,    41,    42,    48,    49,    53,    54,    55,
      56,    58,    60,    64,    38,    38,    49,    31,    38,    49,
      38,     8,     9,    40,    55,    57,    55,    55,     0,    54,
      42,    43,    44,    45,    47,     8,     9,    58,    55,    50,
      55,    57,    58,    38,     4,    55,    58,    33,    34,    19,
      49,    38,    37,    40,    37,    40,    50,    55,    55,    55,
      55,    55,    49,    38,    50,    50,    37,    50,     6,    31,
      50,    50,    55,    55,    55,    55,    58,    55,    61,    26,
      37,    26,    62,    55,    62,    10,    11,    13,    14,    15,
      16,    21,    22,    23,    24,    27,    55,    58,    63,    65,
      66,    67,    68,    72,    73,    74,    75,    37,    27,    38,
      38,    38,    38,    38,    38,    38,    38,    38,    38,    31,
      31,    58,    19,    19,    19,    19,    19,    58,    71,    55,
      59,    59,    71,    25,    33,    69,    59,    59,    50,    31,
      32,    51,    31,    58,    55,    28,    42,    43,    70,    31,
      29,    31,    29,    31,    69,    31,    69
};

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    52,    53,    53,    53,    53,    54,    54,    54,    54,
      54,    55,    55,    55,    55,    55,    55,    55,    55,    55,
      55,    55,    55,    55,    55,    55,    55,    55,    56,    57,
      57,    58,    58,    58,    58,    58,    58,    58,    58,    58,
      58,    58,    59,    61,    60,    60,    62,    62,    63,    63,
      63,    63,    63,    63,    63,    63,    64,    65,    66,    67,
      68,    69,    69,    70,    70,    70,    71,    71,    72,    72,
      73,    73,    74,    74,    75,    75,    75,    75,    75
};

  /* YYR2[YYN] -- Number of symbols on the right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     1,     2,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     3,     3,     3,     5,     4,     3,
       3,     3,     3,     2,     3,     3,     1,     1,     4,     3,
       3,     2,     3,     3,     3,     4,     4,     3,     1,     3,
      10,     1,     3,     0,     8,     7,     0,     2,     1,     1,
       1,     1,     1,     1,     1,     1,     3,     3,     3,     3,
      13,     1,     1,     1,     1,     1,     1,     3,     3,     4,
       3,     3,     3,     3,     1,     1,     3,     3,     3
};


#define yyerrok         (yyerrstatus = 0)
#define yyclearin       (yychar = YYEMPTY)
#define YYEMPTY         (-2)
#define YYEOF           0

#define YYACCEPT        goto yyacceptlab
#define YYABORT         goto yyabortlab
#define YYERROR         goto yyerrorlab


#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)                                  \
do                                                              \
  if (yychar == YYEMPTY)                                        \
    {                                                           \
      yychar = (Token);                                         \
      yylval = (Value);                                         \
      YYPOPSTACK (yylen);                                       \
      yystate = *yyssp;                                         \
      goto yybackup;                                            \
    }                                                           \
  else                                                          \
    {                                                           \
      yyerror (YY_("syntax error: cannot back up")); \
      YYERROR;                                                  \
    }                                                           \
while (0)

/* Error token number */
#define YYTERROR        1
#define YYERRCODE       256



/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)                        \
do {                                            \
  if (yydebug)                                  \
    YYFPRINTF Args;                             \
} while (0)

/* This macro is provided for backward compatibility. */
#ifndef YY_LOCATION_PRINT
# define YY_LOCATION_PRINT(File, Loc) ((void) 0)
#endif


# define YY_SYMBOL_PRINT(Title, Type, Value, Location)                    \
do {                                                                      \
  if (yydebug)                                                            \
    {                                                                     \
      YYFPRINTF (stderr, "%s ", Title);                                   \
      yy_symbol_print (stderr,                                            \
                  Type, Value); \
      YYFPRINTF (stderr, "\n");                                           \
    }                                                                     \
} while (0)


/*----------------------------------------.
| Print this symbol's value on YYOUTPUT.  |
`----------------------------------------*/

static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
{
  FILE *yyo = yyoutput;
  YYUSE (yyo);
  if (!yyvaluep)
    return;
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# endif
  YYUSE (yytype);
}


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

static void
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
{
  YYFPRINTF (yyoutput, "%s %s (",
             yytype < YYNTOKENS ? "token" : "nterm", yytname[yytype]);

  yy_symbol_value_print (yyoutput, yytype, yyvaluep);
  YYFPRINTF (yyoutput, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

static void
yy_stack_print (yytype_int16 *yybottom, yytype_int16 *yytop)
{
  YYFPRINTF (stderr, "Stack now");
  for (; yybottom <= yytop; yybottom++)
    {
      int yybot = *yybottom;
      YYFPRINTF (stderr, " %d", yybot);
    }
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)                            \
do {                                                            \
  if (yydebug)                                                  \
    yy_stack_print ((Bottom), (Top));                           \
} while (0)


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

static void
yy_reduce_print (yytype_int16 *yyssp, YYSTYPE *yyvsp, int yyrule)
{
  unsigned long int yylno = yyrline[yyrule];
  int yynrhs = yyr2[yyrule];
  int yyi;
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
             yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      YYFPRINTF (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr,
                       yystos[yyssp[yyi + 1 - yynrhs]],
                       &(yyvsp[(yyi + 1) - (yynrhs)])
                                              );
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)          \
do {                                    \
  if (yydebug)                          \
    yy_reduce_print (yyssp, yyvsp, Rule); \
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
#ifndef YYINITDEPTH
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
#  if defined __GLIBC__ && defined _STRING_H
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
static YYSIZE_T
yystrlen (const char *yystr)
{
  YYSIZE_T yylen;
  for (yylen = 0; yystr[yylen]; yylen++)
    continue;
  return yylen;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
static char *
yystpcpy (char *yydest, const char *yysrc)
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
      YYSIZE_T yyn = 0;
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

/* Copy into *YYMSG, which is of size *YYMSG_ALLOC, an error message
   about the unexpected token YYTOKEN for the state stack whose top is
   YYSSP.

   Return 0 if *YYMSG was successfully written.  Return 1 if *YYMSG is
   not large enough to hold the message.  In that case, also set
   *YYMSG_ALLOC to the required number of bytes.  Return 2 if the
   required number of bytes is too large to store.  */
static int
yysyntax_error (YYSIZE_T *yymsg_alloc, char **yymsg,
                yytype_int16 *yyssp, int yytoken)
{
  YYSIZE_T yysize0 = yytnamerr (YY_NULLPTR, yytname[yytoken]);
  YYSIZE_T yysize = yysize0;
  enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
  /* Internationalized format string. */
  const char *yyformat = YY_NULLPTR;
  /* Arguments of yyformat. */
  char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
  /* Number of reported tokens (one for the "unexpected", one per
     "expected"). */
  int yycount = 0;

  /* There are many possibilities here to consider:
     - If this state is a consistent state with a default action, then
       the only way this function was invoked is if the default action
       is an error action.  In that case, don't check for expected
       tokens because there are none.
     - The only way there can be no lookahead present (in yychar) is if
       this state is a consistent state with a default action.  Thus,
       detecting the absence of a lookahead is sufficient to determine
       that there is no unexpected or expected token to report.  In that
       case, just report a simple "syntax error".
     - Don't assume there isn't a lookahead just because this state is a
       consistent state with a default action.  There might have been a
       previous inconsistent state, consistent state with a non-default
       action, or user semantic action that manipulated yychar.
     - Of course, the expected token list depends on states to have
       correct lookahead information, and it depends on the parser not
       to perform extra reductions after fetching a lookahead from the
       scanner and before detecting a syntax error.  Thus, state merging
       (from LALR or IELR) and default reductions corrupt the expected
       token list.  However, the list is correct for canonical LR with
       one exception: it will still contain any token that will not be
       accepted due to an error action in a later state.
  */
  if (yytoken != YYEMPTY)
    {
      int yyn = yypact[*yyssp];
      yyarg[yycount++] = yytname[yytoken];
      if (!yypact_value_is_default (yyn))
        {
          /* Start YYX at -YYN if negative to avoid negative indexes in
             YYCHECK.  In other words, skip the first -YYN actions for
             this state because they are default actions.  */
          int yyxbegin = yyn < 0 ? -yyn : 0;
          /* Stay within bounds of both yycheck and yytname.  */
          int yychecklim = YYLAST - yyn + 1;
          int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
          int yyx;

          for (yyx = yyxbegin; yyx < yyxend; ++yyx)
            if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR
                && !yytable_value_is_error (yytable[yyx + yyn]))
              {
                if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
                  {
                    yycount = 1;
                    yysize = yysize0;
                    break;
                  }
                yyarg[yycount++] = yytname[yyx];
                {
                  YYSIZE_T yysize1 = yysize + yytnamerr (YY_NULLPTR, yytname[yyx]);
                  if (! (yysize <= yysize1
                         && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
                    return 2;
                  yysize = yysize1;
                }
              }
        }
    }

  switch (yycount)
    {
# define YYCASE_(N, S)                      \
      case N:                               \
        yyformat = S;                       \
      break
      YYCASE_(0, YY_("syntax error"));
      YYCASE_(1, YY_("syntax error, unexpected %s"));
      YYCASE_(2, YY_("syntax error, unexpected %s, expecting %s"));
      YYCASE_(3, YY_("syntax error, unexpected %s, expecting %s or %s"));
      YYCASE_(4, YY_("syntax error, unexpected %s, expecting %s or %s or %s"));
      YYCASE_(5, YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s"));
# undef YYCASE_
    }

  {
    YYSIZE_T yysize1 = yysize + yystrlen (yyformat);
    if (! (yysize <= yysize1 && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
      return 2;
    yysize = yysize1;
  }

  if (*yymsg_alloc < yysize)
    {
      *yymsg_alloc = 2 * yysize;
      if (! (yysize <= *yymsg_alloc
             && *yymsg_alloc <= YYSTACK_ALLOC_MAXIMUM))
        *yymsg_alloc = YYSTACK_ALLOC_MAXIMUM;
      return 1;
    }

  /* Avoid sprintf, as that infringes on the user's name space.
     Don't have undefined behavior even if the translation
     produced a string with the wrong number of "%s"s.  */
  {
    char *yyp = *yymsg;
    int yyi = 0;
    while ((*yyp = *yyformat) != '\0')
      if (*yyp == '%' && yyformat[1] == 's' && yyi < yycount)
        {
          yyp += yytnamerr (yyp, yyarg[yyi++]);
          yyformat += 2;
        }
      else
        {
          yyp++;
          yyformat++;
        }
  }
  return 0;
}
#endif /* YYERROR_VERBOSE */

/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep)
{
  YYUSE (yyvaluep);
  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  YYUSE (yytype);
  YY_IGNORE_MAYBE_UNINITIALIZED_END
}




/* The lookahead symbol.  */
int yychar;

/* The semantic value of the lookahead symbol.  */
YYSTYPE yylval;
/* Number of syntax errors so far.  */
int yynerrs;


/*----------.
| yyparse.  |
`----------*/

int
yyparse (void)
{
    int yystate;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus;

    /* The stacks and their tools:
       'yyss': related to states.
       'yyvs': related to semantic values.

       Refer to the stacks through separate pointers, to allow yyoverflow
       to reallocate them elsewhere.  */

    /* The state stack.  */
    yytype_int16 yyssa[YYINITDEPTH];
    yytype_int16 *yyss;
    yytype_int16 *yyssp;

    /* The semantic value stack.  */
    YYSTYPE yyvsa[YYINITDEPTH];
    YYSTYPE *yyvs;
    YYSTYPE *yyvsp;

    YYSIZE_T yystacksize;

  int yyn;
  int yyresult;
  /* Lookahead token as an internal (translated) token number.  */
  int yytoken = 0;
  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;

#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  yyssp = yyss = yyssa;
  yyvsp = yyvs = yyvsa;
  yystacksize = YYINITDEPTH;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY; /* Cause a token to be read.  */
  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
        /* Give user a chance to reallocate the stack.  Use copies of
           these so that the &'s don't force the real ones into
           memory.  */
        YYSTYPE *yyvs1 = yyvs;
        yytype_int16 *yyss1 = yyss;

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
        yytype_int16 *yyss1 = yyss;
        union yyalloc *yyptr =
          (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
        if (! yyptr)
          goto yyexhaustedlab;
        YYSTACK_RELOCATE (yyss_alloc, yyss);
        YYSTACK_RELOCATE (yyvs_alloc, yyvs);
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

  if (yystate == YYFINAL)
    YYACCEPT;

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

  /* Do appropriate processing given the current state.  Read a
     lookahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to lookahead token.  */
  yyn = yypact[yystate];
  if (yypact_value_is_default (yyn))
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid lookahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = yylex ();
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
      if (yytable_value_is_error (yyn))
        goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the lookahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

  /* Discard the shifted token.  */
  yychar = YYEMPTY;

  yystate = yyn;
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END

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
     '$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
        case 2:
#line 1250 "glb_parser.y" /* yacc.c:1646  */
    {}
#line 2606 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 3:
#line 1251 "glb_parser.y" /* yacc.c:1646  */
    {}
#line 2612 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 4:
#line 1252 "glb_parser.y" /* yacc.c:1646  */
    {YYABORT;}
#line 2618 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 5:
#line 1253 "glb_parser.y" /* yacc.c:1646  */
    {YYABORT;}
#line 2624 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 6:
#line 1257 "glb_parser.y" /* yacc.c:1646  */
    {}
#line 2630 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 7:
#line 1258 "glb_parser.y" /* yacc.c:1646  */
    {}
#line 2636 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 8:
#line 1259 "glb_parser.y" /* yacc.c:1646  */
    {}
#line 2642 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 9:
#line 1260 "glb_parser.y" /* yacc.c:1646  */
    { glb_copy_buff();  glbReset(); }
#line 2648 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 10:
#line 1261 "glb_parser.y" /* yacc.c:1646  */
    { glbNewDetector(); }
#line 2654 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 11:
#line 1265 "glb_parser.y" /* yacc.c:1646  */
    { (yyval.val) = (yyvsp[0].val);                     }
#line 2660 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 12:
#line 1266 "glb_parser.y" /* yacc.c:1646  */
    { (yyval.val) = (yyvsp[0].nameptr)->value;              }
#line 2666 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 13:
#line 1267 "glb_parser.y" /* yacc.c:1646  */
    { (yyval.val) = (yyvsp[0].tptr)->value.var;          }
#line 2672 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 14:
#line 1268 "glb_parser.y" /* yacc.c:1646  */
    { (yyval.val) = (yyvsp[0].val); (yyvsp[-2].tptr)->value.var = (yyvsp[0].val); }
#line 2678 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 15:
#line 1269 "glb_parser.y" /* yacc.c:1646  */
    {
  if(set_exp((yyvsp[-2].name),(yyvsp[0].val),0)==1) yyerror("Unknown identifier: %s", (yyvsp[-2].name));
  (yyval.val) = (yyvsp[0].val);
  if ((yyvsp[-2].name))  { glb_free((yyvsp[-2].name));  (yyvsp[-2].name)=NULL; }
}
#line 2688 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 16:
#line 1274 "glb_parser.y" /* yacc.c:1646  */
    {
  if(set_fnct((yyvsp[-2].name),(yyvsp[0].nameptr)->sf)==1) yyerror("Unknown identifier: %s", (yyvsp[-2].name));
  if ((yyvsp[-2].name))  { glb_free((yyvsp[-2].name));  (yyvsp[-2].name)=NULL; }
}
#line 2697 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 17:
#line 1278 "glb_parser.y" /* yacc.c:1646  */
    {
  if(set_pair((yyvsp[-4].name),(yyvsp[-2].val),(yyvsp[0].val),0)==1) yyerror("Unknown identifier: %s", (yyvsp[-4].name));
  (yyval.val) = (yyvsp[-2].val);
  if ((yyvsp[-4].name))  { glb_free((yyvsp[-4].name));  (yyvsp[-4].name)=NULL; }
}
#line 2707 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 18:
#line 1283 "glb_parser.y" /* yacc.c:1646  */
    {
  /* added safety in case the function pointer is NULL, which is
     sometimes useful for special functions */
  if((yyvsp[-3].tptr)->value.fnctptr==NULL) yyerror("Improper use of special function %s", (yyvsp[-3].tptr)->name);
  else (yyval.val) = (*((yyvsp[-3].tptr)->value.fnctptr))((yyvsp[-1].val)); }
#line 2717 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 19:
#line 1289 "glb_parser.y" /* yacc.c:1646  */
    { (yyval.val) = (yyvsp[-2].val) + (yyvsp[0].val);      }
#line 2723 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 20:
#line 1290 "glb_parser.y" /* yacc.c:1646  */
    { (yyval.val) = (yyvsp[-2].val) - (yyvsp[0].val);      }
#line 2729 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 21:
#line 1291 "glb_parser.y" /* yacc.c:1646  */
    { (yyval.val) = (yyvsp[-2].val) * (yyvsp[0].val);      }
#line 2735 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 22:
#line 1292 "glb_parser.y" /* yacc.c:1646  */
    { (yyval.val) = (yyvsp[-2].val) / (yyvsp[0].val);      }
#line 2741 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 23:
#line 1293 "glb_parser.y" /* yacc.c:1646  */
    { (yyval.val) = -(yyvsp[0].val);          }
#line 2747 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 24:
#line 1294 "glb_parser.y" /* yacc.c:1646  */
    { (yyval.val) = pow ((yyvsp[-2].val), (yyvsp[0].val)); }
#line 2753 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 25:
#line 1295 "glb_parser.y" /* yacc.c:1646  */
    { (yyval.val) = (yyvsp[-1].val);           }
#line 2759 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 26:
#line 1296 "glb_parser.y" /* yacc.c:1646  */
    { (yyval.val) = 0;            }
#line 2765 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 27:
#line 1297 "glb_parser.y" /* yacc.c:1646  */
    { yyerror("Unknown name: %s", (yyvsp[0].name)); YYERROR; }
#line 2771 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 28:
#line 1301 "glb_parser.y" /* yacc.c:1646  */
    {
  glb_List *ltemp;
  ltemp=thread_list(&glb_list_copy,0,0,(yyvsp[0].tptr)->list);
  if(set_exp_list((yyvsp[-3].name),ltemp,3)==1) yyerror("Unknown identifier");
  (yyval.ptr) = ltemp;
  if ((yyvsp[-3].name))  { glb_free((yyvsp[-3].name));  (yyvsp[-3].name)=NULL; }
}
#line 2783 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 29:
#line 1310 "glb_parser.y" /* yacc.c:1646  */
    {
   glb_List *buf = list_cons(NULL, (yyvsp[-2].val));
   buf = list_cons(buf, (yyvsp[0].val));
   (yyval.ptr)  = buf;
}
#line 2793 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 30:
#line 1315 "glb_parser.y" /* yacc.c:1646  */
    { (yyval.ptr) = list_cons((yyvsp[-2].ptr), (yyvsp[0].val)); }
#line 2799 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 31:
#line 1319 "glb_parser.y" /* yacc.c:1646  */
    {(yyval.ptr)=NULL;}
#line 2805 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 32:
#line 1320 "glb_parser.y" /* yacc.c:1646  */
    {(yyval.ptr)=(yyvsp[-1].ptr); }
#line 2811 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 33:
#line 1321 "glb_parser.y" /* yacc.c:1646  */
    {(yyval.ptr)=list_cons(NULL,(yyvsp[-1].val)); }
#line 2817 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 34:
#line 1322 "glb_parser.y" /* yacc.c:1646  */
    {
  if(set_exp_list((yyvsp[-2].name),(yyvsp[0].ptr),3)==1)  yyerror("Unknown identifier");
  (yyval.ptr) = (yyvsp[0].ptr);
  if ((yyvsp[-2].name))  { glb_free((yyvsp[-2].name));  (yyvsp[-2].name)=NULL; }
}
#line 2827 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 35:
#line 1327 "glb_parser.y" /* yacc.c:1646  */
    {(yyval.ptr) = thread_list((yyvsp[-3].tptr)->value.fnctptr,(yyvsp[-3].tptr)->reverse,(yyvsp[-3].tptr)->destroy,(yyvsp[-1].ptr));}
#line 2833 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 36:
#line 1328 "glb_parser.y" /* yacc.c:1646  */
    {(yyval.ptr) = thread_list((yyvsp[-3].tptr)->value.fnctptr,(yyvsp[-3].tptr)->reverse,(yyvsp[-3].tptr)->destroy,(yyvsp[-1].ptr));}
#line 2839 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 37:
#line 1329 "glb_parser.y" /* yacc.c:1646  */
    {(yyval.ptr) = (*((yyvsp[-2].tptr)->value.lfnctptr))();}
#line 2845 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 38:
#line 1330 "glb_parser.y" /* yacc.c:1646  */
    { (yyval.ptr) = (yyvsp[0].tptr)->list;              }
#line 2851 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 39:
#line 1331 "glb_parser.y" /* yacc.c:1646  */
    { (yyval.ptr) = (yyvsp[0].ptr); (yyvsp[-2].tptr)->list = (yyvsp[0].ptr); }
#line 2857 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 40:
#line 1332 "glb_parser.y" /* yacc.c:1646  */
    {(yyval.ptr)=glb_interpolation((yyvsp[-7].ptr),(yyvsp[-5].ptr),floor((yyvsp[-3].val)),(yyvsp[-1].ptr));}
#line 2863 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 41:
#line 1333 "glb_parser.y" /* yacc.c:1646  */
    {}
#line 2869 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 42:
#line 1337 "glb_parser.y" /* yacc.c:1646  */
    {
  double *buf;
  buf=(double*) glb_malloc(sizeof(double)*2);
  buf[0]=(yyvsp[-2].val);
  buf[1]=(yyvsp[0].val)-1;
  (yyval.dpt)=buf;
}
#line 2881 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 43:
#line 1348 "glb_parser.y" /* yacc.c:1646  */
    { if((yyvsp[-1].nameptr)->value==-1) {(yyvsp[-1].nameptr)->value=step_counter((yyvsp[-3].name)); }
  loc_count=(yyvsp[-1].nameptr)->value;
  glb_free(context);
  context =(char *) strdup((yyvsp[-3].name));
  grp_start(context, loc_count);
  if ((yyvsp[-3].name))  { glb_free((yyvsp[-3].name));  (yyvsp[-3].name)=NULL; }
}
#line 2893 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 44:
#line 1355 "glb_parser.y" /* yacc.c:1646  */
    {
  grp_end(context, loc_count);
}
#line 2901 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 45:
#line 1358 "glb_parser.y" /* yacc.c:1646  */
    {
    yyerror("Redefinition of an automatic variable %s", (yyvsp[-4].nameptr)->name); YYERROR;
    if ((yyvsp[-6].name))  { glb_free((yyvsp[-6].name));  (yyvsp[-6].name)=NULL; }
}
#line 2910 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 48:
#line 1370 "glb_parser.y" /* yacc.c:1646  */
    {}
#line 2916 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 49:
#line 1371 "glb_parser.y" /* yacc.c:1646  */
    {}
#line 2922 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 50:
#line 1372 "glb_parser.y" /* yacc.c:1646  */
    {}
#line 2928 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 51:
#line 1373 "glb_parser.y" /* yacc.c:1646  */
    {}
#line 2934 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 52:
#line 1374 "glb_parser.y" /* yacc.c:1646  */
    {}
#line 2940 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 53:
#line 1375 "glb_parser.y" /* yacc.c:1646  */
    {}
#line 2946 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 54:
#line 1376 "glb_parser.y" /* yacc.c:1646  */
    {}
#line 2952 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 55:
#line 1377 "glb_parser.y" /* yacc.c:1646  */
    {}
#line 2958 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 56:
#line 1381 "glb_parser.y" /* yacc.c:1646  */
    {
//  buff.version=strdup($3);
  if (set_string((yyvsp[-2].name), (yyvsp[0].name)) != 0)
    yyerror("Unknown identifier: %s", (yyvsp[-2].name));
  if ((yyvsp[-2].name))  { glb_free((yyvsp[-2].name));  (yyvsp[-2].name)=NULL; }
  if ((yyvsp[0].name))  { glb_free((yyvsp[0].name));  (yyvsp[0].name)=NULL; }
}
#line 2970 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 57:
#line 1391 "glb_parser.y" /* yacc.c:1646  */
    {
  //load_cross($3,loc_count-1);
  xsc.file_name=strdup((yyvsp[0].name));
  (yyval.name)=(yyvsp[0].name);
  if ((yyvsp[-2].name))  { glb_free((yyvsp[-2].name));  (yyvsp[-2].name)=NULL; }
  if ((yyvsp[0].name))  { glb_free((yyvsp[0].name));  (yyvsp[0].name)=NULL; }
}
#line 2982 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 58:
#line 1401 "glb_parser.y" /* yacc.c:1646  */
    {
  //load_flux($3,loc_count-1,1);
  flt.file_name=strdup((yyvsp[0].name));

  //if(set_exp($1,$3,0)==1) yyerror("Unknown identifier");
  (yyval.name)=(yyvsp[0].name);
  if ((yyvsp[-2].name))  { glb_free((yyvsp[-2].name));  (yyvsp[-2].name)=NULL; }
  if ((yyvsp[0].name))  { glb_free((yyvsp[0].name));  (yyvsp[0].name)=NULL; }
}
#line 2996 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 59:
#line 1413 "glb_parser.y" /* yacc.c:1646  */
    {
  //load_flux($3,loc_count-1,1);
  flt.file_name=strdup((yyvsp[0].name));

  //if(set_exp($1,$3,0)==1) yyerror("Unknown identifier");
  (yyval.name)=(yyvsp[0].name);
  if ((yyvsp[-2].name))  { glb_free((yyvsp[-2].name));  (yyvsp[-2].name)=NULL; }
  if ((yyvsp[0].name))  { glb_free((yyvsp[0].name));  (yyvsp[0].name)=NULL; }
}
#line 3010 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 60:
#line 1426 "glb_parser.y" /* yacc.c:1646  */
    {

  int x[6];
  x[0]=(yyvsp[-10].nameptr)->value - 1 ;
  x[1]=(yyvsp[-8].in);
  x[2]=(yyvsp[-6].in);
  x[3]=(yyvsp[-4].in);
  x[4]=(int) (yyvsp[-2].nameptr)->value -1;
  x[5]=(int) (yyvsp[0].nameptr)->value -1;

  set_channel_data(x,loc_count);
  if ((yyvsp[-12].name))  { glb_free((yyvsp[-12].name));  (yyvsp[-12].name)=NULL; }
}
#line 3028 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 61:
#line 1444 "glb_parser.y" /* yacc.c:1646  */
    {(yyval.nameptr)=(yyvsp[0].nameptr);}
#line 3034 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 62:
#line 1445 "glb_parser.y" /* yacc.c:1646  */
    { yyerror("Unknown name: %s", (yyvsp[0].name)); YYERROR; }
#line 3040 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 63:
#line 1449 "glb_parser.y" /* yacc.c:1646  */
    {(yyval.in)=(yyvsp[0].in);}
#line 3046 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 64:
#line 1450 "glb_parser.y" /* yacc.c:1646  */
    {(yyval.in)=1;}
#line 3052 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 65:
#line 1451 "glb_parser.y" /* yacc.c:1646  */
    {(yyval.in)=-1;}
#line 3058 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 66:
#line 1456 "glb_parser.y" /* yacc.c:1646  */
    {
  glb_List **buf;
  energy_len=1;

  buf=(glb_List**) glb_malloc(sizeof( glb_List* ) );
  buf[0]=(yyvsp[0].ptr);
  (yyval.ptrq)=buf;
}
#line 3071 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 67:
#line 1465 "glb_parser.y" /* yacc.c:1646  */
    {
  glb_List **buf;
  buf=(yyvsp[-2].ptrq);
  energy_len++;

  buf=(glb_List**) glb_realloc((void**) buf , sizeof( glb_List* ) * energy_len);

  buf[energy_len-1]=(yyvsp[0].ptr);
  (yyval.ptrq)=buf;
}
#line 3086 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 68:
#line 1478 "glb_parser.y" /* yacc.c:1646  */
    {
  set_exp_energy("@energy",(yyvsp[0].ptrq));
  if ((yyvsp[-2].name))  { glb_free((yyvsp[-2].name));  (yyvsp[-2].name)=NULL; }
}
#line 3095 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 69:
#line 1482 "glb_parser.y" /* yacc.c:1646  */
    {
  set_exp_energy("@energy",(yyvsp[-1].ptrq)); 
  if ((yyvsp[-3].name))  { glb_free((yyvsp[-3].name));  (yyvsp[-3].name)=NULL; }
}
#line 3104 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 70:
#line 1489 "glb_parser.y" /* yacc.c:1646  */
    {
  glb_List **buf;
  buf=(glb_List**) glb_malloc(sizeof(glb_List*)*2);
  buf[0]=list_cons(NULL,(yyvsp[0].dpt)[0]);
  buf[1]=list_cons(NULL,(yyvsp[0].dpt)[1]);
  glb_free((yyvsp[0].dpt));
  if ((yyvsp[-2].name))  { glb_free((yyvsp[-2].name));  (yyvsp[-2].name)=NULL; }
  (yyval.ptrq)=buf;
}
#line 3118 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 71:
#line 1498 "glb_parser.y" /* yacc.c:1646  */
    {
  glb_List **buf;
  buf=(yyvsp[-2].ptrq);
  buf[0]=list_cons(buf[0],(yyvsp[0].dpt)[0]);
  buf[1]=list_cons(buf[1],(yyvsp[0].dpt)[1]);
  glb_free((yyvsp[0].dpt));
  (yyval.ptrq)=buf;
}
#line 3131 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 72:
#line 1509 "glb_parser.y" /* yacc.c:1646  */
    {
  glb_List **buf;

  buf=(glb_List**) glb_malloc(sizeof(glb_List*)*2);
  buf[0]=list_cons(NULL,(yyvsp[0].dpt)[0]);
  buf[1]=list_cons(NULL,(yyvsp[0].dpt)[1]);
  glb_free((yyvsp[0].dpt));
  if ((yyvsp[-2].name))  { glb_free((yyvsp[-2].name));  (yyvsp[-2].name)=NULL; }
  (yyval.ptrq)=buf;
}
#line 3146 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 73:
#line 1519 "glb_parser.y" /* yacc.c:1646  */
    {
  glb_List **buf;
  buf=(yyvsp[-2].ptrq);
  buf[0]=list_cons(buf[0],(yyvsp[0].dpt)[0]);
  buf[1]=list_cons(buf[1],(yyvsp[0].dpt)[1]);
  glb_free((yyvsp[0].dpt));
  (yyval.ptrq)=buf;
}
#line 3159 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 74:
#line 1530 "glb_parser.y" /* yacc.c:1646  */
    {
  int flag;
  (yyval.ptrq)=(yyvsp[0].ptrq);
  flag=set_exp_list("bgrulescoeff",(yyvsp[0].ptrq)[0],0);
  if(flag==1) yyerror("Invalid coefficient in @background");
  flag=set_exp_list("bgrulechannellist",(yyvsp[0].ptrq)[1],0);
  if(flag==1) yyerror("Invalid channel in @background");
  glb_free((yyvsp[0].ptrq));
}
#line 3173 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 75:
#line 1539 "glb_parser.y" /* yacc.c:1646  */
    {
  int flag;
  (yyval.ptrq)=(yyvsp[0].ptrq);
  flag=set_exp_list("rulescoeff",(yyvsp[0].ptrq)[0],0);
  if(flag==1) yyerror("Invalid coefficient in @signal");
  flag=set_exp_list("rulechannellist",(yyvsp[0].ptrq)[1],0);
  if(flag==1) yyerror("Invalid channel in @signal");
  glb_free((yyvsp[0].ptrq));
}
#line 3187 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 76:
#line 1548 "glb_parser.y" /* yacc.c:1646  */
    {
//JK, 2012-05-17  buff.sys_on_strings[buff.numofrules-1] = strdup($3);
  buff.sys_on_strings[loc_count-1] = strdup((yyvsp[0].name));
  if ((yyvsp[-2].name))  { glb_free((yyvsp[-2].name));  (yyvsp[-2].name)=NULL; }
  if ((yyvsp[0].name))  { glb_free((yyvsp[0].name));  (yyvsp[0].name)=NULL; }
}
#line 3198 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 77:
#line 1554 "glb_parser.y" /* yacc.c:1646  */
    {
//JK, 2012-05-17  buff.sys_off_strings[buff.numofrules-1] = strdup($3);
  buff.sys_off_strings[loc_count-1] = strdup((yyvsp[0].name));
  if ((yyvsp[-2].name))  { glb_free((yyvsp[-2].name));  (yyvsp[-2].name)=NULL; }
  if ((yyvsp[0].name))  { glb_free((yyvsp[0].name));  (yyvsp[0].name)=NULL; }
}
#line 3209 "glb_parser.c" /* yacc.c:1646  */
    break;

  case 78:
#line 1560 "glb_parser.y" /* yacc.c:1646  */
    {
  set_multiex_errors((yyvsp[-2].name), (yyvsp[0].ptrq));
  if ((yyvsp[-2].name))  { glb_free((yyvsp[-2].name));  (yyvsp[-2].name)=NULL; }
}
#line 3218 "glb_parser.c" /* yacc.c:1646  */
    break;


#line 3222 "glb_parser.c" /* yacc.c:1646  */
      default: break;
    }
  /* User semantic actions sometimes alter yychar, and that requires
     that yytoken be updated with the new translation.  We take the
     approach of translating immediately before every use of yytoken.
     One alternative is translating here after every semantic action,
     but that translation would be missed if the semantic action invokes
     YYABORT, YYACCEPT, or YYERROR immediately after altering yychar or
     if it invokes YYBACKUP.  In the case of YYABORT or YYACCEPT, an
     incorrect destructor might then be invoked immediately.  In the
     case of YYERROR or YYBACKUP, subsequent parser actions might lead
     to an incorrect destructor call or verbose syntax error message
     before the lookahead is translated.  */
  YY_SYMBOL_PRINT ("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;

  /* Now 'shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*--------------------------------------.
| yyerrlab -- here on detecting error.  |
`--------------------------------------*/
yyerrlab:
  /* Make sure we have latest lookahead translation.  See comments at
     user semantic actions for why this is necessary.  */
  yytoken = yychar == YYEMPTY ? YYEMPTY : YYTRANSLATE (yychar);

  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (YY_("syntax error"));
#else
# define YYSYNTAX_ERROR yysyntax_error (&yymsg_alloc, &yymsg, \
                                        yyssp, yytoken)
      {
        char const *yymsgp = YY_("syntax error");
        int yysyntax_error_status;
        yysyntax_error_status = YYSYNTAX_ERROR;
        if (yysyntax_error_status == 0)
          yymsgp = yymsg;
        else if (yysyntax_error_status == 1)
          {
            if (yymsg != yymsgbuf)
              YYSTACK_FREE (yymsg);
            yymsg = (char *) YYSTACK_ALLOC (yymsg_alloc);
            if (!yymsg)
              {
                yymsg = yymsgbuf;
                yymsg_alloc = sizeof yymsgbuf;
                yysyntax_error_status = 2;
              }
            else
              {
                yysyntax_error_status = YYSYNTAX_ERROR;
                yymsgp = yymsg;
              }
          }
        yyerror (yymsgp);
        if (yysyntax_error_status == 2)
          goto yyexhaustedlab;
      }
# undef YYSYNTAX_ERROR
#endif
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
         error, discard it.  */

      if (yychar <= YYEOF)
        {
          /* Return failure if at end of input.  */
          if (yychar == YYEOF)
            YYABORT;
        }
      else
        {
          yydestruct ("Error: discarding",
                      yytoken, &yylval);
          yychar = YYEMPTY;
        }
    }

  /* Else will try to reuse lookahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:

  /* Pacify compilers like GCC when the user code never invokes
     YYERROR and the label yyerrorlab therefore never appears in user
     code.  */
  if (/*CONSTCOND*/ 0)
     goto yyerrorlab;

  /* Do not reclaim the symbols of the rule whose action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;      /* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (!yypact_value_is_default (yyn))
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


      yydestruct ("Error: popping",
                  yystos[yystate], yyvsp);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END


  /* Shift the error token.  */
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

#if !defined yyoverflow || YYERROR_VERBOSE
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (yychar != YYEMPTY)
    {
      /* Make sure we have latest lookahead translation.  See comments at
         user semantic actions for why this is necessary.  */
      yytoken = YYTRANSLATE (yychar);
      yydestruct ("Cleanup: discarding lookahead",
                  yytoken, &yylval);
    }
  /* Do not reclaim the symbols of the rule whose action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
                  yystos[*yyssp], yyvsp);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
#if YYERROR_VERBOSE
  if (yymsg != yymsgbuf)
    YYSTACK_FREE (yymsg);
#endif
  return yyresult;
}
#line 1566 "glb_parser.y" /* yacc.c:1906  */


extern glb_symrec *sym_table;



/***************************************************************************
 * Function yyerror                                                        *
 ***************************************************************************
 * Print parser errors including the line number where the error occured   *
 ***************************************************************************/
int yyerror (const char *s, ...)  /* Called by yyparse on error */
{
  va_list args;
  va_start(args, s);

  if(yydebug > 0) fprintf(stderr,"*****************************************\n");
  fprintf (stderr,"%s:%d: error: ",
           glb_file_id, glb_line_num+1);
  vfprintf(stderr, s, args);
  fprintf(stderr, "\n");
  if(yydebug > 0) fprintf(stderr,"*****************************************\n");
  va_end(args);
  return 0;
}

int
yywarn (const char *s)  /* Called by yyparse on warning */
{
  fprintf (stderr,"%s:%d: warning: %s\n",
           glb_file_id, glb_line_num+1, s);
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
  for(i=0;i<n;i++) fprintf(stdout,"\n");
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
    {"interpolation",NULL,0,0},
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


static glb_List *glb_samplingbincenter(void)
{
  int i;
  glb_List *res=NULL;
  glb_smear *test;
  test=glb_smear_alloc();

  if(buff.simbins<0) {glb_error("Cannot compute samplingbincenter. Sampling-binning not set up properly."); return NULL;}


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
    {"samplingbincenter",glb_samplingbincenter},
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
      p=glb_putsym(ptr->name,ptr->type);
      p->value.var=ptr->value.var;
      if(ptr->list!=NULL) p->list=thread_list(&glb_list_copy,0,0,ptr->list);
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
        if(ptr->list!=NULL){ list_free(ptr->list);}
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
  ptr->value.var = GLB_NAN; /* set value to GLB_NAN, which ensures an
                               error if an undefined variable is used  */
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
  ptr->list=NULL;
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

double glbGetAEDLVariable(const char *name)
{
  glb_symrec *ptr;
  ptr=glb_getpresym(name);
  if (!ptr)
    return GLB_NAN;
  return ptr->value.var;
}

void glbClearAEDLVariables()
{
   if(pre_sym_table!=NULL) free_presymtable();
}


/* same for lists */

void glbDefineAEDLList(const char *name, double *list, size_t length)
{
  size_t i;
  glb_symrec *ptr;
  if(name==NULL) return;
  if(name[0]!='%'){ fprintf(stderr,"ERROR: Names of AEDL lists have to start with %%\n");return;}
  ptr=glb_getpresym(name);
  if(ptr==0) ptr = glb_putpresym (name, LVAR);
  for(i=0;i<length;i++) {ptr->list=list_cons(ptr->list,list[i]);
  }

  return;
}




void glb_copy_buff()
{
  /* I am not sure how well this assigment really works */
  buff.names=copy_names(buff.names);
  if (buff.filename)  glb_free(buff.filename);
  buff.filename=strdup(glb_file_id);
  buff_list[exp_count]  = glbAllocExp();
  *buff_list[exp_count] = buff;
  if (buff.parent)
  {
    struct glb_experiment *p = buff.parent;
    glbExpRemoveChild(p, &buff);
    glbExpAddChild(p, buff_list[exp_count]);
  }
  exp_count++;
}


/***************************************************************************
 * Function glbReset                                                       *
 ***************************************************************************
 * Resets the parser's internal data structures to prepare for the parsing *
 * of a new experiment.                                                    *
 ***************************************************************************/
void glbReset()
{
  glb_line_num =  0;
  energy_len   =  1;
  energy_count = -1;
  loc_count    = -1;
  flux_count   = -1;
  glbResetNuisance();
  glbInitExp(&buff);

  if(name_table!=NULL) free_nametable();
  if(sym_table!=NULL) free_symtable();
  name_table =(glb_namerec *) NULL;
  sym_table =(glb_symrec *) NULL;
  init_table ();
}


/***************************************************************************
 * Function glbNewDetector                                                 *
 ***************************************************************************
 * Starts a new detector that inherits everything from the previously      *
 * defined one, including the namespace.                                   *
 ***************************************************************************/
void glbNewDetector()
{
  glb_copy_buff();
  energy_len = 1;
  glbResetNuisance();
  glbInitExpFromParent(&buff, buff_list[exp_count-1]);
     //FIXME FIXME FIXME What if parent is not/incorrectly defined?

  // Remove rule names from namespace (rules are not copied by glbInitExpFromParent)
  glb_namerec **ptr = &name_table;
  while (*ptr)
  {
    if (strcmp((*ptr)->context, "rule") == 0)
    {
      glb_namerec *ptr_old = *ptr;
      *ptr = ptr_old->next;
      glb_free(ptr_old->name);        ptr_old->name    = NULL;
      glb_free(ptr_old->context);     ptr_old->context = NULL;
      ptr_old->next   = NULL;
      ptr_old->sf     = NULL;
      ptr_old->type   = ptr_old->value = -1;
      glb_free(ptr_old);
    }
    else
      ptr = &((*ptr)->next);
  }

}

void glbResetNuisance()
{
  nuis.name        = NULL;
  nuis.error       = GLB_NAN;
  nuis.a           = GLB_NAN;
  nuis.systype = 0;
  nuis.n_energies  = -1;
  nuis.energy_list = NULL;
  nuis.error_list  = NULL;
  nuis.a_list      = NULL;
  nuis.ref_count   = 1;
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
  glbResetNuisance();
  glbInitExp(&buff);
  for(i=0;i < GLB_MAX_EXP; i++)
    if (buff_list[i])
    {
      glbFreeExp(buff_list[i]);
      buff_list[i] = NULL;
    }
  /* this here would be the place to check for unuses variables, but
     for that we need an access counter */
  if(name_table!=NULL) free_nametable();
  if(sym_table!=NULL) free_symtable();

  name_table =(glb_namerec *) NULL;
  sym_table =(glb_symrec *) NULL;
  init_table ();
  glb_reset_flux(&flt);
  glb_reset_xsec(&xsc);

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
  memset(&flt, 0, sizeof(flt));
  glb_reset_flux(&flt);
  memset(&xsc, 0, sizeof(xsc));
  glb_reset_xsec(&xsc);
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

  /* Copy last experiment and reset data structures */
  glb_copy_buff();
  glbReset();

  glb_fclose(yyin);
  glb_free(context);
  glb_free(glb_file_id);

  if(k!=0) return -2;

  k=0;

  if(*counter+exp_count>GLB_MAX_EXP) glb_fatal("Too many experiments!");
  for(i=0;i<exp_count;i++)
    {
      ins[*counter+i] = buff_list[i];
      buff_list[i]    = NULL; /* Remove our pointer to that exp to prevent destruction */
      k              += glbDefaultExp(ins[*counter+i]);
    }
  (*counter)= (*counter) + exp_count;

  /* Reallocate data structures for the minimizer */
  glb_init_minimizer();

  if(k!=0) return -1;

  return 0;

}
