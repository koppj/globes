
/*  A Bison parser, made from glb_parser.y
    by GNU Bison version 1.28  */

#define YYBISON 1  /* Identify Bison output.  */

#define	NUM	257
#define	SFNCT	258
#define	VAR	259
#define	FNCT	260
#define	IDN	261
#define	CROSS	262
#define	FLUXP	263
#define	FLUXM	264
#define	GRP	265
#define	GID	266
#define	FNAME	267
#define	VERS	268
#define	SIGNAL	269
#define	BG	270
#define	GRPOPEN	271
#define	GRPCLOSE	272
#define	PM	273
#define	FLAVOR	274
#define	NOGLOBES	275
#define	CHANNEL	276
#define	RULESEP	277
#define	RULEMULT	278
#define	ENERGY	279
#define	NAME	280
#define	RDF	281
#define	NDEF	282
#define	NEG	283

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


 
#line 804 "glb_parser.y"
typedef union {
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
#include <stdio.h>

#ifndef __cplusplus
#ifndef __STDC__
#define const
#endif
#endif



#define	YYFINAL		152
#define	YYFLAG		-32768
#define	YYNTBASE	45

#define YYTRANSLATE(x) ((unsigned)(x) <= 283 ? yytranslate[x] : 68)

static const char yytranslate[] = {     0,
     2,     2,     2,     2,     2,     2,     2,     2,     2,    41,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,    43,
    44,    36,    35,    29,    34,     2,    37,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,    42,     2,
    30,    31,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,    39,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,    33,     2,    32,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,    40,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     1,     3,     4,     5,     6,
     7,     8,     9,    10,    11,    12,    13,    14,    15,    16,
    17,    18,    19,    20,    21,    22,    23,    24,    25,    26,
    27,    28,    38
};

#if YYDEBUG != 0
static const short yyprhs[] = {     0,
     0,     1,     4,     6,     8,    10,    13,    16,    19,    22,
    25,    28,    31,    35,    38,    41,    44,    46,    48,    50,
    54,    58,    62,    68,    73,    77,    81,    85,    89,    92,
    96,   100,   102,   104,   107,   111,   115,   120,   123,   127,
   131,   135,   139,   141,   144,   145,   154,   162,   163,   165,
   167,   169,   171,   173,   175,   179,   183,   187,   201,   203,
   205,   207,   209,   211,   215,   219,   224,   228,   232,   236,
   240,   242,   244,   245,   249,   250
};

static const short yyrhs[] = {    -1,
    45,    46,     0,    21,     0,    40,     0,    41,     0,    66,
    41,     0,    64,    41,     0,    47,    41,     0,    48,    41,
     0,    51,    41,     0,     1,    41,     0,    57,    41,     0,
    60,    42,    41,     0,    63,    41,     0,    55,    41,     0,
    56,    41,     0,     3,     0,    26,     0,     5,     0,     5,
    30,    47,     0,     7,    30,    47,     0,     7,    30,     4,
     0,     7,    30,    47,    23,    47,     0,     6,    43,    47,
    44,     0,    47,    35,    47,     0,    47,    34,    47,     0,
    47,    36,    47,     0,    47,    37,    47,     0,    34,    47,
     0,    47,    39,    47,     0,    43,    47,    44,     0,    54,
     0,    28,     0,    47,    29,     0,    47,    29,    41,     0,
    48,    47,    29,     0,    48,    47,    29,    41,     0,    48,
    47,     0,    33,    48,    32,     0,    33,    47,    32,     0,
     7,    30,    48,     0,    47,    24,    47,     0,    46,     0,
    50,    46,     0,     0,    12,    43,    26,    44,    52,    17,
    53,    18,     0,    12,    43,    27,    44,    17,    53,    18,
     0,     0,    47,     0,    63,     0,    50,     0,    57,     0,
    55,     0,    56,     0,    14,    30,    13,     0,     8,    30,
    13,     0,     9,    30,    13,     0,    22,    30,    58,    23,
    59,    23,    20,    23,    20,    23,    58,    23,    58,     0,
    26,     0,    28,     0,    19,     0,    35,     0,    34,     0,
    25,    30,    48,     0,    60,    23,    48,     0,    60,    23,
    41,    48,     0,    16,    30,    49,     0,    61,    23,    49,
     0,    15,    30,    49,     0,    62,    23,    49,     0,    61,
     0,    62,     0,     0,    26,    65,    63,     0,     0,    26,
    67,    57,     0
};

#endif

#if YYDEBUG != 0
static const short yyrline[] = { 0,
   852,   853,   854,   855,   858,   859,   860,   861,   862,   863,
   864,   865,   866,   869,   870,   871,   876,   877,   878,   879,
   880,   882,   883,   885,   886,   887,   888,   889,   890,   891,
   892,   893,   894,   900,   901,   902,   909,   915,   921,   922,
   923,   928,   938,   939,   942,   948,   951,   955,   956,   957,
   958,   959,   960,   961,   964,   968,   975,   984,  1003,  1004,
  1008,  1009,  1010,  1014,  1023,  1033,  1047,  1055,  1065,  1074,
  1084,  1097,  1111,  1116,  1120,  1124
};
#endif


#if YYDEBUG != 0 || defined (YYERROR_VERBOSE)

static const char * const yytname[] = {   "$","error","$undefined.","NUM","SFNCT",
"VAR","FNCT","IDN","CROSS","FLUXP","FLUXM","GRP","GID","FNAME","VERS","SIGNAL",
"BG","GRPOPEN","GRPCLOSE","PM","FLAVOR","NOGLOBES","CHANNEL","RULESEP","RULEMULT",
"ENERGY","NAME","RDF","NDEF","','","'='","'>'","'}'","'{'","'-'","'+'","'*'",
"'/'","NEG","'^'","'\\247'","'\\n'","';'","'('","')'","input","line","exp","seq",
"rulepart","expseq","group","@1","ingroup","version","cross","flux","channel",
"name","pm","ene","brule","srule","rule","outrule","@2","outchannel","@3", NULL
};
#endif

static const short yyr1[] = {     0,
    45,    45,    45,    45,    46,    46,    46,    46,    46,    46,
    46,    46,    46,    46,    46,    46,    47,    47,    47,    47,
    47,    47,    47,    47,    47,    47,    47,    47,    47,    47,
    47,    47,    47,    48,    48,    48,    48,    48,    48,    48,
    48,    49,    50,    50,    52,    51,    51,    53,    53,    53,
    53,    53,    53,    53,    54,    55,    56,    57,    58,    58,
    59,    59,    59,    60,    60,    60,    61,    61,    62,    62,
    63,    63,    65,    64,    67,    66
};

static const short yyr2[] = {     0,
     0,     2,     1,     1,     1,     2,     2,     2,     2,     2,
     2,     2,     3,     2,     2,     2,     1,     1,     1,     3,
     3,     3,     5,     4,     3,     3,     3,     3,     2,     3,
     3,     1,     1,     2,     3,     3,     4,     2,     3,     3,
     3,     3,     1,     2,     0,     8,     7,     0,     1,     1,
     1,     1,     1,     1,     3,     3,     3,    13,     1,     1,
     1,     1,     1,     3,     3,     4,     3,     3,     3,     3,
     1,     1,     0,     3,     0,     3
};

static const short yydefact[] = {     1,
     3,     4,     0,     0,    17,    19,     0,     0,     0,     0,
     0,     0,     0,     0,     0,     0,    18,    33,     0,     0,
     5,     0,     2,     0,     0,     0,    32,     0,     0,     0,
     0,    71,    72,     0,     0,     0,    11,     0,     0,     0,
     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
    18,     0,     0,     0,    29,     0,    34,     0,     0,     0,
     0,     0,     8,     9,    38,    10,    15,    16,    12,     0,
     0,     0,     0,    14,     7,     6,    20,     0,    22,    21,
    41,    56,    57,     0,     0,    55,     0,    69,    67,    59,
    60,     0,     0,    64,    74,    76,    40,    39,     0,    31,
    35,    26,    25,    27,    28,    30,    36,     0,    65,    13,
    68,    70,    24,     0,    45,     0,     0,     0,    21,    37,
    66,    23,     0,     0,    42,    61,    63,    62,     0,     0,
    43,    49,     0,     0,    53,    54,    52,    50,     0,     0,
    44,    47,     0,    46,     0,     0,     0,     0,     0,    58,
     0,     0
};

static const short yydefgoto[] = {     3,
   131,    65,    25,    88,   133,    26,   123,   134,    27,    28,
    29,    30,    92,   129,    31,    32,    33,    34,    35,    49,
    36,    50
};

static const short yypact[] = {   -15,
-32768,-32768,   132,   -37,-32768,   -22,   -28,     3,     8,    11,
    17,    35,    38,    52,    56,    71,    78,-32768,   234,   300,
-32768,   300,-32768,   296,   252,   -14,-32768,    63,    65,    66,
   -16,    85,    86,    72,    82,    87,-32768,   300,   300,     6,
   106,   109,    54,   111,   300,   300,    18,   234,    83,   104,
-32768,   310,   268,    99,    95,    81,   101,   300,   300,   300,
   300,   300,-32768,-32768,   116,-32768,-32768,-32768,-32768,    69,
   102,   300,   300,-32768,-32768,-32768,   325,   285,-32768,    27,
   300,-32768,-32768,    92,   105,-32768,    53,-32768,-32768,-32768,
-32768,   133,   319,   300,-32768,-32768,-32768,-32768,   284,-32768,
-32768,    15,    15,    95,    95,    95,   118,   234,   300,-32768,
-32768,-32768,-32768,   300,-32768,   144,   300,    -5,    27,-32768,
   300,   325,   145,   171,   325,-32768,-32768,-32768,   140,   171,
-32768,   296,   210,   146,    63,    65,    66,    72,   147,   150,
-32768,-32768,   148,-32768,   149,   158,    18,   159,    18,-32768,
   170,-32768
};

static const short yypgoto[] = {-32768,
    -2,    -3,   -17,   -25,-32768,-32768,-32768,    58,-32768,  -106,
  -102,   -45,   -76,-32768,-32768,-32768,-32768,   -46,-32768,-32768,
-32768,-32768
};


#define	YYLAST		364


static const short yytable[] = {    24,
    23,    53,    95,    37,    96,     1,    70,    38,     5,    79,
     6,     7,     8,   126,    39,    52,    55,   135,    56,    12,
    89,   136,    81,   135,     2,    71,    66,   136,   127,   128,
    94,    51,    40,    18,    77,    78,    80,    41,    19,    20,
    42,    87,    87,    90,    93,    91,   111,   112,    22,   114,
    60,    61,   109,    62,   102,   103,   104,   105,   106,    43,
    58,    59,    60,    61,    44,    62,    93,    45,    87,    87,
   148,     5,   150,     6,     7,     8,   117,   138,   137,    84,
    85,    46,    12,   138,   137,    47,    58,    59,    60,    61,
   121,    62,   -73,   -73,    51,   119,    18,    13,    14,   -75,
    48,    19,    20,    67,    93,    68,    69,    72,    73,   108,
   122,    22,    74,   125,    58,    59,    60,    61,    82,    62,
   132,    83,    75,    86,   100,    15,   132,    76,    99,    24,
   141,   151,     4,    62,     5,   115,     6,     7,     8,     9,
    10,   101,   110,    11,   107,    12,    13,    14,   116,    58,
    59,    60,    61,    15,    62,   118,    16,    17,   120,    18,
   124,   130,   139,   142,    19,    20,   143,   144,   146,   152,
   145,     4,    21,     5,    22,     6,     7,     8,     9,    10,
   147,   149,    11,     0,    12,    13,    14,   140,   -48,     0,
     0,     0,    15,     0,     0,    16,    17,     0,    18,     0,
     0,     0,     0,    19,    20,     0,     0,     0,     0,     0,
     4,    21,     5,    22,     6,     7,     8,     9,    10,     0,
     0,    11,     0,    12,    13,    14,     0,   -51,     0,     0,
     0,    15,     0,     0,    16,    17,     5,    18,     6,     7,
     8,     0,    19,    20,     0,     0,     0,    12,     0,     0,
    21,     0,    22,     0,     5,     0,     6,     7,    54,    51,
     0,    18,     0,     0,     0,    12,    19,    20,     0,     0,
     5,     0,     6,     7,    54,     0,    22,    51,     0,    18,
     0,    12,     0,     0,     0,    20,     5,    79,     6,     7,
    54,     0,    64,    51,    22,    18,     0,    12,     0,    98,
     0,    20,     5,     0,     6,     7,    54,     0,     0,    51,
    22,    18,     0,    12,     0,     0,     0,    20,    58,    59,
    60,    61,     0,    62,    57,    51,    22,    18,   113,    58,
    59,    60,    61,    20,    62,     0,    63,     0,    57,     0,
     0,    97,    22,    58,    59,    60,    61,    57,    62,     0,
     0,     0,    58,    59,    60,    61,     0,    62,    58,    59,
    60,    61,     0,    62
};

static const short yycheck[] = {     3,
     3,    19,    49,    41,    50,    21,    23,    30,     3,     4,
     5,     6,     7,    19,    43,    19,    20,   124,    22,    14,
    46,   124,    40,   130,    40,    42,    41,   130,    34,    35,
    48,    26,    30,    28,    38,    39,    40,    30,    33,    34,
    30,    45,    46,    26,    48,    28,    72,    73,    43,    23,
    36,    37,    70,    39,    58,    59,    60,    61,    62,    43,
    34,    35,    36,    37,    30,    39,    70,    30,    72,    73,
   147,     3,   149,     5,     6,     7,    24,   124,   124,    26,
    27,    30,    14,   130,   130,    30,    34,    35,    36,    37,
   108,    39,    15,    16,    26,    99,    28,    15,    16,    22,
    30,    33,    34,    41,   108,    41,    41,    23,    23,    41,
   114,    43,    41,   117,    34,    35,    36,    37,    13,    39,
   124,    13,    41,    13,    44,    22,   130,    41,    30,   133,
   133,     0,     1,    39,     3,    44,     5,     6,     7,     8,
     9,    41,    41,    12,    29,    14,    15,    16,    44,    34,
    35,    36,    37,    22,    39,    23,    25,    26,    41,    28,
    17,    17,    23,    18,    33,    34,    20,    18,    20,     0,
    23,     1,    41,     3,    43,     5,     6,     7,     8,     9,
    23,    23,    12,    -1,    14,    15,    16,   130,    18,    -1,
    -1,    -1,    22,    -1,    -1,    25,    26,    -1,    28,    -1,
    -1,    -1,    -1,    33,    34,    -1,    -1,    -1,    -1,    -1,
     1,    41,     3,    43,     5,     6,     7,     8,     9,    -1,
    -1,    12,    -1,    14,    15,    16,    -1,    18,    -1,    -1,
    -1,    22,    -1,    -1,    25,    26,     3,    28,     5,     6,
     7,    -1,    33,    34,    -1,    -1,    -1,    14,    -1,    -1,
    41,    -1,    43,    -1,     3,    -1,     5,     6,     7,    26,
    -1,    28,    -1,    -1,    -1,    14,    33,    34,    -1,    -1,
     3,    -1,     5,     6,     7,    -1,    43,    26,    -1,    28,
    -1,    14,    -1,    -1,    -1,    34,     3,     4,     5,     6,
     7,    -1,    41,    26,    43,    28,    -1,    14,    -1,    32,
    -1,    34,     3,    -1,     5,     6,     7,    -1,    -1,    26,
    43,    28,    -1,    14,    -1,    -1,    -1,    34,    34,    35,
    36,    37,    -1,    39,    29,    26,    43,    28,    44,    34,
    35,    36,    37,    34,    39,    -1,    41,    -1,    29,    -1,
    -1,    32,    43,    34,    35,    36,    37,    29,    39,    -1,
    -1,    -1,    34,    35,    36,    37,    -1,    39,    34,    35,
    36,    37,    -1,    39
};
/* -*-C-*-  Note some compilers choke on comments on `#line' lines.  */
#line 3 "/usr/share/bison.simple"
/* This file comes from bison-1.28.  */

/* Skeleton output parser for bison,
   Copyright (C) 1984, 1989, 1990 Free Software Foundation, Inc.

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

/* This is the parser code that is written into each bison parser
  when the %semantic_parser declaration is not specified in the grammar.
  It was written by Richard Stallman by simplifying the hairy parser
  used when %semantic_parser is specified.  */

#ifndef YYPARSE_RETURN_TYPE
#define YYPARSE_RETURN_TYPE int
#endif


#ifndef YYSTACK_USE_ALLOCA
#ifdef alloca
#define YYSTACK_USE_ALLOCA
#else /* alloca not defined */
#ifdef __GNUC__
#define YYSTACK_USE_ALLOCA
#define alloca __builtin_alloca
#else /* not GNU C.  */
#if (!defined (__STDC__) && defined (sparc)) || defined (__sparc__) || defined (__sparc) || defined (__sgi) || (defined (__sun) && defined (__i386))
#define YYSTACK_USE_ALLOCA
#include <alloca.h>
#else /* not sparc */
/* We think this test detects Watcom and Microsoft C.  */
/* This used to test MSDOS, but that is a bad idea
   since that symbol is in the user namespace.  */
#if (defined (_MSDOS) || defined (_MSDOS_)) && !defined (__TURBOC__)
#if 0 /* No need for malloc.h, which pollutes the namespace;
	 instead, just don't use alloca.  */
#include <malloc.h>
#endif
#else /* not MSDOS, or __TURBOC__ */
#if defined(_AIX)
/* I don't know what this was needed for, but it pollutes the namespace.
   So I turned it off.   rms, 2 May 1997.  */
/* #include <malloc.h>  */
 #pragma alloca
#define YYSTACK_USE_ALLOCA
#else /* not MSDOS, or __TURBOC__, or _AIX */
#if 0
#ifdef __hpux /* haible@ilog.fr says this works for HPUX 9.05 and up,
		 and on HPUX 10.  Eventually we can turn this on.  */
#define YYSTACK_USE_ALLOCA
#define alloca __builtin_alloca
#endif /* __hpux */
#endif
#endif /* not _AIX */
#endif /* not MSDOS, or __TURBOC__ */
#endif /* not sparc */
#endif /* not GNU C */
#endif /* alloca not defined */
#endif /* YYSTACK_USE_ALLOCA not defined */

#ifdef YYSTACK_USE_ALLOCA
#define YYSTACK_ALLOC alloca
#else
#define YYSTACK_ALLOC malloc
#endif

/* Note: there must be only one dollar sign in this file.
   It is replaced by the list of actions, each action
   as one case of the switch.  */

#define yyerrok		(yyerrstatus = 0)
#define yyclearin	(yychar = YYEMPTY)
#define YYEMPTY		-2
#define YYEOF		0
#define YYACCEPT	goto yyacceptlab
#define YYABORT 	goto yyabortlab
#define YYERROR		goto yyerrlab1
/* Like YYERROR except do call yyerror.
   This remains here temporarily to ease the
   transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  */
#define YYFAIL		goto yyerrlab
#define YYRECOVERING()  (!!yyerrstatus)
#define YYBACKUP(token, value) \
do								\
  if (yychar == YYEMPTY && yylen == 1)				\
    { yychar = (token), yylval = (value);			\
      yychar1 = YYTRANSLATE (yychar);				\
      YYPOPSTACK;						\
      goto yybackup;						\
    }								\
  else								\
    { yyerror ("syntax error: cannot back up"); YYERROR; }	\
while (0)

#define YYTERROR	1
#define YYERRCODE	256

#ifndef YYPURE
#define YYLEX		yylex()
#endif

#ifdef YYPURE
#ifdef YYLSP_NEEDED
#ifdef YYLEX_PARAM
#define YYLEX		yylex(&yylval, &yylloc, YYLEX_PARAM)
#else
#define YYLEX		yylex(&yylval, &yylloc)
#endif
#else /* not YYLSP_NEEDED */
#ifdef YYLEX_PARAM
#define YYLEX		yylex(&yylval, YYLEX_PARAM)
#else
#define YYLEX		yylex(&yylval)
#endif
#endif /* not YYLSP_NEEDED */
#endif

/* If nonreentrant, generate the variables here */

#ifndef YYPURE

int	yychar;			/*  the lookahead symbol		*/
YYSTYPE	yylval;			/*  the semantic value of the		*/
				/*  lookahead symbol			*/

#ifdef YYLSP_NEEDED
YYLTYPE yylloc;			/*  location data for the lookahead	*/
				/*  symbol				*/
#endif

int yynerrs;			/*  number of parse errors so far       */
#endif  /* not YYPURE */

#if YYDEBUG != 0
int yydebug;			/*  nonzero means print parse trace	*/
/* Since this is uninitialized, it does not stop multiple parsers
   from coexisting.  */
#endif

/*  YYINITDEPTH indicates the initial size of the parser's stacks	*/

#ifndef	YYINITDEPTH
#define YYINITDEPTH 200
#endif

/*  YYMAXDEPTH is the maximum size the stacks can grow to
    (effective only if the built-in stack extension method is used).  */

#if YYMAXDEPTH == 0
#undef YYMAXDEPTH
#endif

#ifndef YYMAXDEPTH
#define YYMAXDEPTH 10000
#endif

/* Define __yy_memcpy.  Note that the size argument
   should be passed with type unsigned int, because that is what the non-GCC
   definitions require.  With GCC, __builtin_memcpy takes an arg
   of type size_t, but it can handle unsigned int.  */

#if __GNUC__ > 1		/* GNU C and GNU C++ define this.  */
#define __yy_memcpy(TO,FROM,COUNT)	__builtin_memcpy(TO,FROM,COUNT)
#else				/* not GNU C or C++ */
#ifndef __cplusplus

/* This is the most reliable way to avoid incompatibilities
   in available built-in functions on various systems.  */
static void
__yy_memcpy (to, from, count)
     char *to;
     char *from;
     unsigned int count;
{
  register char *f = from;
  register char *t = to;
  register int i = count;

  while (i-- > 0)
    *t++ = *f++;
}

#else /* __cplusplus */

/* This is the most reliable way to avoid incompatibilities
   in available built-in functions on various systems.  */
static void
__yy_memcpy (char *to, char *from, unsigned int count)
{
  register char *t = to;
  register char *f = from;
  register int i = count;

  while (i-- > 0)
    *t++ = *f++;
}

#endif
#endif

#line 222 "/usr/share/bison.simple"

/* The user can define YYPARSE_PARAM as the name of an argument to be passed
   into yyparse.  The argument should have type void *.
   It should actually point to an object.
   Grammar actions can access the variable by casting it
   to the proper pointer type.  */

#ifdef YYPARSE_PARAM
#ifdef __cplusplus
#define YYPARSE_PARAM_ARG void *YYPARSE_PARAM
#define YYPARSE_PARAM_DECL
#else /* not __cplusplus */
#define YYPARSE_PARAM_ARG YYPARSE_PARAM
#define YYPARSE_PARAM_DECL void *YYPARSE_PARAM;
#endif /* not __cplusplus */
#else /* not YYPARSE_PARAM */
#define YYPARSE_PARAM_ARG
#define YYPARSE_PARAM_DECL
#endif /* not YYPARSE_PARAM */

/* Prevent warning if -Wstrict-prototypes.  */
#ifdef __GNUC__
#ifdef YYPARSE_PARAM
YYPARSE_RETURN_TYPE
yyparse (void *);
#else
YYPARSE_RETURN_TYPE
yyparse (void);
#endif
#endif

YYPARSE_RETURN_TYPE
yyparse(YYPARSE_PARAM_ARG)
     YYPARSE_PARAM_DECL
{
  register int yystate;
  register int yyn;
  register short *yyssp;
  register YYSTYPE *yyvsp;
  int yyerrstatus;	/*  number of tokens to shift before error messages enabled */
  int yychar1 = 0;		/*  lookahead token as an internal (translated) token number */

  short	yyssa[YYINITDEPTH];	/*  the state stack			*/
  YYSTYPE yyvsa[YYINITDEPTH];	/*  the semantic value stack		*/

  short *yyss = yyssa;		/*  refer to the stacks thru separate pointers */
  YYSTYPE *yyvs = yyvsa;	/*  to allow yyoverflow to reallocate them elsewhere */

#ifdef YYLSP_NEEDED
  YYLTYPE yylsa[YYINITDEPTH];	/*  the location stack			*/
  YYLTYPE *yyls = yylsa;
  YYLTYPE *yylsp;

#define YYPOPSTACK   (yyvsp--, yyssp--, yylsp--)
#else
#define YYPOPSTACK   (yyvsp--, yyssp--)
#endif

  int yystacksize = YYINITDEPTH;
#ifndef YYSTACK_USE_ALLOCA
  int yyfree_stacks = 0;
#endif

#ifdef YYPURE
  int yychar;
  YYSTYPE yylval;
  int yynerrs;
#ifdef YYLSP_NEEDED
  YYLTYPE yylloc;
#endif
#endif

  YYSTYPE yyval;		/*  the variable used to return		*/
				/*  semantic values from the action	*/
				/*  routines				*/

  int yylen;

#if YYDEBUG != 0
  if (yydebug)
    fprintf(stderr, "Starting parse\n");
#endif

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY;		/* Cause a token to be read.  */

  /* Initialize stack pointers.
     Waste one element of value and location stack
     so that they stay on the same level as the state stack.
     The wasted elements are never initialized.  */

  yyssp = yyss - 1;
  yyvsp = yyvs;
#ifdef YYLSP_NEEDED
  yylsp = yyls;
#endif

/* Push a new state, which is found in  yystate  .  */
/* In all cases, when you get here, the value and location stacks
   have just been pushed. so pushing a state here evens the stacks.  */
yynewstate:

  *++yyssp = yystate;

  if (yyssp >= yyss + yystacksize - 1)
    {
      /* Give user a chance to reallocate the stack */
      /* Use copies of these so that the &'s don't force the real ones into memory. */
      YYSTYPE *yyvs1 = yyvs;
      short *yyss1 = yyss;
#ifdef YYLSP_NEEDED
      YYLTYPE *yyls1 = yyls;
#endif

      /* Get the current used size of the three stacks, in elements.  */
      int size = yyssp - yyss + 1;

#ifdef yyoverflow
      /* Each stack pointer address is followed by the size of
	 the data in use in that stack, in bytes.  */
#ifdef YYLSP_NEEDED
      /* This used to be a conditional around just the two extra args,
	 but that might be undefined if yyoverflow is a macro.  */
      yyoverflow("parser stack overflow",
		 &yyss1, size * sizeof (*yyssp),
		 &yyvs1, size * sizeof (*yyvsp),
		 &yyls1, size * sizeof (*yylsp),
		 &yystacksize);
#else
      yyoverflow("parser stack overflow",
		 &yyss1, size * sizeof (*yyssp),
		 &yyvs1, size * sizeof (*yyvsp),
		 &yystacksize);
#endif

      yyss = yyss1; yyvs = yyvs1;
#ifdef YYLSP_NEEDED
      yyls = yyls1;
#endif
#else /* no yyoverflow */
      /* Extend the stack our own way.  */
      if (yystacksize >= YYMAXDEPTH)
	{
	  yyerror("parser stack overflow");
#ifndef YYSTACK_USE_ALLOCA
	  if (yyfree_stacks)
	    {
	      free (yyss);
	      free (yyvs);
#ifdef YYLSP_NEEDED
	      free (yyls);
#endif
	    }
#endif	    
	  return 2;
	}
      yystacksize *= 2;
      if (yystacksize > YYMAXDEPTH)
	yystacksize = YYMAXDEPTH;
#ifndef YYSTACK_USE_ALLOCA
      yyfree_stacks = 1;
#endif
      yyss = (short *) YYSTACK_ALLOC (yystacksize * sizeof (*yyssp));
      __yy_memcpy ((char *)yyss, (char *)yyss1,
		   size * (unsigned int) sizeof (*yyssp));
      yyvs = (YYSTYPE *) YYSTACK_ALLOC (yystacksize * sizeof (*yyvsp));
      __yy_memcpy ((char *)yyvs, (char *)yyvs1,
		   size * (unsigned int) sizeof (*yyvsp));
#ifdef YYLSP_NEEDED
      yyls = (YYLTYPE *) YYSTACK_ALLOC (yystacksize * sizeof (*yylsp));
      __yy_memcpy ((char *)yyls, (char *)yyls1,
		   size * (unsigned int) sizeof (*yylsp));
#endif
#endif /* no yyoverflow */

      yyssp = yyss + size - 1;
      yyvsp = yyvs + size - 1;
#ifdef YYLSP_NEEDED
      yylsp = yyls + size - 1;
#endif

#if YYDEBUG != 0
      if (yydebug)
	fprintf(stderr, "Stack size increased to %d\n", yystacksize);
#endif

      if (yyssp >= yyss + yystacksize - 1)
	YYABORT;
    }

#if YYDEBUG != 0
  if (yydebug)
    fprintf(stderr, "Entering state %d\n", yystate);
#endif

  goto yybackup;
 yybackup:

/* Do appropriate processing given the current state.  */
/* Read a lookahead token if we need one and don't already have one.  */
/* yyresume: */

  /* First try to decide what to do without reference to lookahead token.  */

  yyn = yypact[yystate];
  if (yyn == YYFLAG)
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* yychar is either YYEMPTY or YYEOF
     or a valid token in external form.  */

  if (yychar == YYEMPTY)
    {
#if YYDEBUG != 0
      if (yydebug)
	fprintf(stderr, "Reading a token: ");
#endif
      yychar = YYLEX;
    }

  /* Convert token to internal form (in yychar1) for indexing tables with */

  if (yychar <= 0)		/* This means end of input. */
    {
      yychar1 = 0;
      yychar = YYEOF;		/* Don't call YYLEX any more */

#if YYDEBUG != 0
      if (yydebug)
	fprintf(stderr, "Now at end of input.\n");
#endif
    }
  else
    {
      yychar1 = YYTRANSLATE(yychar);

#if YYDEBUG != 0
      if (yydebug)
	{
	  fprintf (stderr, "Next token is %d (%s", yychar, yytname[yychar1]);
	  /* Give the individual parser a way to print the precise meaning
	     of a token, for further debugging info.  */
#ifdef YYPRINT
	  YYPRINT (stderr, yychar, yylval);
#endif
	  fprintf (stderr, ")\n");
	}
#endif
    }

  yyn += yychar1;
  if (yyn < 0 || yyn > YYLAST || yycheck[yyn] != yychar1)
    goto yydefault;

  yyn = yytable[yyn];

  /* yyn is what to do for this token type in this state.
     Negative => reduce, -yyn is rule number.
     Positive => shift, yyn is new state.
       New state is final state => don't bother to shift,
       just return success.
     0, or most negative number => error.  */

  if (yyn < 0)
    {
      if (yyn == YYFLAG)
	goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }
  else if (yyn == 0)
    goto yyerrlab;

  if (yyn == YYFINAL)
    YYACCEPT;

  /* Shift the lookahead token.  */

#if YYDEBUG != 0
  if (yydebug)
    fprintf(stderr, "Shifting token %d (%s), ", yychar, yytname[yychar1]);
#endif

  /* Discard the token being shifted unless it is eof.  */
  if (yychar != YYEOF)
    yychar = YYEMPTY;

  *++yyvsp = yylval;
#ifdef YYLSP_NEEDED
  *++yylsp = yylloc;
#endif

  /* count tokens shifted since error; after three, turn off error status.  */
  if (yyerrstatus) yyerrstatus--;

  yystate = yyn;
  goto yynewstate;

/* Do the default action for the current state.  */
yydefault:

  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;

/* Do a reduction.  yyn is the number of a rule to reduce with.  */
yyreduce:
  yylen = yyr2[yyn];
  if (yylen > 0)
    yyval = yyvsp[1-yylen]; /* implement default value of the action */

#if YYDEBUG != 0
  if (yydebug)
    {
      int i;

      fprintf (stderr, "Reducing via rule %d (line %d), ",
	       yyn, yyrline[yyn]);

      /* Print the symbols being reduced, and their result.  */
      for (i = yyprhs[yyn]; yyrhs[i] > 0; i++)
	fprintf (stderr, "%s ", yytname[yyrhs[i]]);
      fprintf (stderr, " -> %s\n", yytname[yyr1[yyn]]);
    }
#endif


  switch (yyn) {

case 3:
#line 854 "glb_parser.y"
{YYABORT;;
    break;}
case 4:
#line 855 "glb_parser.y"
{YYABORT;;
    break;}
case 7:
#line 860 "glb_parser.y"
{ /*showlist($1[0]);showlist($1[1]);*/ ;
    break;}
case 8:
#line 861 "glb_parser.y"
{ /*printf ("\t%.10g\n", $1);*/  ;
    break;}
case 9:
#line 862 "glb_parser.y"
{ /*showlist($1);printf("\n");*/ ;
    break;}
case 10:
#line 863 "glb_parser.y"
{;
    break;}
case 11:
#line 864 "glb_parser.y"
{yyerrok;;
    break;}
case 12:
#line 865 "glb_parser.y"
{;
    break;}
case 13:
#line 866 "glb_parser.y"
{
 set_exp_energy("@energy",yyvsp[-2].ptrq);
;
    break;}
case 14:
#line 869 "glb_parser.y"
{;
    break;}
case 15:
#line 870 "glb_parser.y"
{/* printf("%s\n",$1);*/;
    break;}
case 16:
#line 871 "glb_parser.y"
{/* printf("%s\n",$1);*/;
    break;}
case 17:
#line 876 "glb_parser.y"
{ yyval.val = yyvsp[0].val;                         ;
    break;}
case 18:
#line 877 "glb_parser.y"
{ yyval.val = yyvsp[0].nameptr->value;              ;
    break;}
case 19:
#line 878 "glb_parser.y"
{ yyval.val = yyvsp[0].tptr->value.var;              ;
    break;}
case 20:
#line 879 "glb_parser.y"
{ yyval.val = yyvsp[0].val; yyvsp[-2].tptr->value.var = yyvsp[0].val;     ;
    break;}
case 21:
#line 880 "glb_parser.y"
{ if(set_exp(yyvsp[-2].name,yyvsp[0].val,0)==1) yyerror("Unknown identifier");
yyval.val = yyvsp[0].val; ;
    break;}
case 22:
#line 882 "glb_parser.y"
{ if(set_fnct(yyvsp[-2].name,yyvsp[0].nameptr->sf)==1) yyerror("Unknown identifier");;
    break;}
case 23:
#line 883 "glb_parser.y"
{ if(set_pair(yyvsp[-4].name,yyvsp[-2].val,yyvsp[0].val,0)==1) yyerror("Unknown identifier");
yyval.val = yyvsp[-2].val; ;
    break;}
case 24:
#line 885 "glb_parser.y"
{ yyval.val = (*(yyvsp[-3].tptr->value.fnctptr))(yyvsp[-1].val); ;
    break;}
case 25:
#line 886 "glb_parser.y"
{ yyval.val = yyvsp[-2].val + yyvsp[0].val;                    ;
    break;}
case 26:
#line 887 "glb_parser.y"
{ yyval.val = yyvsp[-2].val - yyvsp[0].val;                    ;
    break;}
case 27:
#line 888 "glb_parser.y"
{ yyval.val = yyvsp[-2].val * yyvsp[0].val;                    ;
    break;}
case 28:
#line 889 "glb_parser.y"
{ yyval.val = yyvsp[-2].val / yyvsp[0].val;                    ;
    break;}
case 29:
#line 890 "glb_parser.y"
{ yyval.val = -yyvsp[0].val;                        ;
    break;}
case 30:
#line 891 "glb_parser.y"
{ yyval.val = pow (yyvsp[-2].val, yyvsp[0].val);               ;
    break;}
case 31:
#line 892 "glb_parser.y"
{ yyval.val = yyvsp[-1].val;                         ;
    break;}
case 32:
#line 893 "glb_parser.y"
{ yyval.val = 0;;
    break;}
case 33:
#line 894 "glb_parser.y"
{yyerror("Unknown name");YYERROR;;
    break;}
case 34:
#line 900 "glb_parser.y"
{yyval.ptr=list_cons(NULL,yyvsp[-1].val); ;
    break;}
case 35:
#line 901 "glb_parser.y"
{yyval.ptr=list_cons(NULL,yyvsp[-2].val); ;
    break;}
case 36:
#line 902 "glb_parser.y"
{
  glb_List *buf;
  buf=yyvsp[-2].ptr;
  buf=list_cons(buf,yyvsp[-1].val);
  yyval.ptr=buf;   
;
    break;}
case 37:
#line 909 "glb_parser.y"
{
  glb_List *buf;
  buf=yyvsp[-3].ptr;
  buf=list_cons(buf,yyvsp[-2].val);
  yyval.ptr=buf;   
;
    break;}
case 38:
#line 915 "glb_parser.y"
{
  glb_List *buf;
  buf=yyvsp[-1].ptr;
  buf=list_cons(buf,yyvsp[0].val);
  yyval.ptr=buf;   
;
    break;}
case 39:
#line 921 "glb_parser.y"
{yyval.ptr=yyvsp[-1].ptr; ;
    break;}
case 40:
#line 922 "glb_parser.y"
{yyval.ptr=list_cons(NULL,yyvsp[-1].val); ;
    break;}
case 41:
#line 923 "glb_parser.y"
{ if(set_exp_list(yyvsp[-2].name,yyvsp[0].ptr,3)==1) 
  yyerror("Unknown identifier");
 yyval.ptr = yyvsp[0].ptr; ;
    break;}
case 42:
#line 928 "glb_parser.y"
{ 

  double *buf;
  buf=(double*) glb_malloc(sizeof(double)*2);
  buf[0]=yyvsp[-2].val;
  buf[1]=yyvsp[0].val-1;
  yyval.dpt=buf;
;
    break;}
case 43:
#line 938 "glb_parser.y"
{;
    break;}
case 44:
#line 939 "glb_parser.y"
{;
    break;}
case 45:
#line 943 "glb_parser.y"
{ if(yyvsp[-1].nameptr->value==-1) {yyvsp[-1].nameptr->value=step_counter(yyvsp[-3].name); }
  loc_count=yyvsp[-1].nameptr->value;
  glb_free(context);
  context =(char *) strdup(yyvsp[-3].name);
;
    break;}
case 46:
#line 948 "glb_parser.y"
{ 
   
grp_end(yyvsp[-7].name);;
    break;}
case 47:
#line 951 "glb_parser.y"
{
  yyerror("Redefinition of an automatic variable");YYERROR;;
    break;}
case 49:
#line 956 "glb_parser.y"
{;
    break;}
case 50:
#line 957 "glb_parser.y"
{;
    break;}
case 51:
#line 958 "glb_parser.y"
{;
    break;}
case 52:
#line 959 "glb_parser.y"
{;
    break;}
case 53:
#line 960 "glb_parser.y"
{;
    break;}
case 54:
#line 961 "glb_parser.y"
{;
    break;}
case 55:
#line 964 "glb_parser.y"
{
  buff.version=strdup(yyvsp[0].name);
;
    break;}
case 56:
#line 968 "glb_parser.y"
{
  //load_cross($3,loc_count-1);
  xsc.file_name=strdup(yyvsp[0].name);
  yyval.name=yyvsp[0].name;
;
    break;}
case 57:
#line 975 "glb_parser.y"
{
  //load_flux($3,loc_count-1,1);
  flt.file_name=strdup(yyvsp[0].name);
  //if(set_exp($1,$3,0)==1) yyerror("Unknown identifier");
  yyval.name=yyvsp[0].name;
;
    break;}
case 58:
#line 985 "glb_parser.y"
{

  int x[6];
  x[0]=yyvsp[-10].nameptr->value - 1 ;
  x[1]=yyvsp[-8].in;
  x[2]=yyvsp[-6].in;
  x[3]=yyvsp[-4].in;
  x[4]=(int) yyvsp[-2].nameptr->value -1; 
  x[5]=(int) yyvsp[0].nameptr->value -1;


  set_channel_data(x,loc_count);



;
    break;}
case 59:
#line 1003 "glb_parser.y"
{yyval.nameptr=yyvsp[0].nameptr;
    break;}
case 60:
#line 1004 "glb_parser.y"
{yyerror("Unknown name");YYERROR;;
    break;}
case 61:
#line 1008 "glb_parser.y"
{yyval.in=yyvsp[0].in;;
    break;}
case 62:
#line 1009 "glb_parser.y"
{yyval.in=1;;
    break;}
case 63:
#line 1010 "glb_parser.y"
{yyval.in=-1;;
    break;}
case 64:
#line 1014 "glb_parser.y"
{
  glb_List **buf;
  energy_len=1;
 
  buf=(glb_List**) glb_malloc(sizeof( glb_List* ) ); 
  buf[0]=yyvsp[0].ptr;
 
  yyval.ptrq=buf;
;
    break;}
case 65:
#line 1024 "glb_parser.y"
{
  glb_List **buf;
  buf=yyvsp[-2].ptrq;
  energy_len++;
  
  buf=(glb_List**) glb_realloc((void**) buf , sizeof( glb_List* ) * energy_len);
 
  buf[energy_len-1]=yyvsp[0].ptr; 
  yyval.ptrq=buf; ;
    break;}
case 66:
#line 1034 "glb_parser.y"
{
  glb_List **buf;
  buf=yyvsp[-3].ptrq;
  energy_len++;

  buf=(glb_List**) glb_realloc((void**) buf , sizeof( glb_List* ) * energy_len);

  buf[energy_len-1]=yyvsp[0].ptr; 
  yyval.ptrq=buf; ;
    break;}
case 67:
#line 1047 "glb_parser.y"
{
  glb_List **buf;
  buf=(glb_List**) glb_malloc(sizeof(glb_List*)*2);
  buf[0]=list_cons(NULL,yyvsp[0].dpt[0]);
  buf[1]=list_cons(NULL,yyvsp[0].dpt[1]);
  glb_free(yyvsp[0].dpt);
  yyval.ptrq=buf;
;
    break;}
case 68:
#line 1055 "glb_parser.y"
{
  glb_List **buf;
  buf=yyvsp[-2].ptrq;
  buf[0]=list_cons(buf[0],yyvsp[0].dpt[0]);
  buf[1]=list_cons(buf[1],yyvsp[0].dpt[1]);
  glb_free(yyvsp[0].dpt);
  yyval.ptrq=buf; 
;
    break;}
case 69:
#line 1065 "glb_parser.y"
{
  glb_List **buf;  
 
  buf=(glb_List**) glb_malloc(sizeof(glb_List*)*2);
  buf[0]=list_cons(NULL,yyvsp[0].dpt[0]);
  buf[1]=list_cons(NULL,yyvsp[0].dpt[1]);
  glb_free(yyvsp[0].dpt);
  yyval.ptrq=buf;
;
    break;}
case 70:
#line 1074 "glb_parser.y"
{
  glb_List **buf;
  buf=yyvsp[-2].ptrq;
  buf[0]=list_cons(buf[0],yyvsp[0].dpt[0]);
  buf[1]=list_cons(buf[1],yyvsp[0].dpt[1]);
  glb_free(yyvsp[0].dpt);
  yyval.ptrq=buf; 
;
    break;}
case 71:
#line 1084 "glb_parser.y"
{

int flag;  
  
  yyval.ptrq=yyvsp[0].ptrq;
  
  flag=set_exp_list("bgrulescoeff",yyvsp[0].ptrq[0],0);
  if(flag==1) yyerror("Unknown identifier");
  flag=set_exp_list("bgrulechannellist",yyvsp[0].ptrq[1],0);
  if(flag==1) yyerror("Unknown identifier");
 
  glb_free(yyvsp[0].ptrq);
;
    break;}
case 72:
#line 1097 "glb_parser.y"
{
  int flag;  
  yyval.ptrq=yyvsp[0].ptrq;
  flag=set_exp_list("rulescoeff",yyvsp[0].ptrq[0],0);
  if(flag==1) yyerror("Unknown identifier");
  flag=set_exp_list("rulechannellist",yyvsp[0].ptrq[1],0);
  if(flag==1) yyerror("Unknown identifier"); 

  glb_free(yyvsp[0].ptrq);

 
;
    break;}
case 73:
#line 1111 "glb_parser.y"
{
 if(yyvsp[0].nameptr->value==-1) {yyvsp[0].nameptr->value=step_counter("rule");printf("loc step\n"); }
  loc_count=yyvsp[0].nameptr->value;
  printf("loc_count %d \n",loc_count);
  YYERROR;;
    break;}
case 74:
#line 1116 "glb_parser.y"
{ YYERROR;;
    break;}
case 75:
#line 1120 "glb_parser.y"
{
 if(yyvsp[0].nameptr->value==-1) {yyvsp[0].nameptr->value=step_counter("channel");printf("loc step\n"); }
  loc_count=yyvsp[0].nameptr->value;
  printf("cha loc_count %d \n",loc_count); YYERROR;;
    break;}
case 76:
#line 1124 "glb_parser.y"
{ YYERROR;;
    break;}
}
   /* the action file gets copied in in place of this dollarsign */
#line 554 "/usr/share/bison.simple"

  yyvsp -= yylen;
  yyssp -= yylen;
#ifdef YYLSP_NEEDED
  yylsp -= yylen;
#endif

#if YYDEBUG != 0
  if (yydebug)
    {
      short *ssp1 = yyss - 1;
      fprintf (stderr, "state stack now");
      while (ssp1 != yyssp)
	fprintf (stderr, " %d", *++ssp1);
      fprintf (stderr, "\n");
    }
#endif

  *++yyvsp = yyval;

#ifdef YYLSP_NEEDED
  yylsp++;
  if (yylen == 0)
    {
      yylsp->first_line = yylloc.first_line;
      yylsp->first_column = yylloc.first_column;
      yylsp->last_line = (yylsp-1)->last_line;
      yylsp->last_column = (yylsp-1)->last_column;
      yylsp->text = 0;
    }
  else
    {
      yylsp->last_line = (yylsp+yylen-1)->last_line;
      yylsp->last_column = (yylsp+yylen-1)->last_column;
    }
#endif

  /* Now "shift" the result of the reduction.
     Determine what state that goes to,
     based on the state we popped back to
     and the rule number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTBASE] + *yyssp;
  if (yystate >= 0 && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTBASE];

  goto yynewstate;

yyerrlab:   /* here on detecting error */

  if (! yyerrstatus)
    /* If not already recovering from an error, report this error.  */
    {
      ++yynerrs;

#ifdef YYERROR_VERBOSE
      yyn = yypact[yystate];

      if (yyn > YYFLAG && yyn < YYLAST)
	{
	  int size = 0;
	  char *msg;
	  int x, count;

	  count = 0;
	  /* Start X at -yyn if nec to avoid negative indexes in yycheck.  */
	  for (x = (yyn < 0 ? -yyn : 0);
	       x < (sizeof(yytname) / sizeof(char *)); x++)
	    if (yycheck[x + yyn] == x)
	      size += strlen(yytname[x]) + 15, count++;
	  msg = (char *) malloc(size + 15);
	  if (msg != 0)
	    {
	      strcpy(msg, "parse error");

	      if (count < 5)
		{
		  count = 0;
		  for (x = (yyn < 0 ? -yyn : 0);
		       x < (sizeof(yytname) / sizeof(char *)); x++)
		    if (yycheck[x + yyn] == x)
		      {
			strcat(msg, count == 0 ? ", expecting `" : " or `");
			strcat(msg, yytname[x]);
			strcat(msg, "'");
			count++;
		      }
		}
	      yyerror(msg);
	      free(msg);
	    }
	  else
	    yyerror ("parse error; also virtual memory exceeded");
	}
      else
#endif /* YYERROR_VERBOSE */
	yyerror("parse error");
    }

  goto yyerrlab1;
yyerrlab1:   /* here on error raised explicitly by an action */

  if (yyerrstatus == 3)
    {
      /* if just tried and failed to reuse lookahead token after an error, discard it.  */

      /* return failure if at end of input */
      if (yychar == YYEOF)
	YYABORT;

#if YYDEBUG != 0
      if (yydebug)
	fprintf(stderr, "Discarding token %d (%s).\n", yychar, yytname[yychar1]);
#endif

      yychar = YYEMPTY;
    }

  /* Else will try to reuse lookahead token
     after shifting the error token.  */

  yyerrstatus = 3;		/* Each real token shifted decrements this */

  goto yyerrhandle;

yyerrdefault:  /* current state does not do anything special for the error token. */

#if 0
  /* This is wrong; only states that explicitly want error tokens
     should shift them.  */
  yyn = yydefact[yystate];  /* If its default is to accept any token, ok.  Otherwise pop it.*/
  if (yyn) goto yydefault;
#endif

yyerrpop:   /* pop the current state because it cannot handle the error token */

  if (yyssp == yyss) YYABORT;
  yyvsp--;
  yystate = *--yyssp;
#ifdef YYLSP_NEEDED
  yylsp--;
#endif

#if YYDEBUG != 0
  if (yydebug)
    {
      short *ssp1 = yyss - 1;
      fprintf (stderr, "Error: state stack now");
      while (ssp1 != yyssp)
	fprintf (stderr, " %d", *++ssp1);
      fprintf (stderr, "\n");
    }
#endif

yyerrhandle:

  yyn = yypact[yystate];
  if (yyn == YYFLAG)
    goto yyerrdefault;

  yyn += YYTERROR;
  if (yyn < 0 || yyn > YYLAST || yycheck[yyn] != YYTERROR)
    goto yyerrdefault;

  yyn = yytable[yyn];
  if (yyn < 0)
    {
      if (yyn == YYFLAG)
	goto yyerrpop;
      yyn = -yyn;
      goto yyreduce;
    }
  else if (yyn == 0)
    goto yyerrpop;

  if (yyn == YYFINAL)
    YYACCEPT;

#if YYDEBUG != 0
  if (yydebug)
    fprintf(stderr, "Shifting error token, ");
#endif

  *++yyvsp = yylval;
#ifdef YYLSP_NEEDED
  *++yylsp = yylloc;
#endif

  yystate = yyn;
  goto yynewstate;

 yyacceptlab:
  /* YYACCEPT comes here.  */
#ifndef YYSTACK_USE_ALLOCA
  if (yyfree_stacks)
    {
      free (yyss);
      free (yyvs);
#ifdef YYLSP_NEEDED
      free (yyls);
#endif
    }
#endif
  return 0;

 yyabortlab:
  /* YYABORT comes here.  */
#ifndef YYSTACK_USE_ALLOCA
  if (yyfree_stacks)
    {
      free (yyss);
      free (yyvs);
#ifdef YYLSP_NEEDED
      free (yyls);
#endif
    }
#endif    
  return 1;
}
#line 1127 "glb_parser.y"


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
