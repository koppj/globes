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
#include <time.h>
#include <string.h>
#include <argp.h>
#include <ctype.h>
#include <ltdl.h>


#include <globes/globes.h>

#if HAVE_CONFIG_H
#  include "config.h"
#endif


#ifdef TEST
#include "glb_modules.h"
#endif /* TEST */

const char *argp_program_version =
"globes "VERSION"\n(C) 2002 - 2005 The GLoBES Team\n"
"This is free software see the source for copying conditions. There is NO\n"
"warranty; not even for MERCHANTABILITY or"
" FITNESS FOR A PARTICULAR PURPOSE.";
const char *argp_program_bug_address =
"<globes@ph.tum.de>";

static void parse_definition(const char *in)
{
  const char *delim="=";
  double val;
  char *token=NULL;
  char *lhs=NULL;
  char *rhs=NULL;
  char *inc=NULL;
  size_t length=0;

  inc=strdup(in);
  if(inc==NULL) return;
  token=strtok(inc,delim);
  if(token!=NULL) 
    {
     lhs=strdup(token);
     length++;
    }
  token=strtok(NULL,delim); 
  if(token!=NULL) 
    {
     rhs=strdup(token);
     length++;
    }
  
  if(length!=2) {
    glb_free(lhs);
    glb_free(rhs);
    fprintf(stderr,"ERROR: Definition is not of form 'DEFINITION=VALUE'\n");
    return;}
  val=atof(rhs);
  glbDefineAEDLVariable(lhs,val);
  glb_free(lhs);
  glb_free(rhs);
  return;
}

/* Program documentation. */
static char doc[] ="Rate computation for GLoBES";

/* A description of the arguments we accept. */
static char args_doc[] = "glb-file";

/* The options we understand. */
static struct argp_option options[] ={
  {"channel",'c',  "NUMBER", OPTION_ARG_OPTIONAL,  
   "Show rates for channel given by NUMBER,\n"
   "   if no argument is given all channels are shown" }, 
  {"rule",'r', "NUMBER", OPTION_ARG_OPTIONAL, 
   "Show rates for a rule given by NUMBER,\n"
   "   if no argument is given all rules are shown" }, 
  {"experiment",'e', "NUMBER",0 ,   
   "For multiple experiments, chose experiment number " },
  {"output",   'o', "FILE", 0,
   "Output to FILE instead of standard output" },
  {"parameters",'p',"PARAMETERS",0,
   "Set of oscillation parameters\n   for which the rates are computed\n"
   "   Format 'th12,th13,th23,delta,dm21,dm31'"},
  {"before",'b',0,0,"Channel rates before smearing (implies -c)"},
  {"after",'a',0,0,"Rates after smearing"},
  {"efficiencies",'f',0,0,"Rates without efficiencies"},
  {"backgrounds",'g',0,0,"Rates without backgrounds"},
  {"spectrum",'s',0,0,"Differential event rates, i.e. spectrum"},
  {"total",'t',0,0,"Total event rates"},
  {"mathematica",'m',0,0,"Output in Mathematica (TM) list format"},
  {"coefficients",'i',0,0,"Rule rates without rule coefficients"},
  {"verbosity",'v',"LEVEL",OPTION_ARG_OPTIONAL,"Verbosity level for warnings"
   " and errors, 0 everthing off, 1 errors on (default), 2 warnings"
   " on, 3 debug info, 4 more debug info"},
  {"pretty-printing",'P',0,0,"pretty-printing - default"},
  {"simple-printing",'S',0,0,"simple-printing"},
  {"Oscillation",'O',0,0,"Standards oscillations"},
  {"No-oscillation",'N',0,0,"Oscillations switched off"},
  {"Define",'D',"DEFINITION",0,"Define AEDL variable"},
  {"Left",'L',"STRING",0,"Left delimiter used in formatting output"},
  {"Right",'R',"STRING",0,"Right delimiter used in"
   " formatting output"},
  {"Middle",'M',"STRING",0,"Middle delimiter used in"
   " formatting output"},
  {"user",'u',0,0,"Output formatting in user defined mode"},
  {"wahrscheinlichkeit",'w',0,0,"probability"},
  { 0 } 
};

/* Used by `main' to communicate with `parse_opt'. */
struct arguments
{
  char* args[1];                /* many arguments*/
  int channel,rule,experiment,user;
  int smearing,eff,bg,spectrum,mathematica,coeff,verbosity,pretty,oscillation,probability;
  char *output_file;
  char* params;
  char *left,*middle,*right;
};


/* Parse a single option. */
static error_t
parse_opt (int key, char *arg, struct argp_state *state)
{
  /* Get the INPUT argument from `argp_parse', which we
     know is a pointer to our arguments structure. */
  struct arguments *arguments = state->input;
  
  switch (key)
    {
    case 'v':
      arguments->verbosity = arg ? atoi(arg) : 1;
      break;
    case 'i':
      arguments->coeff = GLB_WO_COEFF;
      break;
    case 'c':
      arguments->channel = arg ? atoi(arg) : GLB_ALL;
      arguments->rule=-2;
      break;
    case 'o':
      arguments->output_file = arg;
      break;

    case 'u':
      arguments->user = 1;
      break;
    case 'p':
      arguments->params = arg;
      break;
    case 'e':
      arguments->experiment = arg ? atoi(arg) : 0;
      break;
    case 'r':
      arguments->rule = arg ? atoi(arg) : GLB_ALL;
      break;
    case 'f':
      arguments->eff = GLB_WO_EFF;
      break;
    case 'g':
      arguments->bg = GLB_WO_BG;
      break;
    case 'a':
      arguments->smearing = GLB_POST;
      break;
    case 'b':
      arguments->smearing = GLB_PRE;
      break;
    case 's':
      arguments->spectrum = 1;
      break;
    case 't':
      arguments->spectrum = 0;
      break;
    case 'S':
      arguments->pretty = 0;
      break;
    case 'P':
      arguments->pretty = 1;
      break;
    case 'm':
      arguments->mathematica = 1;
      break;
    case 'O':
      arguments->oscillation = 1;
      break;
    case 'N':
      arguments->oscillation = 0;
      break;
    case 'D':
      parse_definition(arg);
      break;
    case 'L':
      arguments->left = arg;
      break;
    case 'M':
      arguments->middle = arg;
      break;
    case 'R':
      arguments->right = arg;
      break;
    case 'w':
      arguments->probability = 1;
      break;
      

      
    case ARGP_KEY_ARG:
      if (state->arg_num > 1)
	/* Too many arguments. */
	argp_usage (state);
     
      arguments->args[state->arg_num]=arg;
      
      break;
      
    case ARGP_KEY_END:
      if (state->arg_num < 1)
	/* Not enough arguments. */
	argp_usage (state);
      break;
      
    default:
      return ARGP_ERR_UNKNOWN;
    }
  return 0;
}

/* Our argp parser. */
static struct argp argp = { options, parse_opt, args_doc, doc };




static void glb_channel_printf_sub_total(FILE *stream,
				       const double *energy,
				       const double **res, size_t l,size_t c)
{
  int i,k;
  double *sum;
  sum=(double *) malloc(sizeof(double)*c);
  if(sum==NULL) 
    {
      fprintf(stderr,"glb: FATAL: Memory allocation failed\n");
      exit(1);
    }
  for(k=0;k<c;k++) sum[k]=0.0;
  for(i=0;i<l;i++)
    for(k=0;k<c;k++)
      sum[k]+=res[k][i];
  
  fprintf(stream,"\n");
  for(k=0;k<c;k++)
    {
      fprintf(stream,"%12.6g\t",sum[k]);
    }
  fprintf(stream,"\n\n");
  free(sum);
}

static void glb_channel_printf_total(FILE *stream,
				       const double *energy,
				       const double **res, size_t l,size_t c)
{
  int i,k;
  double *sum,tsum;
  sum=(double *) malloc(sizeof(double)*c);
  if(sum==NULL) 
    {
      fprintf(stderr,"glb: FATAL: Memory allocation failed\n");
      exit(1);
    }
  for(k=0;k<c;k++) sum[k]=0.0;
  for(i=0;i<l;i++)
    for(k=0;k<c;k++)
      sum[k]+=res[k][i];
  tsum=0.0;
  for(k=0;k<c;k++) tsum += sum[k];
  fprintf(stream,"%12.6g ||\t",tsum);
  for(k=0;k<c;k++)
    {
      fprintf(stream,"%12.6g\t",sum[k]);
    }
  fprintf(stream,"\n");
  free(sum);
}



static void glb_channel_printf_mathematica(FILE *stream,
				       const double *energy,
				       const double **res, size_t l,size_t c)
{
  int i,k;
  fprintf(stream,"\n");
  glbPrintDelimiter(stream,'l');
  for(k=0;k<c;k++)
    {
       glbPrintDelimiter(stream,'l');
      for(i=0;i<l;i++)
	{  
	   glbPrintDelimiter(stream,'l');
	  fprintf(stream,"%f",energy[i]);
	   glbPrintDelimiter(stream,'m');
	  fprintf(stream,"%f",res[k][i]);
	   glbPrintDelimiter(stream,'r');
	  if(i<l-1)  glbPrintDelimiter(stream,'m');
	}
       glbPrintDelimiter(stream,'r');
      if(k<c-1)  glbPrintDelimiter(stream,'m');
     
      
    }
   glbPrintDelimiter(stream,'r');
  fprintf(stream,"\n");
}

static void rule_name(FILE *stream,int exp,int number)
{
  fprintf(stream,"--------- %s ---------\n",glbValueToName(exp,"rule",number));
}

static void channel_name(FILE *stream,int exp,int number,int spec)
{
  size_t c,i;
  
  c=glbGetNumberOfChannels(exp);
  if(spec!=1) fprintf(stream,"\n      \t");
  else fprintf(stream,"\n          \t");

  for(i=0;i<c;i++) 
    {
      if(number!=GLB_ALL) i=number;
      fprintf(stream,"%s\t",glbValueToName(exp,"channel",i));
      glbPrintDelimiter(stream,'m');
      if(number!=GLB_ALL) break;
    }

  fprintf(stream,"\n");
  
  return;
}


// ---------------------------------------------
// -----           Main-Routine            -----
// ---------------------------------------------

int main(int argc, char *argv[])
{  
 

#ifdef TEST
  glb_projection pro;
  size_t layers;
  double *len,*den;
  double a,b;
  double *ll,*dd;
  double (*prf)(const glb_params);
#endif /* TEST */
  char *central_values;
  void *print_buf;
  FILE *stream;
  int i,rv,s;
  struct arguments arguments;
  glb_params oscp;
  double osc[]={0.553574,0.160875,M_PI/4,0.0,0.00007,0.003};
  
  /* This serves to make globes working, even with --disable-shared */
  //  LTDL_SET_PRELOADED_SYMBOLS(); 

  arguments.args[0]="-";
  arguments.channel=GLB_ALL;
  arguments.rule=GLB_ALL;
  arguments.experiment=0;
  arguments.bg=GLB_W_BG;
  arguments.eff=GLB_W_EFF;
  arguments.coeff=GLB_W_COEFF;
  arguments.smearing=GLB_POST;
  arguments.params=NULL;
  arguments.output_file=NULL;
  arguments.spectrum=0;
  arguments.mathematica=0;
  arguments.verbosity=1;
  arguments.pretty=1;
  arguments.oscillation=1;
  arguments.left="";
  arguments.right="\n";
  arguments.middle="\t";
  arguments.user=0;
  arguments.probability=0;
 


  /* Parse our arguments; every option seen by `parse_opt' will
     be reflected in `arguments'. */
  argp_parse (&argp, argc, argv, 0, 0, &arguments);  

  /* Initialize libglobes */
  glbInit(argv[0]);
  glbSetPrintDelimiters(arguments.left,arguments.middle,arguments.right);
  glbSetVerbosityLevel(arguments.verbosity);
  
  /* Redirecting the output stream according to -o=FILE */
  if(arguments.output_file==NULL) stream=stdout;
  else
    {
      stream=fopen(arguments.output_file,"w");
      if(stream==NULL) 
	{ 
	  fprintf(stderr,"%s: FATAL: Could not open file for output\n"
		  ,argv[0]);
	  exit(1);
	}
    }

  /* Lexing the oscillation parameters as given by the environment variable
   * GLB_CENTRAL_VALUES
   */
  central_values=getenv("GLB_CENTRAL_VALUES");
  if(central_values!=NULL)
    {
      rv=sscanf(central_values,"%lf , %lf , %lf , %lf , %lf , %lf",&osc[0],
	     &osc[1],&osc[2],&osc[3],&osc[4],&osc[5]);
      if(rv!=glbGetNumOfOscParams()) 
	{
	  fprintf(stderr,"%s: FATAL: Wrong format for oscillation"
		  " parameters from environment variable GLB_CENTRAL_VALUES\n"
		  ,argv[0]);
	  exit(1);
	}
    }


  /* Lexing the oscillation parameters as given by -p='1,2,3...' */
  if(arguments.params!=NULL)
    {
      rv=sscanf(arguments.params,"%lf , %lf , %lf , %lf , %lf , %lf",&osc[0],
	     &osc[1],&osc[2],&osc[3],&osc[4],&osc[5]);
      if(rv!=glbGetNumOfOscParams()) 
	{
	  fprintf(stderr,"%s: FATAL: Wrong format for oscillation"
		  " parameters\n ",argv[0]);
	  exit(1);
	}
    }
 
#ifdef TEST
  
  fprintf(stderr,"... seems to work\n now to something real ...\n");
  glb_load_prior("prior-template.la",&prf);
  
  
  glbDefineAEDLVariable("REP",1.0);

#endif /* TEST */
  /* Processing the argument file */
  s=glbInitExperiment(arguments.args[0],&glb_experiment_list[0],
		      &glb_num_of_exps);

  
  

  /* Testing for failure */
  if(s<-1) {fprintf(stderr,"%s: FATAL: Unrecoverable parse error\n",
		    argv[0]);exit(1);}

#ifdef TEST
  
  
  fprintf(stdout,"x-sec %lf\n",glbXSection(0,0,0.005,1,-1));
  fprintf(stdout,"x-sec %lf\n",glbXSection(0,0,0.01,1,-1));
#endif /* TEST */


  /* Computing the rates */
  oscp=glbAllocParams();

  oscp=glbDefineParams(oscp,osc[0],
		       osc[1],osc[2],osc[3],osc[4],osc[5]);
  
  if(arguments.oscillation==0) oscp=glbDefineParams(oscp,0,
		       0,0,0,osc[4],osc[5]);
  
  glbSetOscillationParameters(oscp);
  
  glbSetRates();

  
#ifdef TEST

  printf("fluxes %d\n",glbGetNumberOfFluxes(0));
  printf("fluxes %d\n",glbGetNumberOfFluxes(1));
  printf("fluxes %d\n",glbGetNumberOfFluxes(2));

  printf("probs %g\n",glbProfileProbability(0,1,3,1,10.0));
  printf("probs %g\n",glbProfileProbability(1,1,3,1,10.0));
  fprintf(stderr,"glbNameToValue ... %d\n",glbNameToValue(0,"channel","#ch0"));
  fprintf(stderr,"glbValueToName ... %s\n",glbValueToName(0,"flux",0));
#endif /* TEST */

 
  /* Chosing the right outputformat */
  print_buf= glbSetChannelPrintFunction(NULL);
  if(arguments.spectrum==0)
    {
      if(arguments.rule!=-2)
	glbSetChannelPrintFunction(glb_channel_printf_total);
      else
	glbSetChannelPrintFunction(glb_channel_printf_sub_total);
    }
  
  if(arguments.spectrum==1) 
    glbSetChannelPrintFunction(print_buf);
  if(arguments.mathematica==1) 
    {
      glbSetChannelPrintFunction(glb_channel_printf_mathematica);
      glbSetPrintDelimiters("{",",","}");
    }

  if(arguments.user==1) 
    {
      glbSetChannelPrintFunction(glb_channel_printf_mathematica);     
    }
 
  /* Probability-only output */
  if(arguments.probability==1)
    {

      if(arguments.pretty==1) channel_name(stream,arguments.experiment,
					   arguments.channel,
					   arguments.spectrum);
      glbShowChannelProbs(stream,arguments.experiment,arguments.channel,
			  arguments.smearing,arguments.eff,arguments.bg);
      
      /* Cleaning up */
      glbFreeParams(oscp);
      if(stream!=stdout) fclose(stream);
      exit(0);  
    }


  /* Displaying the channel rates */
  if(arguments.rule==-2) 
    {
      if(arguments.pretty==1) channel_name(stream,arguments.experiment,
					   arguments.channel,
					   arguments.spectrum);
      glbShowChannelRates(stream,arguments.experiment,arguments.channel,
		      arguments.smearing,arguments.eff,arguments.bg);
    }
  
  /* Displaying the rule level rates */
  if(arguments.rule!=-2)
    {     
      if(arguments.rule!=GLB_ALL)
	{
	 if(arguments.pretty==1) rule_name(stream,arguments.experiment,
					   arguments.rule);
	   glbPrintDelimiter(stream,'l');
	  glbShowRuleRates(stream,arguments.experiment,arguments.rule,GLB_ALL,
			   arguments.eff,arguments.bg,arguments.coeff,GLB_SIG);
	   glbPrintDelimiter(stream,'m');
	  glbShowRuleRates(stream,arguments.experiment,arguments.rule,GLB_ALL,
			   arguments.eff,arguments.bg,arguments.coeff,GLB_BG);
	   glbPrintDelimiter(stream,'r');
	}
      else
	{
	  
	   glbPrintDelimiter(stream,'l');
	  fprintf(stream,"\n");
	  for(i=0;i<glbGetNumberOfRules(arguments.experiment);i++)
	    {
	      fprintf(stream,"\n");
	      if(arguments.pretty==1) rule_name(stream,arguments.experiment,i);
	      glbPrintDelimiter(stream,'l');
	      fprintf(stream,"\t");
	      glbShowRuleRates(stream,arguments.experiment,i,
			       GLB_ALL,arguments.eff,arguments.bg,
			       arguments.coeff,GLB_SIG);
	       glbPrintDelimiter(stream,'m');
	      glbShowRuleRates(stream,arguments.experiment,i,
			       GLB_ALL,arguments.eff,arguments.bg,
			       arguments.coeff,GLB_BG);
	      glbPrintDelimiter(stream,'r');
	      if(i< glbGetNumberOfRules(arguments.experiment)-1) 
		 glbPrintDelimiter(stream,'m');
	   
	    }
	   glbPrintDelimiter(stream,'r');
	}
    }

  /* Cleaning up */  
  glbFreeParams(oscp);
  if(stream!=stdout) fclose(stream);
  exit(0); 
}


