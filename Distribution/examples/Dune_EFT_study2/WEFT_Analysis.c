#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>




#include <globes/globes.h>   /* GLoBES library */
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include "smeft.h"

#define FLAVS 4
#define ROTS 6

#define EXP_FAR  0
#define EXP_NEAR 1
int number_param = 600;


int smeft_remove(glb_params p){
    
    for(int i = 6; i < number_param; i++){
       
            glbSetOscParams(p,0.0,i);
            
        
    }
    return 0;
}

int get_index(char smeft_param_strings[600][64], int size, char *str) {
    for(int i = 0; i < size; i++) {
        if(strcmp(smeft_param_strings[i], str) == 0)
            return i;
    }
    return -1;
}


int main(int argc, char *argv[]) {


   glbInit(argv[0]);	

  /* Intitialize output */

  
  /* Define my standard oscillation parameters */
    double theta12 = 0.5903;
    double theta13 = 0.1503;
    double theta23 = 0.8674;
    double deltacp = 217*M_PI/180;
    double sdm = 7.39e-5;
    double ldm = 2.525e-3;
    
  
	
   /* Initialize one experiment */
   
    smeft_init_probability_engine_3();
    glbRegisterProbabilityEngine(number_param, &smeft_probability_matrix, &smeft_set_oscillation_parameters, &smeft_get_oscillation_parameters, NULL);
   
    //glbInitExperiment("DUNE_GLoBES_FD.glb",&glb_experiment_list[0],&glb_num_of_exps);
    glbInitExperiment("DUNE_GLoBES_ND.glb",&glb_experiment_list[0],&glb_num_of_exps);
    
    /* Initialize a number of parameter vector(s) */
    glb_params true_values = glbAllocParams();
    
    //glb_params central_values = glbAllocParams();
    //glb_params input_errors = glbAllocParams();
    
    glbDefineParams(true_values,theta12,theta13,theta23,deltacp,sdm,ldm);
    //glbDefineParams(test_values,theta12,theta13,theta23,deltacp,sdm,ldm);
    //glbDefineParams(central_values,theta12,theta13,theta23,deltacp,sdm,ldm);
    
    
    
    //glbDefineParams(input_errors,theta12*0.023,theta13*0.015,theta23*0.022,0,sdm*0.028,ldm*0.013);
    smeft_remove(true_values);
    //smeft_remove(test_values);
    //smeft_remove(input_errors);
    //smeft_remove(central_values);
   //glbSetDensityParams(test_values,1.0,GLB_ALL);
    glbSetOscillationParameters(true_values);
    //glbSetDensityParams(input_errors,0.02,GLB_ALL);
    //glbSetCentralValues(true_values);
    //glbSetInputErrors(input_errors);
    
    //glbSetOscillationParameters(true_values);
    
    
    //glbSetRates();
    
  /*FILE *Signal_rate = fopen("./results/signal.txt", "w");
  FILE *bkg_rate = fopen("./results/bkg.txt", "w");
  
  FILE *Signal_rate2 = fopen("./results/signal_NSI.txt", "w");
  FILE *bkg_rate2 = fopen("./results/bkg_NSI.txt", "w");*/
  
  
     
   /* for (int k=0; k < 600; k++){
    
      fprintf(Flux,"%s\n", smeft_param_strings[k]);
      
    }*/




  /*int index = get_index(smeft_param_strings, 600, "ABS_EPS_CC_L_ud_EE");
  
  printf( "%d \n", index);*/
    
    
  
    
    
    
    
  // snprintf(buffer, sizeof(char) * 150, ".//results//Final_results.txt");
   
  
  
 double min = log10(100);
    double max = log10(100);
    //double min = log10(1);
    //double max = log10(2);
    
    int steps = 1;
   double chi2 =0;
   glb_params test_values = glbAllocParams();
   

  //glbSetOscParams(true_values,1,1);
  glbSetRates();
   //int a =  glbShowRuleRates(Signal_rate,0, 0, GLB_ALL,GLB_W_EFF,GLB_W_BG, GLB_W_COEFF, GLB_SIG);
   //int b =  glbShowRuleRates(bkg_rate,0, 0, GLB_ALL,GLB_W_EFF,GLB_W_BG, GLB_W_COEFF, GLB_BG);
   glbCopyParams(true_values,test_values);
    //int b =  glbShowRuleRates(Flux2,0, 0, GLB_ALL,GLB_W_EFF,GLB_W_BG, GLB_W_COEFF, GLB_SIG);
    
    const char *laptonflavors[] = { "E", "MU", "TAU", "S1", "S2", "S3", "S4", "S5", "S6" };
  const char *upquarks[] = { "u", "c" };
  const char *downquarks[] = { "d", "s", "b" };
  const char *matterfermions[] = { "E", "u", "d" };
  const char *interactions[] = { "L", "R", "S", "P", "T" };
  char name[40];

//FILE *Flux = fopen("./results/CC_P_us_new_resultsBOTH_test.txt", "w"); 
            
    for (int l = 0; l < 1; l++) {
        for (int m = 0; m < 1; m++) {
            for (int i = 0; i < 1; i++) {
                //for (int j = 0; j < 3; j++) {
                int j = i;
                    sprintf(name, "ABS_EPS_CC_%s_%s%s_%s%s", interactions[3], upquarks[l], downquarks[m], laptonflavors[i], laptonflavors[j]);
                    int index = get_index(smeft_param_strings, 600, name);
                    for(int l1 = 0; l1< 1;l1++) {
    	  		 double logval = min + l1 * (max - min) / steps;
           		 //double val = pow(10, logval);
           		 double val = 1e-3;
           		  smeft_remove(test_values);
            	         glbSetOscParams(test_values,val,index);
            
                         chi2 = glbChiSys(test_values,GLB_ALL,GLB_ALL);
                         printf("%s %f %g \n",name,val,chi2);
                         //fprintf(Flux,"%s %f %g \n",name,val,chi2); 
                         /*if (chi2 > 10) {
                           break;
                        }*/
    
    
    		    }
                    
                }
            }
        //}
    }
    
    
    double E = 2.5;
      printf("E=%g, e det. coeff=(%g)\n",
         E,
         glbEFTFluxCoeff(0, 0, 3, 3, 0, E));
          printf("E=%g, mu det. coeff=(%g)\n",
         E,
         glbEFTFluxCoeff(0, 0, 3, 3, 1, E));
    
    
    
/*FILE *Flux = fopen("./results/check_tau_detection.txt", "w"); 
    sprintf(name, "ABS_EPS_CC_%s_%s%s_%s%s", interactions[3], upquarks[0], downquarks[0], laptonflavors[2], laptonflavors[0]);
     //sprintf(name, "ABS_EPS_CC_%s_%s%s_%s%s", interactions[3], upquarks[0], downquarks[0], laptonflavors[0], laptonflavors[0]);
                    int index = get_index(smeft_param_strings, 600, name);
                    for(int l1 = 0; l1<= 10;l1++) {
                  int l1 =0;
    	  		 double logval = min + l1 * (max - min) / steps;
           		 double val = pow(10, logval);
           		  smeft_remove(test_values);
            	         glbSetOscParams(test_values,val,index);
            
                         chi2 = glbChiSys(test_values,GLB_ALL,GLB_ALL);
                         glbSetOscParams(true_values,val,index);
                          glbSetRates();
                          int a2 =  glbShowRuleRates(Signal_rate2,0, 0, GLB_ALL,GLB_W_EFF,GLB_W_BG, GLB_W_COEFF, GLB_SIG);
   int b2 =  glbShowRuleRates(bkg_rate2,0, 0, GLB_ALL,GLB_W_EFF,GLB_W_BG, GLB_W_COEFF, GLB_BG);
                         printf("%s %f %g \n",name,val,chi2);
                         fprintf(Flux,"%s %f %g \n",name,val,chi2); 
                         if (chi2 > 10) {
                           break;
                        }
                    }
                    

*/
   
    	
     
    
   

   
   // fclose(Flux);
    //fclose(Signal_rate);
    //fclose(bkg_rate);
    exit(0);



     
    
}




