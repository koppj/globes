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
    int i;
    for(i = 6; i < number_param; i++){
       
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
   
    glbInitExperiment("DUNE_GLoBES_FD.glb",&glb_experiment_list[0],&glb_num_of_exps);
    glbInitExperiment("DUNE_GLoBES_ND.glb",&glb_experiment_list[0],&glb_num_of_exps);
    
    /* Initialize a number of parameter vector(s) */
    glb_params true_values = glbAllocParams();
    glb_params test_values = glbAllocParams();
    glb_params central_values = glbAllocParams();
    glb_params input_errors = glbAllocParams();
    
    glbDefineParams(true_values,theta12,theta13,theta23,deltacp,sdm,ldm);
    glbDefineParams(test_values,theta12,theta13,theta23,deltacp,sdm,ldm);
    glbDefineParams(central_values,theta12,theta13,theta23,deltacp,sdm,ldm);
    
    
    
    glbDefineParams(input_errors,theta12*0.023,theta13*0.015,theta23*0.022,0,sdm*0.028,ldm*0.013);
    smeft_remove(true_values);
    smeft_remove(test_values);
    smeft_remove(input_errors);
    smeft_remove(central_values);
   glbSetDensityParams(true_values,1.0*0,GLB_ALL);
    glbSetDensityParams(test_values,1.0*0,GLB_ALL);
    glbSetDensityParams(central_values,1.0*0,GLB_ALL);
    glbSetOscillationParameters(true_values);
    glbSetDensityParams(input_errors,0.02*0,GLB_ALL);
    glbSetCentralValues(central_values);
    glbSetInputErrors(input_errors);
    
    glbSetOscillationParameters(true_values);
    glbSetRates();
    glbSetOscillationParameters(test_values);
    char buffer[150];
    // snprintf(buffer, sizeof(char) * 150, ".//results//parameters.txt");
    FILE*Flux = fopen(buffer,"w");
     
   /* for (int k=0; k < 600; k++){
    
      fprintf(Flux,"%s\n", smeft_param_strings[k]);
      
    }*/




  /*int index = get_index(smeft_param_strings, 600, "EPS_NC_R_dd_TAUTAU");
  
  printf( "%d \n", index);*/
    
    
  
    
    
    
    
   snprintf(buffer, sizeof(char) * 150, ".//results//Final_results.txt");
   
  int d =  glbShowRuleRates(Flux, 0, 0, GLB_ALL, GLB_W_COEFF, GLB_W_EFF, GLB_W_BG, GLB_SIG);
  
 /* double min = log10(0.0001);
    double max = log10(0.2);
    int steps = 10;
   double chi2 =0;GLB_W_EFF
  
    for(int i = 6; i< 7;i++){
    	for(int j = 0; j <= steps; j++) {
    	   double logval = min + i * (max - min) / steps;
           double val = pow(10, logval);
           printf("Hello\n");
           glbSetOscParams(test_values,val,i);
            printf("Hello\n");
            chi2 = glbChiSys(test_values,GLB_ALL,GLB_ALL);
             
            
            
            fprintf(Flux,"%d %d %g \n",i,j,chi2); 
    
    
    }
    }
   
     
    
    */

   
    fclose(Flux);
    exit(0);



     
    
}




