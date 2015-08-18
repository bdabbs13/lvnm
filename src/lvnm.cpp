//#include <stdio.h>
//#include <math.h>
//#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <ctime>
#include <cmath>
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>
#include <R.h>
#include <Rmath.h>
#include "lvnm.h"
#include "helper.h"
#include "sbm.h"
#include "mmsbm.h"

//using namespace std;
using std::ifstream;
//using std::ofstream;

//gsl_rng *rng;
//ofstream myfile;


extern "C" {
  
  
  // Function to be called by R
  void sbm(int *iters, int *nn_t, int *kk_t, int *YY,
	   double *betaPrior, double *eta,
	   double *flatTable, int *burn_t, int *thin_t,
	   int *start_t, int *multi_t,double *logLik)
  {
    
    GetRNGstate();

    int start = *start_t;
    int total = *iters, burnIn = *burn_t, thin = *thin_t;

    //    int verbose = 0;
    double qq = 3.2;
    int shift_size = 100, extend_max = 10;


    /*****  INITIALIZATION  *****/
    //  Initializing SBM object
    sbm_t mySBM (*nn_t, *kk_t, YY, betaPrior, eta, *multi_t);
    if(start > 0){
      mySBM.loadTable(start, flatTable);
    }
    
    sbmMCMC(mySBM, start, total, burnIn, thin,
	    shift_size, extend_max, qq, flatTable, logLik);
    
    PutRNGstate();
  }
  


  
  
  //  Function to be called by R
  void mmsbm(int *iters, int *nn_t, int *kk_t, int *YY,
	     double *betaPrior, double *alpha,
	     double *flatTable, int *burn_t, int *thin_t,
	     int *start_t,int *multi_t)
  {
    
    /* Initializing Random Number Generator */
    //  time_t seed = time(NULL);
    //  rng = gsl_rng_alloc(gsl_rng_mt19937); /* start the RNG */
    //  gsl_rng_set(rng, static_cast<unsigned long int>(seed));
    GetRNGstate();
    Rprintf("this is the correct function\n");
    //  Reading in values from pointers
    int total = *iters;
    int burnIn = *burn_t, thin = *thin_t;
    int start = *start_t, multiImpute = *multi_t;
    int dd = *kk_t, nn = *nn_t;
    
    //  Allocating memory for variables
    int yyComplete[nn*nn];
    double BB[dd*dd];
    double PP[nn*dd];
    int sendMat[nn*nn]; int recMat[nn*nn];
    
    //  Performing MCMC Algorithm
    mmsbmMCMC(total,burnIn,thin,YY,nn,dd,alpha,betaPrior,
	      BB,PP,sendMat,recMat,yyComplete,
	      flatTable,start,multiImpute);
    
    //  Sending RNG state back to R
    PutRNGstate();
  }
}



