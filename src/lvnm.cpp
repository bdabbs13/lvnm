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
	    //double *flatTable,
	    double *rBlockMat, int *rBlockMemb,
	    int *burn_t, int *thin_t,
	    int *start_t, int *multi_t,double *logLik,
	    int *extend_max_t, int *shift_t, double *qq_t,
	    double *postMat, int *verbose_t)
   {

      GetRNGstate();

      int start = *start_t, verbose = *verbose_t;
      //  MCMC Control Parameters
      int total = *iters, burnIn = *burn_t, thin = *thin_t;

      //  Convergence Checking Criteria
      double qq = *qq_t;
      int shift_size = *shift_t;
      int extend_max = *extend_max_t;


      /*****  INITIALIZATION  *****/
      //  Initializing SBM object
      CSBM *mySBM = new CSBM(*nn_t, *kk_t, YY, betaPrior, eta, *multi_t);

      //  Loading Previous Chain
      if(start > 0){
	 // 1234
	 //mySBM->loadTable(start, flatTable);

	 /* Public should not be able to see these functions
	 if(start == 1){
	    mySBM->drawPP();
	    mySBM->updatePosteriorMemb(0,postMat);
	 }
	 */
      }

      sbmMCMC(mySBM, start, total, burnIn, thin,
	      shift_size, extend_max, qq, //flatTable,
	      rBlockMat, rBlockMemb,
	      logLik, postMat, verbose);

      delete mySBM;
      PutRNGstate();
   }


   void sbmEMout(int *iter_max_t, int *nn_t, int *kk_t, int *YY,
		 double *eta,
		 double *rPosteriorMemb, double *rBlockMat, int *rBlockMemb,
		 double *threshold_t,
		 double *logLik,int *verbose_t)
   {

      GetRNGstate();

      int verbose = *verbose_t;
      int iter_max = *iter_max_t;
      int multi = 1;
      double threshold = *threshold_t;
      double betaPrior[2] = {0,0};

      /*****  INITIALIZATION  *****/
      //  Initializing SBM object

      CSBM *mySBM = new CSBM(*nn_t, *kk_t, YY, betaPrior, eta, multi);
      mySBM->RLoadSBM(rBlockMat, rBlockMemb);
      mySBM->initPPem(0.1);

      sbmEM(mySBM, iter_max, threshold, rPosteriorMemb, rBlockMat, rBlockMemb,
	    logLik, eta, verbose);

      delete mySBM;
      PutRNGstate();
   }



   //  Function to be called by R
   void mmsbm(int *iters, int *nn_t, int *kk_t, int *YY,
	      double *betaPrior, double *alpha,
	      double *flatTable, int *burn_t, int *thin_t,
	      int *start_t,int *multi_t, double *logLik,
	      int *extend_max_t, int *shift_t, double *qq_t, int *verbose_t)
   {

      GetRNGstate();

      //  Reading in values from pointers
      int total = *iters, verbose = *verbose_t;
      int burnIn = *burn_t, thin = *thin_t;
      int start = *start_t;//, multiImpute = *multi_t;

      mmsbm_t *myMMSBM = new mmsbm_t(*nn_t, *kk_t, YY, betaPrior, alpha, *multi_t);
      if(start > 0){
	 myMMSBM->loadTable(start, flatTable);
      }
      double qq = *qq_t;
      int shift_size = *shift_t;
      int extend_max = *extend_max_t;

      mmsbmMCMC(myMMSBM, start, total, burnIn, thin,
		shift_size, extend_max, qq, flatTable, logLik,verbose);

      //  Sending RNG state back to R
      delete myMMSBM;
      PutRNGstate();
   }

   void getBlockMatMLE(int *nn_t, int *kk_t, int *YY, double *rBlockMat, int *rBlockMemb){

      GetRNGstate();
      double betaPrior[2] = {1,1};
      double *eta = new double[*kk_t];
      int multi = 1;
      std::fill(eta,eta + *kk_t, 1.0);

      /*
	double test1 = DBL_EPSILON;
	double test2 = log(test1/2);
	double test3 = log(MIN_LOG/2);
	Rprintf("eps: %f \n log(eps): %f \n log(min): %f",test1,test2,test3);
      */

      /*****  INITIALIZATION  *****/
      //  Initializing SBM object
      CSBM *mySBM = new CSBM(*nn_t, *kk_t, YY, betaPrior, eta, multi);
      mySBM->RLoadSBM(rBlockMat, rBlockMemb);
	 //      mySBM->RLoadBlockMemb(mmb);
      mySBM->computeBlockMatMLE();
      mySBM->GetBlockMat(rBlockMat);
      delete mySBM;

      PutRNGstate();

   }

   void sbmTesting(int *nn_t, int *kk_t, int *adjMat){
      double beta[2];
      double eta[*kk_t];

      CSBM *mySBM = new CSBM(*nn_t, *kk_t, adjMat, beta, eta, 1);
      mySBM->print(true);

   }
}




