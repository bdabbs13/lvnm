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
//#include "sbm.h"
#include "wsbm.h"
#include "dynamic-sbm.h"
//#include "mmsbm.h"

//using namespace std;
using std::ifstream;
//using std::ofstream;

//gsl_rng *rng;
//ofstream myfile;

template <class T>
void MCMC(T *Mod, int start, int total, int burnIn, int thin,
	  int extend_max, double qq, int verbose){

    int ii, step_count;
    int extend_count = -1;
    std::vector<double> ll_vec;

    // Main MCMC Loop
    int raw_count = start;
    ll_vec.resize(total,0.0);

    bool converged = false;
    do{
    	step_count = 0;
	extend_count++;

	if(extend_count > 0){
	    if(verbose > 1) Rprintf("Extending Chain...\n");
	}

    	while(step_count < total){
    	    // Perform an MCMC Step
    	    Mod->step();
    	    Mod->adapt();

    	    if((raw_count >= burnIn) && ((raw_count - burnIn) % thin == 0)){
    		//  Save the result every thin iterations
    		Mod->write(step_count);
		ll_vec[step_count] = Mod->LogLike();
		// if(verbose > 2) Rprintf("ll: %.4f\n",ll_vec[step_count]);

    		step_count++;
    	    }
	    raw_count++;
    	}

	converged = convergenceCheck(ll_vec,qq);

	if(converged){
	    if(verbose > 0){
		Rprintf("Converged after %d extensions\n",
			extend_count);
	    }
	}
    }while((!(converged)) && (extend_count < extend_max));

    if(!converged){
	if(verbose > 0)
	    Rprintf("Warning: Failed to converge after %d extensions\n",
		    extend_count);
    }
    // if(verbose > 2) Mod->print();

}



extern "C" {


    // // Function to be called by R
    // void sbm(int *iters, int *nn_t, int *kk_t, int *YY,
    // 	     double *betaPrior, double *eta,
    // 	     double *rBlockMat, int *rBlockMemb,
    // 	     int *burn_t, int *thin_t,
    // 	     int *start_t, int *multi_t,double *logLik,
    // 	     int *extend_max_t, int *shift_t, double *qq_t,
    // 	     double *postMat, int *verbose_t)
    // {

    // 	GetRNGstate();

    // 	int verbose = *verbose_t;

    // 	/*****  INITIALIZATION  *****/
    // 	//  Initializing SBM object
    // 	CSBM *mySBM = new CSBM(*nn_t, *kk_t, *multi_t);
    // 	mySBM->loadDataR(YY, betaPrior, eta);
    // 	mySBM->loadStateR(rBlockMat,rBlockMemb,postMat,logLik);

    // 	//  Currently Initialization in R isn't accurate??
    // 	mySBM->initRandom();

    // 	//  Loading Previous Chain
    // 	int start = *start_t;
    // 	if(start > 0){
    // 	    // Currently do Nothing
    // 	}

    // 	int total = *iters, burnIn = *burn_t, thin = *thin_t;

    // 	//  Convergence Checking Criteria
    // 	double qq = *qq_t;
    // 	int shift_size = *shift_t;
    // 	int extend_max = *extend_max_t;

    // 	//  Run MCMC Chain

    // 	MCMC(mySBM, start, total, burnIn, thin,
    // 	     extend_max, qq, verbose);

    // 	// sbmMCMC(mySBM, start, total, burnIn, thin,
    // 	// 	shift_size, extend_max, qq, //flatTable,
    // 	// 	rBlockMat, rBlockMemb,
    // 	// 	logLik, postMat, verbose);

    // 	delete mySBM;
    // 	PutRNGstate();
    // }


    // void sbmEMout(int *iter_max_t, int *nn_t, int *kk_t, int *YY,
    // 		  double *eta,
    // 		  double *rPosteriorMemb, double *rBlockMat, int *rBlockMemb,
    // 		  double *threshold_t,
    // 		  double *logLik,int *verbose_t)
    // {

    // 	GetRNGstate();

    // 	int verbose = *verbose_t;
    // 	int iter_max = *iter_max_t;
    // 	int multi = 1;
    // 	double threshold = *threshold_t;
    // 	double betaPrior[2] = {0,0};

    // 	/*****  INITIALIZATION  *****/
    // 	//  Initializing SBM object


    // 	// CSBM *mySBM = new CSBM(*nn_t, *kk_t, YY, betaPrior, eta, multi);
    // 	CSBM *mySBM = new CSBM(*nn_t, *kk_t,multi);
    // 	mySBM->loadDataR(YY, betaPrior, eta);

    // 	mySBM->loadStateR(rBlockMat,rBlockMemb,rPosteriorMemb,logLik);
    // 	// mySBM->RLoadSBM(rBlockMat, rBlockMemb);
    // 	mySBM->initPPem(0.1);

    // 	sbmEM(mySBM, iter_max, threshold, rPosteriorMemb, rBlockMat, rBlockMemb,
    // 	      logLik, eta, verbose);

    // 	delete mySBM;
    // 	PutRNGstate();
    // }


    //  Weighted Stochastic Block Model Wrapper Function (R)
    void wsbm_R(int *iters, int *nn_t, int *kk_t, int *YY,
		double *rPriorSender, double *rPriorReceiver,
		double *rPriorBlockMat, double *rPriorBlockMemb,
		double *rBlockMat, int *rBlockMemb,
		double *rSenderEffects, double *rReceiverEffects,
		int *burn_t, int *thin_t,
		int *start_t, int *multi_t,double *logLik,
		int *extend_max_t, double *qq_t,
		double *postMat, double *rHours, int *verbose_t)
    {

	GetRNGstate();

	/*****  Loading Data/Initialization  *****/
	int verbose = *verbose_t;

	//  Initializing WSBM Object
	CWSBM *myWSBM = new CWSBM(*nn_t, *kk_t,*multi_t);

	//  Load Data From R
	myWSBM->loadDataR(YY,*rHours,
			  rPriorSender,rPriorReceiver,
			  rPriorBlockMat,rPriorBlockMemb);

	myWSBM->loadStateR(rBlockMat, rBlockMemb,
			   rSenderEffects, rReceiverEffects,
			   postMat, logLik);


	/*****  Setting Up MCMC Sampler  *****/
	//  MCMC Control Parameters
	int total = *iters, burnIn = *burn_t, thin = *thin_t;
	int start = *start_t;

	//  Convergence Checking Criteria
	double qq = *qq_t;
	int extend_max = *extend_max_t;

	//  Run MCMC Chain
	// wsbmMCMC(myWSBM, start, total, burnIn, thin,
	// 	 shift_size, extend_max, qq, verbose);

	MCMC(myWSBM, start, total, burnIn, thin,
	     extend_max, qq, verbose);


	/*****  Cleaning Up  *****/
	delete myWSBM;
	PutRNGstate();
    }




    //  Dynamic Weighted Stochastic Block Model Wrapper Function (R)
    void dynsbm_R(int *nn_t, int *kk_t, int *TT_t, int *ee_t, int *multi_t,
		  // Data Values 6-8
		  int *YY, int *rTimeMap, double *rHours,
		  // Hyper Prior Parameter Storage 9-12
		  double *rHyperSender, double *rHyperReceiver,
		  double *rHyperBlockMat, double *rPriorBlockMemb,
		  // Prior Parameter Storage 13-16
		  double *rPriorSender, double *rPriorReceiver,
		  double *rPriorBlockMat, int *rBlockMemb,
		  // Parameter Storage 17-19
		  double *rSenderEffects, double *rReceiverEffects,
		  double *rBlockEffects,
		  // Auxiliary Statistic Storage 20-21
		  double *rLogLik, double *rPosteriorMemb,
		  // Model Flags 22-23
		  int *verbose_t, int *update_mmb_t,
		  // MCMC Control Parameters
		  int *iters, int *burnin_t, int *thin_t, int *start_t,
		  int *extend_max_t, double *qq_t)
    {

	GetRNGstate();


	/*****  Loading Data/Initialization  *****/
	int verbose = *verbose_t;

	//  Initializing DynSBM object
	CDynSBM *myDynSBM = new CDynSBM(*nn_t, *kk_t, *TT_t, *ee_t,
					*iters,
					*multi_t);
	//  Load Data From R
	myDynSBM->LoadDataR(YY, rTimeMap, rHours,
			    rHyperSender, rHyperReceiver, rHyperBlockMat);

	myDynSBM->LoadStateR(rSenderEffects,rReceiverEffects,rBlockEffects,
			     rBlockMemb,rPosteriorMemb,*update_mmb_t,
			     rPriorSender,rPriorReceiver,rPriorBlockMat,
			     rPriorBlockMemb,rLogLik);

	/*****  Setting Up MCMC Sampler  *****/
	//  Load Control Parameters
	int total = *iters, burnIn = *burnin_t, thin = *thin_t;
	int start = *start_t;

	//  Convergence Checking Criteria
	double qq = *qq_t;
	int extend_max = *extend_max_t;

	//  Run MCMC Chain
	// dynSBMMCMC(myDynSBM,
	//  	   start, total, burnIn, thin,
	//  	   shift_size, extend_max, qq, verbose);
	MCMC(myDynSBM, start, total, burnIn, thin,
	     extend_max, qq, verbose);

	/*****  Cleaning Up  *****/
	delete myDynSBM;
	PutRNGstate();
    }





//     //  Function to be called by R
//     void mmsbm(int *iters, int *nn_t, int *kk_t, int *YY,
// 	       double *betaPrior, double *alpha,
// 	       double *flatTable, int *burn_t, int *thin_t,
// 	       int *start_t,int *multi_t, double *logLik,
// 	       int *extend_max_t, int *shift_t, double *qq_t, int *verbose_t)
//     {

// 	GetRNGstate();

// 	//  Reading in values from pointers
// 	int total = *iters, verbose = *verbose_t;
// 	int burnIn = *burn_t, thin = *thin_t;
// 	int start = *start_t;//, multiImpute = *multi_t;

// 	mmsbm_t *myMMSBM = new mmsbm_t(*nn_t, *kk_t, YY, betaPrior, alpha, *multi_t);
// 	if(start > 0){
// 	    myMMSBM->loadTable(start, flatTable);
// 	}
// 	double qq = *qq_t;
// 	int shift_size = *shift_t;
// 	int extend_max = *extend_max_t;

// 	mmsbmMCMC(myMMSBM, start, total, burnIn, thin,
// 		  shift_size, extend_max, qq, flatTable, logLik,verbose);

// 	//  Sending RNG state back to R
// 	delete myMMSBM;
// 	PutRNGstate();
//     }

//     void getBlockMatMLE(int *nn_t, int *kk_t, int *YY, double *rBlockMat, int *rBlockMemb){

// 	GetRNGstate();
// 	double betaPrior[2] = {1,1};
// 	double *eta = new double[*kk_t];
// 	int multi = 1;
// 	std::fill(eta,eta + *kk_t, 1.0);

// 	/*****  INITIALIZATION  *****/
// 	//  Initializing SBM object
// 	// CSBM *mySBM = new CSBM(*nn_t, *kk_t, YY, betaPrior, eta, multi);
// 	CSBM *mySBM = new CSBM(*nn_t, *kk_t, multi);
// 	mySBM->loadDataR(YY, betaPrior, eta);

// 	mySBM->RLoadSBM(rBlockMat, rBlockMemb);
// 	//      mySBM->RLoadBlockMemb(mmb);
// 	mySBM->computeBlockMatMLE();
// 	mySBM->GetBlockMat(rBlockMat);
// 	delete mySBM;

// 	PutRNGstate();

//     }

//     void sbmTesting(int *nn_t, int *kk_t, int *adjMat){
// 	double beta[2];
// 	double eta[*kk_t];

// 	// CSBM *mySBM = new CSBM(*nn_t, *kk_t, adjMat, beta, eta, 1);
// 	CSBM *mySBM = new CSBM(*nn_t, *kk_t,1);
// 	mySBM->loadDataR(adjMat, beta, eta);
// 	mySBM->print(true);

//     }
}




