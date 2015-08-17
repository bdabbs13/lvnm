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

    int total = *iters;
    int burnIn = *burn_t, thin = *thin_t;
    int start = *start_t, multiImpute = *multi_t;
    int nn = *nn_t, dd = *kk_t;
    double BB[dd*dd], PP[nn*dd];
    int yyComplete[nn*nn], PPint[nn];
    

    sbmMCMC(total, burnIn, thin, YY, nn, dd, eta, betaPrior, logLik,
	    BB, PP, PPint, yyComplete, flatTable, start, multiImpute);

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






  /********************************************
   ********  MCMC CONTROL FUNCTIONS  **********
   *******************************************/

  /**********  SBM  **********/
  
  void sbmInit(int *YY, int nn, int dd, double *eta, double *betaPrior, 
	       double *BB, double *PP, int *PPint, int *yyComplete,
	       double *logLik, double *flatTable, int start, int multiImpute){
    
    int ii,kk;
    //    Rprintf("start: %d\n ",start);
    if(start > 0){
    //      Rprintf("start is greater than 0\n");
      sbmLoadTable(start,nn,dd,BB,PP,PPint,flatTable);
      
      // Creating Initial Imputation if Necessary
      if(multiImpute==1){
	for(ii = 0 ; ii < nn*nn ; ii++){
	  if(yyComplete[ii] < 0){
	    yyComplete[ii] = 0;
	  }else{
	    yyComplete[ii] = YY[ii];
	  }
	}
	sbmImputeMissingValues(nn,dd,YY,yyComplete,BB,PP,PPint);
      }else{
	for(ii = 0 ; ii < nn*nn ; ii++){
	  yyComplete[ii] = YY[ii];
	}	
      }
    }else{
      
      
      //  Initializing Complete YY Matrix
      if(multiImpute==1){
	for(ii = 0 ; ii < nn*nn ; ii++){
	  if(yyComplete[ii] < 0){
	    yyComplete[ii] = 0;
	  }else{
	    yyComplete[ii] = YY[ii];
	  }
	}
      }else{
	for(ii = 0 ; ii < nn*nn ; ii++){
	  yyComplete[ii] = YY[ii];
	}	
      }
      //  Initializing BB Matrix
      double holder[2];
      for(ii = 0 ; ii < dd*dd ; ii++){
	rdirichlet(2,betaPrior,holder);
	BB[ii] = holder[0];
      }
      
      // Initializing PP Matrix
      int ind[dd];
      for(ii = 0 ; ii < nn ; ii++){
	rmultinom(1,eta,dd,ind);
	
	// Setting PP[ii,] Value
	for(kk = 0 ; kk < dd ; kk++){
	  if(ind[kk] == 1){
	    PP[ii*dd + kk] = 1.0;
	    PPint[ii] = kk;
	  }else{
	    PP[ii*dd + kk] = 0.0;
	  }
	}   
      }
    }
  }
  void sbmStep(int *YY, int nn, int dd, double *eta, double *betaPrior,
	       double *BB, double *PP, int *PPint){
    //  Drawing new BB Matrix Values
    sbmDrawBB(nn, dd, YY, PP, BB, betaPrior, PPint);
    //  Drawing new membership vectors
    sbmDrawPP(nn, dd, YY, BB, PP, eta, PPint);
  }
  


  void sbmImputeMissingValues(int nn, int dd, int *YY, int *yyComplete,
			      double* BB, double *PP, int *PPint){
    int rr, cc, index;
    for(rr = 0 ; rr < nn ; rr++){ // row loop
      for(cc = 0 ; cc < nn ; cc++){ // col loop
	if(rr != cc){ // ignore diagonal entries
	  if(YY[rr*nn + cc] < 0){ // missing values coded -1
	    index = PPint[rr] * dd + PPint[cc];
	    yyComplete[rr*nn + cc] = (int) rbinom(1.0,BB[index]);
	  }
	}
      }
    }
  }
  
  void sbmLoadTable(int start, int nn, int dd, double *BB, 
		    double *PP, int *PPint, double *flatTable){


    int ii, jj, offset;
    //  skip first (start - 1) sets of draws
    offset = (start - 1) * (dd * (nn + dd));

    //  Loading BB
    for(ii = 0 ; ii < dd*dd ; ii++){
      BB[ii] = flatTable[offset + ii];
    }

    //  Loading PP
    offset = offset + dd*dd;
    for(ii = 0 ; ii < nn*dd ; ii++){
      PP[ii] = flatTable[offset + ii];
    }
    //  Loading PPint
    for(ii = 0 ; ii < nn ; ii++){
      jj = 0;
      while(PP[ii*dd + jj] != 1.0){
	jj++;
      }
      PPint[ii] = jj;
    }

  }
  

  void sbmMCMC(int total,int burnIn,int thin,
	       int *YY,int nn,int dd,
	       double *eta,double *betaPrior, double *logLik,
	       double *BB, double *PP, int *PPint, int *yyComplete,
	       double *flatTable,int start, int multiImpute){
    
    int ii, converged;
    double qq = 3.2;
    int shift_size = 100;
    int extend_max = 10;
    int extend_count = 0;
    //    int burnIn = 0;
    //    int thin = 1;
    
    sbmInit(YY,nn,dd,eta,betaPrior,BB,PP,PPint,
	    yyComplete,logLik,flatTable,start,multiImpute);
    
    //  MCMC Loop
    converged = 0;
    while((converged != 1) & (extend_count <= extend_max)){
      for(ii = start ; ii < total ; ii++){
	if(multiImpute == 1){
	  sbmImputeMissingValues(nn,dd,YY,yyComplete,BB,PP,PPint);
	}
	sbmStep(yyComplete, nn, dd, eta, betaPrior, BB, PP,PPint);
	sbmRotate(nn,dd,BB,PP,PPint);
	
	if(ii >= burnIn && ((ii - burnIn) % thin == 0)){
	  logLik[(ii - burnIn)/thin] = sbmLogLikYY(nn,dd,YY,BB,PP,PPint);
	  updateFlatTable((ii - burnIn)/thin, nn, dd, BB, PP, flatTable);
	}//    printTableMMSBM(nn, dd, BB, PP, sendMat, recMat);
      }
      //Rprintf("total = %d\n",total);
      converged = convergenceCheck(logLik,(total - burnIn)/thin,qq);
      Rprintf("Converged Status %d, after %d extensions\n",converged,extend_count);
      if(converged != 1){
	extend_count = extend_count + 1;
	shiftFlatTable(shift_size,nn,dd,(total - burnIn)/thin,flatTable);
	start = total - (shift_size*thin);
      }
    }
    if(converged != 1){
      Rprintf("MCMC Failed to Converge\n");
    }
  }
  


  /**********  MMSBM  **********/

  void mmsbmInit(int *YY, int nn, int dd, double *alpha, double *betaPrior,
		 double *BB, double *PP, int *sendMat, int *recMat,
		 double *flatTable, int *yyComplete,int start, int multiImpute){
    int ii,jj;

    if(start > 0){
      mmsbmLoadTable(start, nn, dd, BB, PP, flatTable);
      mmsbmDrawZZ(nn, dd, YY, BB, PP, sendMat, recMat);
      for(ii = 0 ; ii < nn*nn ; ii++){
	yyComplete[ii] = YY[ii];
      }
      if(multiImpute == 1){
	mmsbmImputeMissingValues(nn,dd,YY,yyComplete,
				 BB,sendMat,recMat);
      }else{
	for(ii = 0 ; ii < nn*nn ; ii++){
	  yyComplete[ii] = YY[ii];
	}	
      }
    }else{
      if(multiImpute == 1){
	for(ii = 0 ; ii < nn*nn ; ii++){
	  if(yyComplete[ii] < 0){
	    yyComplete[ii] = 0;
	  }else{
	    yyComplete[ii] = YY[ii];
	  }
	}
      }else{
	for(ii = 0 ; ii < nn*nn ; ii++){
	  yyComplete[ii] = YY[ii];
	}
      }
      
      //  Initializing BB Matrix
      double holder[2];
      for(ii = 0 ; ii < dd*dd ; ii++){
	rdirichlet(2,betaPrior,holder);
	BB[ii] = holder[0];
      }
      
      //  Initializing PP Matrix
      double pDraws[dd];
      for(ii = 0 ; ii < nn ; ii++){
	rdirichlet(dd,alpha,pDraws);
	for(jj = 0 ; jj < dd ; jj++){
	  PP[ii*dd + jj] = pDraws[jj];
	}
      }
    }
  }

  void mmsbmImputeMissingValues(int nn, int dd, int *YY, int *yyComplete,
				double* BB, int *sendMat,int *recMat){
    int ii;
    for(ii = 0 ; ii < nn*nn ; ii++){
      if(ii/nn != ii % nn){
	if(YY[ii] < 0){
	  yyComplete[ii] = (int) rbinom(1.0,BB[sendMat[ii]*dd + recMat[ii]]);
	  //      yyComplete[ii] = gsl_ran_bernoulli(rng,BB[sendMat[ii]*dd + recMat[ii]]);
	}
      }
    }
  }

  
  void mmsbmLoadTable(int start, int nn, int dd, double *BB, double *PP, double *flatTable){
    
    int ii, offset;
    offset = (start - 1) * (dd * (nn + dd));
    
    for(ii = 0 ; ii < dd*dd ; ii++){
      BB[ii] = flatTable[offset + ii];
    }
    
    offset = offset + dd*dd;
    for(ii = 0 ; ii < nn*dd ; ii++){
      PP[ii] = flatTable[offset + ii];
    }
  }
  
  
  
  void mmsbmStep(int *YY, int nn, int dd, double *alpha, double *betaPrior,
		 double *BB, double *PP, int *sendMat, int *recMat){
    
    // Updating BB
    mmsbmDrawBB(nn, dd, YY, sendMat, recMat, BB, betaPrior);
    
    // Updating PP
    mmsbmDrawPP(nn, dd, sendMat, recMat, BB, alpha, PP);
    
    // Updating Sender and Receiver Matrices
    mmsbmDrawZZ(nn, dd, YY, BB, PP, sendMat, recMat);
    
  }
  
  
  void mmsbmMCMC(int total,int burnIn,int thin,int *YY,int nn,int dd,
		 double *alpha,double *betaPrior,
		 double *BB, double *PP, int *sendMat, int *recMat,
		 int *yyComplete,double *flatTable, int start,
		 int multiImpute){
    int ii;
    
    //  If start is 0, random initialization
    mmsbmInit(YY,nn,dd,alpha,betaPrior,BB,PP,
	      sendMat,recMat,flatTable,yyComplete,start,multiImpute);

    //  MCMC Loop
    for(ii = start ; ii < total ; ii++){
      if(multiImpute == 1){
	mmsbmImputeMissingValues(nn,dd,YY,yyComplete,
				 BB,sendMat,recMat);
      }
      mmsbmStep(yyComplete, nn, dd, alpha, betaPrior, BB, PP, sendMat,recMat);

      if(ii >= burnIn && ((ii - burnIn) % thin == 0)){
	updateFlatTable((ii - burnIn)/thin, nn, dd, BB, PP, flatTable);
      }
    }
  }
  
  
  
  /***********************************************
 *********  SBM Sampling Functions  *************
 ***********************************************/

  void sbmDrawBB(int nn, int dd, int *YY, double *PP, double *BB, double *betaPrior, int *PPint){
    
    int dd2 = dd*dd;
    double hit[dd2], miss[dd2];
    int ii;

    // Setting Up Hit and Miss Counts
    for(ii = 0 ; ii < dd2 ; ii++){
      hit[ii] = 0.0;
      miss[ii] = 0.0;
    }

    //  Populating Hit and Miss Counts
    int rr, cc, index;
    for(rr = 0 ; rr < nn ; rr++){
      for(cc = 0 ; cc < nn ; cc++){
	index = PPint[rr] * dd + PPint[cc];
	if(YY[rr * nn + cc] == 1){
	  hit[index] = hit[index] + 1;
	}else if(YY[rr * nn + cc] == 0){
	  miss[index] = miss[index] + 1;
	}
      }
    }

    //  Drawing From Posterior Distribution
    double alpha_new[2];
    double holder[2];
    for(ii = 0 ; ii < dd2 ; ii++){
      alpha_new[0] = betaPrior[0] + hit[ii];
      alpha_new[1] = betaPrior[1] + miss[ii];
      rdirichlet(2,alpha_new,holder);
      BB[ii] = holder[0];
    }
  }

  //  Computes full log-likelihood for sbm model
  double sbmLogLikYY(int nn, int dd, int *YY, double *BB, double *PP, int *PPint){
    
    int ii;
    
    // Reading in BB Matrix Values
    double pMat[nn*nn];
    int rr, cc, index;

    // Creating matrix of tie probabilities
    for(rr = 0 ; rr < nn ; rr++){
      for(cc = 0 ; cc < nn ; cc++){
	index = PPint[rr] * dd + PPint[cc];
	pMat[rr * nn + cc] = BB[index];
      }
    }

    // Summing over all tie probabilities, ignoring diagonal
    double total = 0.0;
    for(ii = 0 ; ii < nn*nn ; ii++){
      //      if((ii / nn) != (ii % nn)){
      if(YY[ii] == 1){
	total = total + log(pMat[ii]);
      }else if(YY[ii] == 0){
	total = total + log(1 - pMat[ii]);
      }
    }
    //    }
    return(total);
  }


  //  Computes log-likelihood affected by single node, ii
  double sbmLogLikYY_ii(int ii, int nn, int dd, int *YY, double *BB, double *PP, int *PPint){
    int jj;
    int PPii = PPint[ii];
    int index;
    double total = 0.0;

    //  sums the (i,j)th and (j,i)th probabilities for all j
    for(jj = 0 ; jj < nn ; jj++){
      if(jj != ii){

	//  (i,j)th term
	index = PPii * dd + PPint[jj];
	if(YY[ii*nn + jj] == 1){
	  total = total + log(BB[index]);
	}else if(YY[ii*nn + jj] == 1){
	  total = total + log(1 - BB[index]);
	}

	//  (j,i)th term
	index = PPint[jj] * dd + PPii;
	if(YY[jj*nn + ii] == 1){
	  total = total + log(BB[index]);
	}else if(YY[jj*nn + ii] == 0){
	  total = total + log(1 - BB[index]);
	}
      }
    }

    return(total);
  }


  void sbmDrawPP(int nn, int dd, int *YY, double *BB, double *PP, double *eta,
		 int *PPint){
    
    int ii,jj,kk;
    double pp[dd], total, pMax;
    int ind[dd];
    for(ii = 0 ; ii < nn ; ii++){
      total = 0.0;
      for(jj = 0 ; jj < dd ; jj++){
	PPint[ii] = jj;
	
	// Calculating Log Posterior Probability
	pp[jj] = sbmLogLikYY_ii(ii, nn, dd, YY, BB, PP,PPint) + log(eta[jj]);
	if(jj == 0){
	  pMax = pp[jj];
	}else{
	  if(pp[jj] > pMax){
	    pMax = pp[jj];
	  }
	}
      }
      
      // Exponentiating Probabilities
      for(jj = 0 ; jj < dd ; jj++){
	pp[jj] = exp(pp[jj] - pMax);

	total = total + pp[jj];
      }

      // Normalizing Probabilities
      for(jj = 0 ; jj < dd ; jj++){
	pp[jj] = pp[jj] / total;
      }
      
      // Drawing from Multinomial
      rmultinom(1,pp,dd,ind);
      
      // Setting PP[ii,] Value
      for(kk = 0 ; kk < dd ; kk++){
	if(ind[kk] == 1){
	  PP[ii*dd + kk] = 1.0;
	  PPint[ii] = kk;
	}else{
	  PP[ii*dd + kk] = 0.0;
	}
      }
      
    }
  }

  void sbmRotate(int nn, int dd, double *BB, double *PP, int *PPint){
    int ii, jj, kk;
    int assigned[dd], replaced[dd];
    for(ii = 0 ; ii < dd ; ii++){
      assigned[ii] = -1;
      replaced[ii] = -1;
    }
    int PPtmp[nn];
    for(ii = 0 ; ii < nn ; ii++){
      PPtmp[ii] = PPint[ii];
    }
    int replacements[nn];
    for(ii = 0 ; ii < nn ; ii++){
      replacements[ii] = -1;
    }
    
    for(ii = 0 ; ii < dd ; ii++){
      if(replacements[ii] == -1){
	assigned[PPtmp[ii]] = ii;
	replaced[ii] = PPtmp[ii];
	for(jj = 0 ; jj < nn ; jj++){
	  if(PPtmp[jj] == PPtmp[ii]){
	    replacements[jj] = ii;
	  }
	}
      }
    }

    int uncount = 0;
    for(ii = 0 ; ii < dd ; ii++){
      if(assigned[ii] == -1){
	uncount++;
      }
    }
    int vacancy, newlabel;
    while(uncount > 0){
      
      //  Setting Vacancy to Minimum Unassigned
      for(ii = 0 ; ii < dd ; ii++){
	if(assigned[ii] == -1){
	  vacancy = ii;
	  break;
	}
      }
      //  Setting newlabel to Minimum unreplaced
      for(ii = 0 ; ii < dd ;  ii++){
	if(replaced[ii] == -1){
	  newlabel = ii;
	  break;
	}
      }
      /*
      Rprintf("uncount: %d ",uncount);
      Rprintf("vacancy: %d ",vacancy);
      Rprintf("newlabel: %d ",newlabel);
      Rprintf("\n");
      RprintIntMat(1, dd, assigned);
      */	   

      assigned[vacancy] = newlabel;
      replaced[newlabel] = vacancy;
      for(jj = 0 ; jj < nn ; jj++){
	if(PPtmp[jj] == vacancy){
	  replacements[jj] = newlabel;
	}
      }
      uncount = 0;
      for(ii = 0 ; ii < dd ; ii++){
	if(assigned[ii] == -1){
	  uncount++;
	}
      }
      
    }

    //  Updating Membership Vectors
    for(ii = 0 ; ii < nn ; ii++){
      PPint[ii] = assigned[PPint[ii]];
      for(kk = 0 ; kk < dd ; kk++){
	if(PPint[ii] == kk){
	  PP[ii*dd + kk] = 1.0;
	}else{
	  PP[ii*dd + kk] = 0.0;
	}
      }
    }
    //  Updating Block Matrix
    for(ii = 0 ; ii < dd ; ii++){
      for(jj = 0 ; jj < dd ; jj++){
	BB[ii*dd + jj] = BB[assigned[ii]*dd + assigned[jj]];
      }
    }
    
  }



  
  /******************************************************
   *************  MMSBM Sampling Functions  *************
   *****************************************************/

  void rdirichlet(int k, double *alpha, double *x){
    int ii;
    double draw[k];
    double total = 0.0;
    
    for(ii=0 ; ii < k ; ii++){
      draw[ii] = rgamma(alpha[ii],1);
      //draw[ii] = gsl_ran_gamma(rng,alpha[ii],1);
      total = total + draw[ii];
    }
    
    for(ii = 0; ii < k; ii++){
      x[ii] = draw[ii] / total;
    }
  }
  
  
  void mmsbmDrawPP(int n, int d, int *sendMat, int *recMat, double *BB, double *alpha, double *PP){
    
    int ii,jj,sub1;
    double alpha_post[d];
    for(ii=0 ; ii < n ; ii++){
      //  Computing alpha for dirichlet posterior
      //  Initializing with Prior
      for(jj=0 ; jj < d ; jj++){
	alpha_post[jj] = alpha[jj];
      }
      
      //  Counting occurences of each group membership
      for(jj=0 ; jj < n ; jj++){
	if(jj != ii){
	  if(sendMat[ii * n + jj] >= 0){
	    sub1 = sendMat[ii * n + jj];
	    alpha_post[sub1] = alpha_post[sub1] + 1 ;
	  }
	  if(recMat[jj * n + ii] >= 0){
	    sub1 = recMat[jj * n + ii];
	    alpha_post[sub1] = alpha_post[sub1] + 1 ;
	  }
	}
      }
      
      //  Drawing from Posterior
      rdirichlet(d,alpha_post,&PP[d*ii]);
    }
  }
  
  void mmsbmDrawBB(int nn, int dd, int *YY, int *sendMat, int *recMat,
		   double *BB, double *betaPrior){
    int dd2 = dd*dd;
    double hit[dd2], miss[dd2];
    int ii;

    //  Initializing Hit and Miss Counts
    for(ii = 0 ; ii < dd2 ; ii++){
      hit[ii] = 0.0;
      miss[ii] = 0.0;
    }

    //  Counting Hits and Misses
    int sender, receiver, index;
    for(ii = 0 ; ii < nn*nn ; ii++){
      if(ii != (ii % nn)){
	sender = sendMat[ii];
	receiver = recMat[ii];
	index = sender * dd + receiver;
	if(YY[ii] == 1){
	  hit[index] = hit[index] + 1;
	}else if(YY[ii] == 0){
	  miss[index] = miss[index] + 1;
	}
      }
    }

    /*
      Rprintf("Hits:\n");
      RprintDoubleMat(dd,dd,hit);
      Rprintf("Misses:\n");
      RprintDoubleMat(dd,dd,miss);
      Rprintf("YY:\n");
      RprintIntMat(nn,nn,YY);
    */

    double alpha_new[2];
    double holder[2];
    //  Drawing from Posterior
    for(ii = 0 ; ii < dd2 ; ii++){
      alpha_new[0] = betaPrior[0] + hit[ii];
      alpha_new[1] = betaPrior[1] + miss[ii];
      rdirichlet(2,alpha_new,holder);
      BB[ii] = holder[0];
    }
  }
  

  void mmsbmDrawZZ(int nn, int dd,int *YY, double *BB, double *PP, int *sendMat, int *recMat){
    
    int ii,jj,rr;
    int dd2 = dd*dd;
    size_t kk = dd2;
    int ind[dd2];
    double pp[dd2];
    double total;
    
    // Setting up probability of no edge
    double negBB[dd2];
    for(ii = 0; ii < dd2 ; ii++){
      negBB[ii] = 1 - BB[ii];
    }

    //  For each (i,j) pair, where order matters
    for(ii = 0 ; ii < nn ; ii++){
      for(jj = 0 ; jj < nn ; jj++){

	//  Calculating Probability of Class Pair rr
	total = 0.0;
	
	//  If YY is 1, use prob of edge
	if(YY[ii*nn + jj] < 0){
	  sendMat[ii*nn + jj] = -1;
	  recMat[ii*nn+jj] = -1;
	}else{
	  if(YY[ii*nn + jj] == 1){
	    for(rr = 0 ; rr < dd2 ; rr++){
	      pp[rr] = PP[ii*dd + (rr % dd)] * PP[jj*dd + rr/dd] * BB[rr];
	      total = total + pp[rr];
	    }
	  }else if(YY[ii*nn + jj] == 0){ //  If YY is 0, use prob of no edge
	    for(rr = 0 ; rr < dd2 ; rr++){
	      pp[rr] = PP[ii*dd + (rr % dd)] * PP[jj*dd + rr/dd] * negBB[rr];
	      total = total + pp[rr];
	    } 
	  }
	  
	  //  Normal1izing Probabilities
	  for(rr = 0 ; rr < dd2 ; rr++){
	    pp[rr] = pp[rr] / total;
	  }
	  
	  rmultinom(1,pp,kk,ind);
	  //      gsl_ran_multinomial(rng,kk,1,pp,ind);
	  
	  //  Setting Sender and Receiver Classes
	  for(rr = 0 ; rr < dd2 ; rr++){
	    if(ind[rr] == 1){
	      sendMat[ii*nn + jj] = rr % dd;
	      recMat[ii*nn + jj]  = rr / dd;
	      break;
	    }
	  }
	}
      }
    }
  }
  

/********************************************
 **********  Output Functions  **************
 *******************************************/

/*  These functions allow easy printing of 
 *  matrix-like objects, given the number of
 *  rows and columns.
 */

  void RprintDoubleMat(int rows, int cols, double *mat){
    int ii, jj;
    for(ii = 0 ; ii < rows ; ii++){
      for(jj = 0 ; jj < cols ; jj++){
	Rprintf("%f ",mat[ii*cols + jj]);
      }
      Rprintf("\n");
    }
  }
  
  void RprintIntMat(int rows, int cols, int *mat){
    int ii, jj;
    for(ii = 0 ; ii < rows ; ii++){
      for(jj = 0 ; jj < cols ; jj++){
	Rprintf("%d ",mat[ii*cols + jj]);
      }
      Rprintf("\n");
    }
  }
  

  /*
  //  Prints all of the matrices involved in an MMSBM step
  void printMMSBM(int nn, int dd, double *BB, double *PP,
		  int *sendMat, int *recMat){
    
    printf("\nBlock Matrix\n");
    printDoubleMat(dd,dd,BB);
    
    printf("\nPP Matrix\n");
    printDoubleMat(nn,dd,PP);
    
    printf("\nSender Matrix:\n");
    printIntMat(nn,nn,sendMat);
    
    printf("Receiver Matrix:\n");
    printIntMat(nn,nn,recMat);
    
  }
  */
  /*
void printTableMMSBM(int nn, int dd, double *BB, double *PP,
		     int *sendMat, int *recMat){
  int ii;
  for(ii = 0 ; ii < dd*dd ; ii++){
    //    printf("%f\t",BB[ii]);
    myfile << BB[ii] << "\t";
  }
  for(ii = 0 ; ii < nn*dd ; ii++){
    //    printf("%f\t",PP[ii]);
    myfile << PP[ii] << "\t";
  }
  myfile << "\n";
}
  */
  
  //  Saves the current MMSBM state to the flatTable
  void updateFlatTable(int iter, int nn, int dd, double *BB, double *PP, 
		       double *flatTable){
    int ii, offset;
    offset = iter * (dd * (nn + dd));
    for(ii = 0 ; ii < dd*dd ; ii++){
      flatTable[offset + ii] = BB[ii];
    }
    offset = offset + dd*dd;
    for(ii = 0 ; ii < nn*dd ; ii++){
      flatTable[offset + ii] = PP[ii];
    }
    
  }

  void shiftFlatTable(int shift_size, int nn, int dd, int total,
		      double *flatTable){
    int ii;
    int flatTotal = total * (dd * (nn + dd));
    int flatShift = shift_size * (dd *(nn + dd));
    for(ii = 0 ; ii < flatTotal - flatShift ; ii++){
      flatTable[ii] = flatTable[ii + flatShift];
    }
    
  }



  void getMeanVar(double *vec, int lower, int upper,
		  double *Mean_t, double *Var_t, double * Len_t){

    
    double Var;
    double Mean = 0.0;
    double Len = upper - lower;
    double SumSq = 0.0;
    int ii;
    for(ii = lower ; ii < upper ; ii++){
      Mean = Mean + vec[ii];
      SumSq = SumSq + vec[ii] * vec[ii];
    }
    Mean = Mean / Len;
    Var = ((SumSq/Len) - (Mean * Mean))*(Len / (Len-1));
    *Mean_t = Mean;
    *Var_t = Var;
    *Len_t = Len;
    
  }
  
  int convergenceCheck(double *logLik, int total, double qq){
    
    double llMean1, llVar1, llLen1;
    double llMean2, llVar2, llLen2;
    double sdTot;
    
    
    //Rprintf("lower1 = %d, upper1 = %d\n",0,total/10);
    getMeanVar(logLik,0,total/10,
	       &llMean1, &llVar1, & llLen1);
    //Rprintf("lower2 = %d, upper2 = %d\n",total/2,total);    
    getMeanVar(logLik,total/2,total,
	       &llMean2, &llVar2, & llLen2);
    sdTot = sqrt((llVar1/llLen1) + (llVar2/llLen2));
    
    //Rprintf("Mean 1 = %f \t Mean 2 = %f \t sd.est = %f\n",
    //llMean1,llMean2,sdTot);
    

    if((llMean2 - llMean1) <= (qq * sdTot)){
      if((llMean1 - llMean2) <= (qq * sdTot)){
	return(1);
      }
    }
    return(0);
      
  }



  /////////////  UNUSED OLD FUNCTIONS  /////////////////

  /*
int *readCSV(const char *file, int *rows, int *cols){
  ifstream in;
  in.open(file);

  in >> *rows >> *cols;
  //  printf("Rows: %d  Columns: %d\n",*rows, *cols);
  int *table;
  int total = ( *rows) * ( *cols);
  table = new int [total];
  int ii;
  for(ii = 0 ; ii < total ; ii++){
    in >> table[ii];
  }
  in.close();  
  return(table);
}
  */


}
