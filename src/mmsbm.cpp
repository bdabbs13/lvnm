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
#include "helper.h"
//using namespace std;
using std::ifstream;
//using std::ofstream;

//gsl_rng *rng;
//ofstream myfile;



/********************************************
 ********  MCMC CONTROL FUNCTIONS  **********
 *******************************************/


/***********  MMSBM **********/

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
  
  

/******************************************************
 *************  MMSBM Sampling Functions  *************
 *****************************************************/

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
  


void updateFlatTable(int iter, int nn, int dd, double *BB, double *PP, double *flatTable){
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

