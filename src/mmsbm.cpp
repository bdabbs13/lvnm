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
 ********  MCMC CONTROL FUNCTION  **********
 *******************************************/

  
void mmsbmMCMC(mmsbm_t myMMSBM, int start, int total, int burnIn, int thin,
	       int shift_size, int extend_max, double qq,
	       double *flatTable, double *logLik){
  
  int ii;
  int flatLength = myMMSBM.dd * (myMMSBM.nn + myMMSBM.dd);
  int flatTotal = (total - burnIn) / thin;

  int converged = 0, extend_count = 0;
  while((converged != 1) & (extend_count <= extend_max)){
    // Main MCMC Loop
    for(ii = start ; ii < total ; ii++){
      myMMSBM.step();
      
      if(ii >= burnIn && ((ii - burnIn) % thin == 0)){
	logLik[(ii - burnIn)/thin] = myMMSBM.LL();
	myMMSBM.updateFlatTable((ii - burnIn)/thin, flatTable);
      }
    }
    converged = convergenceCheck(logLik,(total - burnIn)/thin,qq);
    Rprintf("Converged Status %d, after %d extensions\n",converged,extend_count);
    
    if(converged != 1){
      extend_count = extend_count + 1;
      shiftFlatTable(shift_size,flatLength,flatTotal,flatTable);
      std::copy(logLik + shift_size, logLik + flatTotal,logLik);
      start = total - (shift_size*thin);
    }
  }
  if(converged != 1){
    Rprintf("MCMC Failed to Converge\n");
  }
}





mmsbm_t::mmsbm_t (int nodes, int blocks, int* adjMat, double* betaP, double* alphaP, int mImpute){
  
  int ii, jj;
  nn = nodes; dd = blocks;
  
  // Loading Adjacency Matrix
  YY = new int[nn*nn];
  yyComplete = new int[nn*nn];
  
  std::copy(adjMat,adjMat +(nn*nn),YY);
  multiImpute = mImpute == 1;
  
  if(multiImpute){
    for(ii = 0 ;ii < nn*nn ; ii++){
      if(yyComplete[ii] < 0){
	yyComplete[ii] = 0;
      }else{
	yyComplete[ii] = YY[ii];
      }
    }
  }else{
    std::copy(adjMat,adjMat + (nn*nn),yyComplete);
  }
  
  //  Initializing BB Matrix
  std::copy(betaP,betaP+2,betaPrior);
  BB = new double[dd*dd];
  double holder[2];
  for(ii = 0 ; ii < dd*dd ; ii++){
    rdirichlet(2,betaPrior,holder);
    BB[ii] = holder[0];
  }

  //  Initializing PP Matrix
  PP = new double[nn*dd];
  
  alpha = new double[dd];
  std::copy(alphaP,alphaP+dd,alpha);
  
  double pDraws[dd];
  for(ii = 0 ; ii < nn ; ii++){
    rdirichlet(dd,alpha,pDraws);
    for(jj = 0 ; jj < dd ; jj++){
      PP[ii*dd + jj] = pDraws[jj];
    }
  }
  
  //  Initializing Sender and Receiver Membership Mats
  sendMat = new int[nn*nn];
  recMat = new int[nn*nn];
    
}


void mmsbm_t::loadTable (int start, double *flatTable){
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

void mmsbm_t::step(){
  // Updating Sender and Receiver Matrices
  drawZZ();
  
  // Updated BB
  drawBB();
  
  // Updating PP
  drawPP();
  
  if(multiImpute){
    imputeMissingValues();
  }

}

void mmsbm_t::drawBB(){
  
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


void mmsbm_t::drawPP(){
  
  int ii,jj,sub1;
  double alpha_post[dd];
  for(ii=0 ; ii < nn ; ii++){
    //  Computing alpha for dirichlet posterior
    //  Initializing with Prior
    for(jj=0 ; jj < dd ; jj++){
      alpha_post[jj] = alpha[jj];
    }
      
    //  Counting occurences of each group membership
    for(jj=0 ; jj < nn ; jj++){
      if(jj != ii){
	if(sendMat[ii * nn + jj] >= 0){
	  sub1 = sendMat[ii * nn + jj];
	  alpha_post[sub1] = alpha_post[sub1] + 1 ;
	}
	if(recMat[jj * nn + ii] >= 0){
	  sub1 = recMat[jj * nn + ii];
	  alpha_post[sub1] = alpha_post[sub1] + 1 ;
	}
      }
    }
      
    //  Drawing from Posterior
    rdirichlet(dd,alpha_post,&PP[dd*ii]);
  }
}
  
void mmsbm_t::drawZZ(){
  
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


void mmsbm_t::rotate(){

}
void mmsbm_t::imputeMissingValues(){
  int ii;
  for(ii = 0 ; ii < nn*nn ; ii++){
    if(ii/nn != ii % nn){
      if(YY[ii] < 0){
	yyComplete[ii] = (int) rbinom(1.0,BB[sendMat[ii]*dd + recMat[ii]]);
      }
    }
  }
}

double mmsbm_t::LL(){

  int ii,jj,rr,cc;
  // Reading in BB Matrix Values
  double pMat[nn*nn];
  
  // Creating matrix of tie probabilities
  for(rr = 0 ; rr < nn ; rr++){
    for(cc = 0 ; cc < nn ; cc++){
      pMat[rr*nn + cc] = 0.0;
      for(ii = 0 ; ii < dd ; ii++){
	for(jj = 0 ; jj < dd ; jj++){
	  pMat[rr*nn + cc] += PP[rr*dd + jj] * PP[cc*dd + ii] * BB[ii*dd + jj];
	}
      }
    }
  }
  
  // Summing over all tie probabilities, ignoring diagonal
  double total = 0.0;
  for(ii = 0 ; ii < nn*nn ; ii++){
    if(YY[ii] == 1){
      total = total + log(pMat[ii]);
    }else if(YY[ii] == 0){
      total = total + log(1 - pMat[ii]);
    }
  }
  
  return(total);
}

void mmsbm_t::updateFlatTable(int iter, double* flatTable){
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

void mmsbm_t::print (bool printNetwork){
  Rprintf("dd = %d, nn = %d\n",dd,nn);
  Rprintf("bbPrior = (%.2f, %.2f)\n",betaPrior[0],betaPrior[1]);
  Rprintf("dirichletPrior = ");
  RprintDoubleMat(1,dd,alpha);
  Rprintf("BB = \n"); RprintDoubleMat(dd, dd, BB);
  Rprintf("PI = \n"); RprintDoubleMat(nn,dd,PP);
  if(printNetwork){
    Rprintf("YY = \n"); RprintIntMat(nn,nn,YY);
  }
  

}



