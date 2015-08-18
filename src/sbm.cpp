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
#include "sbm.h"
#include "helper.h"
//using namespace std;
using std::ifstream;
//using std::ofstream;

//gsl_rng *rng;
//ofstream myfile;


/****************************************************************
 *********  FUNCTION DEFINITIONS FOR SBM OBJECT CLASS  ***********
 ****************************************************************/

void sbmMCMC(sbm_t mySBM, int start, int total, int burnIn, int thin,
	     int shift_size, int extend_max, double qq,
	     double *flatTable, double *logLik){
  
  int ii;
  int flatLength = mySBM.dd * (mySBM.nn + mySBM.dd);
  int flatTotal = (total - burnIn) / thin;
  
  //  MCMC Loop
  int converged = 0, extend_count = 0;
  while((converged != 1) & (extend_count <= extend_max)){
    for(ii = start ; ii < total ; ii++){
      mySBM.step();
      
      if(ii >= burnIn && ((ii - burnIn) % thin == 0)){
	logLik[(ii - burnIn)/thin] = mySBM.LL();
	mySBM.updateFlatTable((ii - burnIn)/thin,flatTable);
      }
    }
    
    //  Checking for convergence of log-likelihood chain
    converged = convergenceCheck(logLik,(total - burnIn)/thin,qq);
    Rprintf("Converged Status %d, after %d extensions\n",converged,extend_count);
    
    if(converged != 1){
      extend_count = extend_count + 1;
      shiftFlatTable(shift_size,flatLength,flatTotal,flatTable);
      start = total - (shift_size*thin);
    }
  }
  if(converged != 1){
    Rprintf("MCMC Failed to Converge\n");
  }
}
	     



//  Constructor Function for SBM object
sbm_t::sbm_t (int nodes, int blocks, int *adjMat, 
	      double *betaP, double * etaP, int mImpute){
  int ii, kk;
  // Loading Constants
  nn = nodes; dd = blocks;
  
  // Loading Adjacency Matrix
  YY = new int[nn*nn];
  yyComplete = new int[nn*nn];
  
  std::copy(adjMat,adjMat+(nn*nn),YY);
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
  
  PP = new double[nn*dd];
  PPint = new int[nn];
  
  // Initializing PP Matrix
  eta = new double[dd];
  std::copy(etaP,etaP+dd,eta);
  
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


void sbm_t::loadTable(int start, double *flatTable){
  //  int ii, jj;
  //  sbmLoadTable(start,nn,dd,BB,PP,flatTable);

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

  imputeMissingValues();
}


void sbm_t::step(){
  drawBB();
  drawPP();
  rotate();
  if(multiImpute){
    imputeMissingValues();
  }
}

void sbm_t::drawBB(){
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

void sbm_t::drawPP(){
  int ii,jj,kk;
  double pp[dd], total, pMax;
  int ind[dd];
  for(ii = 0 ; ii < nn ; ii++){
    total = 0.0;
    for(jj = 0 ; jj < dd ; jj++){
      PPint[ii] = jj;
	
      // Calculating Log Posterior Probability
      pp[jj] = nodeLL(ii);
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

void sbm_t::rotate(){
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

void sbm_t::imputeMissingValues(){
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

double sbm_t::LL(){
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

double sbm_t::nodeLL(int ii){
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

void sbm_t::updateFlatTable(int iter, double *flatTable){
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


void sbm_t::print(bool printNetwork){
  Rprintf("dd = %d, nn = %d\n",dd,nn);
  Rprintf("bbPrior = (%.2f, %.2f)\n",betaPrior[0],betaPrior[1]);
  Rprintf("multiPrior = ");
  RprintDoubleMat(1,dd,eta);
  Rprintf("BB = \n"); RprintDoubleMat(dd, dd, BB);
  Rprintf("PI = \n"); RprintIntMat(1,nn,PPint);
  if(printNetwork){
    Rprintf("YY = \n"); RprintIntMat(nn,nn,YY);
  }
  
}






