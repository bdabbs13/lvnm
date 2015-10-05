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

void sbmMCMC(sbm_t *mySBM, int start, int total, int burnIn, int thin,
	     int shift_size, int extend_max, double qq,
	     double *flatTable, double *logLik, int verbose){
  
  int ii;
  int flatLength = mySBM->dd * (mySBM->nn + mySBM->dd);
  int flatTotal = (total - burnIn) / thin;
  
  //  MCMC Loop
  int converged = 0, extend_count = 0;
  while((converged != 1) & (extend_count <= extend_max)){

    // Main MCMC Loop
    for(ii = start ; ii < total ; ii++){
      mySBM->step();
      
      if(ii >= burnIn && ((ii - burnIn) % thin == 0)){
	logLik[(ii - burnIn)/thin] = mySBM->LL();
	mySBM->updateFlatTable((ii - burnIn)/thin,flatTable);
      }
    }
    
    //  Checking for convergence of log-likelihood chain
    converged = convergenceCheck(logLik,(total - burnIn)/thin,qq);
    if(verbose > 0){
      Rprintf("Converged Status %d, after %d extensions\n",converged,extend_count);
    }
    
    if(converged != 1){
      extend_count = extend_count + 1;
      if(extend_count <= extend_max){
	shiftFlatTable(shift_size,flatLength,flatTotal,flatTable);
	std::copy(logLik + shift_size, logLik + flatTotal,logLik);
	start = total - (shift_size*thin);
      }
    }
  }
  if(converged != 1){
    if(verbose > -1){
      Rprintf("Warning: MCMC Failed to Converge\n");
    }
  }else{
    if(verbose > 0){
      Rprintf("MCMC Converged\n");
    }
  }
  
}


void sbmEM(sbm_t *mySBM, int iter_max, double threshold,
	   double *flatTable, double *logLik, double *eta, int verbose){
  int dd = mySBM->dd;
  int iter = 0;
  double *BB_old = new double[dd*dd];
  double delta = threshold + 1.0;
  while((iter < iter_max) & (delta >= threshold)){
    iter = iter + 1;
    mySBM->getBB(BB_old);
    mySBM->logBB();
    //    mySBM->print(1 == 0);

    mySBM->iterEM();
    if(verbose > 2){
      mySBM->print(false);
    }
    delta = mySBM->BBdiff(BB_old);
    
    if(verbose > 1){
      Rprintf("iter - %d delta - %.6f\n",iter,delta);
    }
  }
  if(verbose == 1){
    Rprintf("iter - %d delta - %.6f\n",iter,delta);
  }
  //mySBM->print(1 == 0);
  mySBM->updateFlatTable(0,flatTable);
  mySBM->geteta(eta);
  logLik[0] = mySBM->LL();
  delete[] BB_old;
}


void sbm_t::iterEM (){
  int ss, rr, ii, jj;
  int rowd, rown,bbind;
  double total, pipi;
  // E Step
  getMultinomPosterior();  //  Update Posterior Mean HH
  std::copy(HH,HH+(nn*dd),PP);  //  Copy into PP
  
  // M Step
  colSums(PP,nn,dd,eta);
  for(ii = 0 ; ii < dd ; ii++){
    eta[ii] = eta[ii] / nn;
  }
  
  // Updating BB with the MLE conditional on PP and YY
  for(ss = 0 ; ss < dd ; ss++){
    for(rr = 0 ; rr < dd ; rr++){
      total = 0.0;
      bbind = ss*dd + rr;
      BB[bbind] = 0.0;
      for(ii = 0 ; ii < nn ; ii++){
	rowd = ii*dd;
	rown = ii*nn;
	for(jj = 0 ; jj < nn ; jj++){
	  if((ii != jj) & (YY[ii*nn + jj] >= 0)){
	    pipi = PP[rowd + ss] * PP[jj*dd + rr];
	    BB[bbind] = BB[bbind] + pipi * YY[rown + jj];
	    total = total + pipi;
	  }
	}
      }
      if(total < MIN_TOTAL){
	BB[bbind] = MIN_TOTAL;
      }else{
	BB[bbind] = BB[bbind] / total;
      }
      BB_inv[bbind] = 1.0 - BB[bbind];
    }
  }
  is_BB_logged = false;
}



//  This function updates PP to be the posterior mean given BB and YY
void sbm_t::getMultinomPosterior(){
  
  int ii,jj;
  double total, pMax;

  for(ii = 0 ; ii < nn ; ii++){
    for(jj = 0 ; jj < dd ; jj++){
      PPint[ii] = jj;

      // Calculating Log Posterior Probability
      HH[ii*dd + jj] = nodeLL_long(ii) + log(eta[jj]);
      if(jj == 0){
	pMax = HH[ii*dd + jj];
      }else{
	if(HH[ii*dd + jj] > pMax){
	  pMax = HH[ii*dd + jj];
	}
      }
    }
    
    // Exponentiating Probabilities
    total = 0.0;
    for(jj = 0 ; jj < dd ; jj++){
      HH[ii*dd + jj] = exp(HH[ii*dd + jj] - pMax);
      total = total + HH[jj];
    }
    
    // Normalizing Probabilities
    for(jj = 0 ; jj < dd ; jj++){
     HH[ii*dd + jj] = HH[ii*dd + jj] / total;
    }
  }
}



void sbm_t::getBB(double *BB_out){
  std::copy(BB,BB + dd*dd,BB_out);
}

void sbm_t::geteta(double *eta_out){
  std::copy(eta, eta + dd, eta_out);
}

double sbm_t::BBdiff(double *BB_old){
  int ii;
  double delta = 0.0;
  for(ii = 0 ; ii < dd*dd; ii++){
    delta = delta + ((BB[ii] - BB_old[ii]) * (BB[ii] - BB_old[ii]));
  }
  return(delta);
}

void sbm_t::logBB(){
  int ii;
  if(!is_BB_logged){
    for(ii = 0 ; ii < dd*dd ; ii++){
      BB[ii] = log(BB[ii]);
      BB_inv[ii] = log(BB_inv[ii]);
    }
    is_BB_logged = true;
  }
}

void sbm_t::expBB(){
  int ii;
  if(is_BB_logged){
    for(ii = 0 ; ii < dd*dd ; ii++){
      BB[ii] = exp(BB[ii]);
      BB_inv[ii] = exp(BB_inv[ii]);
    }
    is_BB_logged = false;
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
  BB_inv = new double[dd*dd];
  double holder[2];
  for(ii = 0 ; ii < dd*dd ; ii++){
    rdirichlet(2,betaPrior,holder);
    BB[ii] = holder[0];
    BB_inv[ii] = 1.0 - BB[ii];
  }
  is_BB_logged = false;

  PP = new double[nn*dd];
  PPint = new int[nn];
  HH = new double[nn*dd]();
  
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

sbm_t::~sbm_t (){
  delete[] YY;
  delete[] yyComplete;
  delete[] BB;
  delete[] BB_inv;
  delete[] PP;
  delete[] HH;
  delete[] PPint;
  delete[] eta;
  //  Rprintf("Freeing Allocated Memory for SBM Object\n");
}



void sbm_t::loadTable(int start, double *flatTable){
  //  int ii, jj;
  //  sbmLoadTable(start,nn,dd,BB,PP,flatTable);

  int ii, jj, offset;
  double pMax;
  //  skip first (start - 1) sets of draws
  offset = (start - 1) * (dd * (nn + dd));
  
  //  Loading BB
  for(ii = 0 ; ii < dd*dd ; ii++){
    BB[ii] = flatTable[offset + ii];
    BB_inv[ii] = 1.0 - BB[ii];
  }
  
  //  Loading PP
  offset = offset + dd*dd;
  for(ii = 0 ; ii < nn*dd ; ii++){
    PP[ii] = flatTable[offset + ii];
  }
  //  Loading PPint
  for(ii = 0 ; ii < nn ; ii++){
    for(jj = 0 ; jj < dd; jj++){
      if(jj == 0){
	PPint[ii] = jj;
	pMax = PP[ii*dd + jj];
      }else if(PP[ii*dd + jj] > pMax){
	pMax = PP[ii*dd + jj];
	PPint[ii] = jj;
      }
    }
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
    BB_inv[ii] = 1.0 - BB[ii];
  }
  is_BB_logged = false;
}

void sbm_t::drawPP(){
  int ii,jj,kk;
  double pp[dd], total, pMax;
  int ind[dd];

  logBB();
  for(ii = 0 ; ii < nn ; ii++){
    total = 0.0;
    for(jj = 0 ; jj < dd ; jj++){
      PPint[ii] = jj;
	
      // Calculating Log Posterior Probability
      pp[jj] = nodeLL(ii) + log(eta[jj]);
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
  expBB();
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
      BB_inv[ii*dd + jj] = 1.0 - BB[ii*dd + jj];
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
  int rr, cc, index;
  logBB();
  
  double total = 0.0;
  for(rr = 0 ; rr < nn ; rr++){
    for(cc = 0 ; cc < nn ; cc++){
      index = PPint[rr] * dd + PPint[cc];
      if(YY[rr*nn + cc] == 1){
	total = total + BB[index];
      }else if(YY[rr * nn + cc] == 0){
	total = total + BB_inv[index];
      }
    }
  }
  expBB();
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
	total = total + (BB[index]);
      }else if(YY[ii*nn + jj] == 1){
	total = total + (BB_inv[index]);
      }
      
      //  (j,i)th term
      index = PPint[jj] * dd + PPii;
      if(YY[jj*nn + ii] == 1){
	total = total + (BB[index]);
      }else if(YY[jj*nn + ii] == 0){
	total = total + (BB_inv[index]);
      }
    }
  }
  
  return(total);
}


double sbm_t::nodeLL_long(int ii){

  int jj,kk;
  int PPii = PPint[ii];
  int index;
  double total = 0.0;
  
  //  sums the (i,j)th and (j,i)th probabilities for all j
  for(jj = 0 ; jj < nn ; jj++){
    if(jj != ii){
      
      //  (i,j)th term
      if(YY[ii * nn + jj] == 1){
	for(kk = 0 ; kk < dd ; kk++){
	  index = PPii * dd + kk;
	  total = total + PP[jj*dd + kk] * (BB[index]);
	}
      }else if(YY[ii * nn + jj] == 0){
	for(kk = 0 ; kk < dd ; kk++){
	  index = PPii * dd + kk;
	  total = total + PP[jj*dd + kk] * (BB_inv[index]);
	}
      }
	
      //  (j,i)th term
      if(YY[jj * nn + ii] == 1){
	for(kk = 0 ; kk < dd ; kk++){
	  index = kk * dd + PPii;
	  total = total + PP[jj*dd + kk] * (BB[index]);
	}
      }else if(YY[jj * nn + ii] == 0){
	for(kk = 0 ; kk < dd ; kk++){
	  index = kk * dd + PPii;
	  total = total + PP[jj*dd + kk] * (BB_inv[index]);
	}
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
  Rprintf("PI = \n"); RprintDoubleMat(nn,dd,PP);
  if(printNetwork){
    Rprintf("YY = \n"); RprintIntMat(nn,nn,YY);
  }
  
}






