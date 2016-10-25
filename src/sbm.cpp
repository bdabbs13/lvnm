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


/***********  SBM CLASS CONSTRUCTOR AND DESTRUCTOR  **********/

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
  BB_old = new double[dd*dd];
  double holder[2];
  for(ii = 0 ; ii < dd*dd ; ii++){
    rdirichlet(2,betaPrior,holder);
    BB[ii] = holder[0];
    BB_inv[ii] = 1.0 - BB[ii];
  }
  
  is_BB_logged = false;
  
  hit = new double[dd*dd];
  miss = new double[dd*dd];
  
  
  PP = new double[nn*dd];
  PPint = new int[nn];
  HH = new double[nn*dd]();
  
  // Initializing PP Matrix
  eta = new double[dd];
  std::copy(etaP,etaP+dd,eta);
  
  int ind[dd];
  for(ii = 0 ; ii < nn ; ii++){
    double *etaInit = new double[dd];
    double etaTotal = 0.0;
    for(kk = 0 ; kk < dd ; kk++){
      etaTotal += eta[kk];
    }
    for(kk = 0 ; kk < dd ; kk++){
      etaInit[kk] = eta[kk] / etaTotal;
    }
    //    RprintDoubleMat(1,dd,etaInit);
    rmultinom(1,etaInit,dd,ind);
    //    RprintIntMat(1,dd,ind);
      
    delete[] etaInit;
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
  delete[] BB_old;
  delete[] PP;
  delete[] HH;
  delete[] PPint;
  delete[] eta;
  //  Rprintf("Freeing Allocated Memory for SBM Object\n");
}


/***********  INPUT FUNCTION  ************/
void sbm_t::loadTable(int start, double *flatTable){
  //  int ii, jj;
  //  sbmLoadTable(start,nn,dd,BB,PP,flatTable);
  
  int ii, jj, offset;
  double pMax;
  //  skip first (start - 1) sets of draws
  offset = (start - 1) * (dd * (nn + dd));
  
  //  Loading BB
  for(ii = 0 ; ii < dd*dd ; ii++){
    BB[ii] = logCheck(flatTable[offset + ii]);
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

void sbm_t::loadSBM(double *BBout, int *MMBout){
  loadBB(BBout);
  loadPPint(MMBout);
}

void sbm_t::loadBB(double *BBout){
  int ii;
  std::copy(BBout,BBout + (dd*dd), BB);
  for(ii = 0; ii < dd*dd ; ii++){
    BB_inv[ii] = 1.0 - BB[ii];
  }
}

void sbm_t::loadHH(double *HHout){
  std::copy(HHout,HHout+nn*dd,PP);
}

void sbm_t::loadPPint(int *MMBout){
  int ii;
  for(ii = 0 ; ii < nn ; ii++){
    PPint[ii] = MMBout[ii] - 1;
  }
}

void sbm_t::initPPem(double minVal){
  int ii, jj;
  for(ii = 0 ; ii < nn ; ii++){
    for(jj = 0 ; jj < dd ; jj++){
      if(PPint[ii] == jj){
	PP[ii*dd + jj] = 1.0 - minVal;
      }else{
	PP[ii*dd + jj] = minVal / (dd-1);
      }
    }
  }
}


/******  OUTPUT FUNCTIONS  *******/
void sbm_t::updateFlatTable(int iter, double *flatTable){
  expBB();
  int ii, offset;
  
  // Saving BB
  offset = iter * (dd * (nn + dd));
  for(ii = 0 ; ii < dd*dd ; ii++){
    flatTable[offset + ii] = BB[ii];
  }
  
  // Saving PP
  offset = offset + dd*dd;
  for(ii = 0 ; ii < nn*dd ; ii++){
    flatTable[offset + ii] = PP[ii];
  }
}

void sbm_t::updateBB(int iter, double *BBout){
  expBB();
  std::copy(BB,BB+(dd*dd),BBout + (iter * dd*dd));
}

void sbm_t::updateMMB(int iter, int *MMBout){
  std::copy(PPint,PPint + nn, MMBout + (iter * nn));
}

void sbm_t::updateEta(double *eta_out){
  std::copy(eta, eta + dd, eta_out);
}

void sbm_t::updatePosteriorMat(int iter, double *postMat){
  std::copy(HH,HH +(nn*dd),postMat + (iter*dd*nn));
}





/************************  MCMC FUNCTIONS  **********************/


void sbmMCMC(sbm_t *mySBM, int start, int total, int burnIn, int thin,
	     int shift_size, int extend_max, double qq,  //double *flatTable,
	     double *BBout, int * MMBout,  double *logLik, double *postMat,
	     int verbose){
  
  int ii;
  //int flatLength = mySBM->dd * (mySBM->nn + mySBM->dd);
  int flatTotal = (total - burnIn) / thin;
  
  //  MCMC Loop
  int converged = 0, extend_count = 0;
  while((converged != 1) & (extend_count <= extend_max)){
    
    // Main MCMC Loop
    for(ii = start ; ii < total ; ii++){
      mySBM->step();
      
      if(ii >= burnIn && ((ii - burnIn) % thin == 0)){
	logLik[(ii - burnIn)/thin] = mySBM->LL();
	//mySBM->updateFlatTable((ii - burnIn)/thin,flatTable);
	mySBM->updateBB((ii - burnIn)/thin,BBout);
	mySBM->updateMMB((ii - burnIn)/thin,MMBout);
	mySBM->updatePosteriorMat((ii - burnIn)/thin,postMat);
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
	//shiftFlatTable(shift_size,flatLength,flatTotal,flatTable);
	shiftFlatTable(shift_size,mySBM->dd * mySBM->dd,flatTotal,BBout);
	shiftFlatTable(shift_size,mySBM->nn,flatTotal,MMBout);
	shiftFlatTable(shift_size,mySBM->dd * mySBM->nn,flatTotal,postMat);
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



void sbm_t::step(){
  drawPP(); // logs then exps BB
  drawBB(); // updates BB, but doesn't use them, then sets BB_logged to false
  rotate();
  if(multiImpute){
    imputeMissingValues(); // assumes exp BB
  }
}

void sbm_t::drawBB(){
  int dd2 = dd*dd;
  //  double hit[dd2], miss[dd2];
  int ii;
  

  // Setting Up Hit and Miss Counts
  computeHitMiss();

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
  double pMax;//, total;
  int ind[dd];

  //  logBB();
  for(ii = 0 ; ii < nn ; ii++){
    //total = 0.0;
    for(jj = 0 ; jj < dd ; jj++){
      PPint[ii] = jj;
	
      // Calculating Log Posterior Probability
      HH[ii*dd + jj] = nodeLL(ii) + log(eta[jj]);
      //pp[jj] = nodeLL(ii) + log(eta[jj]);
      if(jj == 0){
	pMax = HH[ii*dd + jj];
      }else{
	if(HH[ii*dd + jj] > pMax){
	  pMax = HH[ii*dd + jj];
	}
      }
    }
      
    // Exponentiating Probabilities
    for(jj = 0 ; jj < dd ; jj++){
      HH[ii*dd + jj] = exp(HH[ii*dd + jj] - pMax);
      //      total = total + HH[ii*dd + jj];
    }

    //    int hhOne = -1;
    // Normalizing Probabilities
    /*
    for(jj = 0 ; jj < dd ; jj++){
      HH[ii*dd + jj] = HH[ii*dd + jj] / total;
    }
    */
    normalizeVec(HH+ii*dd,dd);
    logZeroFix(HH + ii*dd,dd);
    
    // Drawing from Multinomial
    rmultinom(1,HH + ii*dd,dd,ind);
      
    // Setting PP[ii,] Value
    for(kk = 0 ; kk < dd ; kk++){
      if(ind[kk] == 1){
	//	PP[ii*dd + kk] = 1.0;
	PPint[ii] = kk;
	break;
      }
      //else{
      //PP[ii*dd + kk] = 0.0;
      //}
    }
  }
  //  expBB();
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
  expBB();
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





/************************  EM FUNCTIONS  *********************/


void sbmEM(sbm_t *mySBM, int iter_max, double threshold,
	   double *HHout, double *BBout, int *MMBout,
	   double *logLik, double *eta, int verbose){
  int iter = 0;
  double delta = threshold + 1.0;
  while((iter < iter_max) & (delta >= threshold)){
    iter = iter + 1;
    mySBM->saveBB_old();

    if(verbose > 2){
      mySBM->print(false);
    }
    mySBM->iterEM();
    delta = mySBM->BBdiff();
    
    if(verbose > 1){
      Rprintf("iter - %d delta - %.6f\n",iter,delta);
    }
  }
  if(verbose == 1){
    Rprintf("iter - %d delta - %.6f\n",iter,delta);
  }
  //mySBM->print(1 == 0);
  mySBM->updatePosteriorMat(0,HHout);
  mySBM->updateBB(0,BBout);
  mySBM->updateMMB(0,MMBout);
  mySBM->updateEta(eta);
  logLik[0] = mySBM->LL();
  //  delete[] BB_log;
}


void sbm_t::iterEM (){
  int ss, rr, ii, jj;
  int rowd, rown,bbind;
  double total, pipi;


  //  E Step
  getMultinomPosterior();  //  Update Posterior Mean HH 
  //Rprintf("\nHH:\n");
  //RprintDoubleMat(nn,dd,HH);
  std::copy(HH,HH+(nn*dd),PP);  //  Copy into PP
  

  // M Step
  colSums(HH,nn,dd,eta);
  double etaTotal = 0.0;
  for(ii = 0 ; ii < dd ; ii++){
    eta[ii] = eta[ii] / nn;
    if(eta[ii] < MIN_LOG){
      eta[ii] = MIN_LOG;
    }
    etaTotal = etaTotal + eta[ii];
  }
  for(ii = 0 ; ii < dd ; ii++){
    eta[ii] = eta[ii] / etaTotal;
  }
  
  // Updating BB with the MLE conditional on HH and YY
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
	    pipi = HH[rowd + ss] * HH[jj*dd + rr];
	    BB[bbind] = BB[bbind] + pipi * YY[rown + jj];
	    total = total + pipi;
	  }
	}
      }
      if(total < MIN_LOG){
	BB[bbind] = MIN_LOG;
      }else{
	BB[bbind] = BB[bbind] / total;
	BB[bbind] = logCheck(BB[bbind]);
	/*
	if(BB[bbind] > MAX_LOG){
	  BB[bbind] = MAX_LOG;
	}else if(BB[bbind] < MIN_LOG){
	  BB[bbind] = MIN_LOG;
	}
	*/
      }
      BB_inv[bbind] = 1.0 - BB[bbind];
    }
  }
  is_BB_logged = false;
}



//  This function updates PP to be the posterior mean given BB and YY
void sbm_t::getMultinomPosterior(){
  
  int ii,jj;
  double pMax;//, total;

  for(ii = 0 ; ii < nn ; ii++){
    //    Rprintf("\nHH[%d,] = ",ii);
    for(jj = 0 ; jj < dd ; jj++){
      PPint[ii] = jj;

      // Calculating Log Posterior Probability
      HH[ii*dd + jj] = nodeLL_long(ii) + log(eta[jj]);
      //      Rprintf("%f, ",HH[ii*dd + jj]);
      if(jj == 0){
	pMax = HH[ii*dd + jj];
      }else{
	if(HH[ii*dd + jj] > pMax){
	  pMax = HH[ii*dd + jj];
	}
      }
    }
    
    // Exponentiating Probabilities
    //    total = 0.0;
    for(jj = 0 ; jj < dd ; jj++){
      HH[ii*dd + jj] = exp(HH[ii*dd + jj] - pMax);
      //      total = total + HH[ii*dd + jj];
    }
    
    // Normalizing Probabilities
    /*
    for(jj = 0 ; jj < dd ; jj++){
     HH[ii*dd + jj] = HH[ii*dd + jj] / total;
    }
    */
    normalizeVec(HH + ii*dd,dd);
  }
}


/*****************  HELPER FUNCTIONS  *****************/
void sbm_t::computeHitMiss(){
  int dd2 = dd*dd; 
  int rr, cc, index;

  //  Setting values to zero
  std::fill(hit, hit+dd2, 0.0);
  std::fill(miss,miss+dd2,0.0);

  //  Counting hits and misses
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
}

void sbm_t::computeBBmle(){
  computeHitMiss();
  
  int ii,dd2 = dd*dd;
  for(ii = 0 ; ii < dd2 ; ii++){
    if(hit[ii] == 0.0){
      BB[ii] = 0.0;
    }else{
      BB[ii] = hit[ii] / (hit[ii] + miss[ii]);
    }
    BB[ii] = logCheck(BB[ii]);
}
}


/******************  Log-Likelihood Functions  *****************/
double sbm_t::LL(){
  int rr, cc;//, index;
  logBB();
  
  double total = 0.0;
  for(rr = 0 ; rr < nn ; rr++){
    for(cc = 0 ; cc < nn ; cc++){
      total = total + tieLL(YY[rr*nn+cc],PPint[rr],PPint[cc]);
    }
  }
  expBB();
  return(total);
}


double sbm_t::nodeLL(int ii){
  logBB();
  int jj;
  int PPii = PPint[ii];
  //  int index;
  double total = 0.0;
    
  //  sums the (i,j)th and (j,i)th probabilities for all j
  for(jj = 0 ; jj < nn ; jj++){
    if(jj != ii){
      
      //  (i,j)th term
      total = total + tieLL(YY[ii*nn + jj],PPii,PPint[jj]);
      
      //  (j,i)th term
      total = total + tieLL(YY[jj*nn + ii],PPint[jj],PPii);
    }
  }
  
  return(total);
}


double sbm_t::tieLL(int yy, int sendBlock, int recBlock){

  if(yy == 1){
    return(BB[sendBlock * dd + recBlock]);
  }else if(yy == 0){
    return(BB_inv[sendBlock * dd + recBlock]);
  }else{
    return(0.0);
  }
}



double sbm_t::nodeLL_long(int ii){
  logBB();
  int jj;
  int PPii = PPint[ii];
  double total = 0.0;
  
  //  sums the (i,j)th and (j,i)th probabilities for all j
  for(jj = 0 ; jj < nn ; jj++){
    if(jj != ii){
      //  (i,j)th term
      total = total + tieLL_sender(YY[ii*nn + jj],PPii,jj);
      //  (j,i)th term
      total = total + tieLL_receiver(YY[jj*nn + ii],PPii,jj);
    }
  }
  
  return(total);
}


double sbm_t::tieLL_sender(int yy, int sendBlock, int receiver){
  int kk;
  int baseBlock = sendBlock*dd;
  int baseNode = receiver*dd;
  double total = 0.0;
  if(yy == 1){
    for(kk = 0 ; kk < dd ; kk++){
      total = total + PP[baseNode + kk] * BB[baseBlock + kk];
    }
  }else if(yy == 0){
    for(kk = 0 ; kk < dd ; kk++){
      total = total + PP[baseNode + kk] * BB_inv[baseBlock + kk];
    }
  }
  return(total);
}

double sbm_t::tieLL_receiver(int yy, int recBlock, int sender){
  double total = 0.0;
  int kk;
  if(yy == 1){
    for(kk = 0 ; kk < dd ; kk++){
      total = total + PP[sender * dd + kk] * BB[kk * dd + recBlock];
    }
  }
  if(yy == 0){
    for(kk = 0 ; kk < dd ; kk++){
      total = total + PP[sender * dd + kk] * BB_inv[kk * dd + recBlock];
    }
  }
  return(total);
}


/*******************  BB UPDATING FUNCTIONS  ******************/

void sbm_t::getBB(double *BB_out){
  std::copy(BB,BB + dd*dd,BB_out);
}

double sbm_t::BBdiff(){
  expBB();
  int ii;
  double delta = 0.0;
  for(ii = 0 ; ii < dd*dd; ii++){
    delta = delta + ((BB[ii] - BB_old[ii]) * (BB[ii] - BB_old[ii]));
  }
  return(delta);
}

void sbm_t::saveBB_old(){
  expBB();
  std::copy(BB,BB+dd*dd,BB_old);
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






