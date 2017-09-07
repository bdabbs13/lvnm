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
//using namespace std using std::ifstream;
//using std::ofstream;

//gsl_rng *rng;
//ofstream myfile;


/*****************************************************************
 *********  FUNCTION DEFINITIONS FOR SBM OBJECT CLASS  ***********
 *****************************************************************/


//  Funciton for Performing EM Algorithm
void sbmEM(CSBM *mySBM, int iter_max, double threshold,
	   double *rPosteriorMemb, double *rBlockMat, int *rBlockMemb,
	   double *logLik, double *rPriorBlockMat, int verbose){
    int iter = 0;
    double delta = threshold + 1.0;

    while((iter < iter_max) & (delta >= threshold)){
	iter = iter + 1;
	mySBM->saveBlockMatOld();

	if(verbose > 2){
	    mySBM->print(false);
	}
	mySBM->iterEM();
	delta = mySBM->BlockMatdiff();

	if(verbose > 1){
	    Rprintf("iter - %d delta - %.6f\n",iter,delta);
	}
    }
    if(verbose == 1){
	Rprintf("iter - %d delta - %.6f\n",iter,delta);
    }
    //mySBM->print(1 == 0);
    mySBM->write(0);
    mySBM->writeRPriorBlockMemb(rPriorBlockMat);
    logLik[0] = mySBM->LogLike();
}





/***********  SBM CLASS CONSTRUCTOR AND DESTRUCTOR  **********/

//  Constructor Function for SBM object

CSBM::CSBM (int rNodes, int rBlocks, int mImpute) : missingVal(-1){
    int ii, jj, kk;

    // Loading Constants
    aNodes = rNodes; aBlocks = rBlocks;

    // Loading Adjacency Matrix
    aImputeFlag = (mImpute == 1);

    //  Initializing aBlockMat Matrix
    aBlockMat.resize(aBlocks);
    aBlockMatInv.resize(aBlocks);
    aBlockMatOld.resize(aBlocks);
    for(ii = 0; ii < aBlocks ; ii++){
	aBlockMat[ii].resize(aBlocks);
	aBlockMatInv[ii].resize(aBlocks);
	aBlockMatOld[ii].resize(aBlocks,0.0);
    }

    is_BlockMat_logged = false;

    aHitMat.resize(aBlocks);
    aMissMat.resize(aBlocks);
    for(ii = 0 ; ii < aBlocks ; ii++){
	aHitMat[ii].resize(aBlocks,0.0);
	aMissMat[ii].resize(aBlocks,0.0);
    }

    aPosteriorMemb.resize(aNodes);
    aPosteriorMembOld.resize(aNodes);

    aBlockMemb.resize(aNodes,0);
    for(ii = 0 ; ii < aNodes ; ii++){
	aPosteriorMemb[ii].resize(aBlocks);
	aPosteriorMembOld[ii].resize(aBlocks);
    }

}



CSBM::~CSBM (){

}

void CSBM::loadDataR(int *adjMat, double *rPriorBlockMat,
		     double *rPriorBlockMemb){
    int ii, jj;

    aAdjMat.resize(aNodes);
    for(ii = 0; ii < aNodes ; ii++){
	aAdjMat[ii].resize(aNodes);
	for(jj = 0; jj < aNodes; jj++){
	    aAdjMat[ii][jj] = adjMat[ii + jj*aNodes];
	}
    }

    if(aImputeFlag){
	aAdjPartial.resize(aNodes);
	for(ii = 0; ii < aNodes; ii++){
	    aAdjPartial[ii].resize(aNodes);
	    for(jj = 0; jj < aNodes; jj++){
		aAdjPartial[ii][jj] = adjMat[ii + jj*aNodes];
	    }
	}
    }

    //  Loading Block Matrix Prior
    std::copy(rPriorBlockMat,rPriorBlockMat+2,aPriorBlockMat);

    // Loading Block Membership Prior
    aPriorBlockMemb.insert(aPriorBlockMemb.end(),
			   &rPriorBlockMemb[0],&rPriorBlockMemb[aBlocks]);

}

void CSBM::loadStateR(double *rBlockMat, int *rBlockMemb,
		      double *rPosteriorMemb, double *rLogLik){

    //  Setting up R Writer
    aWriter = R_WRITER;

    aRBlockMat = rBlockMat;
    aRBlockMemb = rBlockMemb;
    aRPosteriorMemb = rPosteriorMemb;
    aRLogLike = rLogLik;

    //  Loading Current Memory
    RLoadBlockMat(rBlockMat);
    RLoadBlockMemb(rBlockMemb);
}

void CSBM::RLoadSBM(double *rBlockMat, int *rBlockMemb){
    RLoadBlockMat(rBlockMat);
    RLoadBlockMemb(rBlockMemb);
}


void CSBM::initRandom(){

    int ii, jj, kk;

    double holder[2];
    for(ii = 0; ii < aBlocks ; ii++){
	for(jj = 0 ; jj < aBlocks ; jj++){
	    rdirichlet(2,aPriorBlockMat,holder);
	    aBlockMat[ii][jj] = holder[0];
	    aBlockMatInv[ii][jj] = 1 - holder[0];
	}
    }

    //  Initializing Block Memberships
    double *etaInit = new double[aBlocks];
    double etaTotal = 0.0;
    for(kk = 0 ; kk < aBlocks ; kk++){
	etaTotal += aPriorBlockMemb[kk];
    }
    for(kk = 0 ; kk < aBlocks ; kk++){
	etaInit[kk] = aPriorBlockMemb[kk] / etaTotal;
    }

    int *ind = new int[aBlocks];

    for(ii = 0 ; ii < aNodes ; ii++){
	rmultinom(1,etaInit,aBlocks,ind);

	// Setting aPosteriorMemb[ii,] Value

	for(kk = 0 ; kk < aBlocks ; kk++){
	    if(ind[kk] == 1){
		aPosteriorMemb[ii][kk] = 1.0;
		aBlockMemb[ii] = kk;
	    }else{
		//	    aPosteriorMemb[ii*aBlocks + kk] = 0.0;
		aPosteriorMemb[ii][kk] = 0.0;
	    }
	}
    }

    delete[] etaInit;
    delete[] ind;
}


void CSBM::step(){
    drawBlockMemb();
    drawBlockMat();
    rotate();
    if(aImputeFlag){
	//      Rprintf("Imputing Missing Values...\n");
	imputeMissingValues();
    }
}

void CSBM::write(int iter){
    if(aWriter == R_WRITER){
	writeR(iter);
    }
}

double CSBM::LogLike(){
    int rr, cc;
    logBlockMat();

    double total = 0.0;
    for(rr = 0 ; rr < aNodes ; rr++){
	for(cc = 0 ; cc < aNodes ; cc++){
	    total = total + tieLogLike(aAdjMat[rr][cc],aBlockMemb[rr],aBlockMemb[cc]);
	}
    }
    return(total);
}


/********************  Functions used by EM  ********************/

void CSBM::writeRPriorBlockMemb(double *rPriorBlockMemb){

    std::copy(&aPriorBlockMemb[0], &aPriorBlockMemb[aBlocks],
	      rPriorBlockMemb);
}

void CSBM::initPPem(double minVal){
    int ii, kk;
    for(ii = 0 ; ii < aNodes ; ii++){
	for(kk = 0 ; kk < aBlocks ; kk++){
	    if(aBlockMemb[ii] == kk){
		aPosteriorMemb[ii][kk] = 1.0 - minVal;
	    }else{
		aPosteriorMemb[ii][kk] = minVal / (aBlocks-1);
	    }
	}
    }
}


void CSBM::iterEM (){
    int ss, rr, ii, jj;
    int rowd, rown,bbind;
    double total, pipi;

    //  E Step
    savePosteriorMembOld();
    getMultinomPosterior();  //  Update Posterior Mean aPosteriorMemb

    // M Step
    colSums(aPosteriorMemb,aPriorBlockMemb);
    double priorTotal = 0.0;

    for(ii = 0 ; ii < aBlocks ; ii++){
	aPriorBlockMemb[ii] = aPriorBlockMemb[ii] / aNodes;
	if(aPriorBlockMemb[ii] < MIN_LOG){
	    aPriorBlockMemb[ii] = MIN_LOG;
	}
	priorTotal = priorTotal + aPriorBlockMemb[ii];
    }
    for(ii = 0 ; ii < aBlocks ; ii++){
	aPriorBlockMemb[ii] = aPriorBlockMemb[ii] / priorTotal;
    }

    // Updating BlockMat with the MLE conditional on aPosteriorMemb and YY
    for(ss = 0 ; ss < aBlocks ; ss++){
	for(rr = 0 ; rr < aBlocks ; rr++){
	    total = 0.0;
	    bbind = ss*aBlocks + rr;
	    aBlockMat[ss][rr] = 0.0;
	    //aBlockMat[bbind] = 0.0;
	    for(ii = 0 ; ii < aNodes ; ii++){
		rowd = ii*aBlocks;
		rown = ii*aNodes;
		for(jj = 0 ; jj < aNodes ; jj++){
		    if(!isMissing(aAdjMat[ii][jj])){
			pipi = aPosteriorMemb[ii][ss] * aPosteriorMemb[jj][rr];
			aBlockMat[ss][rr] += pipi * aAdjMat[ii][jj]; // Numerator
			total = total + pipi;  //  Denominator
		    }
		}
	    }
	    if(total < MIN_LOG){
		aBlockMat[ss][rr] = MIN_LOG;
	    }else{
		aBlockMat[ss][rr] = aBlockMat[ss][rr] / total;
		aBlockMat[ss][rr] = logCheck(aBlockMat[ss][rr]);
	    }
	    aBlockMatInv[ss][rr] = 1.0 - aBlockMat[ss][rr];
	}
    }
    is_BlockMat_logged = false;
}

//  This function updates aPosteriorMembOld to be the posterior mean given BlockMat and YY
void CSBM::getMultinomPosterior(){

    int ii,jj;
    double pMax;

    for(ii = 0 ; ii < aNodes ; ii++){
	for(jj = 0 ; jj < aBlocks ; jj++){
	    aBlockMemb[ii] = jj;

	    aPosteriorMemb[ii][jj] = nodeLogLike_long(ii) + log(aPriorBlockMemb[jj]);
	    if(jj == 0){
		pMax = aPosteriorMemb[ii][jj];
	    }else{
		if(aPosteriorMemb[ii][jj] > pMax){
		    pMax = aPosteriorMemb[ii][jj];
		}
	    }
	}

	// Exponentiating Probabilities
	for(jj = 0 ; jj < aBlocks ; jj++){
	    aPosteriorMemb[ii][jj] = exp(aPosteriorMemb[ii][jj] - pMax);
	}

	normalizeVec(aPosteriorMemb[ii]);
    }
}

void CSBM::saveBlockMatOld(){
    expBlockMat();
    //   std::copy(BlockMat,BlockMat+aBlocks*aBlocks,aBlockMatOld);
    int ii, jj;
    for(ii = 0 ; ii < aBlocks ; ii++){
	for(jj = 0 ; jj < aBlocks ; jj++){
	    aBlockMatOld[ii][jj] = aBlockMat[ii][jj];
	}
    }
}

double CSBM::BlockMatdiff(){
    expBlockMat();
    int ii,jj;
    double delta = 0.0;
    for(ii = 0 ; ii < aBlocks; ii++){
	for(jj = 0 ; jj < aBlocks ; jj++){
	    delta = delta + ((aBlockMat[ii][jj] - aBlockMatOld[ii][jj]) * (aBlockMat[ii][jj] - aBlockMatOld[ii][jj]));
	}
    }
    return(delta);
}

void CSBM::computeBlockMatMLE(){
    computeHitMiss();

    double hit, miss;
    int ii,jj;//dd2 = aBlocks*aBlocks;
    for(ii = 0 ; ii < aBlocks ; ii++){
	for(jj = 0 ; jj < aBlocks ; jj++){
	    if(aHitMat[ii][jj] == 0.0){
		aBlockMat[ii][jj] = 0.0;
	    }else{
		hit = aHitMat[ii][jj];
		miss = aMissMat[ii][jj];
		aBlockMat[ii][jj] = hit / (hit + miss);
	    }
	    aBlockMat[ii][jj] = logCheck(aBlockMat[ii][jj]);
	    aBlockMatInv[ii][jj] = 1.0 - aBlockMat[ii][jj];
	}
    }
}

void CSBM::print(bool printNetwork){
    Rprintf("aBlocks = %d, aNodes = %d\n",aBlocks,aNodes);
    Rprintf("bbPrior = (%.2f, %.2f)\n",aPriorBlockMat[0],aPriorBlockMat[1]);
    Rprintf("multiPrior = ");
    printPriorBlockMemb();
    printBlockMat();
    printPosteriorMemb();
    if(printNetwork){
	printAdjacencyMatrix();
    }

}


/**********************************************************************/
/*************************  PRIVATE FUNCTIONS *************************/
/**********************************************************************/

void CSBM::RLoadBlockMat(double *rBlockMat){
    int ii,jj;
    int loadIter = 0;

    for(jj = 0 ; jj < aBlocks ; jj++){
	for(ii = 0; ii < aBlocks ; ii++){
	    aBlockMat[ii][jj] = rBlockMat[loadIter++];
	    aBlockMatInv[ii][jj] = 1.0 - aBlockMat[ii][jj];
	}
    }
}

void CSBM::RLoadBlockMemb(int *rBlockMemb){
    int ii;
    for(ii = 0 ; ii < aNodes ; ii++){
	aBlockMemb[ii] = rBlockMemb[ii] - 1;
    }
}

void CSBM::RLoadPosteriorMemb(double *rPosteriorMemb){
    int ii, kk;
    int loadIter = 0;
    for(kk = 0 ; kk < aBlocks; kk++){
	for(ii = 0 ; ii < aNodes ; ii++){
	    aPosteriorMemb[ii][kk] = rPosteriorMemb[loadIter++];
	}
    }
}


/*****  step Functions  *****/

void CSBM::drawBlockMemb(){
    int ii,jj,kk;
    double pMax;
    int ind[aBlocks];

    for(ii = 0 ; ii < aNodes ; ii++){
	for(jj = 0 ; jj < aBlocks ; jj++){
	    aBlockMemb[ii] = jj;

	    // Calculating Log Posterior Probability
	    aPosteriorMemb[ii][jj] = nodeLogLike(ii) + log(aPriorBlockMemb[jj]);
	    if(jj == 0){
		pMax = aPosteriorMemb[ii][jj];
	    }else{
		if(aPosteriorMemb[ii][jj] > pMax){
		    pMax = aPosteriorMemb[ii][jj];
		}
	    }
	}

	// Exponentiating Probabilities
	for(jj = 0 ; jj < aBlocks ; jj++){
	    aPosteriorMemb[ii][jj] = exp(aPosteriorMemb[ii][jj] - pMax);
	}

	normalizeVec(aPosteriorMemb[ii]);
	logZeroFix(aPosteriorMemb[ii]);

	// Drawing from Multinomial
	rmultinom(1,aPosteriorMemb[ii].data(),aBlocks,ind);

	// Setting aPosteriorMemb[ii,] Value
	for(kk = 0 ; kk < aBlocks ; kk++){
	    if(ind[kk] == 1){
		aBlockMemb[ii] = kk;
		break;
	    }
	}
    }
}


void CSBM::drawBlockMat(){
    int ii,jj;

    // Setting Up AHitMat and Miss Counts
    computeHitMiss();

    //  Drawing From Posterior Distribution
    double alpha_new[2];
    double holder[2];
    for(ii = 0 ; ii < aBlocks ; ii++){
	for(jj = 0 ; jj < aBlocks ; jj++){
	    alpha_new[0] = aPriorBlockMat[0] + aHitMat[ii][jj];
	    alpha_new[1] = aPriorBlockMat[1] + aMissMat[ii][jj];
	    rdirichlet(2,alpha_new,holder);
	    aBlockMat[ii][jj] = holder[0];
	    aBlockMatInv[ii][jj] = 1.0 - aBlockMat[ii][jj];
	}
    }
    is_BlockMat_logged = false;
}



void CSBM::imputeMissingValues(){
    expBlockMat();
    int rr, cc, index;
    double tie_prob;
    for(rr = 0 ; rr < aNodes ; rr++){ // row loop
	for(cc = 0 ; cc < aNodes ; cc++){ // col loop
	    if(rr != cc){ // ignore diagonal entries
		if(isMissing(aAdjPartial[rr][cc])){
		    tie_prob = aBlockMat[aBlockMemb[rr]][aBlockMemb[cc]];
		    aAdjMat[rr][cc] = (int) rbinom(1.0,tie_prob);
		}
	    }
	}
    }
}


void CSBM::rotate(){
    int ii, jj, kk;
    int assigned[aBlocks], replaced[aBlocks];
    for(ii = 0 ; ii < aBlocks ; ii++){
	assigned[ii] = -1;
	replaced[ii] = -1;
    }
    int PPtmp[aNodes];
    for(ii = 0 ; ii < aNodes ; ii++){
	PPtmp[ii] = aBlockMemb[ii];
    }
    int replacements[aNodes];
    for(ii = 0 ; ii < aNodes ; ii++){
	replacements[ii] = -1;
    }

    for(ii = 0 ; ii < aBlocks ; ii++){
	if(replacements[ii] == -1){
	    assigned[PPtmp[ii]] = ii;
	    replaced[ii] = PPtmp[ii];
	    for(jj = 0 ; jj < aNodes ; jj++){
		if(PPtmp[jj] == PPtmp[ii]){
		    replacements[jj] = ii;
		}
	    }
	}
    }

    int uncount = 0;
    for(ii = 0 ; ii < aBlocks ; ii++){
	if(assigned[ii] == -1){
	    uncount++;
	}
    }
    int vacancy, newlabel;
    while(uncount > 0){

	//  Setting Vacancy to Minimum Unassigned
	for(ii = 0 ; ii < aBlocks ; ii++){
	    if(assigned[ii] == -1){
		vacancy = ii;
		break;
	    }
	}
	//  Setting newlabel to Minimum unreplaced
	for(ii = 0 ; ii < aBlocks ;  ii++){
	    if(replaced[ii] == -1){
		newlabel = ii;
		break;
	    }
	}

	assigned[vacancy] = newlabel;
	replaced[newlabel] = vacancy;
	for(jj = 0 ; jj < aNodes ; jj++){
	    if(PPtmp[jj] == vacancy){
		replacements[jj] = newlabel;
	    }
	}
	uncount = 0;
	for(ii = 0 ; ii < aBlocks ; ii++){
	    if(assigned[ii] == -1){
		uncount++;
	    }
	}

    }

    //  Updating Block Matrix
    saveBlockMatOld();  //  Putting current version of BlockMat in aBlockMatOld
    for(ii = 0 ; ii < aBlocks ; ii++){
	for(jj = 0 ; jj < aBlocks ; jj++){
	    aBlockMat[ii][jj] = aBlockMatOld[assigned[ii]][assigned[jj]];
	    aBlockMatInv[ii][jj] = 1.0 - aBlockMat[ii][jj];
	}
    }
}

/******************  Log-Likelihood Functions  *****************/

double CSBM::nodeLogLike(int ii){
    logBlockMat();
    int jj;
    int iiBlockMemb = aBlockMemb[ii];
    double total = 0.0;

    //  sums the (i,j)th and (j,i)th probabilities for all j
    for(jj = 0 ; jj < aNodes ; jj++){
	if(jj != ii){

	    //  (i,j)th term
	    total = total + tieLogLike(aAdjMat[ii][jj],iiBlockMemb,aBlockMemb[jj]);

	    //  (j,i)th term
	    total = total + tieLogLike(aAdjMat[jj][ii],aBlockMemb[jj],iiBlockMemb);
	}
    }

    return(total);
}


double CSBM::tieLogLike(int yy, int sendBlock, int recBlock){

    if(yy == 1){
	return aBlockMat[sendBlock][recBlock];
    }else if(yy == 0){
	return aBlockMatInv[sendBlock][recBlock];
    }else{
	return(0.0);
    }
}



double CSBM::nodeLogLike_long(int ii){
    logBlockMat();
    int jj;
    int iiBlockMemb = aBlockMemb[ii];
    double total = 0.0;

    //  sums the (i,j)th and (j,i)th probabilities for all j
    for(jj = 0 ; jj < aNodes ; jj++){
	if(jj != ii){
	    //  (i,j)th term
	    total = total + tieLogLike_sender(aAdjMat[ii][jj],iiBlockMemb,jj);
	    //  (j,i)th term
	    total = total + tieLogLike_receiver(aAdjMat[jj][ii],iiBlockMemb,jj);
	}
    }

    return(total);
}


// Log-Likelihood assuming Posterior Matrix
double CSBM::tieLogLike_sender(int yy, int sendBlock, int receiver){
    int kk;
    double total = 0.0;
    if(yy == 1){
	for(kk = 0 ; kk < aBlocks ; kk++){
	    total += aPosteriorMembOld[receiver][kk] * aBlockMat[sendBlock][kk];
	}
    }else if(yy == 0){
	for(kk = 0 ; kk < aBlocks ; kk++){
	    total += aPosteriorMembOld[receiver][kk] * aBlockMatInv[sendBlock][kk];
	}
    }
    return(total);
}

double CSBM::tieLogLike_receiver(int yy, int recBlock, int sender){
    double total = 0.0;
    int kk;
    if(yy == 1){
	for(kk = 0 ; kk < aBlocks ; kk++){
	    total += aPosteriorMembOld[sender][kk] * aBlockMat[kk][recBlock];
	}
    }
    if(yy == 0){
	for(kk = 0 ; kk < aBlocks ; kk++){
	    total += aPosteriorMembOld[sender][kk] * aBlockMatInv[kk][recBlock];
	}
    }
    return(total);
}


/**********  R_WRITER  Functions  **********/

void CSBM::writeR(int iter){

    aRLogLike[iter] = LogLike();
    writeRBlockMat(iter);
    writeRBlockMemb(iter);
    writeRPosteriorMemb(iter);
}

void CSBM::writeRBlockMat(int iter){
    expBlockMat();
    int ii, jj;
    int saveIter = iter * aBlocks * aBlocks;
    for(jj = 0 ; jj < aBlocks ; jj++){
	for(ii = 0 ; ii < aBlocks ; ii++){
	    aRBlockMat[saveIter++] = aBlockMat[ii][jj];
	}
    }
}

void CSBM::writeRBlockMemb(int iter){
    int ii;
    int saveIter = iter * aNodes;
    for(ii = 0 ; ii < aNodes ; ii++){
	aRBlockMemb[saveIter++] = aBlockMemb[ii] + 1;
    }
}

void CSBM::writeRPosteriorMemb(int iter){
    int ii, kk;
    int saveIter = iter * aBlocks * aNodes;
    for(kk = 0 ; kk < aBlocks; kk++){
	for(ii = 0 ; ii < aNodes ; ii++){
	    aRPosteriorMemb[saveIter++] =  aPosteriorMemb[ii][kk];
	}
    }

}




/*****************  HELPER FUNCTIONS  *****************/
void CSBM::savePosteriorMembOld(){
    //   std::copy(BlockMat,BlockMat+aBlocks*aBlocks,aBlockMatOld);
    int ii, kk;
    for(ii = 0 ; ii < aNodes ; ii++){
	for(kk = 0 ; kk < aBlocks ; kk++){
	    aPosteriorMembOld[ii][kk] = aPosteriorMemb[ii][kk];
	}
    }
}

void CSBM::computeHitMiss(){
    int ii, rr, cc, index;

    //  Setting values to zero
    for(ii = 0 ; ii < aBlocks ; ii++){
	aHitMat[ii].assign(aBlocks,0.0);
	aMissMat[ii].assign(aBlocks,0.0);
    }

    //  Counting hits and misses
    for(rr = 0 ; rr < aNodes ; rr++){
	for(cc = 0 ; cc < aNodes ; cc++){
	    if(aAdjMat[rr][cc] == 1){
		aHitMat[aBlockMemb[rr]][aBlockMemb[cc]]++;
	    }else if(aAdjMat[rr][cc] == 0){
		aMissMat[aBlockMemb[rr]][aBlockMemb[cc]]++;
	    }
	}
    }
}


/***** BlockMat Updating Functions *****/
void CSBM::logBlockMat(){
    int ii,jj;
    if(!is_BlockMat_logged){
	for(ii = 0 ; ii < aBlocks ; ii++){
	    for(jj = 0 ; jj < aBlocks ; jj++){
		aBlockMat[ii][jj] = log(aBlockMat[ii][jj]);
		aBlockMatInv[ii][jj] = log(aBlockMatInv[ii][jj]);
	    }
	}
	is_BlockMat_logged = true;
    }
}

void CSBM::expBlockMat(){
    int ii,jj;
    if(is_BlockMat_logged){
	for(ii = 0 ; ii < aBlocks ; ii++){
	    for(jj = 0 ; jj < aBlocks ; jj++){
		aBlockMat[ii][jj] = exp(aBlockMat[ii][jj]);
		aBlockMatInv[ii][jj] = exp(aBlockMatInv[ii][jj]);
	    }
	}
	is_BlockMat_logged = false;
    }
}

void CSBM::GetBlockMat(double *rBlockMat){
    //   std::copy(aBlockMat,aBlockMat + aBlocks*aBlocks,BlockMat_out);
    int ii,jj;
    int saveIter = 0;
    for(jj = 0 ; jj < aBlocks ; jj++){
	for(ii = 0 ; ii < aBlocks ; ii++){
	    rBlockMat[saveIter++] = aBlockMat[ii][jj];
	    //BlockMat_out[ii + jj*aBlocks] = aBlockMat[ii][jj];
	}
    }
}


// Static Version
void CSBM::printAdjacencyMatrix(){

    int ii, jj;
    Rprintf("Adjacency Matrix:\n");
    for(ii = 0 ; ii < aNodes; ii++){
	for(jj = 0; jj < aNodes; jj++){
	    Rprintf("%d ",aAdjMat[ii][jj]);
	}
	Rprintf("\n");
    }

}

void CSBM::printBlockMat(){
    int ii,jj;
    Rprintf("Block Matrix:\n");
    for(ii = 0 ; ii < aBlocks; ii++){
	for(jj = 0; jj < aBlocks; jj++){
	    Rprintf("%.2f ",aBlockMat[ii][jj]);
	}
	Rprintf("\n");
    }
}

void CSBM::printPosteriorMemb(){
    int ii, kk;
    Rprintf("Block Memberships:\n");
    for(ii = 0 ; ii < aNodes ; ii++){
	for(kk = 0 ; kk < aBlocks ; kk++){
	    Rprintf("%.2f ",aPosteriorMemb[ii][kk]);
	}
	Rprintf("\n");
    }
}

void CSBM::printBlockMemb(){
    int ii;
    Rprintf("Block Memberships:\n");
    for(ii = 0; ii < aNodes; ii++){
	Rprintf("%d ",aBlockMemb[ii]);
    }
    Rprintf("\n\n");
}

void CSBM::printPriorBlockMemb(){
    int ii;
    Rprintf("Block Membership Prior:\n");
    for(ii = 0 ; ii < aBlocks ; ii++){
	Rprintf("%d ",aPriorBlockMemb[ii]);
    }
}

