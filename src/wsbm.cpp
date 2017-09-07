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
#include "wsbm.h"
#include "helper.h"
//using namespace std;
using std::ifstream;
//using std::ofstream;

//gsl_rng *rng;
//ofstream myfile;




void wsbmMCMC(CWSBM *myWSBM, int start, int total, int burnIn, int thin,
	      int shift_size, int extend_max, double qq, int verbose){

    int ii, step_count;
    int extend_count = -1;
    std::vector<double> ll_vec;

    // Main MCMC Loop
    int raw_count = start;
    int step_total = (total - burnIn) / thin;
    ll_vec.resize(step_total,0.0);

    bool converged = false;
    do{
    	step_count = 0;
	extend_count++;

	if(extend_count > 0){
	    if(verbose > 1) Rprintf("Extending Chain...\n");
	}

    	while(step_count < step_total){
    	    // Perform an MCMC Step
    	    myWSBM->step();
    	    myWSBM->adapt();


    	    if((raw_count >= burnIn) && ((raw_count - burnIn) % thin == 0)){
    		//  Save the result every thin iterations
    		myWSBM->write(step_count);
		ll_vec[step_count] = myWSBM->LogLike();

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
    }
    while((!(converged)) && (extend_count < extend_max));


    if(!converged){
	Rprintf("Warning: Failed to converge after %d extensions\n",
		extend_count);
    }
    if(verbose > 2) myWSBM->print(false);

}






/*****************************************************************
 *********  FUNCTION DEFINITIONS FOR WSBM OBJECT CLASS  ***********
 *****************************************************************/


/***********  WSBM CLASS CONSTRUCTOR AND DESTRUCTOR  **********/
//  Just Initializing Storage Containers
CWSBM::CWSBM (int rNodes, int rBlocks, int mImpute) : missingVal(-1){
    int ii, jj, kk;

    // Rprintf("Creating WSBM Object\n");
    //   Loading Constants
    aNodes = rNodes; aBlocks = rBlocks;

    //  Allocating space for adjacency matrix
    aAdjMat.resize(aNodes);
    for(ii = 0; ii < aNodes ; ii++){
	aAdjMat[ii].resize(aNodes);
    }

    aImputeFlag = (mImpute == 1);
    if(aImputeFlag){
	aAdjPartial.resize(aNodes);
	for(ii = 0; ii < aNodes; ii++){
	    aAdjPartial[ii].resize(aNodes);
	}
    }

    /*****  Allocating space for parameters  *****/

    // Allocating Space for Block Matrices
    aBlockMat.resize(aBlocks);
    //   aBlockMatLog.resize(aBlocks);
    aBlockMatOld.resize(aBlocks);
    for(ii = 0; ii < aBlocks ; ii++){
	aBlockMat[ii].resize(aBlocks,0.0);
	// aBlockMatLog[ii].resize(aBlocks,0.0);
	aBlockMatOld[ii].resize(aBlocks,0.0);
    }

    // Allocating Space for Intermediate Block Matrix Values
    aBlockTieCounts.resize(aBlocks);
    aBlockTieSums.resize(aBlocks);
    for(ii = 0 ; ii < aBlocks ; ii++){
	aBlockTieSums[ii].resize(aBlocks,0.0);
	aBlockTieCounts[ii].resize(aBlocks,0.0);
    }

    // Allocating Space for Sender/Receiver Effects
    aSenderEffects.resize(aNodes,1.0);
    aReceiverEffects.resize(aNodes,1.0);

    // Allocating Space for Intermediate SR Values
    aAdjMatRowSums.resize(aNodes,0.0);
    aAdjMatColSums.resize(aNodes,0.0);

    aBlockSenderEffects.resize(aNodes,0.0);
    aBlockReceiverEffects.resize(aNodes,0.0);



    /********************  Only needed for raw WSBM  *********************/
    // Rprintf("Initializing Priors\n");
    // Default Block Membership Prior
    aPriorBlockMemb.resize(aBlocks,1.0);

    aPriorBlockMat = new std::vector<std::vector<double *> >(aBlocks);
    for(kk = 0 ; kk < aBlocks ; kk++){
	(*aPriorBlockMat)[kk].resize(aBlocks,aLocalPriorBlockMat);
    }

    aPriorSender = new std::vector<double *>(aNodes,aLocalPriorSender);
    aPriorReceiver = new std::vector<double *>(aNodes,aLocalPriorSender);

    //  Allocating space for block membership vectors
    aBlockMemb = new std::vector<int>;
    (*aBlockMemb).resize(aNodes,0);

    aPriorOwner = true;

    //  Allocating space for posterior membership probabilities
    aPosteriorMemb.resize(aNodes);
    //    aPosteriorMembOld.resize(aNodes);
    for(ii = 0 ; ii < aNodes ; ii++){
	aPosteriorMemb[ii].resize(aBlocks);
	//	aPosteriorMembOld[ii].resize(aBlocks);
    }
    // Rprintf("Finished Initializing\n");
}

CWSBM::~CWSBM (){
    if(aPriorOwner){
	// Rprintf("I Own this Prior\n");
	deallocatePriors();
    }
}

void CWSBM::loadDataR(int *adjMat, double rHours,
		      double *rPriorSender, double *rPriorReceiver,
		      double *rPriorBlockMat, double * rPriorBlockMemb){
    int ii, jj, kk;

    LoadAdjacencyMatrix(adjMat);

    // Loading hours
    aHours = rHours;

    /*****  Loading Priors  *****/
    std::copy(rPriorSender,rPriorSender+2,aLocalPriorSender);
    std::copy(rPriorReceiver,rPriorReceiver+2,aLocalPriorReceiver);
    std::copy(rPriorBlockMat,rPriorBlockMat+2,aLocalPriorBlockMat);

    // Loading Prior Block Membership
    aPriorBlockMemb.insert(aPriorBlockMemb.end(),
			   &rPriorBlockMemb[0],&rPriorBlockMemb[aBlocks]);
}


void CWSBM::loadStateR(double *rBlockMat, int *rBlockMemb,
		       double *rSenderEffects, double *rReceiverEffects,
		       double *rPosteriorMemb, double *rLogLik){

    aStepType = FULL_STEP;

    //  Saving Pointers for Writer
    aWriter = R_WRITER;

    aRBlockMat = rBlockMat;
    aRBlockMemb = rBlockMemb;
    aRSenderEffects = rSenderEffects;
    aRReceiverEffects = rReceiverEffects;
    aRPosteriorMemb = rPosteriorMemb;
    aRLogLike = rLogLik;

    //  Loading Initial Values
    RLoadBlockMat(rBlockMat);
    RLoadBlockMemb(rBlockMemb);
    RLoadPosteriorMemb(rPosteriorMemb);
    RLoadSenderEffects(rSenderEffects);
    RLoadReceiverEffects(rReceiverEffects);
}


void CWSBM::loadStateR(double *rBlockMat,
		       double *rSenderEffects, double *rReceiverEffects){

    aStepType = PARTIAL_STEP;

    //  Saving Pointers for Writer
    aWriter = DYN_R_WRITER;

    aRBlockMat = rBlockMat;
    aRSenderEffects = rSenderEffects;
    aRReceiverEffects = rReceiverEffects;

    //  Loading values into C++ Objects
    RLoadBlockMat(rBlockMat);
    RLoadSenderEffects(rSenderEffects);
    RLoadReceiverEffects(rReceiverEffects);
}

void CWSBM::LoadReferences(double dHours,
			   std::vector<double *> *dPriorSender,
			   std::vector<double *> *dPriorReceiver,
			   std::vector<std::vector<double *> > *dPriorBlockMat,
			   std::vector<int> *dBlockMemb){

    // Rprintf("Loading References\n");
    aHours = dHours;

    aPriorOwner = false;
    deallocatePriors();

    aPriorSender = dPriorSender;
    aPriorReceiver = dPriorReceiver;
    aPriorBlockMat = dPriorBlockMat;
    aBlockMemb = dBlockMemb;
}

void CWSBM::LoadAdjacencyMatrix(int *AdjMat){
    int ii, jj;
    // Rprintf("Loading Adjacency Matrix\n");
    for(ii = 0; ii < aNodes ; ii++){
	for(jj = 0; jj < aNodes; jj++){
	    aAdjMat[ii][jj] = AdjMat[ii + jj*aNodes];
	}
    }

    if(aImputeFlag){
	for(ii = 0; ii < aNodes; ii++){
	    for(jj = 0; jj < aNodes; jj++){
		aAdjPartial[ii][jj] = AdjMat[ii*aNodes + jj];
	    }
	}
    }

    computeRowColSums();
}


//  Currently an unused function, but could be useful in the future
void CWSBM::initRandom(){
    int ii, jj, kk;

    /*********  Random Initialization  ***********/
    double holder[2];
    // Initialization from Prior for Block Matrix
    for(ii = 0; ii < aBlocks ; ii++){
	for(jj = 0 ; jj < aBlocks ; jj++){
	    rdirichlet(2,(*aPriorBlockMat)[ii][jj],holder);
	    aBlockMat[ii][jj] = holder[0];
	}
    }

    // Initialization for Block Membership
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
		(*aBlockMemb)[ii] = kk;
	    }else{
		aPosteriorMemb[ii][kk] = 0.0;
	    }
	}
    }

    //  Initializing Sender and Receiver Effects
    /*
      for(ii = 0 ; ii < aNodes ; ii++){
      aSenderEffects[ii] = rgamma(aPriorSender[0],aPriorSender[1]);
      aReceiverEffects[ii] = rgamma(aPriorReceiver[0],aPriorReceiver[1]);
      }
    */

    delete[] etaInit;
    delete[] ind;

}

void CWSBM::step(){

    if(aStepType == FULL_STEP){
	fullStep();
    }else if(aStepType == PARTIAL_STEP){
	partialStep();
    }
}


void CWSBM::write(int iter){
    if(aWriter == R_WRITER){
	writeR(iter);
    }else if(aWriter == DYN_R_WRITER){
	writeDynR(iter);
    }else{
	Rprintf("Error:  Unknown writer provided\n");
    }
}

// bool CWSBM::check(double qq){


//     //  Checking for convergence of log-likelihood chain
//     // return (convergenceCheck(aRLogLike,save_iter,qq) > 0);
// }

double CWSBM::LogLike(){
    int rr, cc;
    double total = 0.0;
    for(rr = 0 ; rr < aNodes ; rr++){
	for(cc = 0 ; cc < aNodes ; cc++){
	    total += tieLogLike(aAdjMat[rr][cc],rr,cc);
	}
    }
    return(total);
}


double CWSBM::nodeLogLike(int ii){
    int jj;
    double total = 0.0;

    //  sums the (i,j)th and (j,i)th probabilities for all j
    for(jj = 0 ; jj < aNodes ; jj++){
	if(jj != ii){
	    //  (i,j)th term
	    total += tieLogLike(aAdjMat[ii][jj],ii,jj);
	    //  (j,i)th term
	    total += tieLogLike(aAdjMat[jj][ii],jj,ii);
	}else{
	    //  Only one term for i = j
	    total += tieLogLike(aAdjMat[ii][ii],ii,ii);
	}
    }

    return(total);
}

void CWSBM::computeBlockMatMLE(){
    computeBlockTieSums();

    int ll,kk;
    for(ll = 0 ; ll < aBlocks ; ll++){
	for(kk = 0 ; kk < aBlocks ; kk++){
	    if(aBlockTieSums[ll][kk] == 0.0){
		aBlockMat[ll][kk] = 0.0;
	    }else{
		aBlockMat[ll][kk] = aBlockTieSums[ll][kk]/aBlockTieCounts[ll][kk];
	    }
	    aBlockMat[ll][kk] = logCheck(aBlockMat[ll][kk]);
	    //	 aBlockMatLog[ll][kk] = log(aBlockMat[ll][kk]);
	}
    }
}


void CWSBM::print(bool printNetwork){

    Rprintf("aBlocks = %d, aNodes = %d, aHours = %.2f\n",aBlocks,aNodes,aHours);
    Rprintf("bbPrior = (%.2f, %.2f)\n",
	    (*aPriorBlockMat)[0][0][0],(*aPriorBlockMat)[0][0][1]);
    Rprintf("multiPrior = ");

    //   printPriorBlockMemb();
    printBlockMat();
    printBlockMemb();
    printSenderEffects();
    printReceiverEffects();
    if(printNetwork){
	printAdjacencyMatrix();
    }

}



/*******************************************************************************/
/*******************************************************************************/
/**************************  PRIVATE FUNCTIONS  ********************************/
/*******************************************************************************/
/*******************************************************************************/


void CWSBM::RLoadBlockMat(double *rBlockMat){
    int ii,jj;
    int loadIter = 0;

    for(jj = 0 ; jj < aBlocks ; jj++){
	for(ii = 0; ii < aBlocks ; ii++){
	    aBlockMat[ii][jj] = rBlockMat[loadIter++];
	    //	 aBlockMatLog[ii][jj] = log(aBlockMat[ii][jj]);
	}
    }

}

void CWSBM::RLoadBlockMemb(int *rBlockMemb){
    int ii;
    for(ii = 0 ; ii < aNodes ; ii++){
	(*aBlockMemb)[ii] = rBlockMemb[ii] - 1;
    }
}


void CWSBM::RLoadPosteriorMemb(double *rPosteriorMemb){
    int ii, kk;
    int loadIter = 0;
    for(kk = 0 ; kk < aBlocks; kk++){
	for(ii = 0 ; ii < aNodes ; ii++){
	    aPosteriorMemb[ii][kk] = rPosteriorMemb[loadIter++];
	}
    }
}


void CWSBM::RLoadSenderEffects(double * rSenderEffects){
    int ii;
    for(ii = 0 ; ii < aNodes ; ii++){
	aSenderEffects[ii] = rSenderEffects[ii];
    }
}

void CWSBM::RLoadReceiverEffects(double * rReceiverEffects){
    int ii;
    for(ii = 0 ; ii < aNodes ; ii++){
	aReceiverEffects[ii] = rReceiverEffects[ii];
    }
}


void CWSBM::deallocatePriors(){
    // Rprintf("Deallocating Priors\n");

    delete aPriorSender;
    delete aPriorReceiver;
    delete aPriorBlockMat;
    delete aBlockMemb;
}


/************************  MCMC FUNCTIONS  **********************/


void CWSBM::fullStep(){

    drawBlockMat();
    drawBlockMemb();
    drawSenderEffects();
    drawReceiverEffects();

    //  Need to check/update this function
    //   rotate();
    if(aImputeFlag){
	Rprintf("Imputing Missing Values...\n");
	imputeMissingValues();
	computeRowColSums();
    }

}

void CWSBM::partialStep(){

    drawBlockMat();
    drawSenderEffects();
    normalizeSenderEffects();

    drawReceiverEffects();
    normalizeReceiverEffects();

    if(aImputeFlag){
	Rprintf("Imputing Missing Values...\n");
	imputeMissingValues();
	computeRowColSums();
    }

}


void CWSBM::drawBlockMat(){

    int ll,kk;

    // Computing Block Tie Sums
    computeBlockTieSums();

    //  Drawing From Posterior Distribution
    double gammaPost[2];
    for(ll = 0 ; ll < aBlocks ; ll++){
	for(kk = 0 ; kk < aBlocks ; kk++){
	    gammaPost[0] = (*aPriorBlockMat)[ll][kk][0] + aBlockTieSums[ll][kk];
	    gammaPost[1] = (*aPriorBlockMat)[ll][kk][1] + aHours*aBlockTieCounts[ll][kk];

	    aBlockMat[ll][kk] = rgamma(gammaPost[0],1/gammaPost[1]);
	    //	 aBlockMatLog[ll][kk] = log(aBlockMat[ll][kk]);

	}
    }
}

void CWSBM::drawBlockMemb(){
    int ii,jj,kk;
    double pMax;
    int ind[aBlocks];

    //   Rprintf("Log-Likelihood = %.2f\n",LogLike());

    for(ii = 0 ; ii < aNodes ; ii++){
	//total = 0.0;
	for(jj = 0 ; jj < aBlocks ; jj++){
	    (*aBlockMemb)[ii] = jj;

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
	logZeroFix(aPosteriorMemb[ii]); // If aPosteriorMemb ~ 0 log(0) = NaN

	// Drawing from Multinomial
	rmultinom(1,aPosteriorMemb[ii].data(),aBlocks,ind);

	//  Getting Block Membership Number
	kk = 0;
	while(ind[kk] != 1){
	    kk++;
	}
	(*aBlockMemb)[ii] = kk;
    }
}

void CWSBM::drawSenderEffects(){
    int ii;
    computeBlockReceiverEffects();
    double gammaPost[2];
    //   Rprintf("Drawing Sender Effects:\n");

    for(ii = 0 ; ii < aNodes ; ii++){
	gammaPost[0] = aAdjMatRowSums[ii] + (*aPriorSender)[ii][0];
	gammaPost[1] = aHours*aBlockReceiverEffects[ii] + (*aPriorSender)[ii][1];
	aSenderEffects[ii] = rgamma(gammaPost[0],1/gammaPost[1]);

	//      Rprintf("i = %d, alpha = %.3f, beta = %.3f, s[i] = %.3f\n",
	//	      ii,gammaPost[0],gammaPost[1],aSenderEffects[ii]);
    }


}


void CWSBM::drawReceiverEffects(){
    int jj;
    computeBlockSenderEffects();
    double gammaPost[2];
    //   Rprintf("Drawing Receiver Effects:\n");

    for(jj = 0 ; jj < aNodes ; jj++){
	gammaPost[0] = aAdjMatColSums[jj] + (*aPriorReceiver)[jj][0];
	gammaPost[1] = aHours*aBlockSenderEffects[jj] + (*aPriorReceiver)[jj][1];
	aReceiverEffects[jj] = rgamma(gammaPost[0],1/gammaPost[1]);
	//      Rprintf("j = %d, alpha = %.3f, beta = %.3f, r[j] = %.3f\n",
	//	      jj,gammaPost[0],gammaPost[1],aReceiverEffects[jj]);
    }
}

void CWSBM::imputeMissingValues(){

    int rr, cc, index;
    double tie_mean;
    for(rr = 0 ; rr < aNodes ; rr++){ // row loop
	for(cc = 0 ; cc < aNodes ; cc++){ // col loop
	    if(rr != cc){ // ignore diagonal entries
		if(isMissing(aAdjPartial[rr][cc])){
		    tie_mean = GetTieMean(rr,cc);
		    aAdjMat[rr][cc] = (int) rpois(tie_mean);
		}
	    }
	}
    }
}



void CWSBM::normalizeSenderEffects(){

    int ii;

    std::vector<double> norm_vec;
    norm_vec.resize(aBlocks,0.0);
    std::vector<int> block_count;
    block_count.resize(aBlocks,0);

    for(ii = 0 ; ii < aNodes ; ii++){
	norm_vec[(*aBlockMemb)[ii]] += aSenderEffects[ii];
	block_count[(*aBlockMemb)[ii]]++;
    }
    for(int kk = 0 ; kk < aBlocks; kk++){
	norm_vec[kk] /= block_count[kk];
    }
    for(ii = 0 ; ii < aNodes ; ii++){
	aSenderEffects[ii] /= norm_vec[(*aBlockMemb)[ii]];
    }

}


void CWSBM::normalizeReceiverEffects(){

    int ii;

    std::vector<double> norm_vec;
    norm_vec.resize(aBlocks,0.0);
    std::vector<int> block_count;
    block_count.resize(aBlocks,0);

    for(ii = 0 ; ii < aNodes ; ii++){
	norm_vec[(*aBlockMemb)[ii]] += aReceiverEffects[ii];
	block_count[(*aBlockMemb)[ii]]++;
    }
    for(int kk = 0 ; kk < aBlocks; kk++){
	norm_vec[kk] /= block_count[kk];
    }
    for(ii = 0 ; ii < aNodes ; ii++){
	aReceiverEffects[ii] /= norm_vec[(*aBlockMemb)[ii]];
    }

}

void CWSBM::rotate(){
    int ii, jj, kk;
    int assigned[aBlocks], replaced[aBlocks];
    for(ii = 0 ; ii < aBlocks ; ii++){
	assigned[ii] = -1;
	replaced[ii] = -1;
    }
    int PPtmp[aNodes];
    for(ii = 0 ; ii < aNodes ; ii++){
	PPtmp[ii] = (*aBlockMemb)[ii];
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

    //  Updating Membership Vectors
    /*
      for(ii = 0 ; ii < aNodes ; ii++){
      (*aBlockMemb)[ii] = assigned[(*aBlockMemb)[ii]];
      for(kk = 0 ; kk < aBlocks ; kk++){
      if((*aBlockMemb)[ii] == kk){
      aPosteriorMemb[ii*aBlocks + kk] = 1.0;
      }else{
      aPosteriorMemb[ii*aBlocks + kk] = 0.0;
      }
      }
      }
    */

    //  Updating Block Matrix
    saveBlockMatOld();  //  Putting current version of BlockMat in aBlockMatOld
    for(ii = 0 ; ii < aBlocks ; ii++){
	for(jj = 0 ; jj < aBlocks ; jj++){
	    aBlockMat[ii][jj] = aBlockMatOld[assigned[ii]][assigned[jj]];
	    //	 aBlockMatLog[ii][jj] = log(aBlockMat[ii][jj]);
	}
    }
}


/******************  Log-Likelihood Functions  *****************/
//  Consider Converting to a version using:
//         aBlockMatLog and additional aSenderEffectLog, aReceiverEffectLog
double CWSBM::tieLogLike(int yy, int ss, int rr){
    if(yy < 0){
	return 0;
    }else{
	double lambda = GetTieMean(ss,rr);
	//      Rprintf("%d, %d : Mean = %.2f;  LL = %.2f\n",
	//	      ss,rr,lambda,yy *log(lambda) - lambda);
	return yy*log(lambda) - lambda;
    }
}

double CWSBM::GetTieMean(int ss, int rr){

    double mm;

    mm = aBlockMat[(*aBlockMemb)[ss]][(*aBlockMemb)[rr]];
    mm *= aSenderEffects[ss];
    mm *= aReceiverEffects[rr];
    mm *= aHours;

    if(mm < MIN_LOG){
	mm = MIN_LOG;
    }
    return mm;
}


/********************  Writing Functions ********************/

void CWSBM::writeR(int iter){

    aRLogLike[iter] = LogLike();
    writeRBlockMat(iter);
    writeRBlockMemb(iter);

    writeRSenderEffects(iter);
    writeRReceiverEffects(iter);

    writeRPosteriorMemb(iter);
}

void CWSBM::writeDynR(int iter){
    writeRSenderEffects(iter);
    writeRReceiverEffects(iter);
    writeRBlockMat(iter);
}



void CWSBM::writeRBlockMat(int iter){
    int ii, jj;
    int saveIter = iter * aBlocks * aBlocks;
    for(jj = 0 ; jj < aBlocks ; jj++){
	for(ii = 0 ; ii < aBlocks ; ii++){
	    aRBlockMat[saveIter++] = aBlockMat[ii][jj];
	}
    }
}

void CWSBM::writeRBlockMemb(int iter){
    int ii;
    int saveIter = iter * aNodes;
    for(ii = 0 ; ii < aNodes ; ii++){
	aRBlockMemb[saveIter++] = (*aBlockMemb)[ii] + 1;
    }
}

void CWSBM::writeRPosteriorMemb(int iter){
    int ii, kk;
    int saveIter = iter * aBlocks * aNodes;
    for(kk = 0 ; kk < aBlocks; kk++){
	for(ii = 0 ; ii < aNodes ; ii++){
	    aRPosteriorMemb[saveIter++] =  aPosteriorMemb[ii][kk];
	}
    }
}

void CWSBM::writeRSenderEffects(int iter){
    int ii;
    int saveIter = iter*aNodes;
    for(ii = 0 ; ii < aNodes ; ii++){
	aRSenderEffects[saveIter++] = aSenderEffects[ii];
    }
}
void CWSBM::writeRReceiverEffects(int iter){
    int ii;
    int saveIter = iter*aNodes;
    for(ii = 0 ; ii < aNodes ; ii++){
	aRReceiverEffects[saveIter++] = aReceiverEffects[ii];
    }

}


/*****************  HELPER FUNCTIONS  *****************/
void CWSBM::computeBlockTieSums(){

    int kk;
    //  Setting values to zero
    for(kk = 0 ; kk < aBlocks ; kk++){
	aBlockTieSums[kk].assign(aBlocks,0.0);
	aBlockTieCounts[kk].assign(aBlocks,0.0);
    }

    int ss, rr, sb, rb;
    double srEffect;
    //  Summing over each tie
    for(ss = 0 ; ss < aNodes ; ss++){
	sb = (*aBlockMemb)[ss];  // Sending Block
	for(rr = 0 ; rr < aNodes ; rr++){
	    rb = (*aBlockMemb)[rr];  // Receiving Block
	    if(!isMissing(aAdjMat[ss][rr])){
		aBlockTieSums[sb][rb] += aAdjMat[ss][rr];
		srEffect = aSenderEffects[ss]*aReceiverEffects[rr];
		aBlockTieCounts[sb][rb] += srEffect;
	    }
	}
    }
}

void CWSBM::computeRowColSums(){
    int ii, jj;
    aAdjMatRowSums.assign(aNodes,0.0);
    aAdjMatColSums.assign(aNodes,0.0);
    for(ii = 0 ; ii < aNodes ; ii++){
	for(jj = 0 ; jj < aNodes ; jj++){
	    if(!isMissing(aAdjMat[ii][jj])){
		aAdjMatRowSums[ii] += aAdjMat[ii][jj];
		aAdjMatColSums[jj] += aAdjMat[ii][jj];
	    }
	}
    }
}


void CWSBM::computeBlockSenderEffects(){
    int ii,jj;
    double blockEffect;
    aBlockSenderEffects.assign(aNodes,0.0);
    for(ii = 0 ; ii < aNodes; ii++){
	for(jj = 0 ; jj < aNodes; jj++){
	    if(!isMissing(aAdjMat[ii][jj])){
		blockEffect = aBlockMat[(*aBlockMemb)[ii]][(*aBlockMemb)[jj]];
		aBlockSenderEffects[jj] += aSenderEffects[ii] * blockEffect;
	    }
	}
    }

}

//Compute summed block and receiver effects for each sender
void CWSBM::computeBlockReceiverEffects(){
    int ii,jj;
    double blockEffect;
    aBlockReceiverEffects.assign(aNodes,0.0);
    for(ii = 0 ; ii < aNodes ; ii++){
	for(jj = 0 ; jj < aNodes ; jj++){
	    if(!isMissing(aAdjMat[ii][jj])){
		blockEffect = aBlockMat[(*aBlockMemb)[ii]][(*aBlockMemb)[jj]];
		aBlockReceiverEffects[ii] += aReceiverEffects[jj] * blockEffect;
	    }
	}
    }
}


/*******************  BlockMat UPDATING FUNCTIONS  ******************/

void CWSBM::GetBlockMat(double *rBlockMat){
    int ii,jj;
    int saveIter = 0;
    for(jj = 0 ; jj < aBlocks ; jj++){
	for(ii = 0 ; ii < aBlocks ; ii++){
	    rBlockMat[saveIter++] = aBlockMat[ii][jj];
	    //BlockMat_out[ii + jj*aBlocks] = aBlockMat[ii][jj];
	}
    }
}


void CWSBM::saveBlockMatOld(){
    int ii, jj;
    for(ii = 0 ; ii < aBlocks ; ii++){
	for(jj = 0 ; jj < aBlocks ; jj++){
	    aBlockMatOld[ii][jj] = aBlockMat[ii][jj];
	}
    }
}


/********************  Printing Functions ********************/
void CWSBM::printAdjacencyMatrix(){

    int ii, jj;
    Rprintf("Adjacency Matrix:\n");
    for(ii = 0 ; ii < aNodes; ii++){
	for(jj = 0; jj < aNodes; jj++){
	    Rprintf("%d ",aAdjMat[ii][jj]);
	}
	Rprintf("\n");
    }

}

void CWSBM::printBlockMat(){
    int ii,jj;
    Rprintf("Block Matrix:\n");
    for(ii = 0 ; ii < aBlocks; ii++){
	for(jj = 0; jj < aBlocks; jj++){
	    Rprintf("%.2f ",aBlockMat[ii][jj]);
	}
	Rprintf("\n");
    }
}

void CWSBM::printPosteriorMemb(){
    int ii, kk;
    Rprintf("Block Memberships:\n");
    for(ii = 0 ; ii < aNodes ; ii++){
	for(kk = 0 ; kk < aBlocks ; kk++){
	    Rprintf("%.2f ",aPosteriorMemb[ii][kk]);
	}
	Rprintf("\n");
    }
}

void CWSBM::printBlockMemb(){
    int ii;
    Rprintf("Block Memberships:\n");
    for(ii = 0; ii < aNodes; ii++){
	if(ii % 10 == 0){
	    Rprintf("\n");
	}
	Rprintf("%d ",(*aBlockMemb)[ii]);
    }
    Rprintf("\n\n");
}

void CWSBM::printPriorBlockMemb(){
    int ii;
    Rprintf("Block Membership Prior:\n");
    for(ii = 0 ; ii < aBlocks ; ii++){
	Rprintf("%.3f ",aPriorBlockMemb[ii]);
    }
    Rprintf("\n");
}

void CWSBM::printSenderEffects(){
    int ii;
    Rprintf("Sender Effects:");
    for(ii = 0 ; ii < aNodes ; ii++){
	if(ii % 10 == 0){
	    Rprintf("\n");
	}
	Rprintf("%.2f ",aSenderEffects[ii]);
    }
    Rprintf("\n");
}

void CWSBM::printReceiverEffects(){
    int ii;
    Rprintf("Receiver Effects:");
    for(ii = 0 ; ii < aNodes ; ii++){
	if(ii % 10 == 0){
	    Rprintf("\n");
	}
	Rprintf("%.2f ",aReceiverEffects[ii]);
    }
    Rprintf("\n");
}

