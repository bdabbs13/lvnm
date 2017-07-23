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
#include "dynamic-sbm.h"
#include "helper.h"
//using namespace std;
using std::ifstream;
//using std::ofstream;

//gsl_rng *rng;
//ofstream myfile;





void dynSBMMCMC(CDynSBM *myDynSBM, int start, int total, int burnIn, int thin,
		int shift_size, int extend_max, double qq, int verbose){

    int ii;
    for(ii = 1 ; ii < total ; ii++){
	myDynSBM->step();
	myDynSBM->adapt();
	
	if(ii >= burnIn){
	    if(((ii - burnIn) % thin == 0)){
		//  Save the result every thin iterations
		int save_iter = (ii - burnIn)/thin;
		myDynSBM->Update(save_iter);
	    }
	}
    }

}



/*****************************************************************
 *********  FUNCTION DEFINITIONS FOR CDynSBM OBJECT CLASS  ***********
 *****************************************************************/




/***********  WSBM CLASS CONSTRUCTOR AND DESTRUCTOR  **********/

//  Constructor Function for WSBM object

CDynSBM::CDynSBM (int rNodes, int rBlocks, int rTimes, int rTimeClasses,
		  int total, int *rTimeMap, double *rHours,
		  int mImpute) : missingVal(-1), covCount(0),
				 mhEpsilon(0.000001), mhSD(1), mhStart(100)
{
    int ii, jj, tt, ll, kk;

    //  Loading Dimensions
    aNodes = rNodes;
    aBlocks = rBlocks;
    aTimes = rTimes;
    aClasses = rTimeClasses;
    aTotal = total;

    //Loading Time Information
    aHours.resize(aClasses);
    for(tt = 0 ; tt < aClasses ; tt++){
	aHours[tt] = rHours[tt];
    }

    aTimeMap.resize(aTimes);
    for(tt = 0 ; tt < aTimes; tt++){
	aTimeMap[tt] = rTimeMap[tt] - 1; // Converting to zero indexing
    }

    //  Initializing WSBM Objects
    aWsbmList.reserve(aTimes);
    for(tt = 0 ; tt < aTimes ; tt++){
	aWsbmList.push_back(CWSBM(aNodes,aBlocks,aHours[aTimeMap[tt]],mImpute));
    }


    //  Initializing Prior Containers
    aPriorSender.resize(aClasses,0);
    aPriorReceiver.resize(aClasses,0);
    for(tt = 0 ; tt < aClasses; tt++){
	aPriorSender[tt] = new std::vector<double *>(aNodes,0);
	aPriorReceiver[tt] = new std::vector<double *>(aNodes,0);
	for(ii = 0 ; ii < aNodes ; ii++){
	    (*(aPriorSender[tt]))[ii] = new double[2](); // () sets value to 0
	    (*(aPriorReceiver[tt]))[ii] = new double[2](); // () sets value to 0
	}
    }

    aPriorBlockMat.resize(aClasses,0);
    for(tt = 0 ; tt < aClasses; tt++){
	aPriorBlockMat[tt] = new std::vector<std::vector<double *> >(aBlocks);
	for(ll = 0 ; ll < aBlocks ; ll++){
	    (*(aPriorBlockMat[tt]))[ll].resize(aBlocks,0);
	    for(kk = 0; kk < aBlocks; kk++){
		(*(aPriorBlockMat[tt]))[ll][kk] = new double[2](); // () sets to 0
	    }
	}
    }

    //  Initializing Adaptive Sampling Containers
    aPriorSenderCov.resize(aClasses,0);
    aPriorReceiverCov.resize(aClasses,0);
    for(tt = 0 ; tt < aClasses; tt++){
	aPriorSenderCov[tt] = new std::vector<double *>(aNodes,0);
	aPriorReceiverCov[tt] = new std::vector<double *>(aNodes,0);
	for(ii = 0 ; ii < aNodes ; ii++){
	    (*(aPriorSenderCov[tt]))[ii] = new double[5](); // () sets value to 0
	    (*(aPriorReceiverCov[tt]))[ii] = new double[5](); // () sets value to 0
	}
    }

    aPriorBlockMatCov.resize(aClasses,0);
    for(tt = 0 ; tt < aClasses; tt++){
	aPriorBlockMatCov[tt] = new std::vector<std::vector<double *> >(aBlocks);
	for(ll = 0 ; ll < aBlocks ; ll++){
	    (*(aPriorBlockMatCov[tt]))[ll].resize(aBlocks,0);
	    for(kk = 0; kk < aBlocks; kk++){
		(*(aPriorBlockMatCov[tt]))[ll][kk] = new double[5](); // () sets to 0
	    }
	}
    }



    aBlockMemb.resize(aNodes,-1);

    aPosteriorMemb.resize(aNodes);
    for(ii = 0 ; ii < aNodes ; ii++){
	aPosteriorMemb[ii].resize(aBlocks,0);
    }

    aPriorBlockMemb.resize(aBlocks,0.0);


}

CDynSBM::~CDynSBM (){
    int tt, ii;
    for(tt = 0 ; tt < aClasses ; tt++){
	for(ii = 0 ; ii < aNodes; ii++){
	    delete (*(aPriorSender[tt]))[ii];
	    delete (*(aPriorReceiver[tt]))[ii];
	}
    }
    int ll, kk;
    for(tt = 0 ; tt < aClasses ; tt++){
	for(ll = 0 ; ll < aBlocks; ll++){
	    for(kk = 0 ; kk < aBlocks; kk++){
		delete (*(aPriorBlockMat[tt]))[ll][kk];
	    }
	}
    }
}


/**********  LOADING FUNCTIONS  *********/

/**  Public Function for Loading from R  **/
void CDynSBM::RLoadDynSBM(int *AdjMat,
			  double *rHyperSender, double *rHyperReceiver,
			  double *rHyperBlockMat,
			  double *rSenderEffects, double *rReceiverEffects,
			  double *rBlockEffects, int *rBlockMemb,
			  double *rPosteriorMemb, int update_mmb,
			  double *rPriorSender, double *rPriorReceiver,
			  double *rPriorBlockMat, double *rPriorBlockMemb,
			  double *rLogLik){

    //  Loading Adjacency Matrices
    LoadAdjacencyMatrices(AdjMat);

    //  Loading HyperPriors
    LoadHyperPriors(rHyperSender, rHyperReceiver, rHyperBlockMat);

    //  Loading Parameters into WSBM Objects
    LoadParameters(rSenderEffects,rReceiverEffects,rBlockEffects,
		   rBlockMemb,rPosteriorMemb,update_mmb);

    //  Passing pointers to priors to WSBM objects
    LoadPriors(rPriorSender,rPriorReceiver,rPriorBlockMat,
	       rPriorBlockMemb);
    LoadLogLike(rLogLik);
    PassReferences();
}



/** Loading Adjacency Matrices into WSBM Objects **/
void CDynSBM::LoadAdjacencyMatrices(int *AdjMat){
    int tt;
    for(tt = 0 ; tt < aTimes ; tt++){
	aWsbmList[tt].LoadAdjacencyMatrix(&AdjMat[tt*aNodes*aNodes]);
    }

}

/** Loading Hyperpriors from R **/
void CDynSBM::LoadHyperPriors(double *rHyperSender, double *rHyperReceiver,
			      double *rHyperBlockMat){
    int ii;
    for(ii = 0 ; ii < 4 ; ii++){
	aHyperSender[ii] = rHyperSender[ii];
	aHyperReceiver[ii] = rHyperReceiver[ii];
	aHyperBlockMat[ii] = rHyperBlockMat[ii];
    }
}

/** Loading Initial Parameters from R and Passing to WSBM Objects **/
void CDynSBM::LoadParameters(double *rSenderEffects, double *rReceiverEffects,
			     double *rBlockEffects, int *rBlockMemb,
			     double *rPosteriorMemb, int r_update_mmb){

    int tt;
    for(tt = 0 ; tt < aTimes; tt++){
	aWsbmList[tt].RLoadWSBM(&rBlockEffects[tt*aTotal*aBlocks*aBlocks],
				&rSenderEffects[tt*aTotal*aNodes],
				&rReceiverEffects[tt*aTotal*aNodes]);
    }

    aRBlockMemb = rBlockMemb;
    aRPosteriorMemb = rPosteriorMemb;
    int ii;
    for(ii = 0 ; ii < aNodes ; ii++){
	aBlockMemb[ii] = rBlockMemb[ii] - 1;
	aPosteriorMemb[ii][aBlockMemb[ii]] = 1.0;
    }

    update_mmb = (r_update_mmb > 0);

}

/** Loading Initial Priors from R and Keeping Pointers **/
void CDynSBM::LoadPriors(double *rPriorSender, double *rPriorReceiver,
			 double *rPriorBlockMat, double *rPriorBlockMemb){
    aRPriorSender = rPriorSender;
    LoadPriorSender();

    aRPriorReceiver = rPriorReceiver;
    LoadPriorReceiver();

    aRPriorBlockMat = rPriorBlockMat;
    LoadPriorBlockMat();

    //  This value never changes
    int ll;
    for(ll = 0 ; ll < aBlocks ; ll++){
	aPriorBlockMemb[ll] = rPriorBlockMemb[ll];
    }

}

void CDynSBM::LoadLogLike(double *rLogLik){
    aRLogLike = rLogLik;
}



void CDynSBM::LoadPriorSender(){
    int tt,kk,ii;
    for(tt = 0 ; tt < aClasses ; tt++){
	// Indexing to correct section of array
	int loadIter = tt*(aTotal*2*aNodes);
	// Reading values in column-major order
	for(kk = 0 ; kk < 2 ; kk++){
	    for(ii = 0 ; ii < aNodes ; ii++){
		(*(aPriorSender[tt]))[ii][kk] = aRPriorSender[loadIter++];
	    }
	}
    }
}

void CDynSBM::LoadPriorReceiver(){
    int tt,kk,ii;
    for(tt = 0 ; tt < aClasses ; tt++){
	// Indexing to correct section of array
	int loadIter = tt*(aTotal*2*aNodes);
	// Reading values in column-major order
	for(kk = 0 ; kk < 2 ; kk++){
	    for(ii = 0 ; ii < aNodes ; ii++){
		(*(aPriorReceiver[tt]))[ii][kk] = aRPriorReceiver[loadIter++];
	    }
	}
    }
}

void CDynSBM::LoadPriorBlockMat(){
    int tt,kk,ll,ii;
    for(tt = 0 ; tt < aClasses ; tt++){
	// Indexing to correct section of array
	int loadIter = tt*(aTotal*2*aBlocks*aBlocks);
	// Reading values in column-major order
	for(ii = 0 ; ii < 2 ; ii++){
	    for(kk = 0 ; kk < aBlocks ; kk++){
		for(ll = 0 ; ll < aBlocks ; ll++){
		    (*(aPriorBlockMat[tt]))[ll][kk][ii] =aRPriorBlockMat[loadIter++];
		}
	    }
	}
    }
}


/** Passing References to R Storage Locations to WSBM Objects **/
void CDynSBM::PassReferences(){
    int tt;
    for(tt = 0 ; tt < aTimes ; tt++){
	aWsbmList[tt].LoadReferences(aPriorSender[aTimeMap[tt]],
				     aPriorReceiver[aTimeMap[tt]],
				     aPriorBlockMat[aTimeMap[tt]],
				     &aBlockMemb);
    }

}

/**********  MCMC SAMPLING FUNCTIONS  *********/

void CDynSBM::step(){
    //   Rprintf(" Priors");
    DrawPriors();
    //   Rprintf(" Params");
    DrawParameters();
}


void CDynSBM::DrawPriors(){
    // Drawing Sender Priors
    DrawPriorSender();

    // Drawing Receiver Priors
    DrawPriorReceiver();

    // Drawing BlockMat Priors
    DrawPriorBlockMat();

}

void CDynSBM::DrawParameters(){
    int tt;
    //   Rprintf("tt: ");
    for(tt = 0 ; tt < aTimes ; tt++){
	//      Rprintf("%d ",tt);
	aWsbmList[tt].partialStep();
    }
    if(update_mmb){
	DrawBlockMemb();
    }
}

void CDynSBM::DrawBlockMemb(){
    int ii,kk;
    double pMax;
    int ind[aBlocks];

    //   Rprintf("Log-Likelihood = %.2f\n",LogLike());

    for(ii = 0 ; ii < aNodes ; ii++){
	for(kk = 0 ; kk < aBlocks ; kk++){
	    aBlockMemb[ii] = kk;

	    // Calculating Log Posterior Probability
	    //	 Rprintf("ii = %d, kk = %d\n",ii,kk);
	    aPosteriorMemb[ii][kk] = nodeLogLike(ii) + log(aPriorBlockMemb[kk]);
	    //	 Rprintf("p.%d.%d = %.4f\n",ii,kk,aPosteriorMemb[ii][kk]);
	    if(kk == 0){
		pMax = aPosteriorMemb[ii][kk];
	    }else{
		if(aPosteriorMemb[ii][kk] > pMax){
		    pMax = aPosteriorMemb[ii][kk];
		}
	    }
	}

	// Exponentiating Probabilities
	for(kk = 0 ; kk < aBlocks ; kk++){
	    aPosteriorMemb[ii][kk] = exp(aPosteriorMemb[ii][kk] - pMax);
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
	aBlockMemb[ii] = kk;
    }

}

double CDynSBM::nodeLogLike(int node){
    int tt;
    double sum = 0.0;
    //   Rprintf("tt: ");
    for(tt = 0 ; tt < aTimes ; tt++){
	//      Rprintf("%d ",tt);
	sum += aWsbmList[tt].nodeLogLike(node);
    }
    //   Rprintf("\n");
    return sum;
}

double CDynSBM::LogLike(){
    int tt;
    double sum = 0.0;
    for(tt = 0 ; tt < aTimes ; tt++){
	sum += aWsbmList[tt].LogLike();
    }
    return sum;
}


void CDynSBM::DrawPriorSender(){
    int tt, ii;
    std::vector<double> pp,qq,rr,ss;

    pp.resize(aClasses);
    qq.resize(aClasses);
    rr.resize(aClasses);
    ss.resize(aClasses);

    double prop_sd = 1.0;

    for(ii = 0 ; ii < aNodes ; ii++){
	//  Setting Prior Values
	pp.assign(aClasses,aHyperSender[0]);
	qq.assign(aClasses,aHyperSender[1]);
	rr.assign(aClasses,aHyperSender[2]);
	ss.assign(aClasses,aHyperSender[3]);

	//  Calculating Posterior Parameters
	for(tt = 0 ; tt < aTimes ; tt++){
	    int t_class = aTimeMap[tt];
	    pp[t_class] *= aWsbmList[tt].GetSenderEffect(ii);
	    qq[t_class] += aWsbmList[tt].GetSenderEffect(ii);
	    rr[t_class]++;
	    ss[t_class]++;
	}
	//      Rprintf("p: %.2f  q: %.2f  r: %.2f  s: %.2f\n",
	//	      pp[0],qq[0],rr[0],ss[0]);

	//  Making Posterior Draws
	for(tt = 0 ; tt < aClasses; tt++){
	    // rGammaPriorStep((*(aPriorSender[tt]))[ii],
	    // 		 pp[tt],qq[tt],rr[tt],ss[tt],prop_sd);
	    if(covCount < mhStart){
		GammaPriorMHDraw((*(aPriorSender[tt]))[ii][0],
				 (*(aPriorSender[tt]))[ii][1],
				 pp[tt],qq[tt],rr[tt],ss[tt],
				 mhSD,mhSD,0.0);
	    }else{
		GammaPriorMHDraw((*(aPriorSender[tt]))[ii], covCount,
				 pp[tt],qq[tt],rr[tt],ss[tt],
				 (*(aPriorSenderCov[tt]))[ii],
				 mhEpsilon);
	    }
	}
    }
}

void CDynSBM::DrawPriorReceiver(){
    int tt, ii;
    std::vector<double> pp,qq,rr,ss;

    pp.resize(aClasses);
    qq.resize(aClasses);
    rr.resize(aClasses);
    ss.resize(aClasses);

    double prop_sd = 1.0;

    for(ii = 0 ; ii < aNodes ; ii++){
	//  Setting Prior Values
	pp.assign(aClasses,aHyperReceiver[0]);
	qq.assign(aClasses,aHyperReceiver[1]);
	rr.assign(aClasses,aHyperReceiver[2]);
	ss.assign(aClasses,aHyperReceiver[3]);

	//  Calculating Posterior Parameters
	for(tt = 0 ; tt < aTimes ; tt++){
	    int t_class = aTimeMap[tt];
	    pp[t_class] *= aWsbmList[tt].GetReceiverEffect(ii);
	    qq[t_class] += aWsbmList[tt].GetReceiverEffect(ii);
	    rr[t_class]++;
	    ss[t_class]++;
	}

	//  Making Posterior Draws
	for(tt = 0 ; tt < aClasses; tt++){
	    // rGammaPriorStep((*(aPriorReceiver[tt]))[ii],
	    // 		 pp[tt],qq[tt],rr[tt],ss[tt],prop_sd);
	    if(covCount < mhStart){
		GammaPriorMHDraw((*(aPriorReceiver[tt]))[ii][0],
				 (*(aPriorReceiver[tt]))[ii][1],
				 pp[tt],qq[tt],rr[tt],ss[tt],
				 mhSD,mhSD,0.0);
	    }else{
		GammaPriorMHDraw((*(aPriorReceiver[tt]))[ii], covCount,
				 pp[tt],qq[tt],rr[tt],ss[tt],
				 (*(aPriorReceiverCov[tt]))[ii],
				 mhEpsilon);
	    }
	}
    }
}

void CDynSBM::DrawPriorBlockMat(){
    int tt, ll, kk;
    std::vector<double> pp,qq,rr,ss;

    pp.resize(aClasses);
    qq.resize(aClasses);
    rr.resize(aClasses);
    ss.resize(aClasses);

    double prop_sd = 1.0;

    for(kk = 0 ; kk < aBlocks ; kk++){
	for(ll = 0 ; ll < aBlocks ; ll++){
	    //  Setting Prior Values
	    pp.assign(aClasses,aHyperBlockMat[0]);
	    qq.assign(aClasses,aHyperBlockMat[1]);
	    rr.assign(aClasses,aHyperBlockMat[2]);
	    ss.assign(aClasses,aHyperBlockMat[3]);

	    //  Calculating Posterior Parameters
	    for(tt = 0 ; tt < aTimes ; tt++){
		int t_class = aTimeMap[tt];
		pp[t_class] *= aWsbmList[tt].GetBlockMatEffect(ll,kk);
		qq[t_class] += aWsbmList[tt].GetBlockMatEffect(ll,kk);
		rr[t_class]++;
		ss[t_class]++;
	    }

	    //  Making Posterior Draws
	    for(tt = 0 ; tt < aClasses; tt++){
		// rGammaPriorStep((*(aPriorBlockMat[tt]))[ll][kk],
		// 		    pp[tt],qq[tt],rr[tt],ss[tt],prop_sd);
		if(covCount < mhStart){
		    GammaPriorMHDraw((*(aPriorBlockMat[tt]))[ll][kk][0],
				     (*(aPriorBlockMat[tt]))[ll][kk][1],
				     pp[tt],qq[tt],rr[tt],ss[tt],
				     mhSD,mhSD,0.0);
		}else{
		    GammaPriorMHDraw((*(aPriorBlockMat[tt]))[ll][kk], covCount,
				     pp[tt],qq[tt],rr[tt],ss[tt],
				     (*(aPriorBlockMatCov[tt]))[ll][kk],
				     mhEpsilon);
		}
	    }
	}
    }
}



/**********  UPDATING FUNCTIONS  *********/

void CDynSBM::Update(int iter){

    UpdatePriorSender(iter);
    UpdatePriorReceiver(iter);
    UpdatePriorBlockMat(iter);

    int tt;
    for(tt = 0 ; tt < aTimes ; tt++){
	aWsbmList[tt].partialUpdate(iter);
    }
    UpdateBlockMemb(iter);
    UpdatePosteriorMemb(iter);
    UpdateLogLike(iter);
}

void CDynSBM::UpdatePriorSender(int iter){
    int tt,kk,ii;
    for(tt = 0 ; tt < aClasses ; tt++){
	// Indexing to correct section of array
	int saveIter = tt*(aTotal*2*aNodes) + iter * 2 * aNodes;
	// Writing values in column-major order
	for(kk = 0 ; kk < 2 ; kk++){
	    for(ii = 0 ; ii < aNodes ; ii++){
		aRPriorSender[saveIter++] = (*(aPriorSender[tt]))[ii][kk];
	    }
	}
    }
}

void CDynSBM::UpdatePriorReceiver(int iter){
    int tt,kk,ii;
    for(tt = 0 ; tt < aClasses ; tt++){
	// Indexing to correct section of array
	int saveIter = tt*(aTotal*2*aNodes) + iter * 2*aNodes;
	// Writing values in column-major order
	for(kk = 0 ; kk < 2 ; kk++){
	    for(ii = 0 ; ii < aNodes ; ii++){
		aRPriorReceiver[saveIter++] = (*(aPriorReceiver[tt]))[ii][kk];
	    }
	}
    }
}

void CDynSBM::UpdatePriorBlockMat(int iter){
    int tt,kk,ll,ii;
    for(tt = 0 ; tt < aClasses ; tt++){
	// Indexing to correct section of array
	int saveIter = tt*(aTotal*2*aBlocks*aBlocks) + iter * 2*aBlocks*aBlocks;
	// Reading values in column-major order
	for(ii = 0 ; ii < 2 ; ii++){
	    for(kk = 0 ; kk < aBlocks ; kk++){
		for(ll = 0 ; ll < aBlocks ; ll++){
		    aRPriorBlockMat[saveIter++] = (*(aPriorBlockMat[tt]))[ll][kk][ii];
		}
	    }
	}
    }
}


void CDynSBM::UpdateBlockMemb(int iter){
    int ii;
    int saveIter = iter * aNodes;
    for(ii = 0 ; ii < aNodes ; ii++){
	aRBlockMemb[saveIter++] = aBlockMemb[ii] + 1;
    }
}


void CDynSBM::UpdatePosteriorMemb(int iter){
    int ii, kk;
    int saveIter = iter * aBlocks * aNodes;
    for(kk = 0 ; kk < aBlocks; kk++){
	for(ii = 0 ; ii < aNodes ; ii++){
	    aRPosteriorMemb[saveIter++] =  aPosteriorMemb[ii][kk];
	}
    }
}

void CDynSBM::UpdateLogLike(int iter){
    int tt;
    int saveIter = iter * aTimes;

    for(tt = 0 ; tt < aTimes; tt++){
	aRLogLike[saveIter++] = aWsbmList[tt].LogLike();
    }

}

void CDynSBM::adapt(){

    covCount++;
    for(int tt = 0 ; tt < aClasses ; tt++){
	for(int ii = 0 ; ii < aNodes ; ii++){
	    updateSingleSSE((*(aPriorSender[tt]))[ii], covCount,
			    (*(aPriorSenderCov[tt]))[ii]);
	    updateSingleSSE((*(aPriorReceiver[tt]))[ii], covCount,
			    (*(aPriorReceiverCov[tt]))[ii]);
	}

	for(int ll = 0 ; ll < aBlocks ; ll++){
	    for(int kk = 0 ; kk < aBlocks ; kk++){
		updateSingleSSE((*(aPriorBlockMat[tt]))[ll][kk], covCount,
				(*(aPriorBlockMatCov[tt]))[ll][kk]);
	    }
	}
    }

}

void CDynSBM::PrintCovariance(int tt, int ii){
    Rprintf("Covariance Count: %d\n",covCount);
    Rprintf("Sender Covariance:\n");
    Rprintf("a_sse = %.4f, b_sse = %.4f, ab_csse = %.4f, a_mean = %.4f, b_mean = %.4f",
	    (*(aPriorSenderCov[tt]))[ii][0],
	    (*(aPriorSenderCov[tt]))[ii][1],
	    (*(aPriorSenderCov[tt]))[ii][2],
	    (*(aPriorSenderCov[tt]))[ii][3],
	    (*(aPriorSenderCov[tt]))[ii][4]);

    Rprintf("\n\n");
    Rprintf("Receiver Covariance:\n");
    Rprintf("a_sse = %.4f, b_sse = %.4f, ab_csse = %.4f, a_mean = %.4f, b_mean = %.4f",
	    (*(aPriorReceiverCov[tt]))[ii][0],
	    (*(aPriorReceiverCov[tt]))[ii][1],
	    (*(aPriorReceiverCov[tt]))[ii][2],
	    (*(aPriorReceiverCov[tt]))[ii][3],
	    (*(aPriorReceiverCov[tt]))[ii][4]);

}



/**********  PRINTING FUNCTIONS  *********/

void CDynSBM::printPriors(){
    printPriorSender();
    printPriorReceiver();
    printPriorBlockMat();
}

void CDynSBM::printPriorSender(){
    int tt, ii;
    Rprintf("Sender Effect Priors:\n");
    for(tt = 0; tt < aClasses ; tt++){
	Rprintf("Class %d:\n",tt);
	for(ii = 0 ; ii < aNodes ; ii++){
	    Rprintf("%.2f, %.2f\n",
		    (*(aPriorSender[tt]))[ii][0],
		    (*(aPriorSender[tt]))[ii][1]);
	}
    }
}

void CDynSBM::printPriorReceiver(){
    int tt, ii;
    Rprintf("Receiver Effect Priors:\n");
    for(tt = 0; tt < aClasses ; tt++){
	Rprintf("Class %d:\n",tt);
	for(ii = 0 ; ii < aNodes ; ii++){
	    Rprintf("%.2f, %.2f\n",
		    (*(aPriorReceiver[tt]))[ii][0],
		    (*(aPriorReceiver[tt]))[ii][1]);
	}
    }
}

void CDynSBM::printPriorBlockMat(){
    int tt, ll, kk;
    Rprintf("BlockMat Effect Priors:\n");
    for(tt = 0; tt < aClasses ; tt++){
	Rprintf("Class %d:\n",tt);
	for(kk = 0 ; kk < aBlocks ; kk++){
	    for(ll = 0 ; ll < aBlocks ; ll++){
		Rprintf("%.2f, %.2f\n",
			(*(aPriorBlockMat[tt]))[ll][kk][0],
			(*(aPriorBlockMat[tt]))[ll][kk][1]);
	    }
	}
    }
}

void CDynSBM::printAllWSBM(bool printNetworks){
    int tt;
    for(tt = 0 ; tt < aTimes ; tt++){
	Rprintf("\nNetwork Number: %d\n\n",tt);
	aWsbmList[tt].print(printNetworks);
    }
}

void CDynSBM::printBlockMemb(){
    int ii;
    Rprintf("Block Memberships:\n");
    for(ii = 0; ii < aNodes; ii++){
	if(ii % 10 == 0){
	    Rprintf("\n");
	}
	Rprintf("%d ",aBlockMemb[ii]);
    }
    Rprintf("\n\n");
}
