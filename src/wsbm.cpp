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


/*****************************************************************
 *********  FUNCTION DEFINITIONS FOR WSBM OBJECT CLASS  ***********
 *****************************************************************/


/***********  WSBM CLASS CONSTRUCTOR AND DESTRUCTOR  **********/

//  Constructor Function for WSBM object

CWSBM::CWSBM (int rNodes, int rBlocks, int *adjMat,
	      double *rPriorSender, double *rPriorReceiver,
	      double *rPriorBlockMat, double * rPriorBlockMemb,
	      double rHours, int mImpute) : missingVal(-1){
   int ii, jj, kk;
   //   missingVal = -1;

   // Loading Constants
   aNodes = rNodes; aBlocks = rBlocks;

   // Loading Adjacency Matrix
   aImputeFlag = (mImpute == 1);

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
	    aAdjPartial[ii][jj] = adjMat[ii*aNodes + jj];
	 }
      }
   }

   //  Initializing aBlockMat Matrix
   double holder[2];
   aPriorBlockMat = new std::vector<std::vector<double *> >(aBlocks);
   for(kk = 0 ; kk < aBlocks ; kk++){
      (*aPriorBlockMat)[kk].resize(aBlocks,rPriorBlockMat);
   }
   //   std::copy(rPriorBlockMat,rPriorBlockMat+2,aPriorBlockMat);

   aBlockMat.resize(aBlocks);
   aBlockMatLog.resize(aBlocks);
   aBlockMatOld.resize(aBlocks);
   for(ii = 0; ii < aBlocks ; ii++){
      aBlockMat[ii].resize(aBlocks);
      aBlockMatLog[ii].resize(aBlocks);
      aBlockMatOld[ii].resize(aBlocks,0.0);
      for(jj = 0 ; jj < aBlocks ; jj++){
	 rdirichlet(2,(*aPriorBlockMat)[ii][jj],holder);
	 aBlockMat[ii][jj] = holder[0];
	 aBlockMatLog[ii][jj] = log(holder[0]);
      }
   }

   is_BlockMat_logged = false;

   //  Resizing vectors to hold block counts and block tie sums
   aBlockTieCounts.resize(aBlocks);
   aBlockTieSums.resize(aBlocks);
   for(ii = 0 ; ii < aBlocks ; ii++){
      aBlockTieSums[ii].resize(aBlocks,0.0);
      aBlockTieCounts[ii].resize(aBlocks,0.0);
   }

   aPosteriorMemb.resize(aNodes);
   aPosteriorMembOld.resize(aNodes);

   aBlockMemb = new std::vector<int>;
   (*aBlockMemb).resize(aNodes,0);
   for(ii = 0 ; ii < aNodes ; ii++){
      aPosteriorMemb[ii].resize(aBlocks);
      aPosteriorMembOld[ii].resize(aBlocks);
   }

   // Initializing aPosteriorMemb Matrix
   //   aPriorBlockMemb.resize(aBlocks,0.0);
   //   for(ii = 0 ; ii < aBlocks ; ii++){
   //      aPriorBlockMemb[ii] = rPriorBlockMemb[ii];
   //   }

   aPriorBlockMemb.insert(aPriorBlockMemb.end(),
   			  &rPriorBlockMemb[0],&rPriorBlockMemb[aBlocks]);

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
      /*
	Rprintf("ii = %d\n",ii);
	RprintDoubleMat(1,aBlocks,etaInit);
	RprintIntMat(1,aBlocks,ind);
      */
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

   //  Reading in Sender and Receiver Effect Priors
   aPriorSender = new std::vector<double *>(aNodes,rPriorSender);
   aPriorReceiver = new std::vector<double *>(aNodes,rPriorReceiver);
   aPriorOwner = true;
   //   std::copy(rPriorSender,rPriorSender+2,aPriorSender);
   //   std::copy(rPriorReceiver,rPriorReceiver+2,aPriorReceiver);

   //  Initializing Sender and Receiver Effects
   aSenderEffects.resize(aNodes,1.0);
   aReceiverEffects.resize(aNodes,1.0);
   /*
   for(ii = 0 ; ii < aNodes ; ii++){
      aSenderEffects[ii] = rgamma(aPriorSender[0],aPriorSender[1]);
      aReceiverEffects[ii] = rgamma(aPriorReceiver[0],aPriorReceiver[1]);
   }
   */
   aAdjMatRowSums.resize(aNodes,0.0);
   aAdjMatColSums.resize(aNodes,0.0);
   computeRowColSums();

   aBlockSenderEffects.resize(aNodes,0.0);
   aBlockReceiverEffects.resize(aNodes,0.0);

   aHours = rHours;

   delete[] etaInit;
   delete[] ind;
}



//  Just Initializing Storage Containers
CWSBM::CWSBM (int rNodes, int rBlocks, double rHours,
	      int mImpute) : missingVal(-1){
   int ii, jj, kk;

   //   Loading Constants
   aNodes = rNodes; aBlocks = rBlocks;

   // Loading Adjacency Matrix
   aImputeFlag = (mImpute == 1);

   aAdjMat.resize(aNodes);
   for(ii = 0; ii < aNodes ; ii++){
      aAdjMat[ii].resize(aNodes);
   }

   if(aImputeFlag){
      aAdjPartial.resize(aNodes);
      for(ii = 0; ii < aNodes; ii++){
   	 aAdjPartial[ii].resize(aNodes);
      }
   }

   //  Initializing aBlockMat Matrix

   aBlockMat.resize(aBlocks);
   aBlockMatLog.resize(aBlocks);
   aBlockMatOld.resize(aBlocks);
   for(ii = 0; ii < aBlocks ; ii++){
      aBlockMat[ii].resize(aBlocks,0.0);
      aBlockMatLog[ii].resize(aBlocks,0.0);
      aBlockMatOld[ii].resize(aBlocks,0.0);
   }

   is_BlockMat_logged = false;

   //  Resizing vectors to hold block counts and block tie sums
   aBlockTieCounts.resize(aBlocks);
   aBlockTieSums.resize(aBlocks);
   for(ii = 0 ; ii < aBlocks ; ii++){
      aBlockTieSums[ii].resize(aBlocks,0.0);
      aBlockTieCounts[ii].resize(aBlocks,0.0);
   }

   ////  These values are handled by DynSBM
   // aPosteriorMemb.resize(aNodes);
   // aPosteriorMembOld.resize(aNodes);


   // //   aBlockMemb.resize(aNodes,0);
   // for(ii = 0 ; ii < aNodes ; ii++){
   //    aPosteriorMemb[ii].resize(aBlocks,0.0);
   //    aPosteriorMembOld[ii].resize(aBlocks,0.0);
   // }

   //  Initializing Sender and Receiver Effects
   aSenderEffects.resize(aNodes,1.0);
   aReceiverEffects.resize(aNodes,1.0);

   aAdjMatRowSums.resize(aNodes,0.0);
   aAdjMatColSums.resize(aNodes,0.0);

   aBlockSenderEffects.resize(aNodes,0.0);
   aBlockReceiverEffects.resize(aNodes,0.0);

   aPriorBlockMemb.resize(aBlocks,1.0 / aBlocks);
   aHours = rHours;
   aPriorOwner = false;
}

// void CWSBM::LoadPriors(double *PriorSender, double *PriorReceiver,
// 		       double *PriorBlockMat, double *PriorBlockMemb){

//    // Loading Sender/Receiver Priors
//    std::copy(PriorSender,PriorSender+2,aPriorSender);
//    std::copy(PriorReceiver,PriorReceiver+2,aPriorReceiver);

//    // Loading Block Matrix Priors
//    std::copy(PriorBlockMat,PriorBlockMat+2,aPriorBlockMat);

//    Loading Block Membership Priors
//    aPriorBlockMemb.insert(aPriorBlockMemb.end(),
//    			  &PriorBlockMemb[0],&PriorBlockMemb[aBlocks]);
// }

void CWSBM::LoadReferences(std::vector<double *> *dPriorSender,
			   std::vector<double *> *dPriorReceiver,
			   std::vector<std::vector<double *> > *dPriorBlockMat,
			   std::vector<int> *dBlockMemb){
   aPriorSender = dPriorSender;
   aPriorReceiver = dPriorReceiver;
   aPriorBlockMat = dPriorBlockMat;
   aBlockMemb = dBlockMemb;
}


void CWSBM::LoadAdjacencyMatrix(int *AdjMat){
   int ii, jj;

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


CWSBM::~CWSBM (){
   if(aPriorOwner){
      delete aPriorSender;
      delete aPriorReceiver;
      delete aPriorBlockMat;
      delete aBlockMemb;
   }
}



void CWSBM::RLoadWSBM(double *rBlockMat, int *rBlockMemb,
		      double *rSenderEffects, double *rReceiverEffects,
		      double *rPosteriorMemb){
   //  Saving Pointers for Future Updates
   aRBlockMat = rBlockMat;
   aRBlockMemb = rBlockMemb;
   aRSenderEffects = rSenderEffects;
   aRReceiverEffects = rReceiverEffects;
   aRPosteriorMemb = rPosteriorMemb;

   RLoadBlockMat(rBlockMat);
   RLoadBlockMemb(rBlockMemb);
   RLoadSenderEffects(rSenderEffects);
   RLoadReceiverEffects(rReceiverEffects);
}


void CWSBM::RLoadWSBM(double *rBlockMat,
		      double *rSenderEffects, double *rReceiverEffects){
   //  Saving Pointers for Future Updates
   aRBlockMat = rBlockMat;
   aRSenderEffects = rSenderEffects;
   aRReceiverEffects = rReceiverEffects;

   //  Loading values into C++ Objects
   RLoadBlockMat(rBlockMat);
   RLoadSenderEffects(rSenderEffects);
   RLoadReceiverEffects(rReceiverEffects);
}


void CWSBM::RLoadBlockMat(double *rBlockMat){
   int ii,jj;
   int loadIter = 0;

   for(jj = 0 ; jj < aBlocks ; jj++){
      for(ii = 0; ii < aBlocks ; ii++){
	 aBlockMat[ii][jj] = rBlockMat[loadIter++];
	 aBlockMatLog[ii][jj] = log(aBlockMat[ii][jj]);
      }
   }

}

void CWSBM::RLoadBlockMemb(int *rBlockMemb){
   int ii;
   for(ii = 0 ; ii < aNodes ; ii++){
      (*aBlockMemb)[ii] = rBlockMemb[ii] - 1;
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

void CWSBM::RLoadPosteriorMemb(double *rPosteriorMemb){
   int ii, kk;
   int loadIter = 0;
   for(kk = 0 ; kk < aBlocks; kk++){
      for(ii = 0 ; ii < aNodes ; ii++){
	 aPosteriorMemb[ii][kk] = rPosteriorMemb[loadIter++];
      }
   }
   //   std::copy(rPosteriorMemb,rPosteriorMemb+aNodes*aBlocks,aPosteriorMemb);
}

/**************  SAVING STEPS TO MATRIX  ***************/
// void CWSBM::updateWSBM(int iter, int *rBlockMemb, double *rBlockMat,
// 		       double *rSenderEffects, double *rReceiverEffects,
// 		       double *rPosteriorMemb){
//    updateBlockMat(iter, rBlockMat);
//    updateBlockMemb(iter, rBlockMemb);

//    updateSenderEffects(iter, rSenderEffects);
//    updateReceiverEffects(iter, rReceiverEffects);

//    updatePosteriorMemb(iter, rPosteriorMemb);
// }

void CWSBM::updateWSBM(int iter){
   updateBlockMat(iter);
   updateBlockMemb(iter);

   updateSenderEffects(iter);
   updateReceiverEffects(iter);

   updatePosteriorMemb(iter);
}

void CWSBM::partialUpdate(int iter){
   updateSenderEffects(iter);
   updateReceiverEffects(iter);
   updateBlockMat(iter);
}



void CWSBM::updateBlockMat(int iter){
   int ii, jj;
   int saveIter = iter * aBlocks * aBlocks;
   for(jj = 0 ; jj < aBlocks ; jj++){
      for(ii = 0 ; ii < aBlocks ; ii++){
	 aRBlockMat[saveIter++] = aBlockMat[ii][jj];
      }
   }
}

void CWSBM::updateBlockMemb(int iter){
   int ii;
   int saveIter = iter * aNodes;
   for(ii = 0 ; ii < aNodes ; ii++){
      aRBlockMemb[saveIter++] = (*aBlockMemb)[ii] + 1;
   }
}

void CWSBM::updateSenderEffects(int iter){
   int ii;
   int saveIter = iter*aNodes;
   for(ii = 0 ; ii < aNodes ; ii++){
      aRSenderEffects[saveIter++] = aSenderEffects[ii];
   }
}
void CWSBM::updateReceiverEffects(int iter){
   int ii;
   int saveIter = iter*aNodes;
   for(ii = 0 ; ii < aNodes ; ii++){
      aRReceiverEffects[saveIter++] = aReceiverEffects[ii];
   }

}


void CWSBM::updatePosteriorMemb(int iter){
   int ii, kk;
   int saveIter = iter * aBlocks * aNodes;
   for(kk = 0 ; kk < aBlocks; kk++){
      for(ii = 0 ; ii < aNodes ; ii++){
	 aRPosteriorMemb[saveIter++] =  aPosteriorMemb[ii][kk];
      }
   }
}





/************************  MCMC FUNCTIONS  **********************/


void wsbmMCMC(CWSBM *myWSBM, int start, int total, int burnIn, int thin,
	      int shift_size, int extend_max, double qq,
	      double *rBlockMat, int * rBlockMemb,
	      double *rSenderEffects, double *rReceiverEffects,
	      double *logLik, double *rPosteriorMemb, int verbose){

   int ii;
   //int flatLength = myWSBM->aBlocks * (myWSBM->aNodes + myWSBM->aBlocks);
   int flatTotal = (total - burnIn) / thin;

   //  MCMC Loop
   int converged = 0, extend_count = 0;
   while((converged != 1) & (extend_count <= extend_max)){

      // Main MCMC Loop
      for(ii = start ; ii < total ; ii++){
	 // Perform an MCMC Step
	 myWSBM->step();

	 if(ii >= burnIn && ((ii - burnIn) % thin == 0)){
	    //  Save the result every thin iterations
	    logLik[(ii - burnIn)/thin] = myWSBM->LogLike();
	    myWSBM->updateWSBM((ii - burnIn)/thin);
	    // rBlockMemb, rBlockMat,
	    // rSenderEffects, rReceiverEffects,
	    // rPosteriorMemb);
	    //myWSBM->print(false);
	 }
      }

      //  Checking for convergence of log-likelihood chain
      converged = convergenceCheck(logLik,(total - burnIn)/thin,qq);
      if(verbose > 0){
	 Rprintf("Converged Status %d, after %d extensions\n",
		 converged,extend_count);
      }

      if(converged != 1){
	 extend_count = extend_count + 1;
	 //  Shift Saved Variables to Extend Chain
	 if(extend_count <= extend_max){
	    shiftFlatTable(shift_size,myWSBM->GetBlocks() * myWSBM->GetBlocks(),
			   flatTotal,rBlockMat);
	    shiftFlatTable(shift_size,myWSBM->GetNodes(),flatTotal,rBlockMemb);

	    shiftFlatTable(shift_size,myWSBM->GetNodes(),flatTotal,
			   rSenderEffects);
	    shiftFlatTable(shift_size,myWSBM->GetNodes(),flatTotal,
			   rReceiverEffects);

	    shiftFlatTable(shift_size,myWSBM->GetBlocks() * myWSBM->GetNodes(),
			   flatTotal,rPosteriorMemb);
	    std::copy(logLik + shift_size, logLik + flatTotal,logLik);
	    start = total - (shift_size*thin);
	 }
      }
   }
   if(converged != 1){
      if(verbose > -1) Rprintf("Warning: MCMC Failed to Converge\n");
   }else{
      if(verbose > 0) Rprintf("MCMC Converged\n");
   }

   if(verbose > 2) myWSBM->print(false);

}



void CWSBM::step(){

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
   drawReceiverEffects();

   if(aImputeFlag){
      Rprintf("Imputing Missing Values...\n");
      imputeMissingValues();
      computeRowColSums();
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
	 aBlockMatLog[ll][kk] = log(aBlockMat[ll][kk]);

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
	 aBlockMatLog[ii][jj] = log(aBlockMat[ii][jj]);
      }
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


//  This function updates aPosteriorMembOld to be the posterior mean given BlockMat and YY
/*  Currently Obsolete for WWSBMs

void CWSBM::getMultinomPosterior(){

   int ii,jj;
   double pMax;//, total;

   for(ii = 0 ; ii < aNodes ; ii++){
      //    Rprintf("\naPosteriorMemb[%d,] = ",ii);
      for(jj = 0 ; jj < aBlocks ; jj++){
	 (*aBlockMemb)[ii] = jj;

	 // Calculating Log Posterior Probability
	 aPosteriorMemb[ii][jj] = nodeLogLike_long(ii) + log(aPriorBlockMemb[jj]);
	 //      Rprintf("%f, ",aPosteriorMemb[ii*aBlocks + jj]);
	 if(jj == 0){
	    pMax = aPosteriorMemb[ii][jj];
	 }else{
	    if(aPosteriorMemb[ii][jj] > pMax){
	       pMax = aPosteriorMemb[ii][jj];
	    }
	 }
      }

      // Exponentiating Probabilities
      //    total = 0.0;
      for(jj = 0 ; jj < aBlocks ; jj++){
	 aPosteriorMemb[ii][jj] = exp(aPosteriorMemb[ii][jj] - pMax);
	 //      total = total + aPosteriorMemb[ii*aBlocks + jj];
      }

      //      normalizeVec(aPosteriorMemb + ii*aBlocks,aBlocks);

      normalizeVec(aPosteriorMemb[ii]);

   }
}
*/


/*****************  HELPER FUNCTIONS  *****************/
void CWSBM::computeBlockTieSums(){

   int kk;
   //  Setting values to zero
   for(kk = 0 ; kk < aBlocks ; kk++){
      aBlockTieSums[kk].assign(aBlocks,0.0);
      aBlockTieCounts[kk].assign(aBlocks,0.0);
   }

   int rr, cc;
   double srEffect;
   //  Summing over each tie
   //   Rprintf("Computing Block Tie Counts:\n");
   for(rr = 0 ; rr < aNodes ; rr++){
      for(cc = 0 ; cc < aNodes ; cc++){
	 if(!isMissing(aAdjMat[rr][cc])){
	    aBlockTieSums[(*aBlockMemb)[rr]][(*aBlockMemb)[cc]] += aAdjMat[rr][cc];
	    srEffect = aSenderEffects[rr]*aReceiverEffects[cc];
	    aBlockTieCounts[(*aBlockMemb)[rr]][(*aBlockMemb)[cc]] += srEffect;
	    //	    Rprintf("%.2f ",aBlockTieCounts[(*aBlockMemb)[rr]][(*aBlockMemb)[cc]]);
	 }
      }
      //      Rprintf("\n");
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
	 aBlockMatLog[ll][kk] = log(aBlockMat[ll][kk]);
      }
   }
}


/******************  Log-Likelihood Functions  *****************/
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
   //   int iiBlockMemb = (*aBlockMemb)[ii];
   //  int index;
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
   double mm = aBlockMat[(*aBlockMemb)[ss]][(*aBlockMemb)[rr]];
   mm *= aSenderEffects[ss];
   mm *= aReceiverEffects[rr];
   mm *= aHours;
   return mm;
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

/*
double CWSBM::BlockMatdiff(){
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
*/

void CWSBM::saveBlockMatOld(){
   int ii, jj;
   for(ii = 0 ; ii < aBlocks ; ii++){
      for(jj = 0 ; jj < aBlocks ; jj++){
	 aBlockMatOld[ii][jj] = aBlockMat[ii][jj];
      }
   }
}


void CWSBM::savePosteriorMembOld(){
   //   std::copy(BlockMat,BlockMat+aBlocks*aBlocks,aBlockMatOld);
   int ii, kk;
   for(ii = 0 ; ii < aNodes ; ii++){
      for(kk = 0 ; kk < aBlocks ; kk++){
	 aPosteriorMembOld[ii][kk] = aPosteriorMemb[ii][kk];
      }
   }
}



//  This function doesn't work with aAdjMat implementation
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

// Static Version
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
//
/*  Pointer Version
    void CWSBM::printAdjacencyMatrix(){

    int nodes = aAdjMat->size();
    int ii, jj;
    for(ii = 0 ; ii < nodes; ii++){
    for(jj = 0; jj < nodes; jj++){
    Rprintf("%d ",(*aAdjMat)[ii][jj]);
    }
    Rprintf("\n");
    }

    }
*/




