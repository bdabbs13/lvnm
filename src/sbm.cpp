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


/*****************************************************************
 *********  FUNCTION DEFINITIONS FOR SBM OBJECT CLASS  ***********
 *****************************************************************/


/***********  SBM CLASS CONSTRUCTOR AND DESTRUCTOR  **********/

//  Constructor Function for SBM object

CSBM::CSBM (int rNodes, int rBlocks, int *adjMat,
	    double *rPriorBlockMat, double * rPriorBlockMemb,
	    int mImpute) : missingVal(-1){
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
	    aAdjPartial[ii][jj] = adjMat[ii + jj*aNodes];
	 }
      }
   }

   //  Initializing aBlockMat Matrix
   double holder[2];
   std::copy(rPriorBlockMat,rPriorBlockMat+2,aPriorBlockMat);
   aBlockMat.resize(aBlocks);
   aBlockMatInv.resize(aBlocks);
   aBlockMatOld.resize(aBlocks);
   for(ii = 0; ii < aBlocks ; ii++){
      aBlockMat[ii].resize(aBlocks);
      aBlockMatInv[ii].resize(aBlocks);
      aBlockMatOld[ii].resize(aBlocks,0.0);
      for(jj = 0 ; jj < aBlocks ; jj++){
	 rdirichlet(2,aPriorBlockMat,holder);
	 aBlockMat[ii][jj] = holder[0];
	 aBlockMatInv[ii][jj] = 1 - holder[0];
      }
   }
   /*
     aBlockMat = new double[aBlocks*aBlocks];
     aBlockMatInv = new double[aBlocks*aBlocks];
     aBlockMatOld = new double[aBlocks*aBlocks];
     for(ii = 0 ; ii < aBlocks*aBlocks ; ii++){
     rdirichlet(2,aPriorBlockMat,holder);
     aBlockMat[ii] = holder[0];
     aBlockMatInv[ii] = 1.0 - aBlockMat[ii];
     }
   */
   is_BlockMat_logged = false;

   aHitMat.resize(aBlocks);
   aMissMat.resize(aBlocks);
   for(ii = 0 ; ii < aBlocks ; ii++){
      aHitMat[ii].resize(aBlocks,0.0);
      aMissMat[ii].resize(aBlocks,0.0);
   }
   //   aHitMat = new double[aBlocks*aBlocks];
   //   aMissMat = new double[aBlocks*aBlocks];


   //aPosteriorMemb = new double[aNodes*aBlocks];
   //aPosteriorMembOld = new double[aNodes*aBlocks]();
   aPosteriorMemb.resize(aNodes);
   aPosteriorMembOld.resize(aNodes);

   //aBlockMemb = new int[aNodes];
   aBlockMemb.resize(aNodes,0);
   for(ii = 0 ; ii < aNodes ; ii++){
      aPosteriorMemb[ii].resize(aBlocks);
      aPosteriorMembOld[ii].resize(aBlocks);
   }

   // Initializing aPosteriorMemb Matrix
   //   aPriorBlockMemb = new double[aBlocks];
   //   std::copy(rPriorBlockMemb,rPriorBlockMemb+aBlocks,aPriorBlockMemb);
   aPriorBlockMemb.insert(aPriorBlockMemb.end(),
			  &rPriorBlockMemb[0],&rPriorBlockMemb[aBlocks]);
   //   aPriorBlockMemb.resize(aBlocks);
   //   for(ii = 0 ; ii < aBlocks ; ii++){
   //      aPriorBlockMemb[ii] = rPriorBlockMemb[ii];
   //   }

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



CSBM::~CSBM (){
   //   delete[] YY;
   //   delete[] yyComplete;
   //   delete[] aBlockMat;
   //   delete[] aBlockMatInv;
   //   delete[] aBlockMatOld;
   //   delete[] aPosteriorMembOld;
   //   delete[] aPosteriorMemb;
   //   delete[] aBlockMemb;
   //   delete[] aPriorBlockMemb;
   //  Rprintf("Freeing Allocated Memory for SBM Object\n");
}



/***********  INPUT FUNCTION  ************/
/******  OBSOLETE
	 void CSBM::loadTable(int start, double *flatTable){
	 //  int ii, jj;
	 //  sbmLoadTable(start,aNodes,aBlocks,aBlockMat,aPosteriorMembOld,flatTable);

	 int ii, jj, offset;
	 double pMax;
	 //  skip first (start - 1) sets of draws
	 offset = (start - 1) * (aBlocks * (aNodes + aBlocks));

	 //  Loading BlockMat
	 for(ii = 0 ; ii < aBlocks*aBlocks ; ii++){
	 aBlockMat[ii] = logCheck(flatTable[offset + ii]);
	 aBlockMatInv[ii] = 1.0 - aBlockMat[ii];
	 }

	 //  Loading aPosteriorMembOld
	 offset = offset + aBlocks*aBlocks;
	 for(ii = 0 ; ii < aNodes*aBlocks ; ii++){
	 aPosteriorMembOld[ii] = flatTable[offset + ii];
	 }
	 //  Loading aBlockMemb
	 for(ii = 0 ; ii < aNodes ; ii++){
	 for(jj = 0 ; jj < aBlocks; jj++){
	 if(jj == 0){
	 aBlockMemb[ii] = jj;
	 pMax = aPosteriorMembOld[ii*aBlocks + jj];
	 }else if(aPosteriorMembOld[ii*aBlocks + jj] > pMax){
	 pMax = aPosteriorMembOld[ii*aBlocks + jj];
	 aBlockMemb[ii] = jj;
	 }
	 }
	 }

	 imputeMissingValues();
	 }
*/
void CSBM::RLoadSBM(double *rBlockMat, int *rBlockMemb){
   RLoadBlockMat(rBlockMat);
   RLoadBlockMemb(rBlockMemb);
}

void CSBM::RLoadBlockMat(double *rBlockMat){
   int ii,jj;
   int loadIter = 0;
   //   std::copy(rBlockMat,rBlockMat + (aBlocks*aBlocks), BB;
   for(jj = 0 ; jj < aBlocks ; jj++){
      for(ii = 0; ii < aBlocks ; ii++){
	 //	 BB[ii][jj] = rBlockMat[ii + jj*aBlocks];
	 aBlockMat[ii][jj] = rBlockMat[loadIter++];
	 aBlockMatInv[ii][jj] = 1.0 - aBlockMat[ii][jj];
      }
   }
   /*
     for(ii = 0; ii < aBlocks*aBlocks ; ii++){
     aBlockMatInv[ii] = 1.0 - aBlockMat[ii];
     }
   */
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
   //   std::copy(rPosteriorMemb,rPosteriorMemb+aNodes*aBlocks,aPosteriorMemb);
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


/******  OUTPUT FUNCTIONS  *******/
/**** OBSOLETE
      void CSBM::updateFlatTable(int iter, double *flatTable){
      expBlockMat();
      int ii, offset;

      // Saving BlockMat
      offset = iter * (aBlocks * (aNodes + aBlocks));
      for(ii = 0 ; ii < aBlocks*aBlocks ; ii++){
      flatTable[offset + ii] = aBlockMat[ii];
      }

      // Saving aPosteriorMembOld
      offset = offset + aBlocks*aBlocks;
      for(ii = 0 ; ii < aNodes*aBlocks ; ii++){
      flatTable[offset + ii] = aPosteriorMembOld[ii];
      }
      }
*****/
void CSBM::updateSBM(int iter, int *rBlockMemb, double *rBlockMat,
		     double *rPosteriorMemb){
   updateBlockMat(iter, rBlockMat);
   updateBlockMemb(iter, rBlockMemb);
   updatePosteriorMemb(iter, rPosteriorMemb);
}

void CSBM::updateSBM(int iter, int *rBlockMemb, double *rBlockMat,
		     double *rPosteriorMemb, double *rPriorBlockMat){
   updateBlockMat(iter, rBlockMat);
   updateBlockMemb(iter, rBlockMemb);
   updatePosteriorMemb(iter, rPosteriorMemb);
   updatePriorBlockMemb(rPriorBlockMat);
}

void CSBM::updateBlockMat(int iter, double *rBlockMat){
   expBlockMat();
   //   std::copy(aBlockMat,aBlockMat+(aBlocks*aBlocks),rBlockMat + (iter * aBlocks*aBlocks));
   int ii, jj;
   int saveIter = iter * aBlocks * aBlocks;
   for(jj = 0 ; jj < aBlocks ; jj++){
      for(ii = 0 ; ii < aBlocks ; ii++){
	 rBlockMat[saveIter++] = aBlockMat[ii][jj];
      }
   }

   /*
     int offset = iter * aBlocks * aBlocks;
     for(jj = 0 ; jj < aBlocks ; jj++){
     for(ii = 0 ; ii < aBlocks ; ii++){
     rBlockMat[ii + jj*aBlocks + offset] = aBlockMat[ii][jj];
     }
     }*/

}

void CSBM::updateBlockMemb(int iter, int *rBlockMemb){
   int ii;
   int saveIter = iter * aNodes;
   for(ii = 0 ; ii < aNodes ; ii++){
      rBlockMemb[saveIter++] = aBlockMemb[ii] + 1;
   }
   //   std::copy(aBlockMemb.begin(),aBlockMemb.data() + aNodes,
   //	     rBlockMemb + (iter * aNodes));
}

void CSBM::updatePriorBlockMemb(double *rPriorBlockMemb){

   std::copy(&aPriorBlockMemb[0], &aPriorBlockMemb[aBlocks],
	     rPriorBlockMemb);
   //   std::copy(aPriorBlockMemb, aPriorBlockMemb + aBlocks, rPriorBlockMat);
}

void CSBM::updatePosteriorMemb(int iter, double *rPosteriorMemb){
   int ii, kk;
   int saveIter = iter * aBlocks * aNodes;
   for(kk = 0 ; kk < aBlocks; kk++){
      for(ii = 0 ; ii < aNodes ; ii++){
	 rPosteriorMemb[saveIter++] =  aPosteriorMemb[ii][kk];
	 //	 Rprintf("loadIter = %d, aPosteriorMemb = %.4f\n",loadIter,
	 //aPosteriorMemb[ii][kk]);
      }
   }
   //   std::copy(aPosteriorMemb,aPosteriorMemb +(aNodes*aBlocks),
   //	     rPosteriorMemb + (iter*aBlocks*aNodes));
}





/************************  MCMC FUNCTIONS  **********************/


void sbmMCMC(CSBM *mySBM, int start, int total, int burnIn, int thin,
	     int shift_size, int extend_max, double qq,  //double *flatTable,
	     double *rBlockMat, int * rBlockMemb,  double *logLik,
	     double *rPosteriorMemb, int verbose){

   int ii;
   //int flatLength = mySBM->aBlocks * (mySBM->aNodes + mySBM->aBlocks);
   int flatTotal = (total - burnIn) / thin;

   //  MCMC Loop
   int converged = 0, extend_count = 0;
   while((converged != 1) & (extend_count <= extend_max)){

      // Main MCMC Loop
      for(ii = start ; ii < total ; ii++){
	 mySBM->step();

	 if(ii >= burnIn && ((ii - burnIn) % thin == 0)){
	    logLik[(ii - burnIn)/thin] = mySBM->LogLike();
	    //mySBM->updateFlatTable((ii - burnIn)/thin,flatTable);
	    mySBM->updateSBM((ii - burnIn)/thin, rBlockMemb, rBlockMat, rPosteriorMemb);
	    //mySBM->updateBlockMat((ii - burnIn)/thin,rBlockMat);
	    //mySBM->updateBlockMemb((ii - burnIn)/thin,rBlockMemb);
	    //mySBM->updatePosteriorMemb((ii - burnIn)/thin,rPosteriorMemb);
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

	    shiftFlatTable(shift_size,mySBM->GetBlocks() * mySBM->GetBlocks(),flatTotal,rBlockMat);
	    shiftFlatTable(shift_size,mySBM->GetNodes(),flatTotal,rBlockMemb);
	    shiftFlatTable(shift_size,mySBM->GetBlocks() * mySBM->GetNodes(),flatTotal,rPosteriorMemb);
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



void CSBM::step(){
   drawBlockMemb(); // logs then exps BlockMat
   drawBlockMat(); // updates BlockMat, but doesn't use them, then sets BlockMat_logged to false
   rotate();
   if(aImputeFlag){
      //      Rprintf("Imputing Missing Values...\n");
      imputeMissingValues(); // assumes exp BlockMat
   }
}

void CSBM::drawBlockMat(){
   //  int dd2 = aBlocks*aBlocks;
   //  double aHitMat[dd2], aMissMat[dd2];
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
   /*
     for(ii = 0 ; ii < dd2 ; ii++){
     alpha_new[0] = aPriorBlockMat[0] + aHitMat[ii];
     alpha_new[1] = aPriorBlockMat[1] + aMissMat[ii];
     rdirichlet(2,alpha_new,holder);
     aBlockMat[ii] = holder[0];
     aBlockMatInv[ii] = 1.0 - aBlockMat[ii];
     }
   */
   is_BlockMat_logged = false;
}

void CSBM::drawBlockMemb(){
   int ii,jj,kk;
   double pMax;
   int ind[aBlocks];

   for(ii = 0 ; ii < aNodes ; ii++){
      //total = 0.0;
      for(jj = 0 ; jj < aBlocks ; jj++){
	 aBlockMemb[ii] = jj;

	 // Calculating Log Posterior Probability
	 aPosteriorMemb[ii][jj] = nodeLogLike(ii) + log(aPriorBlockMemb[jj]);
	 //pp[jj] = nodeLogLike(ii) + log(aPriorBlockMemb[jj]);
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
	 //      total = total + aPosteriorMemb[ii*aBlocks + jj];
      }

      //    int hhOne = -1;
      // Normalizing Probabilities
      /*
	for(jj = 0 ; jj < aBlocks ; jj++){
	aPosteriorMemb[ii*aBlocks + jj] = aPosteriorMemb[ii*aBlocks + jj] / total;
	}
      */
      //      normalizeVec(aPosteriorMemb+ii*aBlocks,aBlocks);
      //      logZeroFix(aPosteriorMemb + ii*aBlocks,aBlocks);
      normalizeVec(aPosteriorMemb[ii]);
      logZeroFix(aPosteriorMemb[ii]);

      // Drawing from Multinomial
      rmultinom(1,aPosteriorMemb[ii].data(),aBlocks,ind);

      // Setting aPosteriorMemb[ii,] Value
      for(kk = 0 ; kk < aBlocks ; kk++){
	 if(ind[kk] == 1){
	    //	aPosteriorMembOld[ii*aBlocks + kk] = 1.0;
	    aBlockMemb[ii] = kk;
	    break;
	 }
	 //else{
	 //aPosteriorMembOld[ii*aBlocks + kk] = 0.0;
	 //}
      }
   }
   //  expBlockMat();
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
      /*
	Rprintf("uncount: %d ",uncount);
	Rprintf("vacancy: %d ",vacancy);
	Rprintf("newlabel: %d ",newlabel);
	Rprintf("\n");
	RprintIntMat(1, aBlocks, assigned);
      */

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
      aBlockMemb[ii] = assigned[aBlockMemb[ii]];
      for(kk = 0 ; kk < aBlocks ; kk++){
	 if(aBlockMemb[ii] == kk){
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
	 aBlockMatInv[ii][jj] = 1.0 - aBlockMat[ii][jj];
	 //	 aBlockMat[ii*aBlocks + jj] = aBlockMat[assigned[ii]*aBlocks + assigned[jj]];
	 //aBlockMatInv[ii*aBlocks + jj] = 1.0 - aBlockMat[ii*aBlocks + jj];
      }
   }
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



/************************  EM FUNCTIONS  *********************/


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
   mySBM->updateSBM(0,rBlockMemb,rBlockMat,rPosteriorMemb,rPriorBlockMat);
   /*
     mySBM->updatePosteriorMemb(0,rPosteriorMemb);
     mySBM->updateBlockMat(0,rBlockMat);
     mySBM->updateBlockMemb(0,rBlockMemb);
     mySBM->updatePriorBlockMemb(rPriorBlockMat);
   */
   logLik[0] = mySBM->LogLike();
   //  delete[] BlockMat_log;
}


void CSBM::iterEM (){
   int ss, rr, ii, jj;
   int rowd, rown,bbind;
   double total, pipi;


   //  E Step
   //   std::copy(aPosteriorMemb,aPosteriorMemb+(aNodes*aBlocks),aPosteriorMembOld);
   savePosteriorMembOld();
   getMultinomPosterior();  //  Update Posterior Mean aPosteriorMemb
   //Rprintf("\naPosteriorMemb:\n");
   //RprintDoubleMat(aNodes,aBlocks,aPosteriorMemb);


   // M Step
   //   colSums(aPosteriorMemb,aNodes,aBlocks,aPriorBlockMemb);
   colSums(aPosteriorMemb,aPriorBlockMemb);
   double priorTotal = 0.0;
   //   RprintDoubleMat(1,aBlocks,aPriorBlockMemb);

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
   //   RprintDoubleMat(1,aBlocks,aPriorBlockMemb);

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
	       //	       if((ii != jj) & (YY[ii*aNodes + jj] >= 0)){
	       if(!isMissing(aAdjMat[ii][jj])){
		  pipi = aPosteriorMemb[ii][ss] * aPosteriorMemb[jj][rr];  // How much to contribute
		  //aBlockMat[bbind] = aBlockMat[bbind] + pipi * YY[rown + jj];  // numerator contribution
		  aBlockMat[ss][rr] += pipi * aAdjMat[ii][jj];  // numerator contribution
		  total = total + pipi;  //  Denominator contribution
	       }
	    }
	 }
	 if(total < MIN_LOG){
	    aBlockMat[ss][rr] = MIN_LOG;
	 }else{
	    aBlockMat[ss][rr] = aBlockMat[ss][rr] / total;
	    aBlockMat[ss][rr] = logCheck(aBlockMat[ss][rr]);
	    /*
	      if(aBlockMat[bbind] > MAX_LOG){
	      aBlockMat[bbind] = MAX_LOG;
	      }else if(aBlockMat[bbind] < MIN_LOG){
	      aBlockMat[bbind] = MIN_LOG;
	      }
	    */
	 }
	 aBlockMatInv[ss][rr] = 1.0 - aBlockMat[ss][rr];
      }
   }
   is_BlockMat_logged = false;
}



//  This function updates aPosteriorMembOld to be the posterior mean given BlockMat and YY
void CSBM::getMultinomPosterior(){

   int ii,jj;
   double pMax;//, total;

   for(ii = 0 ; ii < aNodes ; ii++){
      //    Rprintf("\naPosteriorMemb[%d,] = ",ii);
      for(jj = 0 ; jj < aBlocks ; jj++){
	 aBlockMemb[ii] = jj;

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
      /*      Rprintf("ori:  ");
      for(jj = 0 ; jj < aBlocks; jj++){
	 Rprintf("%.2f ",aPosteriorMemb[ii][jj]);
      }
      */
      normalizeVec(aPosteriorMemb[ii]);
      /*
      Rprintf("  new:  ");
      for(jj = 0 ; jj < aBlocks; jj++){
	 Rprintf("%.2f ",aPosteriorMemb[ii][jj]);
      }
      Rprintf("\n");
      */
   }
}


/*****************  HELPER FUNCTIONS  *****************/
void CSBM::computeHitMiss(){
   //   int dd2 = aBlocks*aBlocks;
   int ii, rr, cc, index;

   //  Setting values to zero
   for(ii = 0 ; ii < aBlocks ; ii++){
      aHitMat[ii].assign(aBlocks,0.0);
      aMissMat[ii].assign(aBlocks,0.0);
   }
   //   std::fill(aHitMat, aHitMat+dd2, 0.0);
   //   std::fill(aMissMat,aMissMat+dd2,0.0);

   //  Counting hits and misses
   for(rr = 0 ; rr < aNodes ; rr++){
      for(cc = 0 ; cc < aNodes ; cc++){
	 //	 index = aBlockMemb[rr] * aBlocks + aBlockMemb[cc];
	 //if(YY[rr * aNodes + cc] == 1){
	 if(aAdjMat[rr][cc] == 1){
	    aHitMat[aBlockMemb[rr]][aBlockMemb[cc]]++;
	    //aHitMat[index] = aHitMat[index] + 1;
	    //}else if(YY[rr * aNodes + cc] == 0){
	 }else if(aAdjMat[rr][cc] == 0){
	    aMissMat[aBlockMemb[rr]][aBlockMemb[cc]]++;
	    //aMissMat[index] = aMissMat[index] + 1;
	 }
      }
   }
}

void CSBM::computeBlockMatMLE(){
   computeHitMiss();

   int ii,jj;//dd2 = aBlocks*aBlocks;
   for(ii = 0 ; ii < aBlocks ; ii++){
      for(jj = 0 ; jj < aBlocks ; jj++){
	 if(aHitMat[ii][jj] == 0.0){
	    aBlockMat[ii][jj] = 0.0;
	 }else{
	    aBlockMat[ii][jj] = aHitMat[ii][jj] / (aHitMat[ii][jj] + aMissMat[ii][jj]);
	 }
	 aBlockMat[ii][jj] = logCheck(aBlockMat[ii][jj]);
	 aBlockMatInv[ii][jj] = 1.0 - aBlockMat[ii][jj];
      }
   }
}


/******************  Log-Likelihood Functions  *****************/
double CSBM::LogLike(){
   int rr, cc;//, index;
   logBlockMat();

   double total = 0.0;
   for(rr = 0 ; rr < aNodes ; rr++){
      for(cc = 0 ; cc < aNodes ; cc++){
	 //	 total = total + tieLogLike(YY[rr*aNodes+cc],aBlockMemb[rr],aBlockMemb[cc]);
	 total = total + tieLogLike(aAdjMat[rr][cc],aBlockMemb[rr],aBlockMemb[cc]);
      }
   }
   //   expBlockMat();  123456
   return(total);
}


double CSBM::nodeLogLike(int ii){
   logBlockMat();
   int jj;
   int iiBlockMemb = aBlockMemb[ii];
   //  int index;
   double total = 0.0;

   //  sums the (i,j)th and (j,i)th probabilities for all j
   for(jj = 0 ; jj < aNodes ; jj++){
      if(jj != ii){

	 //  (i,j)th term
	 //	 total = total + tieLogLike(YY[ii*aNodes + jj],iiBlockMemb,aBlockMemb[jj]);
	 total = total + tieLogLike(aAdjMat[ii][jj],iiBlockMemb,aBlockMemb[jj]);

	 //  (j,i)th term
	 //total = total + tieLogLike(YY[jj*aNodes + ii],aBlockMemb[jj],iiBlockMemb);
	 total = total + tieLogLike(aAdjMat[jj][ii],aBlockMemb[jj],iiBlockMemb);
      }
   }

   return(total);
}


double CSBM::tieLogLike(int yy, int sendBlock, int recBlock){

   if(yy == 1){
      return aBlockMat[sendBlock][recBlock];//sendBlock * aBlocks + recBlock];
   }else if(yy == 0){
      return aBlockMatInv[sendBlock][recBlock];// * aBlocks + recBlock];
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
	 //total = total + tieLogLike_sender(YY[ii*aNodes + jj],iiBlockMemb,jj);
	 total = total + tieLogLike_sender(aAdjMat[ii][jj],iiBlockMemb,jj);
	 //  (j,i)th term
	 //total = total + tieLogLike_receiver(YY[jj*aNodes + ii],iiBlockMemb,jj);
	 total = total + tieLogLike_receiver(aAdjMat[jj][ii],iiBlockMemb,jj);
      }
   }

   return(total);
}


// Log-Likelihood assuming Posterior Matrix
double CSBM::tieLogLike_sender(int yy, int sendBlock, int receiver){
   int kk;
   //   int baseBlock = sendBlock*aBlocks;
   //   int baseNode = receiver*aBlocks;
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


/*******************  BlockMat UPDATING FUNCTIONS  ******************/

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

void CSBM::savePosteriorMembOld(){
   //   std::copy(BlockMat,BlockMat+aBlocks*aBlocks,aBlockMatOld);
   int ii, kk;
   for(ii = 0 ; ii < aNodes ; ii++){
      for(kk = 0 ; kk < aBlocks ; kk++){
	 aPosteriorMembOld[ii][kk] = aPosteriorMemb[ii][kk];
      }
   }
}



//  This function doesn't work with aAdjMat implementation
void CSBM::print(bool printNetwork){
   Rprintf("aBlocks = %d, aNodes = %d\n",aBlocks,aNodes);
   Rprintf("bbPrior = (%.2f, %.2f)\n",aPriorBlockMat[0],aPriorBlockMat[1]);
   Rprintf("multiPrior = ");
   printPriorBlockMemb();
   //   RprintDoubleMat(1,aBlocks,aPriorBlockMemb);
   //   Rprintf("aBlockMat = \n"); RprintDoubleMat(aBlocks, aBlocks, aBlockMat);
   printBlockMat();
   printPosteriorMemb();
   //   printPosteriorMemb(aNodes,aBlocks,aPosteriorMemb);
   if(printNetwork){
      //      Rprintf("YY = \n"); RprintIntMat(aNodes,aNodes,YY);
      printAdjacencyMatrix();
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
//
/*  Pointer Version
    void CSBM::printAdjacencyMatrix(){

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




