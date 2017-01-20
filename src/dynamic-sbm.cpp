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


/*****************************************************************
 *********  FUNCTION DEFINITIONS FOR WSBM OBJECT CLASS  ***********
 *****************************************************************/




/***********  WSBM CLASS CONSTRUCTOR AND DESTRUCTOR  **********/

//  Constructor Function for WSBM object

CDynSBM::CDynSBM (int rNodes, int rBlocks, int rTimes, int rTimeClasses,
		  int total, int *rTimeMap, double *rHours,
		  int mImpute) : missingVal(-1)
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

void CDynSBM::LoadAdjacencyMatrices(int *AdjMat){
   int tt;
   for(tt = 0 ; tt < aTimes ; tt++){
      aWsbmList[tt].LoadAdjacencyMatrix(&AdjMat[tt*aNodes*aNodes]);
   }

}

void CDynSBM::LoadHyperPriors(double *rHyperSender, double *rHyperReceiver,
			      double *rHyperBlockMat){
   int ii;
   for(ii = 0 ; ii < 4 ; ii++){
      aHyperSender[ii] = rHyperSender[ii];
      aHyperReceiver[ii] = rHyperReceiver[ii];
      aHyperBlockMat[ii] = rHyperBlockMat[ii];
   }
}

void CDynSBM::LoadParameters(double *rSenderEffects, double *rReceiverEffects,
			     double *rBlockEffects, int *rBlockMemb){

   int tt;
   for(tt = 0 ; tt < aTimes; tt++){
      aWsbmList[tt].RLoadWSBM(&rBlockEffects[tt*aTotal*aBlocks*aBlocks],
			      &rSenderEffects[tt*aTotal*aNodes],
			      &rReceiverEffects[tt*aTotal*aNodes]);
   }

   aRBlockMemb = rBlockMemb;
   int ii;
   for(ii = 0 ; ii < aNodes ; ii++){
      aBlockMemb[ii] = rBlockMemb[ii] - 1;
      //      Rprintf("ii = %d, rBlockMemb = %d, aBlockMemb = %d\n",
      //ii,rBlockMemb[ii],aBlockMemb[ii]);
      aPosteriorMemb[ii][aBlockMemb[ii]] = 1.0;
   }


}

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
   DrawPriors();
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
   for(tt = 0 ; tt < aTimes ; tt++){
      aWsbmList[tt].partialStep();
   }
   DrawBlockMemb();
}

void CDynSBM::DrawBlockMemb(){
   int ii,jj,kk;
   double pMax;
   int ind[aBlocks];

   //   Rprintf("Log-Likelihood = %.2f\n",LogLike());

   for(ii = 0 ; ii < aNodes ; ii++){
      for(jj = 0 ; jj < aBlocks ; jj++){
	 aBlockMemb[ii] = jj;

	 // Calculating Log Posterior Probability
	 //	 Rprintf("ii = %d, jj = %d",ii,jj);
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
      aBlockMemb[ii] = kk;
   }

}

double CDynSBM::nodeLogLike(int node){
   int tt;
   double sum = 0.0;
   for(tt = 0 ; tt < aTimes ; tt++){
      sum += aWsbmList[tt].nodeLogLike(node);
   }
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

      //  Making Posterior Draws
      for(tt = 0 ; tt < aClasses; tt++){
	 rGammaPriorStep((*(aPriorSender[tt]))[ii],
			 pp[tt],qq[tt],rr[tt],ss[tt],prop_sd);
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
	 rGammaPriorStep((*(aPriorReceiver[tt]))[ii],
			 pp[tt],qq[tt],rr[tt],ss[tt],prop_sd);
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
	    rGammaPriorStep((*(aPriorBlockMat[tt]))[ll][kk],
			 pp[tt],qq[tt],rr[tt],ss[tt],prop_sd);
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
