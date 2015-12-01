/*****  lsm.cpp
 *****  
 *****/

//#include <stdio.h>
//#include <math.h>
//#include <stdlib.h>
//#include <iostream>
//#include <fstream>
//#include <ctime>
#include <cmath>
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>
#include <R.h>
#include <Rmath.h>
#include "lsm.h"
#include "mmsbm.h"
//using namespace std;
//using std::ifstream;
using std::isfinite;
//using std::ofstream;


extern "C" {
  /*********************************************************
   ************  Externally Available Functions ************
   ********************************************************/
  void lsm(int *iters, int *nn_t, int *kk_t, int *YY,
	   double *mhControl, double *alphaPrior, double *zzPrior,
	   double *initPrior, double *flatTable, 
	   int *burn_t, int *thin_t, int *start_t,
	   int *multiImpute_t){
    
    
    GetRNGstate();

    int total = *iters;
    int burnIn = *burn_t, thin = *thin_t;
    int nn = *nn_t, dd = *kk_t, start = *start_t;
    int multiImpute = *multiImpute_t;
    int yyComplete[nn*nn];
    double dMat[nn*nn];
    double alpha[1], ZZ[nn*dd];

    lsmMCMC(total,burnIn,thin,YY,nn,dd,
	    mhControl,initPrior,alphaPrior,zzPrior,
	    alpha, ZZ, dMat, yyComplete, flatTable, start,multiImpute);

    PutRNGstate();
  }
  

    /*****************************************************
   *************  MCMC CONTROL FUNCTIONS  **************
   ****************************************************/
  
  /**********  SBM  **********/
  
  void lsmInit(int *YY, int nn, int dd,
	       double *mhControl, double *initPrior,
	       double *alpha, double *ZZ,
	       int *yyComplete,double *dMat,
	       double *flatTable,
	       int start, int multiImpute){

    int ii;
    
    //  Initializing Complete Matrix
    if(start > 0){
      lsmLoadTable(start,nn,dd,alpha,ZZ,flatTable);
      if(multiImpute == 1){
	lsmImputeMissingValues(nn,dd,YY,yyComplete,alpha,ZZ,dMat);
      }else{
	for(ii = 0 ; ii < nn*nn ; ii++){
	  yyComplete[ii] = YY[ii];
	}	
      }
    }else{
      
      if(multiImpute == 1){
	for(ii = 0 ; ii < nn*nn ; ii++){
	  if(YY[ii] > 0 ){
	    yyComplete[ii] = YY[ii];
	  }else{
	    yyComplete[ii] = 0;
	  }
	}
      }else{	
	for(ii = 0 ; ii < nn*nn ; ii++){
	  yyComplete[ii] = YY[ii];
	}	
      }

      double ll = INFINITY;
      ZZ[0] = 0.0; ZZ[1] = 0.0; ZZ[2] = 0.0;
      //  Try Initial Values until a finite log-likelihood is achieved
      while(!isfinite(ll)){
	
	alpha[0] = rnorm(0.0,initPrior[0]);
	for(ii = 4 ; ii < nn*dd ; ii++){
	  ZZ[ii] = rnorm(0.0,initPrior[1]);
	}
	
	//  Enforcing Positive Restrictions
	if(ZZ[4] < 0){
	  ZZ[4] = -1.0 * ZZ[4];
	}
	if(ZZ[3] < 0){
	  ZZ[3] = -1.0 * ZZ[3];
	}
	
	
	ll = lsmLogLik(nn,dd,YY,alpha,ZZ,dMat);
      }
      
    }
  }
  void lsmStep(int *YY, int nn, int dd, 
	       double *mhControl, double *alphaPrior, double *zzPrior,
	       double *alpha, double *ZZ, double *dMat){
    
    //  Drawing a new alpha value;
    lsmDrawAlpha(nn, dd, YY, alpha, ZZ, alphaPrior, mhControl, dMat);

    //  Drawing new ZZ vectors;
    lsmDrawZZ(nn, dd, YY, mhControl, zzPrior, alpha, ZZ, dMat);
    
  }
  
  void lsmImputeMissingValues(int nn, int dd,
			      int *YY, int *yyComplete,
			      double *alpha, double *ZZ, double *dMat){
    
    int ii;
    double prob;
    //  Updating Distance Matrix
    distMat(nn,dd,ZZ,dMat);
    
    for(ii = 0 ; ii < nn*nn ; ii++){
      if(YY[ii] < 0){
	if(ii / nn != ii % nn){
	  prob = logitInverse(alpha[0] - dMat[ii]);
	  yyComplete[ii] = (int) rbinom(1.0,prob);
	}
      }
    }
  }
  
  void lsmLoadTable(int start, int nn, int dd,
		    double *alpha, double *ZZ,
		    double *flatTable);
  
  void lsmMCMC(int total, int burnIn, int thin,
	       int *YY,int nn,int dd,
	       double *mhControl, double *initPrior,
	       double *alphaPrior, double *zzPrior,
	       double *alpha, double *ZZ, double *dMat,
	       int *yyComplete,double *flatTable,int start,
	       int multiImpute){

    int ii;

    lsmInit(YY,nn,dd,mhControl,initPrior,alpha,ZZ,
	    yyComplete,dMat,flatTable,start,multiImpute);

    for(ii = start ; ii < total ; ii++){
      if(multiImpute == 1){
	lsmImputeMissingValues(nn,dd,YY,yyComplete,alpha,ZZ,dMat);
      }
      lsmStep(yyComplete,nn,dd,mhControl,alphaPrior,zzPrior,alpha, ZZ, dMat);
      if((ii >= burnIn) && ((ii - burnIn) % thin ==0)){
	lsmUpdateFlatTable((ii - burnIn)/thin, nn, dd, 
			   YY, alpha, ZZ, dMat, flatTable);
      }
    }
  }
  
  

    
  /*****************************************************
   ***************  SAMPLING FUNCTIONS  ****************
   ****************************************************/
  
  /********  GENERIC  ********/
  
  //void rdirichlet(int k, double *alpha, double *x);
  
  double logitInverse(double x){
    return 1.0/(1.0 + exp(-1.0 * x));
  }

  void distMat(int nn, int dd, double *ZZ, double *dMat){
    int ii,jj,kk,ind;
    double tmp;
    for(ii = 0 ; ii < nn ; ii++){
      for(jj = 0 ; jj < ii ; jj++){
	ind = ii*nn + jj;
	dMat[ind] = 0.0;
	for(kk = 0 ; kk < dd ; kk++){
	  tmp = ZZ[ii*dd + kk] - ZZ[jj*dd + kk];
	  dMat[ind] = dMat[ind] + tmp*tmp;
	}
	dMat[jj*nn + ii] = dMat[ind];
      }
    }
    
  }


  double lsmLogPriorAlpha(double *alpha, double *alphaPrior){
    return dnorm(alpha[0],0.0,alphaPrior[0],1);
  }


  double lsmLogPriorZZ(int dd, double *zz,double zzPrior){
    int ii;
    double total = 0.0;
    for(ii = 0 ; ii < dd ; ii++){
      total = total + dnorm(zz[ii],0.0,zzPrior,1);
    }
    return total;
  }
  
  
  /**********  LSM  **********/
  
  void lsmDrawAlpha(int nn, int dd, int *YY,
		    double *alpha, double *ZZ, 
		    double *alphaPrior,double *mhControl, double *dMat){
    
    double llOld, llNew, lpOld, lpNew;
    double alphaNew = rnorm(alpha[0],mhControl[0]);
    //Rprintf("cur = %f;  prop = %f\n",alpha[0],alphaNew);
    
    llOld = lsmLogLik(nn,dd,YY,alpha,ZZ,dMat);
    llNew = lsmLogLik(nn,dd,YY,&alphaNew,ZZ,dMat);
    
    lpOld = lsmLogPriorAlpha(alpha,alphaPrior);
    lpNew = lsmLogPriorAlpha(&alphaNew,alphaPrior);
    
    double logRR = llNew - llOld + lpNew - lpOld;
    double draw = runif(0.0,1.0);
    //Rprintf("unif = %f;  logRR = %f\n",draw,logRR);
    if(log(draw) < logRR){
      alpha[0] = alphaNew;
    }
  }
  
  double lsmLogLik(int nn, int dd, int *YY, 
		   double *alpha, double *ZZ, double *dMat){
    
    distMat(nn,dd,ZZ,dMat);
    //    RprintDoubleMat(nn,nn,dMat);
    int ii,jj;
    double total = 0.0,tmp;
    for(ii = 0 ; ii < nn ; ii++){
      for(jj = 0 ; jj < ii ; jj++){
	tmp = logitInverse(alpha[0] - dMat[ii*nn + jj]);
	if(YY[ii*nn + jj] == 1){
	  total = total + log(tmp);
	}else if(YY[ii*nn + jj] == 0){
	  total = total + log(1.0 - tmp);
	}
	
	if(YY[jj*nn + ii] == 1){
	  total = total + log(tmp);
	}else if(YY[jj*nn + ii] == 0){
	  total = total + log(1.0 - tmp);
	}
      }
    }
    //    Rprintf("ll = %f\n",total);
    return total;
  }


  double lsmLogLikii(int nn, int dd, int *YY, int ii,
		   double *alpha, double *ZZ, double *dMat){
    
    distMat(nn,dd,ZZ,dMat);
    //    RprintDoubleMat(nn,nn,dMat);
    int jj;
    double total = 0.0,tmp;
    for(jj = 0 ; jj < nn ; jj++){
      if(jj != ii){
	tmp = logitInverse(alpha[0] - dMat[ii*nn + jj]);
	if(YY[ii*nn + jj] == 1){
	  total = total + log(tmp);
	}else if(YY[ii*nn + jj] == 0){
	  total = total + log(1.0 - tmp);
	}
	
	if(YY[jj*nn + ii] == 1){
	  total = total + log(tmp);
	}else if(YY[jj*nn + ii] == 0){
	  total = total + log(1.0 - tmp);
	}
      }
    }
    //    Rprintf("ll = %f\n",total);
    return total;
  }

  
  void lsmDrawZZ(int nn, int dd, int *YY,
		 double *mhControl, double *zzPrior,
		 double *alpha, double *ZZ, double *dMat){
    
    int ii,kk;
    double llOld, llNew, lpOld, lpNew;
    double draw, logRR;
    double ZZnew[nn*dd];
    for(ii = 0 ; ii < nn*dd ; ii++){
      ZZnew[ii] = ZZ[ii];
    }
    
    for(kk = 0 ; kk < dd ; kk++){
      //      llOld = lsmLogLik(nn,dd,YY,alpha,ZZ,dMat);

      for(ii = 1 ; ii < nn ; ii++){
	if(ii != 1 || kk != 0){
	  llOld = lsmLogLikii(nn,dd,YY,ii,alpha,ZZ,dMat);
	  lpOld = lsmLogPriorZZ(dd, &ZZ[ii*dd],*zzPrior);

	  ZZnew[ii*dd + kk] = rnorm(ZZ[ii*dd + kk],mhControl[1]);
	  if((kk == 0 && ii == 2) || (kk == 1 && ii == 1)){
	    if(ZZnew[ii*dd + kk] < 0){
			ZZnew[ii*dd + kk] = -1.0 * ZZnew[ii*dd + kk];
	    }
	  }
	  //	  llNew = lsmLogLik(nn,dd,YY,alpha,ZZnew,dMat);
	  llNew = lsmLogLikii(nn,dd,YY,ii,alpha,ZZnew,dMat);
	  lpNew = lsmLogPriorZZ(dd,&ZZnew[ii*dd],*zzPrior);
	  
	  //	  logRR = llNew - llOld + lpNew - lpOld;
	  logRR = llNew - llOld + lpNew - lpOld;
	  //	  	  Rprintf("logRR = %f;  logR2 = %f\n",logRR,logR2);

	  draw = runif(0.0,1.0);
	  if(log(draw) < logRR){
	    ZZ[ii*dd + kk] = ZZnew[ii*dd + kk];
	    llOld = llNew;
	  }else{
	    ZZnew[ii*dd + kk] = ZZ[ii*dd + kk];
	  }
	}
      }
    }
    
    return;
  }
  
  /*****************************************************
   ****************  OUTPUT FUNCTIONS  *****************
   ****************************************************/
  
  /*  
  //  FUNCTIONS FOR PRINTING IN R
  void RprintDoubleMat(int rows, int cols, double *mat){
    int ii, jj;
    for(ii = 0 ; ii < rows ; ii++){
      for(jj = 0 ; jj < cols ; jj++){
	Rprintf("%f ",mat[ii*cols + jj]);
      }
      Rprintf("\n");
    }
  }
  
  void RprintIntMat(int rows, int cols, int *mat){
    int ii, jj;
    for(ii = 0 ; ii < rows ; ii++){
      for(jj = 0 ; jj < cols ; jj++){
	Rprintf("%d ",mat[ii*cols + jj]);
      }
      Rprintf("\n");
    }
  }
  */
  //  UPDATING FLAT TABLE
  void lsmUpdateFlatTable(int iter, int nn, int dd, int *YY,
		       double *alpha, double *ZZ, double *dMat,
		       double *flatTable){

    int ii, offset;
    offset = iter * ((nn * dd) + 2);
    
    flatTable[offset] = lsmLogLik(nn,dd,YY,alpha,ZZ,dMat);
    offset++;

    flatTable[offset] = alpha[0];
    offset++;

    for(ii = 0 ; ii < nn*dd ; ii++){
      flatTable[offset + ii] = ZZ[ii];
    }
    
  }

  void lsmLoadTable(int start, int nn, int dd,
			  double *alpha, double *ZZ,
			  double *flatTable){
    int ii, offset;
    offset = (start - 1) * ((nn*dd) +2) + 1;
    
    alpha[0] = flatTable[offset];
    offset++;
    
    for(ii = 0 ; ii < nn*dd ; ii++){
      ZZ[ii] = flatTable[offset + ii]; 
    }
    
  }
  
}




