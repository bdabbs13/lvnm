/*****  lsm.h
 *****  Header File for cidR.cpp
 *****/

#ifndef LSM_BD
#define LSM_BD

extern "C" {
  /*********************************************************
   ************  Externally Available Functions ************
   ********************************************************/
  void lsm(int *iters, int *nn_t, int *kk_t, int *YY,
	   double *mhControl, double *alphaPrior, double *zzPrior,
	   double *initPrior, double *flatTable,
	   int *burn_t, int *thin, int *start_t,
	   int *multiImpute_t);


    /*****************************************************
   *************  MCMC CONTROL FUNCTIONS  **************
   ****************************************************/
  
  /**********  SBM  **********/
  
  void lsmInit(int *YY, int nn, int dd,
	       double *mhControl, double *initPrior,
	       double *alpha, double *ZZ,
	       int *yyComplete, double *dMat,
	       double *flatTable,
	       int start, int multiImpute);
  
  void lsmStep(int *YY, int nn, int dd, 
	       double *mhControl, double *alphaPrior, double *zzPrior,
	       double *alpha, double *ZZ, double *dMat);
  
  void lsmImputeMissingValues(int nn, int dd,
			      int *YY, int *yyComplete,
			      double *alpha, double *ZZ, double *dMat);
  
  void lsmLoadTable(int start, int nn, int dd,
		    double *alpha, double *ZZ,
		    double *flatTable);
  
  void lsmMCMC(int total,int burnIn, int thin,
	       int *YY,int nn,int dd,
	       double *mhControl, double *initPrior,
	       double *alphaPrior, double *zzPrior,
	       double *alpha, double *ZZ, double *dMat,
	       int *yyComplete,double *flatTable,int start,
	       int multiImpute);
  
  

    
  /*****************************************************
   ***************  SAMPLING FUNCTIONS  ****************
   ****************************************************/
  
  /********  GENERIC  ********/
  
  //void rdirichlet(int k, double *alpha, double *x);

  double logitInverse(double x);
  void distMat(int nn, int dd, double *ZZ, double *dMat);

  
  /**********  SBM  **********/
  
  void lsmDrawAlpha(int nn, int dd, int *YY,
		    double *alpha, double *ZZ, 
		    double *alphaPrior, double *mhControl, double *dMat);
  
  double lsmLogLik(int nn, int dd, int *YY, 
		   double *alpha, double *ZZ, double *dMat);
  
  void lsmDrawZZ(int nn, int dd, int *YY,
		 double *mhControl, double *zzPrior,
		 double *alpha, double *ZZ, double *dMat);
  
  /*****************************************************
   ****************  OUTPUT FUNCTIONS  *****************
   ****************************************************/
  
  
  //  FUNCTIONS FOR PRINTING IN R
  /*
  void RprintDoubleMat(int rows, int cols, double *mat);
  
  void RprintIntMat(int rows, int cols, int *mat);
  */
  //  UPDATING FLAT TABLE
  void lsmUpdateFlatTable(int iter, int nn, int dd, int *YY,
		       double *alpha, double *ZZ, double *dMat,
		       double *flatTable);
  
}
#endif

