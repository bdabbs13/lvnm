/*****  cid-R.h
 *****  Header File for cidR.cpp
 *****/

#ifndef CID_R
#define CID_R

extern "C" {
  /*********************************************************
   ************  Externally Available Functions ************
   ********************************************************/
  void sbm(int *iters, int *nn_t, int *kk_t, int *YY,
	   double *betaPrior, double *eta,
	   double *flatTable, int *burn_t, int *thin_t,
	   int *start_t, int *multi_t,double *logLik);
  
  void mmsbm(int *iters, int *nn_t, int *kk_t, int *YY,
	     double *betaPrior, double *alpha,
	     double *flatTable, int *burn_t, int *thin_t,
	     int *start_t, int *multi_t);
  
  
  /*****************************************************
   *************  MCMC CONTROL FUNCTIONS  **************
   ****************************************************/
  
  /**********  SBM  **********/
  
  void sbmInit(int *YY, int nn, int dd,
	       double *eta, double *betaPrior,
	       double *BB, double *PP, int *PPint,
	       int *yyComplete, double *logLik, double *flatTable,
	       int start, int multiImpute);
  
  void sbmStep(int *YY, int nn, int dd, 
	       double *eta, double *betaPrior,
	       double *BB, double *PP, int *PPint);
  
  void sbmImputeMissingValues(int nn, int dd,
			      int *YY, int *yyComplete,
			      double* BB, double *PP, int *PPint);
  
  void sbmLoadTable(int start, int nn, int dd,
		    double *BB, double *PP, double *flatTable);
  
  void sbmMCMC(int total,int burnIn, int thin,
	       int *YY,int nn,int dd,
	       double *eta,double *betaPrior, double *logLik,
	       double *BB,double *PP, int *PPint, 
	       int *yyComplete,double *flatTable,
	       int start, int multiImpute);
  
  

  /**********  MMSBM  **********/
  
  void mmsbmInit(int *YY, int nn, int dd, 
		 double *alpha, double *betaPrior,
		 double *BB, double *PP,
		 int *sendMat, int *recMat, double *flatTable,
		 int *yyComplete, int start, int multiImpute);
  
  void mmsbmStep(int *YY, int nn, int dd,
		 double *alpha, double *betaPrior,
		 double *BB, double *PP,
		 int *sendMat, int *recMat);
  
void mmsbmImputeMissingValues(int nn, int dd, int *YY,
			      int *yyComplete, double* BB,
			      int *sendMat,int *recMat);
  
  void mmsbmLoadTable(int start, int nn, int dd,
		      double *BB, double *PP, 
		      double *flatTable);
  
  void mmsbmMCMC(int total,int burnIn, int thin,
		 int *YY,int nn,int dd,
		 double *alpha,double *betaPrior,
		 double *BB, double *PP,
		 int *sendMat, int *recMat,
		 int *yyComplete, double *flatTable,
		 int start, int multiImpute);
  
  
  
  /*****************************************************
   ***************  SAMPLING FUNCTIONS  ****************
   ****************************************************/
  
  /********  GENERIC  ********/
  
  void rdirichlet(int k, double *alpha, double *x);
  
  
  /**********  SBM  **********/
  
  void sbmDrawBB(int nn, int dd, int *YY,
		 double *PP, double *BB, 
		 double *betaPrior, int *PPint);
  
  double sbmLogLikYY(int nn, int dd, int *YY, 
		     double *BB, double *PP, int *PPint);
  
  void sbmDrawPP(int nn, int dd, int *YY,
		 double *BB, double *PP,
		 double *eta, int *PPint);
  void sbmRotate(int nn, int dd, double *BB, 
		 double *PP, int *PPint);

  
  /**********  MMSBM  **********/
  
  void mmsbmDrawPP(int n, int d,
		   int *sendMat, int *recMat,
		   double *BB, double *alpha, double *PP);
  
  void mmsbmDrawZZ(int nn, int dd,int *YY,
		   double *BB, double *PP,
		   int *sendMat, int *recMat);
  
  void mmsbmDrawBB(int nn, int dd, int *YY,
		   int *sendMat, int *recMat,
		   double *BB, double *betaPrior);
  
  
  /*****************************************************
   ****************  OUTPUT FUNCTIONS  *****************
   ****************************************************/
  
  
  //  FUNCTIONS FOR PRINTING IN R
  void RprintDoubleMat(int rows, int cols, double *mat);
  
  void RprintIntMat(int rows, int cols, int *mat);
  
  //  UPDATING FLAT TABLE
  void updateFlatTable(int iter, int nn, int dd,
		       double *BB, double *PP,
		       double *flatTable);

  void shiftFlatTable(int shift_size, int nn, int dd, int total,
		      double *flatTable);

  //  HELPER FUNCTIONS 
  void getMeanVar(double *vec, int lower, int upper,
		  double *Mean_t, double *Var_t, double * Len_t);
  
  int convergenceCheck(double *logLik, int total, double qq);

  
}
#endif

