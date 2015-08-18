/*****  mmsbm.h
 *****  Header File for cidR.cpp
 *****/

#ifndef MMSBM_R
#define MMSBM_R

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
  
void updateFlatTable(int iter, int nn, int dd, double *BB, 
		     double *PP, double *flatTable);
  
#endif

