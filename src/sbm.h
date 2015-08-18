/*****  sbm.h
 *****  Header File for cidR.cpp
 *****/

#ifndef SBM_R
#define SBM_R


//  Definition for SBM class
class sbm_t {
  int *YY;
  int *yyComplete;
  double *BB;
  double *PP;
  int *PPint;
  double betaPrior[2];
  double *eta;
  bool multiImpute;
 public:
  int nn, dd;
  sbm_t (int, int, int*, double*, double*, int);
  void loadTable (int, double*);
  void step ();
  void drawBB();
  void drawPP();
  void rotate();
  void imputeMissingValues();
  double LL();
  double nodeLL(int);
  void updateFlatTable(int, double*);
  void print (bool);
};

//  Function for performing the MCMC algorithm
void sbmMCMC(sbm_t mySBM, int start, int total, int burnIn, int thin,
	     int shift_size, int extend_max, double qq,
	     double *flatTable, double *logLik);


  

#endif

