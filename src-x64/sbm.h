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
  double *BB_inv;
  double *PP;
  double *HH;
  int *PPint;
  double betaPrior[2];
  double *eta;
  bool multiImpute;
  bool is_BB_logged;
 public:
  int nn, dd;
  sbm_t (int, int, int*, double*, double*, int);
  ~sbm_t();
  void loadTable (int, double*);
  void step ();
  void drawBB();
  void drawPP();
  void rotate();
  void imputeMissingValues();
  double LL();
  double nodeLL(int);
  double nodeLL_long(int);
  void getMultinomPosterior();
  void iterEM();
  void updateFlatTable(int, double*);
  void updatePosteriorMat(int, double*);
  void print (bool);
  void getBB(double *BB_out);
  void geteta(double *eta_out);
  double BBdiff(double *BB_old);
  void logBB();
  void expBB();
};

//  Function for performing the MCMC algorithm
void sbmMCMC(sbm_t *mySBM, int start, int total, int burnIn, int thin,
	     int shift_size, int extend_max, double qq,
	     double *flatTable, double *logLik, double *postMat, int verbose);

void sbmEM(sbm_t *mySBM, int iter_max, double threshold,
	   double *flatTable, double *logLik, double *eta,
	   int verbose);

  

#endif

