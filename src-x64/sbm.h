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
  double *BB_old;
  double *hit;
  double *miss;
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

  // I/O Functions
  void loadTable (int, double*);
  void loadSBM(double*,int*);
  void loadBB(double *);
  void loadHH(double *);
  void loadPPint(int *);
  void initPPem(double);
  void updateFlatTable(int, double*);
  void updateBB(int, double*);
  void updateMMB(int, int*);
  void updateEta(double *eta_out);
  void updatePosteriorMat(int, double*);
  
  //  MCMC Functions
  void step ();
  void drawBB();
  void drawPP();
  void rotate();
  void imputeMissingValues();
  
  // EM Functions
  void iterEM();
  void getMultinomPosterior();
  
  // Helper Functions
  void computeHitMiss ();
  void computeBBmle();
  
  // Log-Likelihood Functions
  double LL();
  double nodeLL(int);
  double tieLL(int, int, int);
  double nodeLL_long(int);
  double tieLL_sender(int, int, int);
  double tieLL_receiver(int, int, int);

  // BB management
  void getBB(double *BB_out);
  void saveBB_old();
  double BBdiff();
  void logBB();
  void expBB();
  
  void print (bool);
};

//  Function for performing the MCMC algorithm
void sbmMCMC(sbm_t *mySBM, int start, int total, int burnIn, int thin,
	     int shift_size, int extend_max, double qq, //double *flatTable,
	     double *BBout, int *MMBout,
	     double *logLik, double *postMat, int verbose);

void sbmEM(sbm_t *mySBM, int iter_max, double threshold,
	   double *flatTable, double *BBout, int *MMBout,
	   double *logLik, double *eta,
	   int verbose);

  

#endif

