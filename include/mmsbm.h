/*****  mmsbm.h
 *****  Header File for cidR.cpp
 *****/

#ifndef MMSBM_R
#define MMSBM_R


class mmsbm_t {
  int *YY;
  int *yyComplete;
  double *BB;
  double *PP;
  int *sendMat;
  int *recMat;
  double betaPrior[2];
  double *alpha;
  bool multiImpute;
 public:
  int nn, dd;
  mmsbm_t (int, int, int*, double*, double*, int);
  ~mmsbm_t ();
  void loadTable (int, double*);
  void step ();
  void drawBB();
  void drawPP();
  void drawZZ();
  void rotate();
  void imputeMissingValues();
  double LL();
  void updateFlatTable(int, double*);
  void print (bool);
};



/**********  MMSBM  **********/
  
void mmsbmMCMC(mmsbm_t *myMMSBM, int start, int total, int burnIn, int thin,
	       int shift_size, int extend_max, double qq,
	       double *flatTable, double *logLik, int verbose);
  
#endif

