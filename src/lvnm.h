/*****  lvnm.h
 *****  Header File for lvnm.cpp:  Defines external functions
 *****/

#ifndef CHANGE_R
#define CHANGE_R


extern "C" {
  
  //  Fits a stochastic block model
  void sbm(int *iters, int *nn_t, int *kk_t, int *YY,
	   double *betaPrior, double *eta,
	   double *flatTable, int *burn_t, int *thin_t,
	   int *start_t, int *multi_t,double *logLik);

  //  Fits a mixed membership stochastic block model
  void mmsbm(int *iters, int *nn_t, int *kk_t, int *YY,
	     double *betaPrior, double *alpha,
	     double *flatTable, int *burn_t, int *thin_t,
	     int *start_t,int *multi_t, double *logLik);
  
  
}
#endif

