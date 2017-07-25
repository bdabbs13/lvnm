/*****  lvnm.h
 *****  Header File for lvnm.cpp:  Defines external functions
 *****/

#ifndef CHANGE_R
#define CHANGE_R


extern "C" {

    //  Fits a stochastic block model
    void sbm(int *iters, int *nn_t, int *kk_t, int *YY,
	     double *betaPrior, double *eta, //double *flatTable,
	     double *BBout, int *MMBout,
	     int *burn_t, int *thin_t,
	     int *start_t, int *multi_t,double *logLik,
	     int *extend_max_t, int *shift_t, double *qq_t,
	     double *postMat, int *verbose_t);

    //  Performs EM algorithm and returns final PI and BB
    void sbmEMout(int *iter_max_t, int *nn_t, int *kk_t, int *YY,
		  double *eta,
		  double *HHout, double *BBout, int *MMBout,
		  double *threshold_t,
		  double *logLik,int *verbose_t);

    //  Runs MCMC Algorithm for Weighted SBM
    void wsbm(int *iters, int *nn_t, int *kk_t, int *YY,
	      double *betaPrior, double *eta,
	      double *rPriorSender, double *rPriorReceiver,
	      double *rBlockMat, int *rBlockMemb,
	      double *rSenderEffects, double *rReceiverEffects,
	      int *burn_t, int *thin_t,
	      int *start_t, int *multi_t,double *logLik,
	      int *extend_max_t, int *shift_t, double *qq_t,
	      double *postMat, double *rHours, int *verbose_t);


    //  Fits a mixed membership stochastic block model
    void mmsbm(int *iters, int *nn_t, int *kk_t, int *YY,
	       double *betaPrior, double *alpha,
	       double *flatTable, int *burn_t, int *thin_t,
	       int *start_t,int *multi_t, double *logLik,
	       int *extend_max_t, int *shift_t, double *qq_t, int *verbose_t);

    void getBBmle(int *nn_t, int *kk_t, int *YY, double *BB, int *mmb);

    /* void dynsbm(int *iters, int *burnin_t, int *thin_t, int *start_t, */
    /* 		int *extend_max_t, int *shift_t, int *qq_t, */
		
    /* 		int *nn_t, int *YY, int *kk_t, int *multi_t, */
    /* 		int *TT_t, int *ee_t, int *rTimeMap, double *rHours, */

    /* 		double *rHyperSender, double *rHyperReceiver, */
    /* 		double *rHyperBlockMat, double *rPriorBlockMemb, */

    /* 		double *rPriorSender, double *rPriorReceiver, */
    /* 		double *rPriorBlockMat, int *rBlockMemb, */

    /* 		double *rSenderEffects, double *rReceiverEffects, */
    /* 		double *rBlockEffects, */

    /* 		double *rLogLik, double *rPosteriorMemb, int *verbose_t, int *update_mmb_t); */

    void dynsbm(int *nn_t, int *kk_t, int *TT_t, int *ee_t, int *multi_t,
		// Data Values
		int *YY, int *rTimeMap, double *rHours,
		// Hyper Prior Parameter Storage
		double *rHyperSender, double *rHyperReceiver,
		double *rHyperBlockMat, double *rPriorBlockMemb,
		// Prior Parameter Storage
		double *rPriorSender, double *rPriorReceiver,
		double *rPriorBlockMat, int *rBlockMemb,
		// Parameter Storage
		double *rSenderEffects, double *rReceiverEffects,
		double *rBlockEffects,
		// Auxiliary Statistic Storage
		double *rLogLik, double *rPosteriorMemb,
		// Model Flags
		int *verbose_t, int *update_mmb_t,
		// MCMC Control Parameters
		int *iters, int *burnin_t, int *thin_t, int *start_t,
		int *extend_max_t, int *shift_t, int *qq_t);

}
#endif

