/*****  wsbm.h
 *****  Header File for cidR.cpp
 *****/

#ifndef WSBM_R
#define WSBM_R
#include <vector>

//  Definition for WSBM class
class CWSBM {
 public:
   CWSBM (int rNodes, int aBlocks, int *adjMat, double *rPriorBlockMat,
	 double *rPriorBlockMemb, int mImpute);
   ~CWSBM();

   // Loading Functions
   void RLoadWSBM(double *rBlockMat,int *rBlockMemb);

   //  MCMC Functions
   void step();
   void imputeMissingValues();  //Need to be updated

   // EM Functions
   //   void iterEM();
   //   void getMultinomPosterior();
   //   double BlockMatdiff();


   // Log-Likleihood Function
   double LogLike();


   // I/O Functions
   // Saving Functions
   void updateWSBM(int iter, int *rBlockMemb, double *rBlockMat,
		   double *rPosteriorMemb);
   //   void updateWSBM(int iter, int *rBlockMemb, double *rBlockMat,
   //		  double *rPosteriorMemb, double *rEta);

   //  Retrieval Functions
   int GetNodes() const { return aNodes;}
   int GetBlocks() const { return aBlocks;}
   void GetBlockMat(double *rBlockMat);


   // Helper Functions
   void computeBlockMatMLE();
   void print(bool);




 private:
   int aNodes;
   int aBlocks;
   int missingVal;

   /***** Observed Network *****/
   std::vector<std::vector<int> > aAdjMat;
   std::vector<std::vector<int> > aAdjPartial;

   /***** Parameters *****/
   //  Block Matrix
   std::vector<std::vector<double> > aBlockMat;
   std::vector<std::vector<double> > aBlockMatLog;
   std::vector<std::vector<double> > aBlockMatOld;

   //  Block Membership Vector
   std::vector<int> aBlockMemb;

   //  Posterior Block Membership Probability Matrix
   std::vector<std::vector<double> > aPosteriorMemb;
   std::vector<std::vector<double> > aPosteriorMembOld;

   /***** Priors *****/
   std::vector<double> aPriorBlockMemb;
   double aPriorBlockMat[2];

   /*****  Simplified Calculations *****/
   std::vector<std::vector<int> > aBlockTieCounts;
   std::vector<std::vector<double> > aBlockTieSums;
   //   std::vector<std::vector<double> > aHitMat;
   //   std::vector<std::vector<double> > aMissMat;

   /***** Flags *****/
   bool aImputeFlag;
   bool is_BlockMat_logged;

   // MCMC Functions
   void drawBlockMemb();
   void drawBlockMat();
   void rotate();


   // Internal Loading Functions
   void RLoadBlockMat(double *rBlockMat);
   void RLoadBlockMemb(int *rBlockMemb);
   void RLoadPosteriorMemb(double *rPosteriorMemb);

   // Internal Updating Functions
   void updateBlockMat(int iter , double *rBlockMat);
   void updateBlockMemb(int iter, int *rBlockMemb); // updateMMB
   void updatePriorBlockMemb(double *rPriorBlockMat);
   void updatePosteriorMemb(int iter, double *rPosteriorMemb);

   // Helper Functions
   void savePosteriorMembOld();
   void saveBlockMatOld();
   void computeBlockTieSums();
   bool isMissing(int val) { return (val == missingVal);}

   // Log-Likelihood Functions
   double nodeLogLike(int ii);
   double tieLogLike(int yy, int sendBlock, int recBlock);
   double GetTieMean(int ss, int rr);

   // Block matrix Functions
   //   void logBlockMat();
   //   void expBlockMat();

   // Debugging Functions
   void printAdjacencyMatrix();
   void printBlockMat();
   void printPosteriorMemb();
   void printBlockMemb();
   void printPriorBlockMemb();




};

//  Function for performing the MCMC algorithm
void wsbmMCMC(CWSBM *myWSBM, int start, int total, int burnIn, int thin,
	     int shift_size, int extend_max, double qq, //double *flatTable,
	     double *rBlockMat, int *rBlockMemb,
	     double *logLik, double *rPosteriorMemb, int verbose);

//void printAdjacencyMatrix();

#endif

