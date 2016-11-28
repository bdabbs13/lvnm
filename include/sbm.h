/*****  sbm.h
 *****  Header File for cidR.cpp
 *****/

#ifndef SBM_R
#define SBM_R
#include <vector>

//  Definition for SBM class
class CSBM {
 public:
   CSBM (int rNodes, int aBlocks, int *adjMat, double *rPriorBlockMat,
	 double *rPriorBlockMemb, int mImpute);
   ~CSBM();

   // Loading Functions
   void RLoadSBM(double *rBlockMat,int *rBlockMemb);
   void initPPem(double);

   //  MCMC Functions
   void step();
   void imputeMissingValues();  //Need to be updated

   // EM Functions
   void iterEM();
   void getMultinomPosterior();
   void saveBlockMatOld();
   double BlockMatdiff();


   // Log-Likleihood Function
   double LogLike();


   // I/O Functions
   // Saving Functions
   void updateSBM(int iter, int *rBlockMemb, double *rBlockMat,
		  double *rPosteriorMemb);
   void updateSBM(int iter, int *rBlockMemb, double *rBlockMat,
		  double *rPosteriorMemb, double *rEta);

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
   std::vector<std::vector<int> > aAdjMat;
   std::vector<std::vector<int> > aAdjPartial;

   std::vector<std::vector<double> > aBlockMat;
   std::vector<std::vector<double> > aBlockMatInv;
   std::vector<std::vector<double> > aBlockMatOld;

   std::vector<std::vector<double> > aHitMat;
   std::vector<std::vector<double> > aMissMat;

   std::vector<std::vector<double> > aPosteriorMemb;
   std::vector<std::vector<double> > aPosteriorMembOld;

   std::vector<int> aBlockMemb;
   std::vector<double> aPriorBlockMemb;
   double aPriorBlockMat[2];

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
   void computeHitMiss();
   bool isMissing(int val) { return (val == missingVal);}

   // Log-Likelihood Functions
   double nodeLogLike(int ii);
   double tieLogLike(int yy, int sendBlock, int recBlock);

   // Log-Likelihood Marginalized over Memberships
   // Only used by EM
   double nodeLogLike_long(int ii);
   double tieLogLike_sender(int yy, int sendBlock, int receiver);
   double tieLogLike_receiver(int, int recBlock, int sender);

   // Block matrix Functions
   void logBlockMat();
   void expBlockMat();

   // Debugging Functions
   void printAdjacencyMatrix();
   void printBlockMat();
   void printPosteriorMemb();
   void printBlockMemb();
   void printPriorBlockMemb();




};

//  Function for performing the MCMC algorithm
void sbmMCMC(CSBM *mySBM, int start, int total, int burnIn, int thin,
	     int shift_size, int extend_max, double qq, //double *flatTable,
	     double *rBlockMat, int *rBlockMemb,
	     double *logLik, double *rPosteriorMemb, int verbose);

void sbmEM(CSBM *mySBM, int iter_max, double threshold,
	   double *flatTable, double *rBlockMat, int *rBlockMemb,
	   double *logLik, double *rPriorBlockMat,
	   int verbose);

void printAdjacencyMatrix();

#endif

