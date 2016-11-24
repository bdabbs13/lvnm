/*****  sbm.h
 *****  Header File for cidR.cpp
 *****/

#ifndef SBM_R
#define SBM_R
#include <vector>

//  Definition for SBM class
class CSBM {
 public:
   CSBM (int rNodes, int aBlocks, int *adjMat, double *betaP, double *etaP,
	 int mImpute);
   ~CSBM();

   //  MCMC Functions
   void step();
   void drawBlockMemb();
   void drawBlockMat();
   void rotate();
   void imputeMissingValues();

   // EM Functions
   void iterEM();
   void getMultinomPosterior();

   // Log-Likleihood Functions
   double LogLike();


   // I/O Functions
   void RLoadSBM(double *rBlockMat,int *rBlockMemb);
   void initPPem(double);

   void updateSBM(int iter, int *rBlockMemb, double *rBlockMat,
		  double *rPosteriorMemb);
   void updateSBM(int iter, int *rBlockMemb, double *rBlockMat,
		  double *rPosteriorMemb, double *rEta);

   int GetNodes() const { return aNodes;}
   int GetBlocks() const { return aBlocks;}

   // Helper Functions
   void computeHitMiss();
   void computeBlockMatMLE();

   // Log-Likelihood Functions
   double nodeLogLike(int ii);
   double tieLogLike(int yy, int sendBlock, int recBlock);

   // Log-Likelihood Marginalized over Memberships
   // Only used by EM
   double nodeLogLike_long(int ii);
   double tieLogLike_sender(int yy, int sendBlock, int receiver);
   double tieLogLike_receiver(int, int recBlock, int sender);

   // Block matrix Functions
   void getBlockMat(double *rBlockMat);
   void saveBlockMatOld();
   double BlockMatdiff();
   void logBlockMat();
   void expBlockMat();

   void savePosteriorMembOld();

   void print (bool);
   void printAdjacencyMatrix();
   void printBlockMat();
   void printPosteriorMemb();
   void printBlockMemb();


 private:
   int aNodes;
   int aBlocks;

   std::vector<std::vector<int> > aAdjMat;
   std::vector<std::vector<int> > aAdjPartial;

   std::vector<std::vector<double> > aBlockMat;
   std::vector<std::vector<double> > aBlockMatInv;
   std::vector<std::vector<double> > aBlockMatOld;

   std::vector<std::vector<double> > aHitMat;
   std::vector<std::vector<double> > aMissMat;

   std::vector<std::vector<double> > aPosteriorMemb;
   std::vector<std::vector<double> > aPosteriorMembOld;

   //int *aBlockMemb;
   std::vector<int> aBlockMemb;
   double betaPrior[2];
   double *eta;
   //   std::vector<double> eta;
   bool multiImpute;
   bool imputeFlag;
   bool is_BlockMat_logged;



   // Internal Loading Functions
   void RLoadBlockMat(double *rBlockMat);
   void RLoadBlockMemb(int *rBlockMemb);
   void RLoadPosteriorMemb(double *rPosteriorMemb);

   // Internal Updating Functions
   void updateBlockMat(int iter , double *rBlockMat);
   void updateBlockMemb(int iter, int *rBlockMemb); // updateMMB
   void updateEta(double *rEta);
   void updatePosteriorMemb(int iter, double *rPosteriorMemb);



};

//  Function for performing the MCMC algorithm
void sbmMCMC(CSBM *mySBM, int start, int total, int burnIn, int thin,
	     int shift_size, int extend_max, double qq, //double *flatTable,
	     double *rBlockMat, int *rBlockMemb,
	     double *logLik, double *rPosteriorMemb, int verbose);

void sbmEM(CSBM *mySBM, int iter_max, double threshold,
	   double *flatTable, double *rBlockMat, int *rBlockMemb,
	   double *logLik, double *rEta,
	   int verbose);

void printAdjacencyMatrix();

#endif

