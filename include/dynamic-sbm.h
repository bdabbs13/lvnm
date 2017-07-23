/*****  wsbm.h
 *****  Header File for cidR.cpp
 *****/

#ifndef DYNSBM_R
#define DYNSBM_R
#include <vector>



//  Definition for WSBM class
class CDynSBM {
 public:
   CDynSBM (int rNodes, int rBlocks, int rTimes, int rTimeClasses,
	    int total, int *rTimeMap, double *rHours,
	    int mImpute);
   ~CDynSBM();


   /***** Time Hash Functions *****/
   //   int GetEquivalentTime(int t);

   // Loading Functions
   void RLoadDynSBM(int *AdjMat,
		    double *rHyperSender, double *rHyperReceiver,
		    double *rHyperBlockMat,
		    double *rSenderEffects, double *rReceiverEffects,
		    double *rBlockEffects, int *rBlockMemb,
		    double *rPosteriorMemb, int update_mmb,
		    double *rPriorSender, double *rPriorReceiver,
		    double *rPriorBlockMat, double *rPriorBlockMemb,
		    double *rLogLik);

   //  MCMC Functions
   void step();
   void DrawPriors();
   void DrawParameters();
   void DrawBlockMemb();
   //void imputeMissingValues();  //Need to be updated

   // Log-Likleihood Function
   double LogLike();

   /*****  I/O Functions  *****/

   // Saving Functions
   void Update(int iter);
   void adapt();

   //  Retrieval Functions
   int GetNodes() const { return aNodes;}
   int GetBlocks() const { return aBlocks;}
   int GetTimes() const {return aTimes;}
   int GetTimeClasses() const {return aClasses;}
   void GetBlockMat(double *rBlockMat);

   // Helper Functions
   void computeBlockMatMLE();
   void print(bool);
   void printPriors();
   void printAllWSBM(bool printNetworks);
   void printBlockMemb();
   void PrintCovariance(int tt, int ii);


 private:
   int aNodes;
   int aBlocks;
   int aTimes;
   int aClasses;
   int missingVal;
   int aTotal;
   bool update_mmb;
   double mhEpsilon;
   double mhSD;
   double mhStart;

   /***** Observed Network *****/
   std::vector<CWSBM> aWsbmList;

   /***** Time Hash Functions *****/
   std::vector<double> aHours;
   std::vector<int> aTimeMap;

   //  Block Membership Vector
   std::vector<int> aBlockMemb; // aNodes
   std::vector<std::vector<double> > aPosteriorMemb; // aNodes x aBlocks

   //  Priors
   std::vector<double> aPriorBlockMemb;
   std::vector<std::vector<double *> *> aPriorSender;
   std::vector<std::vector<double *> *> aPriorReceiver;
   std::vector< std::vector<std::vector<double *> > *> aPriorBlockMat;

   //  Adaptive Sampling Characteristics
   int covCount;
   std::vector<std::vector<double *> *> aPriorSenderCov;
   std::vector<std::vector<double *> *> aPriorReceiverCov;
   std::vector< std::vector<std::vector<double *> > *> aPriorBlockMatCov;

   //  Hyperpriors
   double aHyperSender[4];
   double aHyperReceiver[4];
   double aHyperBlockMat[4];

   /***** R Storage Pointers *****/
   double *aRPriorSender;
   double *aRPriorReceiver;
   double *aRPriorBlockMat;
   int *aRBlockMemb;
   double *aRPosteriorMemb;
   double *aRLogLike;

   /***** Flags *****/
   bool aImputeFlag;
   bool is_BlockMat_logged;


   // Loading Functions
   void LoadAdjacencyMatrices(int *AdjMat);
   void LoadHyperPriors(double *rHyperSender, double *rHyperReceiver,
   			double *rHyperBlockMat);
   void LoadParameters(double *rSenderEffects, double *rReceiverEffects,
   		       double *rBlockEffects, int *rBlockMemb,
   		       double *rPosteriorMemb, int update_mmb);
   void LoadPriors(double *rPriorSender, double *rPriorReceiver,
   		   double *rPriorBlockMat, double *rPriorBlockMemb);
   void LoadLogLike(double *rLogLik);
   void PassReferences();

   // MCMC Functions
   void DrawPriorSender();
   void DrawPriorReceiver();
   void DrawPriorBlockMat();

   /* void drawBlockMemb(); */
   /* void drawBlockMat(); */
   /* void drawReceiverEffects(); */
   /* void rotate(); */


   // Internal Loading Functions
   void LoadPriorSender();
   void LoadPriorReceiver();
   void LoadPriorBlockMat();

   /* void RLoadBlockMat(double *rBlockMat); */
   /* void RLoadBlockMemb(int *rBlockMemb); */
   /* void RLoadPosteriorMemb(double *rPosteriorMemb); */
   /* void RLoadSenderEffects(double *rSenderEffects); */
   /* void RLoadReceiverEffects(double *rReceiverEffects); */

   // Internal Updating Functions
   void UpdatePriorSender(int iter);
   void UpdatePriorReceiver(int iter);
   void UpdatePriorBlockMat(int iter);
   void UpdateBlockMemb(int iter);
   void UpdatePosteriorMemb(int iter);
   void UpdateLogLike(int iter);

   /*
   void updateBlockMat(int iter , double *rBlockMat);
   void updateBlockMemb(int iter, int *rBlockMemb); // updateMMB
   void updatePriorBlockMemb(double *rPriorBlockMat);
   void updatePosteriorMemb(int iter, double *rPosteriorMemb);
   void updateSenderEffects(int iter, double *rSenderEffects);
   void updateReceiverEffects(int iter, double *rReceiverEffects);
   */

   // Helper Functions
   void savePosteriorMembOld();
   void saveBlockMatOld();
   void computeBlockTieSums();
   void computeRowColSums();
   void computeBlockSenderEffects();
   void computeBlockReceiverEffects();
   bool isMissing(int val) { return (val == missingVal);}

   // Log-Likelihood Functions
   double nodeLogLike(int ii);
   double tieLogLike(int yy, int sendBlock, int recBlock);
   double GetTieMean(int ss, int rr);

   // Block matrix Functions
   //   void logBlockMat();
   //   void expBlockMat();

   // Debugging Functions
   void printPriorSender();
   void printPriorReceiver();
   void printPriorBlockMat();

   /* void printAdjacencyMatrix(); */
   /* void printBlockMat(); */
   /* void printPosteriorMemb(); */
   /* void printBlockMemb(); */
   /* void printPriorBlockMemb(); */
   /* void printSenderEffects(); */
   /* void printReceiverEffects(); */




};

//  Function for performing the MCMC algorithm
void dynSBMMCMC(CDynSBM *myDynSBM, int start, int total, int burnIn, int thin,
	   int shift_size, int extend_max, double qq, int verbose);

#endif

