/*****  wsbm.h
 *****  Header File for cidR.cpp
 *****/

#ifndef DYNSBM_R
#define DYNSBM_R
#include <vector>



//  Definition for WSBM class
class CDynSBM {
 public:
   /* CDynSBM (int rNodes, int rBlocks, int rTimes, int rTimeClasses, */
   /* 	    int total, int *rTimeMap, double *rHours, */
   /* 	    int mImpute); */
   CDynSBM (int rNodes, int rBlocks, int rTimes, int rTimeClasses,
	    int total, int mImpute);
   ~CDynSBM();


   /***** Time Hash Functions *****/
   //   int GetEquivalentTime(int t);

   // Loading Functions
   void LoadDataR(int *AdjMat, int *rTimeMap, double *rHours,
		  double *rHyperSender, double *rHyperReceiver,
		  double *rHyperBlockMat);

   void LoadStateR(double *rSenderEffects, double *rReceiverEffects,
		   double *rBlockEffects, int *rBlockMemb,
		   double *rPosteriorMemb, int update_mmb,
		   double *rPriorSender, double *rPriorReceiver,
		   double *rPriorBlockMat, double *rPriorBlockMemb,
		   double *rLogLik);

   //  Retrieval Functions
   int GetNodes() const { return aNodes;}
   int GetBlocks() const { return aBlocks;}
   int GetTimes() const {return aTimes;}
   int GetTimeClasses() const {return aClasses;}


   /*****  MCMC Interface Functions  *****/
   //  MCMC Functions
   void step();
   void adapt();

   // Writing Functions
   void write(int iter);

   // Log-Likleihood Function
   double LogLike();


   // Printing Functions (Primarily for Debugging)
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

   enum WriterType { R_WRITER, TEXT_WRITER};
   WriterType aWriter;

   /***** Observed Network *****/
   std::vector<CWSBM *> aWsbmList;

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
   void LoadTime(int *rTimeMap, double *rHours);
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
   void DrawPriors();
   void DrawParameters();
   void DrawBlockMemb();
   //void imputeMissingValues();  //Need to be updated

   void DrawPriorSender();
   void DrawPriorReceiver();
   void DrawPriorBlockMat();

   // Writing Functions
   void writeR(int iter);

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
   void writeRPriorSender(int iter);
   void writeRPriorReceiver(int iter);
   void writeRPriorBlockMat(int iter);
   void writeRBlockMemb(int iter);
   void writeRPosteriorMemb(int iter);
   void writeRLogLike(int iter);

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

