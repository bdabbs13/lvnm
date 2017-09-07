/*****  wsbm.h
 *****  Header File for cidR.cpp
 *****/

#ifndef WSBM_R
#define WSBM_R
#include <vector>

//  Definition for WSBM class
class CWSBM {
 public:
    CWSBM (int rNodes, int rBlocks, int mImpute);
    ~CWSBM();

    // Loading Functions
    void loadDataR(int *adjMat, double rHours,
		   double *rPriorSender, double *rPriorReceiver,
		   double *rPriorBlockMat, double * rPriorBlockMemb);
    // Load State with Block Membership
    void loadStateR(double *rBlockMat,int *rBlockMemb,
		    double *rSenderEffects, double *rReceiverEffects,
		    double *rPosteriorMemb, double *rLogLik);

    /***  Functions Reading from DynSBM  ***/
    // Load State without Block Membership
    void loadStateR(double *rBlockMat,
		    double *rSenderEffects, double *rReceiverEffects);

    // Passing References to Priors from dynSBM
    void LoadReferences(double hours,
			std::vector<double *> *dPriorSender,
			std::vector<double *> *dPriorReceiver,
			std::vector<std::vector<double *> > *dPriorBlockMat,
			std::vector<int> *dBlockMemb);

    void LoadAdjacencyMatrix(int *AdjMat);

    void initRandom();

    /*****  MCMC Interface Functions  *****/
    //  MCMC Functions
    void step();
    void adapt() { };

    // Writing Function
    void write(int iter);

    // Log-Likleihood Functions
    double LogLike();
    double nodeLogLike(int ii);


    //  Retrieval Functions
    int GetNodes() const { return aNodes;}
    int GetBlocks() const { return aBlocks;}
    void GetBlockMat(double *rBlockMat);

    //  Retrieval Functions Used by DynSBM
    double GetSenderEffect(int node) const {return aSenderEffects[node];}
    double GetReceiverEffect(int node) const {return aReceiverEffects[node];}
    double GetBlockMatEffect(int ll, int kk) const {return aBlockMat[ll][kk];}

    // Helper and Printing Functions
    void computeBlockMatMLE();
    void print(bool);


 private:
    /**********  Parameters and Containers **********/

    // Dimensional Parameters
    int aNodes;
    int aBlocks;

    // Flags
    int missingVal;  // Value Indicating Missing
    bool aImputeFlag;

    /***** Observed Network *****/
    std::vector<std::vector<int> > aAdjMat; // aNodes x aNodes
    std::vector<std::vector<int> > aAdjPartial; // aNodes x aNodes

    //  Values Only Needed to Be Computed Once
    std::vector<int> aAdjMatRowSums; // aNodes
    std::vector<int> aAdjMatColSums; // aNodes

    /***  Time Parameter ***/
    double aHours;

    //  Block Membership Params
    std::vector<int> *aBlockMemb; // aNodes
    std::vector<std::vector<double> > aPosteriorMemb; // aNodes x aBlocks

    /*****  Parameter Containers  *****/
    std::vector<std::vector<double> > aBlockMat; // aBlocks x aBlocks
    //  Used in rotate (May be obsolete?)
    std::vector<std::vector<double> > aBlockMatOld; // aBlocks x aBlocks

    //  Sender and Receiver Effects
    std::vector<double> aSenderEffects; // aNodes
    std::vector<double> aReceiverEffects; // aNodes


    /*****  Prior Containers  *****/
    std::vector<double> aPriorBlockMemb;
    std::vector<double *> *aPriorSender;
    std::vector<double *> *aPriorReceiver;
    std::vector<std::vector<double *> > *aPriorBlockMat;

    bool aPriorOwner;
    double aLocalPriorBlockMat[2];
    double aLocalPriorSender[2];
    double aLocalPriorReceiver[2];


    //  Step Type Information
    enum StepType {FULL_STEP,PARTIAL_STEP};
    StepType aStepType;

    //  Writer Information
    enum WriterType {R_WRITER, DYN_R_WRITER};
    WriterType aWriter;

    /*****  Pointers for R_WRITER and DYN_R_WRITER   *****/
    double *aRBlockMat;
    int *aRBlockMemb;
    double *aRSenderEffects;
    double *aRReceiverEffects;
    double *aRPosteriorMemb;
    double *aRLogLike;


    /*****  Simplified Calculations *****/
    std::vector<std::vector<double> > aBlockTieCounts;
    std::vector<std::vector<double> > aBlockTieSums;
    std::vector<double> aBlockSenderEffects;
    std::vector<double> aBlockReceiverEffects;


    /**********  Internal Functions  **********/

    // loadStateR Functions
    void RLoadBlockMat(double *rBlockMat);
    void RLoadBlockMemb(int *rBlockMemb);
    void RLoadPosteriorMemb(double *rPosteriorMemb);
    void RLoadSenderEffects(double *rSenderEffects);
    void RLoadReceiverEffects(double *rReceiverEffects);
    void deallocatePriors();


    // MCMC Functions
    void fullStep();      // Step Updating All Parameters
    void partialStep();  // Step Without Updating Block Membership

    // step Functions
    void drawBlockMat();
    void drawBlockMemb();
    void drawSenderEffects();
    void drawReceiverEffects();

    void imputeMissingValues();  // Need to be updated


    // Normalizing Sender/Receiver Effects Within Blocks
    void normalizeReceiverEffects();
    void normalizeSenderEffects();
    // Simple Label-Switching Solution
    void rotate();

    // Log-Likelihood Functions
    double tieLogLike(int yy, int sendBlock, int recBlock);
    double GetTieMean(int ss, int rr);



    // Writing Functions for R
    void writeR(int iter);
    void writeDynR(int iter);

    // Internal Updating Functions
    void writeRBlockMat(int iter);
    void writeRBlockMemb(int iter);
    void writeRPosteriorMemb(int iter);
    void writeRSenderEffects(int iter);
    void writeRReceiverEffects(int iter);


    // Helper Functions
    void saveBlockMatOld();
    void computeBlockTieSums();
    void computeRowColSums();
    void computeBlockSenderEffects();
    void computeBlockReceiverEffects();
    bool isMissing(int val) { return (val == missingVal);}


    // Printing Functions
    void printAdjacencyMatrix();
    void printBlockMat();
    void printPosteriorMemb();
    void printBlockMemb();
    void printPriorBlockMemb();
    void printSenderEffects();
    void printReceiverEffects();




};

//  Function for performing the MCMC algorithm

/* void wsbmMCMC(CWSBM *myWSBM, int start, int total, int burnIn, int thin, */
/* 	      int shift_size, int extend_max, double qq, int verbose); */


/* void wsbmMCMC(CWSBM *myWSBM, int start, int total, int burnIn, int thin, */
/* 	      int shift_size, int extend_max, double qq, //double *flatTable, */
/* 	      double *rBlockMat, int *rBlockMemb, */
/* 	      double *rSenderEffects, double *rReceiverEffects, */
/* 	      double *logLik, double *rPosteriorMemb, int verbose); */

//void printAdjacencyMatrix();

#endif

