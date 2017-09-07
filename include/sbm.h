/*****  sbm.h
 *****  Header File for cidR.cpp
 *****/

#ifndef SBM_R
#define SBM_R
#include <vector>

//  Definition for SBM class
class CSBM {
 public:
    CSBM (int rNodes, int aBlocks,int mImpute);
    ~CSBM();

    // Loading Functions
    void loadDataR(int *adjMat, double *rPriorBlockMat,double *rPriorBlockMemb);
    void loadStateR(double *rBlockMat, int *rBlockMemb,
		    double *rPosteriorMemb, double *rLogLik);
    void RLoadSBM(double *rBlockMat,int *rBlockMemb);


    //  Functions used by MCMC

    // Initialization Functions
    void initRandom();

    //  MCMC Functions
    void step();
    void adapt() { };

    //  Writing
    void write(int iter);
    //  Log Likelihood
    double LogLike();


    //  Functions used by EM
    void writeRPriorBlockMemb(double *rPriorBlockMat);
    void initPPem(double);

    void iterEM();
    void getMultinomPosterior();
    void saveBlockMatOld();
    double BlockMatdiff();

    //  Retrieval Functions
    int GetNodes() const { return aNodes;}
    int GetBlocks() const { return aBlocks;}
    void GetBlockMat(double *rBlockMat);

    // Helper Functions
    void computeBlockMatMLE();
    void print(bool);


 private:
    /**********  Parameters and Containers  **********/

    // Dimensional Paramters
    int aNodes;
    int aBlocks;

    // Flags
    int missingVal;
    bool aImputeFlag;
    bool is_BlockMat_logged;

    /***** Observed Network *****/
    std::vector<std::vector<int> > aAdjMat;
    std::vector<std::vector<int> > aAdjPartial;


    /*** Block Membership Containers ***/
    std::vector<int> aBlockMemb;

    std::vector<std::vector<double> > aPosteriorMemb;
    std::vector<std::vector<double> > aPosteriorMembOld;

    /***  Block Matrix Containers  ***/
    std::vector<std::vector<double> > aBlockMat;
    std::vector<std::vector<double> > aBlockMatInv;
    std::vector<std::vector<double> > aBlockMatOld;

    /***  Intermediate Containers  ***/
    std::vector<std::vector<double> > aHitMat;
    std::vector<std::vector<double> > aMissMat;

    /***** Prior Containers *****/
    std::vector<double> aPriorBlockMemb;
    double aPriorBlockMat[2];


    //  Writer Information
    enum WriterType {R_WRITER};
    WriterType aWriter;

    /*****  Pointers for R_WRITER  *****/
    double *aRBlockMat;
    int *aRBlockMemb;
    double *aRPosteriorMemb;
    double *aRLogLike;

    /**********  Internal Functions  **********/

    // loadStateR
    void RLoadBlockMat(double *rBlockMat);
    void RLoadBlockMemb(int *rBlockMemb);
    void RLoadPosteriorMemb(double *rPosteriorMemb);

    // step Functions
    void drawBlockMemb();
    void drawBlockMat();
    void imputeMissingValues();  //Need to be updated

    //  Simple Label-Switching Solution
    void rotate();


    // Log-Likelihood Functions
    double nodeLogLike(int ii);
    double tieLogLike(int yy, int sendBlock, int recBlock);

    // Log-Likelihood Marginalized over Memberships
    // Only used by EM
    double nodeLogLike_long(int ii);
    double tieLogLike_sender(int yy, int sendBlock, int receiver);
    double tieLogLike_receiver(int, int recBlock, int sender);

    // Writing Functions for R
    void writeR(int iter);

    // Individual writeR Functions
    void writeRBlockMat(int iter);
    void writeRBlockMemb(int iter);
    void writeRPosteriorMemb(int iter);

    // Helper Functions
    void savePosteriorMembOld();
    void computeHitMiss();
    bool isMissing(int val) { return (val == missingVal);}

    // Block matrix Functions
    void logBlockMat();
    void expBlockMat();

    // Printing Functions
    void printAdjacencyMatrix();
    void printBlockMat();
    void printPosteriorMemb();
    void printBlockMemb();
    void printPriorBlockMemb();

};

//  Function for performing the EM algorithm
void sbmEM(CSBM *mySBM, int iter_max, double threshold,
	   double *flatTable, double *rBlockMat, int *rBlockMemb,
	   double *logLik, double *rPriorBlockMat,
	   int verbose);

// void printAdjacencyMatrix();

#endif

